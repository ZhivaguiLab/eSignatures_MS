import os
import sys
import glob
import numpy as np
import pandas as pd

# Threading limits have to be set BEFORE importing numpy or joblib.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_v, '1')
os.environ.setdefault('LOKY_MAX_CPU_COUNT', '2')


# SigProfilerAssignment hardcodes n_jobs=64 and passes it straight to joblib.
# This overrides LOKY_MAX_CPU_COUNT and will crash the node by running out 
# of file handles (Errno 23). We patch joblib.Parallel right here before SPA 
# imports it to force a strict worker cap.
import joblib

_MAX_WORKERS = int(os.environ.get('SPA_MAX_WORKERS', '4'))
_OriginalParallel = joblib.Parallel

class _CappedParallel(_OriginalParallel):
    def __init__(self, *args, **kwargs):
        if 'n_jobs' in kwargs and kwargs['n_jobs'] not in (None, 1):
            kwargs['n_jobs'] = _MAX_WORKERS
        elif args and isinstance(args[0], int) and args[0] not in (0, 1):
            args = (_MAX_WORKERS,) + args[1:]
        kwargs.setdefault('batch_size', 16)
        super().__init__(*args, **kwargs)

joblib.Parallel = _CappedParallel
joblib.parallel.Parallel = _CappedParallel
print(f"joblib.Parallel successfully capped at {_MAX_WORKERS} workers")

from SigProfilerAssignment import Analyzer as Analyze

# Making sure the correct MutationType column
BASE_DIR    = '/tscc/lustre/restricted/alexandrov-ddn/users/z8jiang/signature_stress_test/Maria_Z_eSignature'
ESS_PATH    = os.path.join(BASE_DIR, 'eSS_ref49_hClust_relabelled_to_old_numbering.txt')
COSMIC_PATH = os.path.join(BASE_DIR, 'COSMIC_v3.6_SBS_GRCh38.txt')
SCRATCH     = os.path.join(BASE_DIR, 'null_run')
RESULT_DIR  = os.path.join(BASE_DIR, 'null_cosines')

MUT_CLASSES  = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
CANONICAL_96 = [f + '[' + s + ']' + t for s in MUT_CLASSES for f in 'ACGT' for t in 'ACGT']
CLASS_BLOCKS = {s: np.array([i for i, m in enumerate(CANONICAL_96) if m[2:5] == s]) 
                for s in MUT_CLASSES}

def canonicalize(in_path, out_path, label):
    """
    Cleans up a 96-channel matrix and writes it out in the exact order SPA expects.
    Prints any formatting repairs it makes so we can track input issues in the logs.
    """
    raw = pd.read_csv(in_path, sep='\t', dtype=str)
    first = raw.columns[0]

    print(f"[{label}] Reading {in_path} | Shape: {raw.shape}")

    # Fix 1: Strip whitespace
    stripped = raw[first].astype(str).str.strip()
    n_ws = int((stripped != raw[first].astype(str)).sum())
    if n_ws:
        print(f"[{label}] Fixed {n_ws} row labels with extra whitespace")
    raw[first] = stripped

    # Fix 2: Drop blank rows
    blank = raw[first].isin(['', 'nan', 'None'])
    if blank.any():
        print(f"[{label}] Dropped {int(blank.sum())} blank rows")
        raw = raw[~blank]

    # Fix 3: Handle duplicate row labels
    dup_mask = raw[first].duplicated(keep=False)
    if dup_mask.any():
        offenders = sorted(set(raw.loc[dup_mask, first]))
        print(f"[{label}] Found duplicated labels: {offenders}. Keeping the first occurrence.")
        raw = raw.drop_duplicates(subset=first, keep='first')

    # Duplicate columns will crash the NNLS fit later, so catch it now
    if raw.columns.duplicated().any():
        bad = sorted(set(raw.columns[raw.columns.duplicated()]))
        raise ValueError(f"[{label}] Duplicate columns found in {in_path}: {bad}")

    df = raw.set_index(first)
    df.index.name = 'MutationType'
    df = df.apply(pd.to_numeric, errors='coerce')

    missing = [m for m in CANONICAL_96 if m not in df.index]
    extra   = [m for m in df.index if m not in CANONICAL_96]
    
    if extra:
        print(f"[{label}] Dropped {len(extra)} non-standard rows (e.g., {extra[:3]})")
    if missing:
        raise ValueError(f"[{label}] Missing {len(missing)} canonical channels (e.g., {missing[:3]}). Fix the source file.")

    df = df.reindex(CANONICAL_96)

    n_nan = int(df.isna().values.sum())
    if n_nan:
        print(f"[{label}] Replaced {n_nan} unreadable values with 0.0")
        df = df.fillna(0.0)

    df.to_csv(out_path, sep='\t')
    return df

def permute_library(matrix, rng):
    """
    Applies a single random permutation across all profiles in the matrix.
    This shuffles the channels while preserving all pairwise similarities 
    between the signatures in the library.
    """
    n_channels, _ = matrix.shape
    perm = np.arange(n_channels)
    
    for idx in CLASS_BLOCKS.values():
        perm[idx] = rng.permutation(idx)
        
    return matrix[perm, :]

def read_cosines(output_dir):
    """
    Pulls the non-negative least squares reconstruction cosines.
    We read the De_Novo_map file instead of Samples_Stats to bypass SPA's 
    biological sparsity penalties (which fail when total mutations = 1).
    """
    hits = glob.glob(os.path.join(output_dir, '**', 'De_Novo_map_to_COSMIC_*.csv'), recursive=True)
    if not hits:
        raise FileNotFoundError(f"Could not find the De_Novo_map CSV in {output_dir}")
        
    path = hits[0]
    mapping = pd.read_csv(path)
    
    # Standardize column names
    mapping.columns = mapping.columns.str.strip()
    out = mapping[['De novo extracted', 'Cosine Similarity']].copy()
    out = out.rename(columns={'De novo extracted': 'cosmic_signature', 'Cosine Similarity': 'cosine_similarity'})
    
    # Clean up the prefix so it matches our target names exactly
    out['cosmic_signature'] = out['cosmic_signature'].astype(str).str.replace('Signature 96-', '').str.strip()
    return out

def check_replicate_spread():
    """
    Sanity check: compares finished replicates to ensure the permutation is 
    actually working. If a target never moves across shuffles, it's just 
    riding on class composition. Run via: python null_library_fit.py --check
    """
    files = sorted(glob.glob(os.path.join(RESULT_DIR, 'null_*.csv')))
    if len(files) < 2:
        print(f"Need at least 2 finished replicates to compare, but only found {len(files)}")
        return

    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    print(f"Loaded {len(files)} replicates\n")

    wide = df.pivot_table(index='cosmic_signature', columns='replicate', values='cosine_similarity', aggfunc='first')

    spread = wide.std(axis=1)
    n_unique = wide.nunique(axis=1)
    n_reps = wide.shape[1]

    frozen = spread[spread < 1e-9].index.tolist()
    print(f"Targets identical across all {n_reps} replicates: {len(frozen)} / {len(wide)}")
    if frozen:
        print("  " + ", ".join(frozen[:25]))
        if len(frozen) > 25:
            print(f"  ...and {len(frozen) - 25} more")

    print(f"\nTargets with only 1 or 2 distinct values: {int((n_unique <= 2).sum())} / {len(wide)}")
    print(f"Median distinct values per target: {n_unique.median():.1f} (out of {n_reps})\n")

    print("--- 10 targets that moved the LEAST ---")
    print(spread.sort_values().head(10).round(5).to_string())
    print("\n--- 10 targets that moved the MOST ---")
    print(spread.sort_values(ascending=False).head(10).round(5).to_string())

def main(replicate):
    os.makedirs(RESULT_DIR, exist_ok=True)
    tag = 'observed' if replicate == 0 else f'null_{replicate:04d}'

    # Check if this replicate already ran successfully so array restarts only fill gaps
    done_file = os.path.join(RESULT_DIR, f'{tag}.csv')
    if os.path.exists(done_file) and os.path.getsize(done_file) > 0:
        print(f"{tag} already exists and is non-empty. Skipping.")
        return

    run_dir = os.path.join(SCRATCH, tag)
    os.makedirs(run_dir, exist_ok=True)

    cosmic_path = os.path.join(run_dir, 'cosmic_canonical.SBS96.txt')
    canonicalize(COSMIC_PATH, cosmic_path, 'COSMIC')

    ess_path = os.path.join(run_dir, 'ess_canonical.SBS96.txt')
    ess = canonicalize(ESS_PATH, ess_path, 'eSS')

    if replicate == 0:
        library = ess.values.astype(float)
    else:
        rng = np.random.default_rng(seed=replicate)
        library = permute_library(ess.values.astype(float), rng)

    library_path = os.path.join(run_dir, 'library.SBS96.txt')
    pd.DataFrame(library, index=ess.index, columns=ess.columns).to_csv(library_path, sep='\t')

    out_dir = os.path.join(run_dir, 'fit')

    # We use the COSMIC file as both the samples and the signatures to fit.
    # The 49-profile eSS library (real or shuffled) serves as the database.
    Analyze.decompose_fit(
        samples=cosmic_path,
        output=out_dir,
        input_type="matrix",
        signatures=cosmic_path,
        signature_database=library_path,
        new_signature_thresh_hold=0.0, # Force decomposition. Default 0.8 rejects too many targets
        collapse_to_SBS96=False,
        make_plots=False,
        export_probabilities=False,
        verbose=False,
        connected_sigs=False
    )

    cosines = read_cosines(out_dir)
    cosines['replicate'] = replicate
    cosines.to_csv(done_file, index=False)
    print(f"Success: Wrote {len(cosines)} cosines for {tag}")

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        check_replicate_spread()
    else:
        rep_num = int(sys.argv[1]) if len(sys.argv) > 1 else 0
        main(rep_num)