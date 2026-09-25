"""
utils/naming.py

Shared sample-name parsing utilities used across the eSignatures pipeline.
All scripts should import from here rather than maintaining local copies.

Functions
---------
load_acronym_mapping      - Load compound → acronym TSV into a dict
apply_acronym_mapping     - Map a sample prefix using that dict (superset version)
extract_sample_prefix     - Strip trailing numeric replicate suffix
remove_replicate_suffix   - Strip trailing numeric suffix from a compound token
apply_celegans_collapsing - Collapse C. elegans sample prefixes to treatment level
parse_sample_into_species_model_compound
                          - Top-level parser: sample name → (species_model, compound)
"""

import os
import pandas as pd


# ---------------------------------------------------------------------------
# Acronym mapping
# ---------------------------------------------------------------------------

def load_acronym_mapping(abbreviation_file_path):
    """
    Load compound → acronym mapping from a TSV file.

    The file must have 'compound' and 'acronym' columns.
    Stores entries under both original capitalisation and lowercase keys,
    and also with underscores replaced by spaces, to support flexible lookup.

    Parameters
    ----------
    abbreviation_file_path : str or None
        Path to the TSV abbreviation file.

    Returns
    -------
    dict
        {compound_string: acronym_string, ...}
        Empty dict if file is missing or malformed.
    """
    if not abbreviation_file_path or not os.path.exists(abbreviation_file_path):
        print(f"Acronym mapping file not found: {abbreviation_file_path}")
        return {}

    try:
        df = pd.read_csv(abbreviation_file_path, sep='\t')

        if 'compound' not in df.columns or 'acronym' not in df.columns:
            print(f"Warning: expected 'compound' and 'acronym' columns in "
                  f"{abbreviation_file_path}")
            return {}

        mapping = {}
        for _, row in df.iterrows():
            compound = str(row['compound']).strip()
            acronym  = str(row['acronym']).strip()

            mapping[compound]                          = acronym
            mapping[compound.lower()]                  = acronym
            mapping[compound.replace('_', ' ')]        = acronym
            mapping[compound.replace('_', ' ').lower()] = acronym

        print(f"Loaded {len(df)} compound-to-acronym mappings "
              f"from {abbreviation_file_path}")
        return mapping

    except Exception as e:
        print(f"Error loading acronym mapping file: {e}")
        return {}


def apply_acronym_mapping(prefix, acronym_mapping):
    """
    Apply acronym mapping to a sample prefix using exact matches only.

    This is the superset version: it handles all special cases present across
    the summary and matrix scripts, including organoid models and combined
    compound exposures.

    Parameters
    ----------
    prefix : str
        Sample prefix string, e.g. "Mouse_Liver_Aflatoxin_B1".
    acronym_mapping : dict
        Mapping returned by load_acronym_mapping().

    Returns
    -------
    str
        Mapped prefix if an exact match is found, otherwise original prefix.
    """
    parts = prefix.split("_")

    # --- Determine prefix_model and compound token ---
    # Prefixes where the model occupies the first THREE parts
    three_part_model_prefixes = [
        "Mouse_Lymph_Node_9,10-dimethyl-1,2-benzanthracene",
        "Mouse_Esophagus_HPV_4-Nitroquinoline_N-oxide",
        "Mouse_Liver_HBV_Aflatoxin_B1",
        "Mouse_Bone_Marrow_5-aza-4-thio-2-deoxycytidine",
        "Mouse_Lymph_Node_5-aza-4-thio-2-deoxycytidine",
    ]

    if prefix in three_part_model_prefixes:
        prefix_model = "_".join(parts[0:3])
        compound     = "_".join(parts[3:])

    elif prefix == "Human_Breast_organoids_GammaRay":
        prefix_model = "Human_Breast\norganoids"
        compound     = parts[-1]

    elif prefix == "Human_Colon_organoids_GammaRay":
        prefix_model = "Human_Colon\norganoids"
        compound     = parts[-1]

    else:
        prefix_model = "_".join(parts[0:2])
        compound     = "_".join(parts[2:])

    if not acronym_mapping:
        return prefix

    # --- Standard exact match ---
    if compound in acronym_mapping:
        new_prefix = f"{prefix_model}_{acronym_mapping[compound]}"
        print(f"  Mapping: '{prefix}' → '{new_prefix}'")
        return new_prefix

    # --- Combined compound: Sodium arsenite + simulated solar radiation ---
    if compound == "Sodium_arsenite_simulated_solar_radiation":
        sa  = acronym_mapping.get("Sodium_arsenite", "Sodium_arsenite")
        ssr = acronym_mapping.get("Simulated_solar_radiation", "Simulated_solar_radiation")
        new_prefix = f"{prefix_model}_{sa}_{ssr}"
        print(f"  Mapping: '{prefix}' → '{new_prefix}'")
        return new_prefix

    # --- Mouse intestinal organoids ---
    if prefix in ["Mouse_m_Intestinal_organoids_High_fat_diet",
                  "Mouse_m_Intestinal_organoids_Normal_fat_diet"]:
        experiment = "_".join(parts[4:])
        new_prefix = f"Mouse_m_Intestinal_organoids_{experiment}"
        print(f"  Mapping: '{prefix}' → '{new_prefix}'")
        return new_prefix

    # --- Human intestinal organoids ---
    if prefix in ["Human_h_Intestinal_organoids_Colibactin",
                  "Human_h_Intestinal_organoids_5-Fluorouracil"]:
        compound_part = parts[-1]
        new_prefix = f"Human_h_Intestinal_organoids_{compound_part}"
        print(f"  Mapping: '{prefix}' → '{new_prefix}'")
        return new_prefix

    # No match — return unchanged
    return prefix


# ---------------------------------------------------------------------------
# Prefix / suffix helpers
# ---------------------------------------------------------------------------

def extract_sample_prefix(sample_name):
    """
    Strip a trailing numeric replicate index from a sample name.

    e.g. "Mouse_Liver_4NQO_3"  → "Mouse_Liver_4NQO"
         "Mouse_Liver_4NQO"    → "Mouse_Liver_4NQO"

    Parameters
    ----------
    sample_name : str

    Returns
    -------
    str
    """
    parts = sample_name.split('_')
    if parts[-1].isnumeric():
        return '_'.join(parts[:-1])
    return sample_name


def remove_replicate_suffix(name):
    """
    Strip a trailing numeric suffix from a compound token.

    e.g. "4NQO_3" → "4NQO"
         "4NQO"   → "4NQO"

    Parameters
    ----------
    name : str

    Returns
    -------
    str
    """
    parts = name.split('_')
    if len(parts) > 1 and parts[-1].isnumeric():
        return '_'.join(parts[:-1])
    return name


# ---------------------------------------------------------------------------
# C. elegans collapsing
# ---------------------------------------------------------------------------

def apply_celegans_collapsing(prefix):
    """
    Collapse a C. elegans sample prefix to treatment level.

    Rules (applied in order):
    1. Prefixes containing '_CX-5461_UVA' (except the NO-treatment control)
       are collapsed to 'celegans_<strain>_CX-5461_UVA'.
    2. Prefixes that do NOT end in '_NO' and do NOT contain '_CX-5461_'
       are collapsed to 'celegans_<treatment>' (third token only).
    3. All other C. elegans prefixes are collapsed to 'celegans_<strain>'
       (first two tokens).

    Parameters
    ----------
    prefix : str
        Must start with 'celegans_'.

    Returns
    -------
    str
        Collapsed prefix, or original prefix unchanged if not C. elegans.
    """
    if not prefix.startswith('celegans_'):
        return prefix

    if '_CX-5461_UVA' in prefix and prefix != 'celegans_CX-5461_NO':
        parts = prefix.split('_CX-5461_UVA')
        if len(parts) == 2:
            return parts[0] + '_CX-5461_UVA'

    if not prefix.endswith('_NO') and '_CX-5461_' not in prefix:
        parts = prefix.split('_')
        if len(parts) >= 3:
            return f'celegans_{parts[2]}'

    parts = prefix.split('_')
    return "_".join(parts[0:2])


# ---------------------------------------------------------------------------
# Top-level parser
# ---------------------------------------------------------------------------

def parse_sample_into_species_model_compound(sample_name, acronym_mapping=None,
                                              debug=False):
    """
    Parse a standardised sample name into (species_model, compound).

    Parameters
    ----------
    sample_name : str
        Standardised sample name from a clustering output file.
    acronym_mapping : dict, optional
        Mapping returned by load_acronym_mapping(). Applied after prefix
        extraction.
    debug : bool, optional
        If True, print step-by-step parsing details for every sample.
        Default False.

    Returns
    -------
    tuple[str, str]
        (species_model, compound) e.g. ("Mouse_MEF", "4NQO")
    None
        If the sample should be excluded (currently unused, reserved for
        future filtering rules).
    """
    def _dbg(msg):
        if debug:
            print(f"  [parse] {msg}")

    _dbg(f"input: {sample_name}")

    # --- C. elegans ---
    if sample_name.startswith('celegans_'):
        prefix           = extract_sample_prefix(sample_name)
        collapsed_prefix = apply_celegans_collapsing(prefix)
        _dbg(f"celegans collapsed: {collapsed_prefix}")

        if acronym_mapping:
            collapsed_prefix = apply_acronym_mapping(collapsed_prefix, acronym_mapping)
            _dbg(f"celegans after acronym: {collapsed_prefix}")

        parts = collapsed_prefix.split('_', 1)
        compound = remove_replicate_suffix(parts[1]) if len(parts) >= 2 else "Unknown"
        _dbg(f"result: ('C.elegans', '{compound}')")
        return "C.elegans", compound

    # --- Extract prefix (strips trailing replicate index) ---
    prefix = extract_sample_prefix(sample_name)
    _dbg(f"prefix: {prefix}")

    # --- Mouse intestinal organoids (before acronym mapping) ---
    if 'Mouse_m_Intestinal_organoids' in prefix:
        if 'High_fat_diet' in prefix or 'HFD' in prefix:
            return 'Mouse_m_Intestinal_organoids', 'High_fat_diet'
        if 'Normal_fat_diet' in prefix or 'NFD' in prefix:
            return 'Mouse_m_Intestinal_organoids', 'Normal_fat_diet'
        return 'Mouse_m_Intestinal_organoids', 'Unknown'

    if 'Human_h_Intestinal_organoids' in prefix:
        if 'Colibactin' in prefix:
            return 'Human_h_Intestinal_organoids', 'Colibactin'
        if '5-Fluorouracil' in prefix:
            return 'Human_h_Intestinal_organoids', '5-Fluorouracil'
        return 'Human_h_Intestinal_organoids', 'Unknown'

    # --- Apply acronym mapping ---
    if acronym_mapping:
        prefix = apply_acronym_mapping(prefix, acronym_mapping)
        _dbg(f"after acronym mapping: {prefix}")

    # --- Breast / Colon organoids (display label contains newline) ---
    if 'Breast_organoids' in prefix:
        parts = prefix.split('_')
        return 'Human_Breast\norganoids', parts[-1]

    if 'Colon_organoids' in prefix:
        parts = prefix.split('_')
        return 'Human_Colon\norganoids', parts[-1]

    # --- General case ---
    parts = prefix.split('_')

    if len(parts) >= 4 and any(kw in prefix for kw in
                                ['Lymph_Node', 'Bone_Marrow',
                                 'Breast_organoids', 'Colon_organoids']):
        # Three-token model (e.g. Mouse_Lymph_Node)
        species_model = '_'.join(parts[:3])
        raw_compound  = '_'.join(parts[3:]) if len(parts) > 3 else "Unknown"

    elif len(parts) >= 3:
        # Standard two-token model (e.g. Mouse_Liver)
        species_model = '_'.join(parts[:2])
        raw_compound  = '_'.join(parts[2:])

    else:
        species_model = '_'.join(parts[:2]) if len(parts) >= 2 else parts[0]
        raw_compound  = "Unknown"

    compound = remove_replicate_suffix(raw_compound)
    _dbg(f"result: ('{species_model}', '{compound}')")
    return species_model, compound
