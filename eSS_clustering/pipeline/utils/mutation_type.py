"""
utils/mutation_type.py

Mutation-type-specific configuration used across the eSignatures pipeline.
All type-dependent constants live here so each script imports from one place
rather than containing if/else chains.

Usage
-----
    from utils.mutation_type import get_config

    cfg = get_config("SBS")
    cfg = get_config("DBS")
    cfg = get_config("ID")
"""

from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# Artifact signature lists (COSMIC v3.x)
# ---------------------------------------------------------------------------

# SBS artifacts — known sequencing / extraction artefacts (COSMIC v3.4)
_SBS_ARTIFACTS = [
    'SBS27', 'SBS43', 'SBS45', 'SBS46', 'SBS47', 'SBS48', 'SBS49', 'SBS50',
    'SBS51', 'SBS52', 'SBS53', 'SBS54', 'SBS55', 'SBS56', 'SBS57', 'SBS58',
    'SBS59', 'SBS60', 'SBS95',
]

# DBS: COSMIC v3.3.1 does not define accepted artefact DBS signatures;
# extend this list if/when COSMIC publishes them.
_DBS_ARTIFACTS: list = []

# ID: COSMIC v3.3.1 artefact indel signatures
_ID_ARTIFACTS = ['ID9']


# ---------------------------------------------------------------------------
# Data class
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MutationTypeConfig:
    """
    All mutation-type-specific constants in one place.

    Attributes
    ----------
    name : str
        Short name used in CLI arguments and output path namespacing.
        One of "SBS", "DBS", "ID".
    n_contexts : int
        Number of mutational contexts (96 / 78 / 83).
    plot_function : str
        Name of the sigProfilerPlotting function to call.
    plot_matrix_type : str
        Second positional argument passed to that function.
    pdf_prefix : str
        Prefix of the PDF file that sigProfilerPlotting produces.
        e.g. "SBS_96_plots" → sigPlt produces "SBS_96_plots_<project>.pdf"
    artifact_signatures : list[str]
        COSMIC signature IDs to flag as artefacts in static heatmap.
    has_prenormalized_files : bool
        True  → load raw counts AND normalized files separately (SBS; the
                normalized files are written by preprocessing).
        False → load raw counts only; pipeline normalises internally (DBS/ID).
    cluster_prefix : str
        Label prefix in output IDs: "eSS", "eDS", "eIS".
    default_custom_thresholds : dict[str, float]
        Per-cluster cosine-distance thresholds ({sample-name substring:
        distance}) applied when --custom_thresholds is not given. A main
        cluster containing a matching sample is re-split at this tighter
        distance. Empty for every type: custom thresholds are off by default.
    custom_threshold_presets : dict[str, dict[str, float]]
        Named sets of custom thresholds that can be passed by name to
        --custom_thresholds (or run_pipeline.sh), e.g. the SBS 'aai-split'
        used for testing.
    """
    name:                  str
    n_contexts:            int
    plot_function:         str
    plot_matrix_type:      str
    pdf_prefix:            str
    artifact_signatures:   list
    has_prenormalized_files: bool
    cluster_prefix:        str   # Label prefix in output IDs: "eSS", "eDS", "eIS"
    default_custom_thresholds: dict = field(default_factory=dict)
    custom_threshold_presets: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Pre-built configs
# ---------------------------------------------------------------------------

_CONFIGS = {
    "SBS": MutationTypeConfig(
        name                  = "SBS",
        n_contexts            = 96,
        plot_function         = "plotSBS",
        plot_matrix_type      = "96",
        pdf_prefix            = "SBS_96_plots",
        artifact_signatures   = _SBS_ARTIFACTS,
        has_prenormalized_files = True,
        cluster_prefix        = "eSS",
        # Optional, for testing: re-split the cluster containing the
        # Aristolochic acid I / Dibenzo[a,l]pyrene samples at distance 0.095
        # (the "manual AAI/DBP split"). Not applied by default.
        custom_threshold_presets = {
            "aai-split": {
                "Aristolochic_acid_I": 0.095,
                "Dibenzo[a,l]pyrene":  0.095,
            },
        },
    ),
    "DBS": MutationTypeConfig(
        name                  = "DBS",
        n_contexts            = 78,
        plot_function         = "plotDBS",
        plot_matrix_type      = "78",
        pdf_prefix            = "DBS_78_plots",
        artifact_signatures   = _DBS_ARTIFACTS,
        has_prenormalized_files = False,
        cluster_prefix        = "eDS",
    ),
    "ID": MutationTypeConfig(
        name                  = "ID",
        n_contexts            = 83,
        plot_function         = "plotID",
        plot_matrix_type      = "83",
        pdf_prefix            = "ID_83_plots",
        artifact_signatures   = _ID_ARTIFACTS,
        has_prenormalized_files = False,
        cluster_prefix        = "eIS",
    ),
}


def get_config(mutation_type: str) -> MutationTypeConfig:
    """
    Return the MutationTypeConfig for the given mutation type string.

    Parameters
    ----------
    mutation_type : str
        Case-insensitive.  One of "SBS", "DBS", "ID".

    Returns
    -------
    MutationTypeConfig

    Raises
    ------
    ValueError
        If the mutation type is not recognised.
    """
    key = mutation_type.upper()
    if key not in _CONFIGS:
        raise ValueError(
            f"Unknown mutation type '{mutation_type}'. "
            f"Must be one of: {', '.join(_CONFIGS.keys())}"
        )
    return _CONFIGS[key]


def call_plot_function(cfg: MutationTypeConfig, matrix_path: str,
                       output_dir: str, project_name: str) -> None:
    """
    Call the appropriate sigProfilerPlotting function for this mutation type.

    Parameters
    ----------
    cfg          : MutationTypeConfig
    matrix_path  : str  Path to the normalised profile matrix TSV.
    output_dir   : str  Directory where the PDF will be written.
    project_name : str  Project label (used in the output PDF filename).
    """
    import sigProfilerPlotting as sigPlt

    fn = getattr(sigPlt, cfg.plot_function)
    fn(matrix_path, output_dir, project_name,
       cfg.plot_matrix_type, percentage=True)
