"""
utils/min_mutations.py

The minimum total mutation count per sample (preprocessing cutoff), resolved
from config/preprocessing.yaml and the optional --min_mutations setting:

    setting            cutoffs used                                   run name
    -----------------  ---------------------------------------------  -----------------
    None / 'default'   <TYPE>: min_mutations (SBS: 307, every species)  min307
    a number, e.g. 250 that number for every species                  min250
    'per-species'      <TYPE>: min_mutations_per_species               min-per-species

min_mutations in the config may be a single number (every species) or a
{species: number} mapping. Cutoffs are returned as {species: number}.
"""

# Species the pipeline recognises (perform_clustering.SPECIES_PATTERNS).
SPECIES = ("mouse", "human", "celegans", "chicken", "rat")

PER_SPECIES = "per-species"


def _expand(value, source):
    if value is None:
        return {}
    if isinstance(value, dict):
        return {sp: _number(v, f"{source}: {sp}") for sp, v in value.items()}
    return {sp: _number(value, source) for sp in SPECIES}


def _number(value, source):
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{source}: expected a number, got {value!r}") from None
    if v < 0:
        raise ValueError(f"{source}: must be >= 0, got {value!r}")
    return int(v) if v.is_integer() else v


def resolve(type_config, setting=None):
    """Return {species: cutoff} ({} for no cutoff) for a type's config section."""
    type_config = type_config or {}
    if setting is None or str(setting).strip().lower() in ("", "default"):
        return _expand(type_config.get("min_mutations"), "min_mutations")
    setting = str(setting).strip().lower()
    if setting == PER_SPECIES:
        per_species = type_config.get("min_mutations_per_species")
        if not per_species:
            raise ValueError("'per-species' needs a min_mutations_per_species section "
                             "for this mutation type in config/preprocessing.yaml.")
        return _expand(per_species, "min_mutations_per_species")
    return _expand(setting, "--min_mutations (a number or 'per-species')")


def label(cutoffs):
    """Run-name label: 'nomin', 'min<N>' (one cutoff) or 'min-per-species'."""
    values = set(cutoffs.values())
    if not values:
        return "nomin"
    if len(values) == 1:
        return f"min{values.pop():g}"
    return "min-per-species"


def describe(cutoffs):
    """Readable form for logs, e.g. '307 for every species'."""
    values = set(cutoffs.values())
    if not values:
        return "none"
    if len(values) == 1:
        return f"{values.pop():g} for every species"
    return ", ".join(f"{sp} {v:g}" for sp, v in cutoffs.items())
