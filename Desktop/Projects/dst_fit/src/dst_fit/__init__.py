"""dst-fit: constitutive identification for dense and shear-thickening suspensions.

Units N1 to N6 of the spec exist: the TRIOS reader, the model library and
fitter, the identifiability and specification verdicts, the whole-series joint
fit, and the cross-curve consistency test that N4, N6 and G4 all apply.

What running those on real data established, and the reason the identifiability
half of this package is the half that matters: **phi and phi_m enter Wyart-Cates
only through phi/phi_J, so a single flow curve approaching jamming determines
the ratio and neither number.** Three published or in-house phi-series all return
a phi_m that tracks the fitted phi. Use ``WYART_CATES_KNOWN_PHI`` where the
volume fraction is known independently; see ``README.md``, gates N4 to G4.

``dst_fit.io_literature`` reads the Edinburgh DataShare deposits and carries
their published constants. It is not imported here because it is only useful
with data that does not ship with the package.
"""

from dst_fit.fitting import (
    DivergenceFit,
    FitResult,
    fit_divergence,
    fit_model,
    synthetic_curve,
)
from dst_fit.identifiability import (
    Assessment,
    ConsistencyCheck,
    assess,
    consistency_across_curves,
    pairwise_degeneracies,
    runs_test,
    transition_coverage,
    weighted_mean,
)
from dst_fit.io_trios import (
    FlowSweep,
    check_stress_consistency,
    read_flow_sweep,
    read_series,
)
from dst_fit.models import MODELS, Model, get_model, thickening_branch

__all__ = [
    "MODELS",
    "Assessment",
    "ConsistencyCheck",
    "DivergenceFit",
    "FitResult",
    "FlowSweep",
    "Model",
    "assess",
    "check_stress_consistency",
    "consistency_across_curves",
    "fit_divergence",
    "fit_model",
    "get_model",
    "pairwise_degeneracies",
    "read_flow_sweep",
    "read_series",
    "runs_test",
    "synthetic_curve",
    "thickening_branch",
    "transition_coverage",
    "weighted_mean",
]

# Single source of truth is pyproject.toml; tests/test_packaging.py fails if the
# two drift. They were 0.0.1 and 0.0.4 simultaneously before that test existed.
__version__ = "0.1.0"
