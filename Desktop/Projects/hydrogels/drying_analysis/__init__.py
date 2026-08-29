"""
drying_analysis -- camera-first measurement layer for the PEG drying-fracture
program (crack / peel / curl mode map).

Each module implements one roadmap observable from hydrogel_MOC.md and is
validated against synthetic ground truth in `validate.py`. This closes the
AUDIT_camera.md S1 gap: the prior align/analyze/annotate/visualize scripts
measure surface-wrinkle morphology, not the drying-fracture observables.

Modules
-------
calibration        : px/mm from a fiducial of known length (no hardcoded scale).
swelling_kinetics  : D, k_m, Bi from gravimetric m(t) curves (MOC 5A).
indentation_modulus: E from gravity-sphere indentation + finite-thickness (MOC 5).
curl_radius        : curl radius rho from a side-view arc (circle fit).
crack_metrics      : polygon area A_p and crack spacing Lambda from a top view.
toughness          : pure-shear Gamma, fractocohesive l_fc, ell = Gamma/E.
kymograph          : side-view kymograph + crack-vs-peel ordering (MOC 5B).
"""
from . import (
    calibration,
    crack_metrics,
    curl_radius,
    indentation_modulus,
    kymograph,
    swelling_kinetics,
    toughness,
)

__all__ = [
    "calibration",
    "crack_metrics",
    "curl_radius",
    "indentation_modulus",
    "kymograph",
    "swelling_kinetics",
    "toughness",
]
