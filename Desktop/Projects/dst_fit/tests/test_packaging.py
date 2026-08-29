"""Does this package actually install and run for someone who is not us?

Every check here corresponds to a real defect found on 2026-08-29, while getting
the package ready to ship. None of them would have been caught by the numerical
tests, because we run from a source tree with every dependency already present.

The one that mattered: ``scipy`` was used on the main code path and absent from
``dependencies``. ``pip install dst-fit`` would have succeeded, ``import dst_fit``
would have succeeded, and the first call to ``fit_model`` would have raised
``ModuleNotFoundError``.
"""

from __future__ import annotations

import ast
import sys
import tomllib
from pathlib import Path

import pytest

import dst_fit

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
SRC = ROOT / "src" / "dst_fit"

# Modules in the standard library, or the package itself, need no declaration.
_LOCAL = {"dst_fit"}


def project() -> dict:
    with PYPROJECT.open("rb") as fh:
        return tomllib.load(fh)["project"]


def declared_dependencies() -> set[str]:
    out = set()
    for spec in project().get("dependencies", []):
        name = spec.split(">")[0].split("<")[0].split("=")[0].split("[")[0]
        out.add(name.strip().lower().replace("-", "_"))
    return out


def imported_top_level_modules() -> set[str]:
    """Every top-level module imported anywhere under ``src/dst_fit``."""
    found = set()
    for path in SRC.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    found.add(alias.name.split(".")[0])
            elif isinstance(node, ast.ImportFrom):
                if node.level == 0 and node.module:
                    found.add(node.module.split(".")[0])
    return found


def test_every_imported_third_party_module_is_declared():
    """The scipy defect, guarded.

    Anything imported by the package and not in the standard library has to be in
    ``dependencies``, or an installed copy raises at first use.
    """
    imported = imported_top_level_modules()
    third_party = {
        m for m in imported
        if m not in sys.stdlib_module_names and m not in _LOCAL and not m.startswith("_")
    }
    missing = sorted(third_party - declared_dependencies())
    assert not missing, (
        f"imported by src/dst_fit but not in pyproject dependencies: {missing}. "
        f"An installed copy would raise ModuleNotFoundError at first use."
    )


def test_scipy_specifically_is_declared():
    """Named explicitly, because it is the one that actually shipped broken."""
    assert "scipy" in declared_dependencies()


def test_version_matches_pyproject():
    """These read 0.0.1 and 0.0.4 at the same time before this test existed."""
    assert dst_fit.__version__ == project()["version"]


def test_license_file_exists_and_is_referenced():
    """pyproject claimed MIT for two days with no LICENSE file in the tree."""
    declared = project().get("license")
    assert isinstance(declared, dict) and "file" in declared, (
        "declare the licence as a file so the text ships with the wheel"
    )
    licence = ROOT / declared["file"]
    assert licence.is_file(), f"{declared['file']} is referenced but missing"
    text = licence.read_text(encoding="utf-8")
    assert "MIT License" in text
    assert "Copyright (c)" in text


def test_readme_referenced_by_pyproject_exists():
    assert (ROOT / project()["readme"]).is_file()


def test_package_exposes_a_usable_public_api():
    """__all__ must not name anything that is not importable."""
    for name in dst_fit.__all__:
        assert hasattr(dst_fit, name), f"dst_fit.__all__ names {name!r}, which is absent"


@pytest.mark.parametrize("module", ["fitting", "identifiability", "models", "profiling",
                                    "joint", "transient", "io_trios", "io_literature"])
def test_every_module_imports_without_the_literature_data(module):
    """Importing must never depend on the gitignored deposits being present."""
    __import__(f"dst_fit.{module}")
