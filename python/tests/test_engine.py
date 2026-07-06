"""Smoke tests for the PyO3 boundary (stages._engine).

These tests prove that the Rust extension built into the wheel:

1. Imports without errors on every CPython we build for.
2. Exposes the M0 surface (version() + the vle-thermo smoke path).
3. Evaluates the vle-thermo dependency end-to-end across the FFI boundary.

CI runs this file via cibuildwheel's `test-command = "pytest {package}/tests"`
on every (OS, arch) combination, so a missing or broken binding fails the
release pipeline before it can publish.
"""

import re

# Importing `stages` exercises the Rust shared object via stages._engine. If
# the wheel was built without the `python` feature, or the abi3 target
# mismatched, this import fails and pytest aborts before the asserts run.
import stages


SEMVER_RE = re.compile(r"^\d+\.\d+\.\d+(?:[-.+][\w.-]+)?$")


def test_version_is_semver_string() -> None:
    """`version()` and `__version__` are semver-shaped strings."""
    assert SEMVER_RE.match(stages.version())
    assert SEMVER_RE.match(stages.__version__)


def test_smoke_bubble_temperature_is_physical() -> None:
    """The vle-thermo smoke path evaluates and returns a plausible bubble T [K].

    Deliberately a wide window — this asserts the dependency works across the
    FFI boundary, not a validated value (methanol/water with a bare cubic EOS
    is only approximate; the pinned literature point arrives at M1).
    """
    t = stages.smoke_bubble_temperature()
    assert isinstance(t, float)
    assert 280.0 < t < 400.0, f"methanol/water bubble T = {t} K outside [280, 400)"
