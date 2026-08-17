"""
Tests for the backwards-compatible module aliases at the old flat import paths.

When the package was reorganised into ``core``/``io``/``report`` subpackages,
every module kept a shim at its previous top-level path so that existing user
code importing e.g. ``teloclip.analysis`` keeps working. The rest of the test
suite imports from the new paths, so without these tests the shims would go
entirely unexercised and could rot silently.

Each shim is checked for two things: that every name it advertises in
``__all__`` actually imports, and that the object it hands back is the *same*
object as the one in the new module rather than a copy. Identity matters
because ``isinstance`` checks and dataclass equality across the two import
paths would otherwise fail in confusing ways.
"""

import importlib

import pytest

# (legacy module path, current module path) for every module that moved.
SHIMS = [
    ('teloclip.overhang', 'teloclip.core.overhang'),
    ('teloclip.analysis', 'teloclip.core.analysis'),
    ('teloclip.extension', 'teloclip.core.extension'),
    ('teloclip.motifs', 'teloclip.core.motifs'),
    ('teloclip.seqops', 'teloclip.core.seqops'),
    ('teloclip.streaming_analysis', 'teloclip.core.streaming_analysis'),
    ('teloclip.samops', 'teloclip.io.sam'),
    ('teloclip.streaming_io', 'teloclip.io.streaming'),
    ('teloclip.extract_io', 'teloclip.io.extract'),
    ('teloclip.reporting', 'teloclip.report.text'),
    ('teloclip.alignment_layout', 'teloclip.report.layout'),
    ('teloclip.html_report', 'teloclip.report.html'),
]


@pytest.mark.parametrize('legacy,current', SHIMS, ids=[s[0] for s in SHIMS])
def test_shim_module_imports(legacy, current):
    """
    Both the legacy and the current module path import cleanly.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : str
        Path the module now lives at.
    """
    assert importlib.import_module(legacy) is not None
    assert importlib.import_module(current) is not None


@pytest.mark.parametrize('legacy,current', SHIMS, ids=[s[0] for s in SHIMS])
def test_shim_reexports_are_the_same_objects(legacy, current):
    """
    Every name in the shim's ``__all__`` is the object from the new module.

    Re-exporting a copy rather than the original would break ``isinstance``
    checks and dataclass equality for callers mixing the two import paths.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : str
        Path the module now lives at.
    """
    old = importlib.import_module(legacy)
    new = importlib.import_module(current)

    assert old.__all__, f'{legacy} re-exports nothing'

    for name in old.__all__:
        assert hasattr(old, name), f'{legacy}.{name} is advertised but missing'
        assert hasattr(new, name), f'{current}.{name} is missing'
        assert getattr(old, name) is getattr(new, name), (
            f'{legacy}.{name} is not the same object as {current}.{name}'
        )


@pytest.mark.parametrize('legacy,current', SHIMS, ids=[s[0] for s in SHIMS])
def test_shim_covers_the_public_surface(legacy, current):
    """
    The shim re-exports every public name defined by the new module.

    Guards against a name being added to a moved module and quietly not
    reaching users who still import from the legacy path.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : str
        Path the module now lives at.
    """
    old = importlib.import_module(legacy)
    new = importlib.import_module(current)

    # Only names the new module actually defines, not what it imported. An
    # imported symbol's __module__ points at where it was defined.
    defined = {
        name
        for name in dir(new)
        if not name.startswith('_')
        and getattr(getattr(new, name), '__module__', None) == current
    }

    missing = defined - set(old.__all__)
    assert not missing, f'{legacy} does not re-export {sorted(missing)}'
