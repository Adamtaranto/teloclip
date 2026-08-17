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

# (legacy module path, tuple of modules it now fronts) for every module that
# moved. Most map one-to-one; teloclip.html_report was split across four.
SHIMS = [
    ('teloclip.overhang', ('teloclip.core.overhang',)),
    ('teloclip.analysis', ('teloclip.core.analysis',)),
    ('teloclip.extension', ('teloclip.core.extension',)),
    ('teloclip.motifs', ('teloclip.core.motifs',)),
    ('teloclip.seqops', ('teloclip.core.seqops',)),
    ('teloclip.streaming_analysis', ('teloclip.core.streaming_analysis',)),
    ('teloclip.samops', ('teloclip.io.sam',)),
    ('teloclip.streaming_io', ('teloclip.io.streaming',)),
    ('teloclip.extract_io', ('teloclip.io.extract',)),
    ('teloclip.reporting', ('teloclip.report.text',)),
    ('teloclip.alignment_layout', ('teloclip.report.layout',)),
    ('teloclip.html_report', ('teloclip.report.html', 'teloclip.report.panels')),
]

IDS = [s[0] for s in SHIMS]

# The split of teloclip.html_report also created teloclip.report.charts and
# teloclip.report.css. Those are new internal seams rather than names the
# legacy module ever exposed, so the shim is not expected to front them.


@pytest.mark.parametrize('legacy,current', SHIMS, ids=IDS)
def test_shim_module_imports(legacy, current):
    """
    Both the legacy and the current module paths import cleanly.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : tuple of str
        Paths the module's contents now live at.
    """
    assert importlib.import_module(legacy) is not None
    for module in current:
        assert importlib.import_module(module) is not None


@pytest.mark.parametrize('legacy,current', SHIMS, ids=IDS)
def test_shim_reexports_are_the_same_objects(legacy, current):
    """
    Every name in the shim's ``__all__`` is the object from the new module.

    Re-exporting a copy rather than the original would break ``isinstance``
    checks and dataclass equality for callers mixing the two import paths.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : tuple of str
        Paths the module's contents now live at.
    """
    old = importlib.import_module(legacy)
    new_modules = [importlib.import_module(m) for m in current]

    assert old.__all__, f'{legacy} re-exports nothing'

    for name in old.__all__:
        assert hasattr(old, name), f'{legacy}.{name} is advertised but missing'

        sources = [m for m in new_modules if hasattr(m, name)]
        assert sources, f'{name} is missing from all of {current}'
        assert any(getattr(old, name) is getattr(m, name) for m in sources), (
            f'{legacy}.{name} is not the same object as the one it fronts'
        )


@pytest.mark.parametrize('legacy,current', SHIMS, ids=IDS)
def test_shim_covers_the_public_surface(legacy, current):
    """
    The shim re-exports every public name defined by the new modules.

    Guards against a name being added to a moved module and quietly not
    reaching users who still import from the legacy path.

    Parameters
    ----------
    legacy : str
        Old top-level import path.
    current : tuple of str
        Paths the module's contents now live at.
    """
    old = importlib.import_module(legacy)

    defined = set()
    for module_name in current:
        module = importlib.import_module(module_name)
        # Only names the module actually defines, not what it imported. An
        # imported symbol's __module__ points at where it was defined.
        defined |= {
            name
            for name in dir(module)
            if not name.startswith('_')
            and getattr(getattr(module, name), '__module__', None) == module_name
        }

    missing = defined - set(old.__all__)
    assert not missing, f'{legacy} does not re-export {sorted(missing)}'
