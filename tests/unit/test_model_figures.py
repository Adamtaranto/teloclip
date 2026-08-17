"""
Smoke tests for the coverage-decay documentation figures.

These verify that the plotting code runs headless and writes non-trivial
PNGs; visual regression is out of scope (the committed docs images are
reviewed by eye when regenerated).
"""

import pytest

pytest.importorskip('numpy')
matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')  # headless before any pyplot import

from teloclip.model.figures import (  # noqa: E402
    NATURE_RC,
    plot_decay_families,
    plot_decay_with_ci,
)


def _assert_png(path):
    """
    Check that a written figure exists and looks like a real PNG.

    Parameters
    ----------
    path : pathlib.Path
        Path the figure was written to.
    """
    assert path.exists()
    data = path.read_bytes()
    assert data[:8] == b'\x89PNG\r\n\x1a\n'
    assert len(data) > 10_000  # an empty axes canvas is far smaller


def test_plot_decay_families(tmp_path):
    """The two-panel analytic figure renders to a valid PNG."""
    out = tmp_path / 'families.png'
    plot_decay_families(str(out))
    _assert_png(out)


def test_plot_decay_with_ci(tmp_path):
    """The bootstrap-band figure renders to a valid PNG (small replicate count)."""
    out = tmp_path / 'ci.png'
    plot_decay_with_ci(str(out), depth=15.0, n_reps=8, seed=1)
    _assert_png(out)


def test_generator_script(tmp_path, monkeypatch):
    """The docs regeneration script writes both figures where asked."""
    import importlib.util
    from pathlib import Path

    import teloclip.model.figures as figures

    script = (
        Path(__file__).resolve().parents[2] / 'scripts' / 'generate_model_figures.py'
    )
    spec = importlib.util.spec_from_file_location('generate_model_figures', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    # Keep the smoke test fast: shrink the CI figure to a few replicates.
    monkeypatch.setattr(
        figures,
        'plot_decay_with_ci',
        lambda out_path: plot_decay_with_ci(out_path, depth=10.0, n_reps=4, seed=1),
    )
    module.main(['--out-dir', str(tmp_path)])
    _assert_png(tmp_path / 'coverage-decay-analytic.png')
    _assert_png(tmp_path / 'coverage-decay-bootstrap.png')


def test_nature_rc_is_headless_safe():
    """The style dict never sets an interactive backend or huge canvas."""
    assert 'backend' not in NATURE_RC
    assert NATURE_RC['savefig.dpi'] == 300
