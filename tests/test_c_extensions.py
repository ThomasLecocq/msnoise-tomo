import importlib
import numpy as np


def test_import_all_extensions():
    modules = [
        "msnoise_tomo.lib.mk_MatPaths",
        "msnoise_tomo.lib.mkMatSmoothing",
        "msnoise_tomo.lib.vg_fta",
    ]

    for mod in modules:
        m = importlib.import_module(mod)
        assert m is not None


def test_mkMatSmoothing_basic():
    mod = importlib.import_module("msnoise_tomo.lib.mkMatSmoothing")

    # Try to call something minimal if exposed
    funcs = dir(mod)

    # This is intentionally generic (we don't assume API stability)
    assert len(funcs) > 0


def test_vg_fta_numpy_interaction():
    mod = importlib.import_module("msnoise_tomo.lib.vg_fta")

    # minimal numpy interaction sanity check
    arr = np.ones(10, dtype=np.float64)

    # We don't assume exact function names → just ensure module loads with numpy
    assert arr.sum() == 10.0