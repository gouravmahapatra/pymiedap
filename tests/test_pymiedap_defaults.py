import importlib
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import numpy as np


def _install_stub_modules():
    native_module_names = [
        "module_mie",
        "module_mieshell",
        "module_readmie",
        "module_dap",
        "module_geos",
    ]

    for name in native_module_names:
        sys.modules.setdefault(name, types.ModuleType(name))

    matplotlib_module = sys.modules.setdefault("matplotlib", types.ModuleType("matplotlib"))
    pyplot_module = sys.modules.setdefault("matplotlib.pyplot", types.ModuleType("matplotlib.pyplot"))
    pyplot_module.ioff = getattr(pyplot_module, "ioff", lambda: None)
    pyplot_module.ion = getattr(pyplot_module, "ion", lambda: None)
    matplotlib_module.pyplot = pyplot_module

    pil_module = sys.modules.setdefault("PIL", types.ModuleType("PIL"))
    image_module = sys.modules.setdefault("PIL.Image", types.ModuleType("PIL.Image"))
    pil_module.Image = image_module


def load_pymiedap_module():
    sys.modules.pop("pymiedap.pymiedap", None)
    try:
        return importlib.import_module("pymiedap.pymiedap")
    except ModuleNotFoundError:
        _install_stub_modules()
        sys.modules.pop("pymiedap.pymiedap", None)
        return importlib.import_module("pymiedap.pymiedap")


class ModelDefaultsTest(unittest.TestCase):
    def test_wavelength_update_propagates_to_default_layers(self):
        pmd = load_pymiedap_module()

        model = pmd.Model()
        new_wavelengths = np.array([0.55, 0.65, 0.75])
        model.wvl_list = new_wavelengths

        np.testing.assert_allclose(model.wvl_list, new_wavelengths)
        np.testing.assert_allclose(model.rindex_gas.shape, new_wavelengths.shape)

        for layer in vars(model.layers).values():
            self.assertEqual(len(layer.tau), len(new_wavelengths))
            self.assertEqual(len(layer.tau_g), len(new_wavelengths))
            self.assertEqual(len(layer.tau_ray), len(new_wavelengths))
            self.assertEqual(len(layer.aerosols.nr), len(new_wavelengths))
            self.assertEqual(len(layer.aerosols.ni), len(new_wavelengths))

    def test_surface_and_albedo_defaults_stay_in_sync(self):
        pmd = load_pymiedap_module()

        model = pmd.Model()
        np.testing.assert_allclose(model.surface, np.diag([0.0, 0.0, 0.0, 0.0]))
        self.assertEqual(model.asurf, 0)

        model.asurf = 0.3
        self.assertAlmostEqual(model.surface[0, 0], 0.3)

        new_surface = np.diag([0.8, 0.1, 0.0, 0.0])
        model.surface = new_surface
        self.assertAlmostEqual(model.asurf, 0.8)
        np.testing.assert_allclose(model.surface, new_surface)


class GeometryUtilitiesTest(unittest.TestCase):
    def test_calc_azimuth_and_rotate_stokes_have_expected_simple_values(self):
        pmd = load_pymiedap_module()

        phase = np.array([90.0])
        sza = np.array([45.0])
        emission = np.array([45.0])

        azimuth = pmd.calc_azimuth(phase, sza, emission)
        np.testing.assert_allclose(azimuth, np.array([180.0]))

        q_rot, u_rot = pmd.rotate_stokes(
            np.array([1.0]),
            np.array([0.0]),
            np.array([np.pi / 4.0]),
        )
        np.testing.assert_allclose(q_rot, np.array([0.0]), atol=1e-12)
        np.testing.assert_allclose(u_rot, np.array([-1.0]), atol=1e-12)


class NativeIntegrationTest(unittest.TestCase):
    def test_exopy_imports_with_native_modules(self):
        pmd = load_pymiedap_module()

        native_modules_ready = all(
            hasattr(module, "__file__")
            for module in (pmd.mie, pmd.mieshell, pmd.readmie, pmd.dap, pmd.geos)
        )
        if not native_modules_ready:
            self.skipTest("Native extensions are not available for exopy import testing.")

        import pymiedap.exopy as exopy

        self.assertTrue(hasattr(exopy, "new_body"))
        self.assertTrue(hasattr(exopy, "run_simulation"))

    def test_default_model_computes_and_reads_single_geometry(self):
        pmd = load_pymiedap_module()

        native_modules_ready = all(
            hasattr(module, "__file__")
            for module in (pmd.mie, pmd.mieshell, pmd.readmie, pmd.dap, pmd.geos)
        )
        if not native_modules_ready:
            self.skipTest("Native extensions are not available for end-to-end testing.")

        model = pmd.Model()
        output_name = "default_model"

        with tempfile.TemporaryDirectory() as tmpdir:
            pmd.compute_model(
                model,
                force=True,
                rename=True,
                output_name=output_name,
                path_input=tmpdir,
            )

            self.assertEqual(len(model.name), 1)
            self.assertTrue(os.path.exists(model.name[0]))
            self.assertEqual(Path(model.name[0]).parent, Path(tmpdir))

            phase = np.array([30.0])
            sza = np.array([20.0])
            emission = np.array([10.0])
            intensity, q, u, v = pmd.read_dap_output(phase, sza, emission, model.name[0])

            self.assertEqual(intensity.shape, (1,))
            self.assertEqual(q.shape, (1,))
            self.assertEqual(u.shape, (1,))
            self.assertEqual(v.shape, (1,))
            self.assertTrue(np.isfinite(intensity[0]))


if __name__ == "__main__":
    unittest.main()
