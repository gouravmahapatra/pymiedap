# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

from __future__ import annotations

import importlib.machinery
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


PROJECT_ROOT = Path(__file__).resolve().parent


class F2PYBuildExt(build_ext):
    """Build Fortran extensions through ``numpy.f2py`` without numpy.distutils."""

    def run(self):
        for extension in self.extensions:
            self.build_extension(extension)

    def build_extension(self, ext):
        build_dir = Path(self.build_temp) / ext.name
        build_dir.mkdir(parents=True, exist_ok=True)

        command = [
            sys.executable,
            "-m",
            "numpy.f2py",
            "-c",
            "-m",
            ext.name,
            *[str((PROJECT_ROOT / source).resolve()) for source in ext.sources],
        ]
        subprocess.run(command, cwd=build_dir, check=True)

        candidates = []
        for suffix in importlib.machinery.EXTENSION_SUFFIXES:
            candidates.extend(build_dir.glob(f"{ext.name}*{suffix}"))

        if not candidates:
            raise RuntimeError(f"Failed to build extension {ext.name!r} with f2py.")

        built_extension = max(candidates, key=lambda path: path.stat().st_mtime)
        target_path = Path(self.get_ext_fullpath(ext.name))
        target_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built_extension, target_path)


EXTENSIONS = [
    Extension(
        name="module_mie",
        sources=[
            "mie_sig.pyf",
            "mie_source/anbn.f",
            "mie_source/devel.f",
            "mie_source/fichid.f",
            "mie_source/gauleg.f",
            "mie_source/pitau.f",
            "mie_source/rminmax.f",
            "mie_source/scatmat.f",
            "mie_source/sizedis.f",
            "mie_source/writsc.f",
        ],
    ),
    Extension(
        name="module_readmie",
        sources=[
            "readmie_sig.pyf",
            "readmie_source/read_mie_output.f",
            "readmie_source/file2coefs.f",
        ],
    ),
    Extension(
        name="module_mieshell",
        sources=[
            "mieshell_sig.pyf",
            "mieshell_source/bhcoat.f",
            "mieshell_source/gauleg.f",
            "mieshell_source/pitau.f",
            "mieshell_source/scatmat.f",
            "mieshell_source/sizedis.f",
        ],
    ),
    Extension(
        name="module_dap",
        sources=[
            "dap_sig.pyf",
            "dap_source/adding.f",
            "dap_source/addlay.f",
            "dap_source/addsm.f",
            "dap_source/assign.f",
            "dap_source/bmolecules.f",
            "dap_source/brack.f",
            "dap_source/bstart.f",
            "dap_source/double.f",
            "dap_source/expbmu.f",
            "dap_source/fillup.f",
            "dap_source/gauleg.f",
            "dap_source/init.f",
            "dap_source/layer0.f",
            "dap_source/layerm.f",
            "dap_source/ldiapr.f",
            "dap_source/newfou.f",
            "dap_source/nobot.f",
            "dap_source/notop.f",
            "dap_source/ord1m.f",
            "dap_source/ord2m.f",
            "dap_source/prod.f",
            "dap_source/rdiapr.f",
            "dap_source/renorm.f",
            "dap_source/scalzm.f",
            "dap_source/setfou.f",
            "dap_source/setmu.f",
            "dap_source/setzm.f",
            "dap_source/star.f",
            "dap_source/top2bot.f",
            "dap_source/trace.f",
            "dap_source/transf.f",
            "dap_source/tstar.f",
        ],
    ),
    Extension(
        name="module_geos",
        sources=[
            "geos_sig.pyf",
            "geos_source/bracks.f",
            "geos_source/getgeos.f",
            "geos_source/rdfous.f",
            "geos_source/read_dap.f",
            "geos_source/spline.f",
            "geos_source/splint.f",
            "geos_source/wrout.f",
        ],
    ),
]


README = PROJECT_ROOT.joinpath("README.md").read_text(encoding="utf-8")


setup(
    name="PyMieDAP",
    version="0.1.1",
    description="Radiative transfer and polarized scattering tools for planetary atmospheres",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Loic Rossi, Javier Berzosa Molina, Daphne Stam",
    author_email="l.c.g.rossi@tudelft.nl",
    maintainer="Loic Rossi",
    maintainer_email="l.c.g.rossi@tudelft.nl",
    license="GPL/CeCILL",
    python_requires=">=3.9",
    packages=find_packages(include=["pymiedap", "pymiedap.*"]),
    install_requires=[
        "numpy>=1.23",
        "scipy",
        "matplotlib",
        "Pillow",
    ],
    extras_require={
        "test": ["pytest"],
        "dev": ["pytest", "jupyter", "ipython"],
    },
    ext_modules=EXTENSIONS,
    cmdclass={"build_ext": F2PYBuildExt},
    include_package_data=True,
    zip_safe=False,
)
