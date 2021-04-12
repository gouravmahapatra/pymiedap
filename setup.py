# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# File setup.py

import setuptools
from numpy.distutils.core import setup, Extension
import os
import os.path as osp

module_mie = Extension(name = 'module_mie',
                       sources=['./src/mie_sig.pyf',
                                './src/mie_source/anbn.f',
                                './src/mie_source/devel.f',
                                './src/mie_source/fichid.f',
                                './src/mie_source/gauleg.f',
                                './src/mie_source/pitau.f',
                                './src/mie_source/rminmax.f',
                                './src/mie_source/scatmat.f',
                                './src/mie_source/sizedis.f',
                                './src/mie_source/writsc.f',
                                ])

module_readmie = Extension(name = 'module_readmie',
                            sources=['./src/readmie_sig.pyf',
                                     './src/readmie_source/read_mie_output.f',
                                     './src/readmie_source/file2coefs.f', ])

module_mieshell = Extension(name = 'module_mieshell',
                            sources=['./src/mieshell_sig.pyf',
                                     './src/mieshell_source/bhcoat.f',
                                     './src/mieshell_source/gauleg.f',
                                     './src/mieshell_source/pitau.f',
                                     './src/mieshell_source/scatmat.f',
                                     './src/mieshell_source/sizedis.f',
                                     ])

module_dap = Extension(name = 'module_dap',
                       sources=['./src/dap_sig.pyf',
                                './src/dap_source/adding.f',
                                './src/dap_source/addlay.f',
                                './src/dap_source/addsm.f',
                                './src/dap_source/assign.f',
                                './src/dap_source/bmolecules.f',
                                './src/dap_source/brack.f',
                                './src/dap_source/bstart.f',
                                './src/dap_source/double.f',
                                './src/dap_source/expbmu.f',
                                './src/dap_source/fillup.f',
                                './src/dap_source/gauleg.f',
                                './src/dap_source/init.f',
                                './src/dap_source/layer0.f',
                                './src/dap_source/layerm.f',
                                './src/dap_source/ldiapr.f',
                                './src/dap_source/newfou.f',
                                './src/dap_source/nobot.f',
                                './src/dap_source/notop.f',
                                './src/dap_source/ord1m.f',
                                './src/dap_source/ord2m.f',
                                './src/dap_source/prod.f',
                                './src/dap_source/rdiapr.f',
                                './src/dap_source/renorm.f',
                                './src/dap_source/scalzm.f',
                                './src/dap_source/setfou.f',
                                './src/dap_source/setmu.f',
                                './src/dap_source/setzm.f',
                                './src/dap_source/star.f',
                                './src/dap_source/top2bot.f',
                                './src/dap_source/trace.f',
                                './src/dap_source/transf.f',
                                './src/dap_source/tstar.f',
                                ])

module_geos = Extension(name = 'module_geos',
                        sources=['./src/geos_sig.pyf',
                                 './src/geos_source/bracks.f',
                                 './src/geos_source/getgeos.f',
                                 './src/geos_source/rdfous.f',
                                 './src/geos_source/read_dap.f',
                                 './src/geos_source/spline.f',
                                 './src/geos_source/splint.f',
                                 './src/geos_source/wrout.f',
                                 ])

longdescp = '''Python Mie DAP package: radiative transfer model and tools to
            simulate flux and polarization of planetary atmospheres.
            Exopy also allows for simulation with complex orbits and Star-Planet-Moon systems'''

setup(name='PyMieDAP',
      version='0.1',
      description='Python Mie DAP program',
      url='',
      author='Loic Rossi, Javier Berzosa Molina, Daphne Stam',
      author_email='l.c.g.rossi@tudelft.nl',
      maintainer ='Loic Rossi',
      maintainer_email='l.c.g.rossi@tudelft.nl',
      license='GPL/CeCILL',
      packages=['pymiedap','pymiedap.exopy_source'],
      ext_modules=[module_mie, module_readmie, module_mieshell, module_dap, module_geos],
      zip_safe=False,
      long_description = longdescp,
      )
