# This file is part of PyMieDAP, released under GNU General Public License.
# See license.md or http://gitlab.com/loic.cg.rossi/pymiedap for details.

# File setup.py

import setuptools
from numpy.distutils.core import setup, Extension
import os
import os.path as osp

module_mie = Extension(name = 'module_mie',
                       sources=['./mie_source/anbn.f',
                                './mie_source/devel.f',
                                './mie_source/fichid.f',
                                './mie_source/gauleg.f',
                                './mie_source/pitau.f',
                                './mie_source/rminmax.f',
                                './mie_source/scatmat.f',
                                './mie_source/sizedis.f',
                                './mie_source/writsc.f',
                                ])

module_readmie = Extension(name = 'module_readmie',
                            sources=['./readmie_source/read_mie_output.f',
                                     './readmie_source/file2coefs.f', ])

module_mieshell = Extension(name = 'module_mieshell',
                            sources=['./mieshell_source/bhcoat.f',
                                     './mieshell_source/gauleg.f',
                                     './mieshell_source/pitau.f',
                                     './mieshell_source/scatmat.f',
                                     './mieshell_source/sizedis.f',
                                     ])

module_dap = Extension(name = 'module_dap',
                       sources=['./dap_source/adding.f',
                                './dap_source/addingascii.f',
                                './dap_source/addlay.f',
                                './dap_source/bmolecules.f',
                                './dap_source/brack.f',
                                './dap_source/bstart.f',
                                './dap_source/double.f',
                                './dap_source/expbmu.f',
                                './dap_source/fillup.f',
                                './dap_source/gauleg.f',
                                './dap_source/init.f',
                                './dap_source/layer0.f',
                                './dap_source/layerm.f',
                                './dap_source/newfouascii.f',
                                './dap_source/newfou.f',
                                './dap_source/newfouderiv.f',
                                './dap_source/RGPD3P.f',
                                './dap_source/nobot.f',
                                './dap_source/notop.f',
                                './dap_source/ord1m.f',
                                './dap_source/ord2m.f',
                                './dap_source/prod.f',
                                './dap_source/renorm.f',
                                './dap_source/scalzm.f',
                                './dap_source/setfou.f',
                                './dap_source/setfouascii.f',
                                './dap_source/setmu.f',
                                './dap_source/setmuascii.f',
                                './dap_source/setzm.f',
                                './dap_source/star.f',
                                './dap_source/top2bot.f',
                                './dap_source/trace.f',
                                './dap_source/transf.f',
                                './dap_source/tstar.f'
                                ])

module_geos = Extension(name = 'module_geos',
                        sources=['./geos_source/bracks.f',
                                 './geos_source/getgeos.f',
                                 './geos_source/rdfousascii.f',
                                 './geos_source/rdfous.f',
                                 './geos_source/rdfousderiv.f',
                                 './geos_source/read_dapascii.f',
                                 './geos_source/spline.f',
                                 './geos_source/splint.f',
                                 './geos_source/read_dap.f',
                                 './geos_source/RGPD3P.f',
                                 './geos_source/RGBI3P.f',
                                 './geos_source/RGLCTN.f',
                                 './geos_source/RGPLNL.f',
                                 './geos_source/wrout.f',
                                 './geos_source/ascii2hdf5.f',
                                 './geos_source/ascii2hdf5deriv.f',
                                 './geos_source/hdf5deriv2ascii.f',
                                 './geos_source/hdf5deriv2hdf5.f',
                                 './geos_source/hdf52ascii.f',
                                 './geos_source/hdf52hdf5deriv.f',
                                 './geos_source/gauleg.f'
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
