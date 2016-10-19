# File setup.py


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('',parent_package,top_path)

    config.add_extension(name = 'module_mie',
                         sources=['./mie_sig.pyf',
                                  './mie_source/anbn.f',
                                  './mie_source/devel.f',
                                  './mie_source/fichid.f',
                                  './mie_source/gauleg.f',
                                  './mie_source/getin.f',
                                  './mie_source/pitau.f',
                                  './mie_source/rdwavfile.f',
                                  './mie_source/rminmax.f',
                                  './mie_source/scatmat.f',
                                  './mie_source/sizedis.f',
                                  './mie_source/writsc.f',
                                  ])

    config.add_extension(name = 'module_readmie',
                         sources=['./readmie_sig.pyf',
                                  './readmie_source/read_mie_output.f', ])

    config.add_extension(name = 'module_mieshell',
                         sources=['./mieshell_sig.pyf',
                                  './mieshell_source/bhcoat.f',
                                  './mieshell_source/gauleg.f',
                                  './mieshell_source/pitau.f',
                                  './mieshell_source/scatmat.f',
                                  './mieshell_source/sizedis.f',
                                  ])


    config.add_extension(name = 'module_dap',
                         sources=['./dap_sig.pyf',
                                  './dap_source/adding.f',
                                  './dap_source/addlay.f',
                                  './dap_source/addsm.f',
                                  './dap_source/assign.f',
                                  './dap_source/baerosols.f',
                                  './dap_source/bmolecules.f',
                                  './dap_source/brack.f',
                                  './dap_source/bstart.f',
                                  './dap_source/double.f',
                                  './dap_source/endfou.f',
                                  './dap_source/expbmu.f',
                                  './dap_source/fillup.f',
                                  './dap_source/gauleg.f',
                                  './dap_source/init.f',
                                  './dap_source/layer0.f',
                                  './dap_source/layerm.f',
                                  './dap_source/ldiapr.f',
                                  './dap_source/newfou.f',
                                  './dap_source/nobot.f',
                                  './dap_source/notop.f',
                                  './dap_source/ord1m.f',
                                  './dap_source/ord2m.f',
                                  './dap_source/prod.f',
                                  './dap_source/rdiapr.f',
                                  './dap_source/renorm.f',
                                  './dap_source/scalzm.f',
                                  './dap_source/setfou.f',
                                  './dap_source/setmu.f',
                                  './dap_source/setzm.f',
                                  './dap_source/star.f',
                                  './dap_source/top2bot.f',
                                  './dap_source/trace.f',
                                  './dap_source/transf.f',
                                  './dap_source/tstar.f',
                                  ])

    config.add_extension(name = 'module_geos',
                         sources=['./geos_sig.pyf',
                                  './geos_source/bracks.f',
                                  './geos_source/getgeos.f',
                                  './geos_source/rdfous.f',
                                  './geos_source/read_dap.f',
                                  './geos_source/spline.f',
                                  './geos_source/splint.f',
                                  './geos_source/wrout.f',
                                  ])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

setup(name='PyMieDAP', version='0.1', description='Python Mie DAP program',
      url='', author='Loic Rossi, Daphne Stam',
      author_email='l.c.g.rossi@tudelft.nl', license='GPL/CeCILL',
      packages=['pymiedap'], zip_safe=False,
      setup(**configuration(top_path='').todict())
