#!/bin/sh

#This little script recompiles the fortran codes related to the radiative
# transfer code with f2py. For a quick switch between 32 and 64 bits arch

#make miecode
#f2py --fcompiler='gfortran' -c ./mie_source/*.f -m module_mie

#make mieshell
#f2py --fcompiler='gfortran' -c ./mieshell_source/*.f -m module_mieshell

#make readmie
#f2py -c ./readmie_source/read_mie_output.f -m module_readmie

#make dap
#f2py --fcompiler='gfortran' -c ./dap_source/*.f -m module_dap

#make read dap
#f2py --fcompiler='gfortran' -c ./geos_source/bracks.f ./geos_source/rdfous.f ./geos_source/read_dap.f ./geos_source/getgeos.f ./geos_source/spline.f ./geos_source/splint.f -m module_geos
#f2py -c ./geos_source/bracks.f ./geos_source/rdfous.f ./geos_source/read_dap.f ./geos_source/getgeos.f ./geos_source/spline.f ./geos_source/splint.f -m module_geos

#SIGNATURE FILES
#make miecode
f2py ./mie_source/*.f -m module_mie -h mie_sig.pyf

#make mieshell
f2py ./mieshell_source/*.f -m module_mieshell -h mieshell_fig.pyf

#make readmie
f2py ./readmie_source/read_mie_output.f -m module_readmie -h readmie_sig.pyf

#make dap
f2py ./dap_source/*.f -m module_dap -h dap_sig.pyf

#make read dap
f2py ./geos_source/bracks.f ./geos_source/rdfous.f ./geos_source/read_dap.f ./geos_source/getgeos.f ./geos_source/spline.f ./geos_source/splint.f -m module_geos -h geos_sig.pyf
