# Makefile for the Venus dap-program with Fourier coefficients
#
# Author: D. M. Stam
# Date: July 2014
# 

FC= gfortran
FFLAGS=  -g -O -Wall

OBJECTS = main.o \
          bmolecules.o \
	  baerosols.o \
          brack.o \
          trace.o \
          gauleg.o \
          adding.o \
          setmu.o \
          setfou.o \
          setzm.o \
          scalzm.o \
          transf.o \
          renorm.o \
          expbmu.o \
          layer0.o \
          layerm.o \
          ord1m.o \
          assign.o \
          ord2m.o \
          fillup.o \
          bstart.o \
          double.o \
          star.o \
          prod.o \
          addsm.o \
          rdiapr.o \
          ldiapr.o \
          addlay.o \
          notop.o \
          nobot.o \
          tstar.o \
          top2bot.o \
          newfou.o 

dapx: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o dapx

main.o: main.f max_incl
brack.o: brack.f 
gauleg.o: gauleg.f
bmolecules.o: bmolecules.f max_incl
baerosols.o: baerosols.f max_incl
adding.o: adding.f max_incl
setmu.o: setmu.f max_incl
setfou.o: setfou.f max_incl
setzm.o: setzm.f max_incl
scalzm.o: scalzm.f max_incl
transf.o: transf.f max_incl
renorm.o: renorm.f max_incl
expbmu.o: expbmu.f max_incl
layer0.o: layerm.f max_incl
layerm.o: layerm.f max_incl
ord1m.o: ord1m.f max_incl
assign.o: assign.f max_incl
ord2m.o: ord2m.f max_incl
fillup.o: fillup.f max_incl
bstart.o: bstart.f max_incl
double.o: double.f max_incl
star.o: star.f max_incl
prod.o: prod.f max_incl
addsm.o: addsm.f max_incl
traces.o: traces.f max_incl
trace.o: trace.f max_incl
rdiapr.o: rdiapr.f max_incl
ldiapr.o: ldiapr.f max_incl
addlay.o: addlay.f max_incl
notop.o: notop.f max_incl
nobot.o: nobot.f max_incl
tstar.o: tstar.f max_incl
top2bot.o: top2bot.f max_incl
newfou.o: newfou.f max_incl

clean: 
	rm -f dapx
	rm -f *.[oc]
