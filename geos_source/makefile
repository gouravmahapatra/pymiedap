# Makefile for the Fourier geometries program
#
# Author: D. M. Stam
# Date: July 2013
# 

FC= gfortran
FFLAGS=  -g -O -fbounds-check

OBJECTS = main.o \
          bracks.o \
          splint.o \
          rdfous.o \
          getgeos.o \
          spline.o \
          wrout.o

geosx: $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o geosx

main.o: main.f max_incl
rdfous.o: rdfous.f max_incl
getgeos.o: getgeos.f max_incl
wrout.o: wrout.f max_incl

clean: 
	rm -f geosx
	rm -f *.[oc]
