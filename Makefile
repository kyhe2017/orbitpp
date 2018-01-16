#---------------------------------------------------------------
#  Makefile to produce output on Unix workstations
#
#  orbitpp
#
#  make
#---------------------------------------------------------------

os:= $(shell uname -s)
arch:= $(shell uname -m)
host:= $(shell hostname)

ORBITPP_OBJECTS:= \
      orbitpp_private.o \
      geometry.o \
      trajectory.o \
      diagloss.o \
      Qloss.o \
      distribution.o \
      aux_distribution.o \
      statistics.o \
      aux_statistics.o \
      momenta.o \
      csdata.o \
      mod_interfaces.o \
      orbitpp_sub.o \
      orbitpp_sub2.o \
      orbitpp.o

ORBITPP_SRCS:= $(ORBITPP_OBJECTS:.o=.f90)

F90:= ifort

FFLAGS:=-r8 -O2 # -fpe0 -ftrapuv

#---------------------------------------------------------------

all: orbitpp
	@echo -e "\E[1;31m"
	@echo -e ">> `date '+%a %d-%h-%y %r'`  `uname -ms` $(LOGNAME)"
	@echo -e "\E[0m"

orbitpp: $(ORBITPP_OBJECTS)
	$(F90) $(ORBITPP_OBJECTS) -o $@

#---------------------------------------------------------------
.SUFFIXES: .o .f90

.f90.o: $(ORBITPP_SRCS)
	$(F90) $(FFLAGS) -c $<

#---------------------------------------------------------------
clean:
	rm -f orbitpp *.o *.mod

