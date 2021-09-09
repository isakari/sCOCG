TARGET = ../bin/a.out

FC = ifort
#USE_OPENMP = 1
# OPT = 1
DBG = 1
OBJS = module_globals.o \
	module_random_numbers.o \
	module_scocg.o \
	main.o
GFORTRAN = gfortran
IFORT = ifort

ifeq ($(FC), $(GFORTRAN))
FDEBUGS = -g -fbounds-check -Wuninitialized -fbacktrace -Wunused
FFLAGS = -cpp
endif
ifeq ($(FC), $(IFORT))
FDEBUGS = -check all -warn all -gen_interfaces -fpe0 -ftrapuv -traceback -g
FFLAGS = -fpp -DIFORT
endif
ifdef DBG
FFLAGS += $(FDEBUGS) -DDBG
endif

ifdef OPT
ifeq ($(FC), $(GFORTRAN))
FFLAGS += -O3
LDFLAGS += -O3
endif
ifeq ($(FC), $(IFORT))
FFLAGS += -O3 -ip
LDFLAGS += -O3 -ip
endif
endif

ifeq ($(FC), $(GFORTRAN))
LIBS += -llapack -lblas
endif
ifeq ($(FC), $(IFORT))
LIBS += -mkl
endif

LDFLAGS = $(LIBS)
INCLUDES =

ifdef USE_OPENMP
ifeq ($(FC), $(GFORTRAN))
FFLAGS += -fopenmp
LDFLAGS += -fopenmp
endif
ifeq ($(FC), $(IFORT))
FFLAGS += -qopenmp
LDFLAGS += -qopenmp
endif
endif

.SUFFIXES:
.SUFFIXES: .o .f90 .f

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

%.mod: %.f90 %.o
	@:

$(TARGET): $(OBJS)
	$(FC) $(OBJS) $(LDFLAGS) -o $@

.PHONY: clean
clean:
	-rm -f *.o *.mod *.d *~ $(TARGET) *genmod*
