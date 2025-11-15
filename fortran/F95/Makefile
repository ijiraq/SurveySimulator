# Things that are needed fort the CLI (User interface)
U	= modelutils GiMeObj 
OBJU	= $(addsuffix .o,${U})
SRCU	= $(addsuffix .f95,${U})
FPPU	= $(addsuffix .fpp,${U})
MODU	= $(addsuffix .mod,${U})

# Things that should end up with python module access
M   = common_data datadec poly_lib effut xvutils numutils getsur ioutils rot elemutils debug surveysub
# surveysub
OBJM	= $(addsuffix .o,${M})
SRCM	= $(addsuffix .f95,${M})
FPPM	= $(addsuffix .fpp,${M})
MODM	= $(addsuffix .mod,${M})

# Things that should be in the library
# A	= common_data debug poly_dec datadec ioutils poly_lib effut elemutils rot xvutils getsur numutils
# 
A	= common_data poly_dec poly_lib debug datadec effut ioutils elemutils rot xvutils  getsur numutils surveysub
OBJA	= $(addsuffix .o,${A})
SRCA	= $(addsuffix .f95,${A})
FPPA	= $(addsuffix .fpp,${A})
MODA	= $(addsuffix .mod,${A})

# SSIM	= SurveySubsF95
SSIM	= SurveySubs
GOBJ	= GiMeObjF95

FC	= gfortran
FPP	= gfortran -E -x f95-cpp-input -fPIC -DPYTHON
# FPP_F90FLAGS = -x f95-cpp-input -fPIC


FFLAGS	= -O3 -fPIC -x f95-cpp-input
FFLAGP	= -O3
CFLAGS  = 
CPPFLAGS = 
LIBTOOL = ar src
F2PYFLAGS = --build-dir build --fcompiler=gfortran

KIND_MAP = kind_map

TOPTARGETS := clean test
# use the GIMEOBJ from the command line or the default value of InnerHotModel
GIMEOBJ ?= InnerHotModel

$(TOPTARGETS): $(SUBDIRS)

.PHONY: $(SUBDIRS) 

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

all:	Driver $(SUBDIRS)

test: 
	$(MAKE) -s -C ../tests test LANGUAGE=F95

%.fpp :	%.f95
	$(FPP) $< -o $@

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean: rmtmp
	rm -f _$(SSIM)*.so _$(GOBJ)*.so  _$(IOBJ)*.so _$(ROBJ)*.so Driver GiMeObj.f95
	rm -rf build

rmtmp:
	rm -f *~ *.o *.mod *.fpp f90wrap*f90 SimulDetect.dat SimulTrack.dat
	rm -f $(IOBJ).py $(ROBJ).py $(SSIM).py $(GOBJ).py _$(SSIM)*.so
	rm -f *.log LOG .f2py_f2cmap

.PHONY: link

link:
	\rm -f GiMeObj.f95
	\rm -f GiMeObj.o
	ln -s $(GIMEOBJ).f95 GiMeObj.f95

Driver: link Driver.f95 $(OBJA) $(OBJU) Makefile GiMeObj.f95
	$(FC) $(FFLAGP) -o Driver $(OBJA) $(OBJU) Driver.f95

magerr: link magerr.f95 $(OBJU) ${SRCM} Makefile
	$(FC) $(FFLAGP) -o magerr $(OBJU) ${SRCM} magerr.f95
	rm GiMeObj.f95

Modules: _$(SSIM).so Makefile
	echo "Modules have been built"

${SSIM}: ${FPPM} $(OBJA) $(OBJM) Makefile
	f90wrap -m ${SSIM} ${FPPM} -k ${KIND_MAP} -v
	mkdir -p build
	ln -s -f $(CURDIR)/*.mod build/
	f2py-f90wrap $(F2PYFLAGS) -c -m _${SSIM} $(OBJA) $(SRCM) f90wrap*.f90
	\rm f90wrap*.f90 *.fpp *.o *.mod 

