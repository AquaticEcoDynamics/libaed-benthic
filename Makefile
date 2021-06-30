###############################################################################
#                                                                             #
# Makefile to build the aed water quality library                             #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2021 -  The University of Western Australia               #
#                                                                             #
#   GLM is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   GLM is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

srcdir=src
incdir=include

ifeq ($(F90),)
  F90=gfortran
endif

ifeq ($(SINGLE),true)
  TARGET=lib/libaed-benthic_s.a
  objdir=obj_s
  moddir=mod_s
else
  TARGET=lib/libaed-benthic.a
  objdir=obj
  moddir=mod
endif

INCLUDES=-I../libaed-water/${incdir}  -I../libaed-water/${moddir}

ifeq ($(F90),ifort)
  INCLUDES+=-I/opt/intel/include
  DEBUG_FFLAGS=-g -traceback
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -warn all -module ${moddir} -static-intel -mp1 -stand f08 -warn nounused $(DEFINES) $(INCLUDES)
  FFLAGS+=-module ../libaed-water/mod
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-check all -check noarg_temp_created
  endif
  ifeq ($(SINGLE),true)
    FFLAGS+=-real-size 32
  else
    FFLAGS+=-real-size 64
  endif
else ifeq ($(F90),flang)
  DEBUG_FFLAGS=-g
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -module ${moddir} $(DEFINES) $(INCLUDES)
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-Mbounds
  endif
  FFLAGS+=-r8
else
  DEBUG_FFLAGS=-g -fbacktrace
  OPT_FFLAGS=-O3
  FFLAGS=-fPIC -Wall -J ${moddir} -ffree-line-length-none -std=f2008 $(DEFINES) $(INCLUDES)
  FFLAGS+=-fall-intrinsics -Wno-unused -Wno-unused-dummy-argument -fno-range-check -Wno-integer-division
  ifeq ($(WITH_CHECKS),true)
    FFLAGS+=-fcheck=all
  endif
  FFLAGS+=-fdefault-real-8 -fdefault-double-8
endif

ifeq ($(DEBUG),true)
  DEBUG_CFLAGS=-g
  OPT_CFLAGS=
  OPT_FFLAGS=
else
  DEBUG_FFLAGS=
  DEBUG_CFLAGS=
  # OPT_CFLAGS=-O4 -Ofast -frounding-math
  OPT_CFLAGS=-O3
  # OPT_CFLAGS=
  # OPT_FFLAGS=
endif

ifeq ($(SINGLE),true)
  FFLAGS += -DSINGLE=1
endif


FFLAGS+=$(DEBUG_FFLAGS) $(OPT_FFLAGS)

OBJS=${objdir}/aed_bivalve.o \
     ${objdir}/aed_habitat_benthic.o \
     ${objdir}/aed_macroalgae.o \
     ${objdir}/aed_macroalgae2.o \
     ${objdir}/aed_macrophyte.o \
     ${objdir}/aed_benthic.o

all: $(TARGET)

lib:
	@mkdir lib

${moddir}:
	@mkdir ${moddir}

${objdir}:
	@mkdir ${objdir}

${TARGET}: ${objdir} ${moddir} ${OBJS} lib
	ar rv $@ ${OBJS}
	ranlib $@

clean: ${objdir}
	@touch ${objdir}/1.o 1.i90
	@/bin/rm ${objdir}/*.o *.i90

distclean: clean
	@touch lib mod_s mod
	@/bin/rm -rf lib
	@/bin/rm -rf obj obj_s
	@/bin/rm -rf mod mod_s

${objdir}/%.o: ${srcdir}/%.F90 ../libaed-water/include/aed.h
	$(F90) $(FFLAGS) -g -c $< -o $@
