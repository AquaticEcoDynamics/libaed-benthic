###############################################################################
#                                                                             #
# Makefile to build libaed-benthic                                            #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2024 -  The University of Western Australia               #
#                                                                             #
#   AED is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   AED is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

LIBAEDBEN=aed-benthic
OUTLIB=lib$(LIBAEDBEN)

INCLUDES=-I../libaed-water/${incdir}  -I../libaed-water/${moddir}

include ../libaed-water/make_defs.inc

OBJS=${objdir}/aed_bivalve.o \
     ${objdir}/aed_habitat_benthic.o \
     ${objdir}/aed_habitat_ruppia.o \
     ${objdir}/aed_habitat_chara.o \
     ${objdir}/aed_habitat_galaxiid.o \
     ${objdir}/aed_habitat_seagrass.o \
     ${objdir}/aed_macroalgae.o \
     ${objdir}/aed_macroalgae2.o \
     ${objdir}/aed_macrophyte.o \
     ${objdir}/aed_benthic.o

include ../libaed-water/make_rules.inc
