!###############################################################################
!#                                                                             #
!# aed_benthic.F90                                                             #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2013 - 2022 -  The University of Western Australia               #
!#                                                                             #
!#   AED is free software: you can redistribute it and/or modify               #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED is distributed in the hope that it will be useful,                    #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Created May 2020                                                            #
!#                                                                             #
!###############################################################################

#include "aed.h"


!###############################################################################
MODULE aed_benthic
!-------------------------------------------------------------------------------
   USE aed_core, ONLY : aed_model_data_t

   USE aed_bivalve
   USE aed_macrophyte
   USE aed_macroalgae
   USE aed_macroalgae2
   USE aed_habitat_benthic

   IMPLICIT NONE

   !#---------------------------------------------------------------------------

   PRIVATE   !# By default make everything private

   PUBLIC aed_new_ben_model, aed_print_ben_version

   !#---------------------------------------------------------------------------

CONTAINS
!===============================================================================


!###############################################################################
FUNCTION aed_new_ben_model(modelname) RESULT(model)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*),INTENT(in) :: modelname
!
!LOCALS
   CLASS (aed_model_data_t),POINTER :: model
   CHARACTER(len=4) :: prefix
!
!-------------------------------------------------------------------------------
!BEGIN
   NULLIFY(model)

   SELECT CASE (modelname)
      CASE ('aed_bivalve');         prefix = 'BIV'; ALLOCATE(aed_bivalve_data_t::model)
      CASE ('aed_macroalgae');      prefix = 'MAG'; ALLOCATE(aed_macroalgae_data_t::model)
      CASE ('aed_macroalgae2');     prefix = 'MA2'; ALLOCATE(aed_macroalgae2_data_t::model)
      CASE ('aed_macrophyte');      prefix = 'MAC'; ALLOCATE(aed_macrophyte_data_t::model)
      CASE ('aed_habitat_benthic'); prefix = 'HAB'; ALLOCATE(aed_habitat_benthic_data_t::model)
   END SELECT

   IF (ASSOCIATED(model)) THEN
      model%aed_model_name = modelname
      model%aed_model_prefix = prefix
   ENDIF
END FUNCTION aed_new_ben_model
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
SUBROUTINE aed_print_ben_version
!-------------------------------------------------------------------------------
!BEGIN
   print*,"    libaed-ben version ", TRIM(AED_VERSION)
#ifdef __INTEL_COMPILER
   print*,"    libaed built using intel fortran version ", __INTEL_COMPILER
#else
# ifdef __PGI
   print*,"    libaed built using pgfortran version ", __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
#  ifdef __GNUC__
    print*,"    libaed built using gfortran version ", __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
#  else
#   ifdef __clang__
     print*,"    libaed built using flang version ", __clang_major__, '.', __clang_minor__, '.', __clang_patchlevel__
#   else
     print*,"    libaed built using unknown fortran version "
#   endif
#  endif
# endif
#endif
END SUBROUTINE aed_print_ben_version
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!===============================================================================
END MODULE aed_benthic
