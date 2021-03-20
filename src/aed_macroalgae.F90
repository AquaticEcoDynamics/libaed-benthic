!###############################################################################
!#                                                                             #
!# aed_macroalgae.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2017 - 2020 -  The University of Western Australia               #
!#                                                                             #
!#   AED2+ is free software: you can redistribute it and/or modify             #
!#   it under the terms of the GNU General Public License as published by      #
!#   the Free Software Foundation, either version 3 of the License, or         #
!#   (at your option) any later version.                                       #
!#                                                                             #
!#   AED2+ is distributed in the hope that it will be useful,                  #
!#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
!#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
!#   GNU General Public License for more details.                              #
!#                                                                             #
!#   You should have received a copy of the GNU General Public License         #
!#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
!#   For use in a commercial seeting please contact the authors.               #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Originally created August 2017 by Matthew Hipsey, UWA                       #
!# Follow updates @ https://github.com/AquaticEcoDynamics/libaed2-plus         #
!#                                                                             #
!###############################################################################

#include "aed+.h"

MODULE aed_macroalgae
!-------------------------------------------------------------------------------
!  aed_macroalgae --- macroalgae biogeochemical model
!-------------------------------------------------------------------------------
   USE aed_core
   USE aed_util,ONLY : find_free_lun, &
                        exp_integral, &
                        aed_bio_temp_function, &
                        fTemp_function,fSal_function, &
                        water_viscosity, in_zone_set
   USE aed_bio_utils

   IMPLICIT NONE

   PRIVATE   ! By default make everything private

   PUBLIC aed_macroalgae_data_t


   TYPE,extends(aed_model_data_t) :: aed_macroalgae_data_t
      !# Variable identifiers
      INTEGER,ALLOCATABLE :: id_p(:),  id_pben(:)
      INTEGER,ALLOCATABLE :: id_in(:), id_inben(:)
      INTEGER,ALLOCATABLE :: id_ip(:), id_ipben(:)
      INTEGER,ALLOCATABLE :: id_rho(:),id_vvel(:)
      INTEGER,ALLOCATABLE :: id_N2P(:),id_C2N(:), id_C2P(:)
      INTEGER,ALLOCATABLE :: id_N2P_ben(:),id_C2N_ben(:), id_C2P_ben(:)
      INTEGER,ALLOCATABLE :: id_fT(:), id_fI(:), id_fNit(:), &
                             id_fPho(:), id_fSal(:)
      INTEGER,ALLOCATABLE :: id_fT_ben(:), id_fI_ben(:), id_fNit_ben(:), &
                             id_fPho_ben(:), id_fSal_ben(:)
      INTEGER :: id_Pexr, id_Pmor, id_Pupttarget(1:2)
      INTEGER :: id_Nexr, id_Nmor, id_Nupttarget(1:4)
      INTEGER :: id_Cexr, id_Cmor, id_Cupttarget
      INTEGER :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      INTEGER :: id_DOupttarget
      INTEGER :: id_par, id_I_0, id_extc, id_taub, id_sedzone
      INTEGER :: id_tem, id_sal, id_dz, id_dens
      INTEGER :: id_GPP, id_NMP, id_PPR, id_NPR, id_dPAR, id_dHGT
      INTEGER :: id_TMALG, id_dEXTC, id_TIN, id_TIP, id_MPB, id_d_MPB, id_d_BPP
      INTEGER :: id_NUP, id_PUP, id_CUP
      INTEGER :: id_mhsi
      INTEGER :: id_mag_ben, id_min_ben, id_mip_ben
      INTEGER :: id_nup_ben, id_pup_ben, id_gpp_ben, id_rsp_ben, id_nmp_ben
      INTEGER :: id_slough_trig, id_tem_avg, id_tau_avg, id_par_avg, id_slg_ben

      !# Model parameters and options
      TYPE(phyto_data),DIMENSION(:),ALLOCATABLE :: malgs
      INTEGER  :: num_malgae
      LOGICAL  :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL  :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL  :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL  :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER  :: nnup, npup
      INTEGER  :: n_zones

      AED_REAL, ALLOCATABLE :: active_zones(:)
      AED_REAL :: min_rho,max_rho
      AED_REAL :: slough_stress
      LOGICAL  :: simMalgFeedback
      INTEGER  :: simSloughing
      INTEGER  :: simMalgHSI
      INTEGER  :: simCGM

     CONTAINS
         PROCEDURE :: define            => aed_define_macroalgae
         PROCEDURE :: initialize        => aed_initialize_macroalgae
         PROCEDURE :: calculate         => aed_calculate_macroalgae
         PROCEDURE :: calculate_benthic => aed_calculate_benthic_macroalgae
         PROCEDURE :: mobility          => aed_mobility_macroalgae
         PROCEDURE :: light_extinction  => aed_light_extinction_macroalgae
         PROCEDURE :: bio_drag          => aed_bio_drag_macroalgae
        !PROCEDURE :: delete            => aed_delete_macroalgae

   END TYPE

   ! MODULE GLOBALS
   AED_REAL, PARAMETER :: DTday       = (15./60.)/24.  ! 15 minutes
   AED_REAL, PARAMETER :: StrAvgTime  = 2.0/24.0       !  2 hours
   AED_REAL, PARAMETER :: LgtAvgTime  = 0.5            ! 12 hours
   AED_REAL, PARAMETER :: TempAvgTime = 1.0            !  1 day
   AED_REAL :: dtlim = 900
   AED_REAL :: macroHgt, macroExt = zero_
   LOGICAL  :: extra_diag = .false.
   INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE aed_macroalgae_load_params(data, dbase, count, list, settling, resuspension, tau_0)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: settling(*)
   AED_REAL,INTENT(in)         :: resuspension(*)
   AED_REAL,INTENT(in)         :: tau_0(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(phyto_nml_data) :: pd(MAX_PHYTO_TYPES)
   NAMELIST /malgae_data/ pd
!-------------------------------------------------------------------------------
!BEGIN
    tfil = find_free_lun()
    open(tfil,file=dbase, status='OLD', iostat=status)
    IF (status /= 0) STOP 'Cannot open malgae_data namelist file: ' !,dbase
    read(tfil,nml=malgae_data,iostat=status)
    close(tfil)
    IF (status /= 0) STOP 'Error reading namelist malgae_data'

    data%simCGM = 0

    data%num_malgae = count
    ALLOCATE(data%malgs(count))
    ALLOCATE(data%id_p(count)) ; data%id_p(:) = 0
    ALLOCATE(data%id_in(count)) ; data%id_in(:) = 0
    ALLOCATE(data%id_ip(count)) ; data%id_ip(:) = 0
    ALLOCATE(data%id_rho(count)) ; data%id_rho(:) = 0
    ALLOCATE(data%id_pben(count)) ; data%id_p(:) = 0
    ALLOCATE(data%id_inben(count)) ; data%id_inben(:) = 0
    ALLOCATE(data%id_ipben(count)) ; data%id_ipben(:) = 0
    IF (diag_level>9) THEN
       ALLOCATE(data%id_c2p(count)) ; data%id_c2p(:) = 0
       ALLOCATE(data%id_n2p(count)) ; data%id_n2p(:) = 0
       ALLOCATE(data%id_c2n(count)) ; data%id_c2n(:) = 0
       ALLOCATE(data%id_fT(count)) ; data%id_fT(:) = 0
       ALLOCATE(data%id_fI(count)) ; data%id_fI(:) = 0
       ALLOCATE(data%id_fNit(count)) ; data%id_fNit(:) = 0
       ALLOCATE(data%id_fPho(count)) ; data%id_fPho(:) = 0
       ALLOCATE(data%id_fSal(count)) ; data%id_fSal(:) = 0
       ALLOCATE(data%id_vvel(count)) ; data%id_vvel(:) = 0
       ALLOCATE(data%id_fT_ben(count)) ; data%id_fT_ben(:) = 0
       ALLOCATE(data%id_fI_ben(count)) ; data%id_fI_ben(:) = 0
       ALLOCATE(data%id_fNit_ben(count)) ; data%id_fNit_ben(:) = 0
       ALLOCATE(data%id_fPho_ben(count)) ; data%id_fPho_ben(:) = 0
       ALLOCATE(data%id_fSal_ben(count)) ; data%id_fSal_ben(:) = 0
       ALLOCATE(data%id_c2p_ben(count)) ; data%id_c2p_ben(:) = 0
       ALLOCATE(data%id_n2p_ben(count)) ; data%id_n2p_ben(:) = 0
       ALLOCATE(data%id_c2n_ben(count)) ; data%id_c2n_ben(:) = 0
    ENDIF

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%malgs(i)%p_name       = pd(list(i))%p_name
       data%malgs(i)%p0           = pd(list(i))%p0
       data%malgs(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%malgs(i)%settling     = settling(i)
       data%malgs(i)%resuspension = resuspension(i)
       data%malgs(i)%tau_0        = tau_0(i)
       data%malgs(i)%Xcc          = pd(list(i))%Xcc
       data%malgs(i)%R_growth     = pd(list(i))%R_growth/secs_per_day
       data%malgs(i)%fT_Method    = pd(list(i))%fT_Method
       data%malgs(i)%theta_growth = pd(list(i))%theta_growth
       data%malgs(i)%T_std        = pd(list(i))%T_std
       data%malgs(i)%T_opt        = pd(list(i))%T_opt
       data%malgs(i)%T_max        = pd(list(i))%T_max
       data%malgs(i)%lightModel   = pd(list(i))%lightModel
       data%malgs(i)%I_K          = pd(list(i))%I_K
       data%malgs(i)%I_S          = pd(list(i))%I_S
       data%malgs(i)%KePHY        = pd(list(i))%KePHY
       data%malgs(i)%f_pr         = pd(list(i))%f_pr
       data%malgs(i)%R_resp       = pd(list(i))%R_resp/secs_per_day
       data%malgs(i)%theta_resp   = pd(list(i))%theta_resp
       data%malgs(i)%k_fres       = pd(list(i))%k_fres
       data%malgs(i)%k_fdom       = pd(list(i))%k_fdom
       data%malgs(i)%salTol       = pd(list(i))%salTol
       data%malgs(i)%S_bep        = pd(list(i))%S_bep
       data%malgs(i)%S_maxsp      = pd(list(i))%S_maxsp
       data%malgs(i)%S_opt        = pd(list(i))%S_opt
       data%malgs(i)%simDINUptake = pd(list(i))%simDINUptake
       data%malgs(i)%simDONUptake = pd(list(i))%simDONUptake
       data%malgs(i)%simNFixation = pd(list(i))%simNFixation
       data%malgs(i)%simINDynamics= pd(list(i))%simINDynamics
       data%malgs(i)%N_o          = pd(list(i))%N_o
       data%malgs(i)%K_N          = pd(list(i))%K_N
       data%malgs(i)%X_ncon       = pd(list(i))%X_ncon
       data%malgs(i)%X_nmin       = pd(list(i))%X_nmin
       data%malgs(i)%X_nmax       = pd(list(i))%X_nmax
       data%malgs(i)%R_nuptake    = pd(list(i))%R_nuptake/secs_per_day
       data%malgs(i)%k_nfix       = pd(list(i))%k_nfix
       data%malgs(i)%R_nfix       = pd(list(i))%R_nfix/secs_per_day
       data%malgs(i)%simDIPUptake = pd(list(i))%simDIPUptake
       data%malgs(i)%simIPDynamics= pd(list(i))%simIPDynamics
       data%malgs(i)%P_0          = pd(list(i))%P_0
       data%malgs(i)%K_P          = pd(list(i))%K_P
       data%malgs(i)%X_pcon       = pd(list(i))%X_pcon
       data%malgs(i)%X_pmin       = pd(list(i))%X_pmin
       data%malgs(i)%X_pmax       = pd(list(i))%X_pmax
       data%malgs(i)%R_puptake    = pd(list(i))%R_puptake/secs_per_day
       data%malgs(i)%simSiUptake  = pd(list(i))%simSiUptake
       data%malgs(i)%Si_0         = pd(list(i))%Si_0
       data%malgs(i)%K_Si         = pd(list(i))%K_Si
       data%malgs(i)%X_sicon      = pd(list(i))%X_sicon

       data%malgs(i)%c1           = 0.0124/60.   ! From Chung et al (2014)
       data%malgs(i)%c3           = 0.0230/60.   !  "
       data%malgs(i)%f1           = 0.675        ! Ross and Sharples (2007)
       data%malgs(i)%f2           = 0.750        !  "
       data%malgs(i)%d_phy        = 1e-5

       ! Register group as a state variable
       data%id_p(i) = aed_define_variable(                                    &
                              TRIM(data%malgs(i)%p_name),                      &
                              'mmol/m**3',                                     &
                              'macroalgae '//TRIM(data%malgs(i)%p_name),       &
                              pd(list(i))%p0,                                  &
                             !pd(list(i))%p_initial,                           &
                              minimum=pd(list(i))%p0,                          &
                              mobility = data%malgs(i)%w_p)

        IF(TRIM(data%malgs(i)%p_name)== 'cgm') THEN
          data%simCGM = i
          IF(data%simSloughing>0) data%simSloughing = 2    ! If sloughing is requested force to 2 for now
        ENDIF

       ! Register rho (internal density) as a state variable, if required
       IF (data%malgs(i)%settling == _MOB_STOKES_) THEN
           data%id_rho(i) = aed_define_variable(                              &
                              TRIM(data%malgs(i)%p_name)//'_rho',              &
                              'kg/m**3',                                       &
                             'macroalgae '//TRIM(data%malgs(i)%p_name)//'_rho',&
                              (data%min_rho+data%max_rho)/2.,                  &
                              minimum=data%min_rho,                            &
                              mobility = data%malgs(i)%w_p)
       ENDIF

       ! Register internal nitrogen group as a state variable, if required
       IF (data%malgs(i)%simINDynamics /= 0) THEN
          IF(data%malgs(i)%simINDynamics == 1)THEN
            minNut = data%malgs(i)%p0*data%malgs(i)%X_ncon
          ELSE
            minNut = data%malgs(i)%p0*data%malgs(i)%X_nmin
          ENDIF
          ! Register IN group as a state variable
          data%id_in(i) = aed_define_variable(                                &
                              TRIM(data%malgs(i)%p_name)//'_IN',               &
                              'mmol/m**3',                                     &
                              'macroalgae '//TRIM(data%malgs(i)%p_name)//'_IN',&
                              pd(list(i))%p0*data%malgs(i)%X_ncon,             &
                             !pd(list(i))%p_initial*data%malgs(i)%X_ncon,      &
                              minimum=minNut,                                  &
                              mobility = data%malgs(i)%w_p)
       ENDIF

       ! Register internal phosphorus group as a state variable, if required
       IF (data%malgs(i)%simIPDynamics /= 0) THEN
          minNut = data%malgs(i)%p0*data%malgs(i)%X_pmin
          ! Register IP group as a state variable
          data%id_ip(i) = aed_define_variable(                                &
                              TRIM(data%malgs(i)%p_name)//'_IP',               &
                              'mmol/m**3',                                     &
                              'macroalgae '//TRIM(data%malgs(i)%p_name)//'_IP',&
                              !pd(list(i))%p0*data%malgs(i)%X_pcon,             &
                              MAX(pd(list(i))%p0,pd(list(i))%p_initial)*data%malgs(i)%X_pcon,&
                              minimum=minNut,                                  &
                              mobility = data%malgs(i)%w_p)
       ENDIF

       ! Group specific diagnostic variables
       IF (diag_level>9) THEN
          data%id_fI(i)   = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_fI', '-', 'fI (0-1)')
          data%id_fNit(i) = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_fNit', '-', 'fNit (0-1)')
          data%id_fPho(i) = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_fPho', '-', 'fPho (0-1)')
          data%id_fT(i)   = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_fT', '-', 'fT (>0)')
          data%id_fSal(i) = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_fSal', '-', 'fSal (>1)')
          data%id_c2p(i)  = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2p', '-', 'C:P')
          data%id_c2n(i)  = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2n', '-', 'C:N')
          data%id_n2p(i)  = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_n2p', '-', 'N:P')

          ! Register vertical velocity diagnostic, where relevant
          IF (data%malgs(i)%settling == _MOB_STOKES_ ) THEN
            data%id_vvel(i) = aed_define_diag_variable( TRIM(data%malgs(i)%p_name)//'_vvel', 'm/s', 'vertical velocity')
          ENDIF
       ENDIF

       ! Allocate benthic variables if the group has attached/benthic growth form
       IF (data%malgs(i)%settling == _MOB_ATTACHED_) THEN
         data%id_pben(i) = aed_define_sheet_variable(                         &
                          TRIM(data%malgs(i)%p_name)//'_ben',                  &
                          'mmolC/m**2',                                        &
                          'macroalgae '//TRIM(data%malgs(i)%p_name),           &
                          pd(list(i))%p_initial,                               &
                          minimum=pd(list(i))%p0,                              &
                          maximum=1e8 )

         IF (data%malgs(i)%simIPDynamics /= 0) THEN
            minNut = data%malgs(i)%p0*data%malgs(i)%X_pmin
            ! Register IP group as a state variable
            data%id_ipben(i) = aed_define_sheet_variable(                     &
                          TRIM(data%malgs(i)%p_name)//'_IP_ben',               &
                          'mmol/m**2',                                         &
                          'macroalgae '//TRIM(data%malgs(i)%p_name)//'_IP_ben',&
                          pd(list(i))%p_initial*data%malgs(i)%X_pcon,          &
                          minimum=minNut)
         ENDIF
         IF (data%malgs(i)%simINDynamics /= 0) THEN
           IF(data%malgs(i)%simINDynamics == 1)THEN
            minNut = data%malgs(i)%p0*data%malgs(i)%X_ncon
           ELSE
            minNut = data%malgs(i)%p0*data%malgs(i)%X_nmin
           ENDIF
           ! Register IN group as a state variable
           data%id_inben(i) = aed_define_sheet_variable(                      &
                          TRIM(data%malgs(i)%p_name)//'_IN_ben',               &
                          'mmol/m**2',                                         &
                          'macroalgae '//TRIM(data%malgs(i)%p_name)//'_IN_ben',&
                          pd(list(i))%p_initial*data%malgs(i)%X_ncon,          &
                          minimum=minNut)
         ENDIF
         ! Group specific diagnostics for the benthic variables
         IF (diag_level>9) THEN
          data%id_fI_ben(i)   = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fI_ben', '-', 'fI (0-1)')
          data%id_fNit_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fNit_ben', '-', 'fNit (0-1)')
          data%id_fPho_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fPho_ben', '-', 'fPho (0-1)')
          data%id_fT_ben(i)   = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fT_ben', '-', 'fT (>0)')
          data%id_fSal_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fSal_ben', '-', 'fSal (-)')
          data%id_c2p_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2p_ben', '-', 'C:P')
          data%id_c2n_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2n_ben', '-', 'C:N')
          data%id_n2p_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_n2p_ben', '-', 'N:P')
         ENDIF

       ENDIF
    ENDDO
END SUBROUTINE aed_macroalgae_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_macroalgae(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macroalgae biogeochemical model
!
!  Here, the aed_macroalgae namelist is read and the variables exported
!  by the model are registered with AED2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER            :: status, i

   !  %% NAMELIST VARS
   INTEGER            :: num_malgae=0
   INTEGER            :: the_malgae(MAX_PHYTO_TYPES)=0
   INTEGER            :: settling(MAX_PHYTO_TYPES)= _MOB_CONST_
   AED_REAL           :: resuspension(MAX_PHYTO_TYPES)=0.
   AED_REAL           :: tau_0(MAX_PHYTO_TYPES)=0.1
   CHARACTER(len=64)  :: p_excretion_target_variable=''
   CHARACTER(len=64)  :: p_mortality_target_variable=''
   CHARACTER(len=64)  :: p1_uptake_target_variable=''
   CHARACTER(len=64)  :: p2_uptake_target_variable=''
   CHARACTER(len=64)  :: n_excretion_target_variable=''
   CHARACTER(len=64)  :: n_mortality_target_variable=''
   CHARACTER(len=64)  :: n1_uptake_target_variable=''
   CHARACTER(len=64)  :: n2_uptake_target_variable=''
   CHARACTER(len=64)  :: n3_uptake_target_variable=''
   CHARACTER(len=64)  :: n4_uptake_target_variable=''
   CHARACTER(len=64)  :: c_excretion_target_variable=''
   CHARACTER(len=64)  :: c_mortality_target_variable=''
   CHARACTER(len=64)  :: c_uptake_target_variable=''
   CHARACTER(len=64)  :: do_uptake_target_variable=''
   CHARACTER(len=128) :: dbase='aed_malgae_pars.nml'
   AED_REAL           :: zerolimitfudgefactor = 15.*60.
   AED_REAL           :: min_rho = 900.
   AED_REAL           :: max_rho = 1200.
   AED_REAL           :: slough_stress = -1.0
   INTEGER            :: simMalgHSI = 0
   INTEGER            :: simSloughing = 0
   INTEGER            :: n_zones = 0
   INTEGER            :: active_zones(1000) = 0
   LOGICAL            :: simMalgFeedback = .true.
   LOGICAL            :: extra_debug = .false.
   !  %% END NAMELIST VARS

   NAMELIST /aed_macroalgae/ num_malgae, the_malgae, settling, resuspension,  &
                    p_excretion_target_variable,p_mortality_target_variable,   &
                     p1_uptake_target_variable, p2_uptake_target_variable,     &
                    n_excretion_target_variable,n_mortality_target_variable,   &
                     n1_uptake_target_variable,n2_uptake_target_variable,      &
                     n3_uptake_target_variable,n4_uptake_target_variable,      &
                    c_excretion_target_variable,c_mortality_target_variable,   &
                      c_uptake_target_variable, do_uptake_target_variable,     &
                    dbase, zerolimitfudgefactor,slough_stress,simSloughing,    &
                     simMalgHSI, n_zones, active_zones, simMalgFeedback,       &
                    extra_debug, extra_diag, diag_level, tau_0, dtlim
!-----------------------------------------------------------------------
!BEGIN

   print *,"        aed_macroalgae initialization"

   dtlim = zerolimitfudgefactor

   ! Read the namelist, and set module parameters
   read(namlst,nml=aed_macroalgae,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_macroalgae'
   IF( extra_debug )  extra_diag = .true.       ! legacy use of extra_debug
   IF( extra_debug )  diag_level = 10           ! legacy use of extra_debug
   IF( extra_diag )   diag_level = 10           ! legacy use of extra_debug
   data%min_rho = min_rho ; data%max_rho = max_rho
   data%simMalgHSI = simMalgHSI
   data%n_zones = n_zones
   IF( n_zones>0 ) THEN
     ALLOCATE(data%active_zones(n_zones))
     DO i=1,n_zones
       data%active_zones(i) = active_zones(i)
     ENDDO
   ENDIF

   data%simSloughing = simSloughing
   data%slough_stress = slough_stress

   data%simMalgFeedback = simMalgFeedback
   PRINT *,'          NOTE - macroalgae feedbacks to water column properties: ',simMalgFeedback

   ! Store parameter values in a local malgae strcutured type
   ! Note: all rates must be provided in values per day,
   !       but are converted in here to rates per second
   CALL aed_macroalgae_load_params(data,dbase,num_malgae,the_malgae,settling,resuspension,tau_0)

   CALL aed_bio_temp_function(data%num_malgae,             &
                               data%malgs%theta_growth,     &
                               data%malgs%T_std,            &
                               data%malgs%T_opt,            &
                               data%malgs%T_max,            &
                               data%malgs%aTn,              &
                               data%malgs%bTn,              &
                               data%malgs%kTn,              &
                               data%malgs%p_name)

   ! Register link to nutrient pools, if variable names are provided in namelist
   data%do_Pexc = p_excretion_target_variable .NE. ''
   IF (data%do_Pexc) THEN
     data%id_Pexr = aed_locate_variable( p_excretion_target_variable)
   ENDIF
   data%do_Nexc = n_excretion_target_variable .NE. ''
   IF (data%do_Nexc) THEN
     data%id_Nexr = aed_locate_variable( n_excretion_target_variable)
   ENDIF
   data%do_Cexc = c_excretion_target_variable .NE. ''
   IF (data%do_Cexc) THEN
     data%id_Cexr = aed_locate_variable( c_excretion_target_variable)
   ENDIF
   data%do_Pmort = p_mortality_target_variable .NE. ''
   IF (data%do_Pmort) THEN
     data%id_Pmor = aed_locate_variable( p_mortality_target_variable)
   ENDIF
   data%do_Nmort = n_mortality_target_variable .NE. ''
   IF (data%do_Nmort) THEN
     data%id_Nmor = aed_locate_variable( n_mortality_target_variable)
   ENDIF
   data%do_Cmort = c_mortality_target_variable .NE. ''
   IF (data%do_Cmort) THEN
     data%id_Cmor = aed_locate_variable( c_mortality_target_variable)
   ENDIF

   data%npup = 0
   IF (p1_uptake_target_variable .NE. '') data%npup = 1
   IF (p2_uptake_target_variable .NE. '') data%npup = 2
   data%do_Puptake = .FALSE.
   IF (data%npup>0) data%do_Puptake=.TRUE.
   IF (data%do_Puptake) THEN
     IF (data%npup>0) data%id_Pupttarget(ifrp) = aed_locate_variable(p1_uptake_target_variable) !; ifrp=1  ! CAB - now constants in bio utils
     IF (data%npup>1) data%id_Pupttarget(idop) = aed_locate_variable(p2_uptake_target_variable) !; idop=2  ! CAB - now constants in bio utils
   ENDIF
   data%nnup = 0
   IF (n1_uptake_target_variable .NE. '') data%nnup = 1
   IF (n2_uptake_target_variable .NE. '') data%nnup = 2
   IF (n3_uptake_target_variable .NE. '') data%nnup = 3
   IF (n4_uptake_target_variable .NE. '') data%nnup = 4
   data%do_Nuptake = .false.
   IF (data%nnup>0) data%do_Nuptake=.true.
   IF (data%do_Nuptake) THEN
     IF (data%nnup>0) data%id_Nupttarget(ino3) = aed_locate_variable( n1_uptake_target_variable) !; ino3=1  ! CAB - now constants in bio utils
     IF (data%nnup>1) data%id_Nupttarget(inh4) = aed_locate_variable( n2_uptake_target_variable) !; inh4=2  ! CAB - now constants in bio utils
     IF (data%nnup>2) data%id_Nupttarget(idon) = aed_locate_variable( n3_uptake_target_variable) !; idon=3  ! CAB - now constants in bio utils
     IF (data%nnup>3) data%id_Nupttarget(in2)  = aed_locate_variable( n4_uptake_target_variable) !; in2 =4  ! CAB - now constants in bio utils
   ENDIF
   data%do_Cuptake = c_uptake_target_variable .NE. ''
   IF (data%do_Cuptake) THEN
     data%id_Cupttarget = aed_locate_variable( c_uptake_target_variable)
   ENDIF
   data%do_DOuptake = do_uptake_target_variable .NE. ''
   IF (data%do_DOuptake) THEN
     data%id_DOupttarget = aed_locate_variable( do_uptake_target_variable)
   ENDIF

   ! Register diagnostic variables
   IF (diag_level>0) THEN
     data%id_TMALG   = aed_define_diag_variable('tmalg','g DW/m**2', 'MAG: total macroalgal biomass')
     data%id_TIN     = aed_define_diag_variable('in','mmol/m**3', 'MAG: total macroalgal nitrogen')
     data%id_TIP     = aed_define_diag_variable('ip','mmol/m**3', 'MAG: total macroalgal phosphorus')
     data%id_mag_ben = aed_define_sheet_diag_variable('mag_ben','mmol/m**2/d', 'BEN MAG: total C biomass')
     data%id_min_ben = aed_define_sheet_diag_variable('in_ben','mmol/m**2/d', 'BEN MAG: total N biomass')
     data%id_mip_ben = aed_define_sheet_diag_variable('ip_ben','mmol/m**2/d', 'BEN MAG: total P biomass')
   ENDIF
   IF (diag_level>1) THEN
     data%id_GPP     = aed_define_diag_variable('gpp','mmol/m**3/d', 'MAG: macroalgal gross primary production')
     data%id_PUP     = aed_define_diag_variable('pup','mmol/m**3/d', 'MAG: macroalgal phosphorous uptake')
     data%id_NUP     = aed_define_diag_variable('nup','mmol/m**3/d','MAG: macroalgal nitrogen uptake')
     data%id_NMP     = aed_define_diag_variable('nmp','mmol/m**3/d',  'net macroalgal production')
     data%id_gpp_ben = aed_define_sheet_diag_variable('gpp_ben','/d', 'BEN MAG: macroalgal gross primary production')
     data%id_rsp_ben = aed_define_sheet_diag_variable('rsp_ben','/d', 'BEN MAG: macroalgal respiration')
     data%id_pup_ben = aed_define_sheet_diag_variable('pup_ben','mmol/m**2/d', 'BEN MAG: macroalgal phosphorous uptake')
     data%id_nup_ben = aed_define_sheet_diag_variable('nup_ben','mmol/m**2/d', 'BEN MAG: macroalgal nitrogen uptake')
     data%id_nmp_ben = aed_define_sheet_diag_variable('nmp_ben','mmol/m**2/d', 'BEN MAG: net macroalgal production')
     data%id_slg_ben = aed_define_sheet_diag_variable('slg_ben','mmol/m**2/d', 'BEN MAG: sloughing rate')
   ENDIF

   IF ( simMalgHSI>0 ) &
     data%id_mhsi  = aed_define_sheet_diag_variable('HSI','-', 'MAG: macroalgae habitat suitability')

   IF ( data%simCGM >0 ) THEN
     data%id_slough_trig = aed_define_sheet_diag_variable('cgm_sltg','-', 'MAG: cgm slough trigger')
     data%id_tem_avg = aed_define_sheet_diag_variable('cgm_tavg','-', 'MAG: cgm temperature average')
     data%id_par_avg = aed_define_sheet_diag_variable('cgm_lavg','-', 'MAG: cgm light average')
     data%id_tau_avg = aed_define_sheet_diag_variable('cgm_savg','-', 'MAG: cgm stress average')
   ENDIF

   IF (diag_level>9) THEN
     data%id_dEXTC = aed_define_diag_variable('extc','/m','MAG: extinction due to macroalgae')
     data%id_dPAR  = aed_define_sheet_diag_variable('parc','%', 'MAG: PAR reaching the bottom canopy')
     data%id_dHGT  = aed_define_sheet_diag_variable('hgtc','m', 'MAG: height of macroalgal canopy')
   ENDIF

   ! Register the required environmental dependencies
   data%id_dz      = aed_locate_global('layer_ht')
   data%id_tem     = aed_locate_global('temperature')
   data%id_sal     = aed_locate_global('salinity')
   data%id_dens    = aed_locate_global('density')
   data%id_extc    = aed_locate_global('extc_coef')
   data%id_par     = aed_locate_global('par')
   data%id_I_0     = aed_locate_global_sheet('par_sf')
   data%id_taub    = aed_locate_global_sheet('taub')
   data%id_sedzone = aed_locate_global_sheet('sed_zone')

END SUBROUTINE aed_define_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_macroalgae(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to initialize bottom diagnostics, and other checks
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mag_i
   AED_REAL :: matz, malg_density
!-------------------------------------------------------------------------------
!BEGIN

   ! Check to ensure this zone is colonisable
   matz = _STATE_VAR_S_(data%id_sedzone)
   IF ( .NOT. in_zone_set(matz, data%active_zones) ) THEN
     DO mag_i=1,data%num_malgae
       IF (data%malgs(mag_i)%settling == _MOB_ATTACHED_) THEN
         _STATE_VAR_S_(data%id_pben(mag_i)) = zero_
         IF (data%malgs(mag_i)%simINDynamics > 0) _STATE_VAR_S_(data%id_inben(mag_i)) = zero_
         IF (data%malgs(mag_i)%simIPDynamics > 0) _STATE_VAR_S_(data%id_ipben(mag_i)) = zero_
       ENDIF
     ENDDO
     RETURN
   ENDIF
   IF ( data%simCGM >0 ) THEN
     _DIAG_VAR_S_(data%id_tem_avg) = 4.     !_STATE_VAR_S_(data%id_tem)
     _DIAG_VAR_S_(data%id_tau_avg) = 0.001  !_STATE_VAR_S_(data%id_taub)
     _DIAG_VAR_S_(data%id_par_avg) = 0.1    !_STATE_VAR_S_(data%id_par)
   ENDIF

END SUBROUTINE aed_initialize_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_macroalgae(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Right hand sides of macroalgae biogeochemical model
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mag_i,c
   AED_REAL :: temp, par, Io, salinity, extc, dz, salt
   AED_REAL :: pup, no3up,nh4up, cup, rsiup
   AED_REAL :: phy, tphy, tin, tip, tchla
   AED_REAL :: INi, IPi
   AED_REAL :: primprod(data%num_malgae), exudation(data%num_malgae)
   AED_REAL :: a_nfix(data%num_malgae), respiration(data%num_malgae)
   AED_REAL :: cuptake(data%num_malgae), cexcretion(data%num_malgae), cmortality(data%num_malgae)
   AED_REAL :: nuptake(data%num_malgae,1:4), nexcretion(data%num_malgae), nmortality(data%num_malgae)
   AED_REAL :: puptake(data%num_malgae,1:2), pexcretion(data%num_malgae), pmortality(data%num_malgae)
   AED_REAL :: siuptake(data%num_malgae), siexcretion(data%num_malgae), simortality(data%num_malgae)
   AED_REAL :: fT, fNit, fPho, fSil, fI, fXl, fSal, PNf
   AED_REAL :: upTot,net_cuptake
   AED_REAL :: flux, available
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Retrieve current environmental conditions.
   temp     = _STATE_VAR_(data%id_tem)    ! local temperature
   salinity = _STATE_VAR_(data%id_sal)    ! local salinity
   par      = _STATE_VAR_(data%id_par)    ! local photosynthetically active radn
   Io       = _STATE_VAR_S_(data%id_I_0)  ! surface short wave radiation

   ! Retrieve current (local) state variable values.
   cup = zero_ ; no3up = zero_ ; nh4up = zero_ ; pup = zero_ ; rsiup = zero_
   IF (data%do_Nuptake) THEN
       no3up = _STATE_VAR_(data%id_Nupttarget(1))
       nh4up = _STATE_VAR_(data%id_Nupttarget(2))
   ENDIF
   IF (data%do_Puptake)  pup   = _STATE_VAR_(data%id_Pupttarget(1))
   IF (data%do_Cuptake)  cup   = _STATE_VAR_(data%id_Cupttarget)
   IF (data%do_Siuptake) rsiup = _STATE_VAR_(data%id_Siupttarget)

   ! Initialise cumualtive and working vars
   tphy = zero_ ; tchla= zero_ ; tin  = zero_ ; tip  = zero_
   INi  = zero_ ; IPi  = zero_
   IF (diag_level>0) THEN
     _DIAG_VAR_(data%id_TIN) = zero_
     _DIAG_VAR_(data%id_TIP) = zero_
   ENDIF
   IF (diag_level>1) THEN
     _DIAG_VAR_(data%id_GPP) = zero_
     _DIAG_VAR_(data%id_PUP) = zero_
     _DIAG_VAR_(data%id_NUP) = zero_
     _DIAG_VAR_(data%id_NMP) = zero_
   ENDIF

   !-- Loop through the floating macroalgal groups
   DO mag_i=1,data%num_malgae

      a_nfix(mag_i)      = zero_
      primprod(mag_i)    = zero_
      exudation(mag_i)   = zero_
      respiration(mag_i) = zero_

      cuptake(mag_i)     = zero_
      cexcretion(mag_i)  = zero_
      cmortality(mag_i)  = zero_
      nuptake(mag_i,:)   = zero_
      nexcretion(mag_i)  = zero_
      nmortality(mag_i)  = zero_
      puptake(mag_i,:)   = zero_
      pexcretion(mag_i)  = zero_
      pmortality(mag_i)  = zero_

      ! Retrieve this macroalgae group
      phy = _STATE_VAR_(data%id_p(mag_i))

      fI = zero_; fNit = zero_; fPho = zero_; fSil = one_; fSal = one_; fXl = one_


      IF ( mag_i==data%simCGM ) THEN

        ! Don't grow CGM slough
        primprod(mag_i) = zero_

        ! Respiration and general metabolic loss
        respiration(mag_i) = bio_respiration(data%malgs(mag_i)%R_resp,data%malgs(mag_i)%theta_resp,temp)

        ! photo-exudation
        exudation(mag_i) = zero_

      ELSE   ! Main functions for floating and growing malg

        ! Get the temperature limitation function
        fT = fTemp_function(data%malgs(mag_i)%fT_Method,    &
                            data%malgs(mag_i)%T_max,        &
                            data%malgs(mag_i)%T_std,        &
                            data%malgs(mag_i)%theta_growth, &
                            data%malgs(mag_i)%aTn,          &
                            data%malgs(mag_i)%bTn,          &
                            data%malgs(mag_i)%kTn,temp)

        !fSal = fSal_function(salinity, minS, Smin, Smax, maxS )
        fSal = one_ ! fSal_function(salinity, 25., 30., 45., 80. )
        IF( TRIM(data%malgs(mag_i)%p_name) == 'ulva') THEN
          salt = salinity
          IF( salt<=5. ) THEN
            fSal = zero_
          ELSE IF ( salt>5. .AND. salt<=18.  ) THEN
            fSal = 0. + ( (salt-5.)/(18.-5.) )
          ELSE IF ( salt>18. .AND. salt<=40. ) THEN
            fSal = one_
          ELSE IF ( salt>40. .AND. salt<=65. ) THEN
            fSal = 1. - ( (salt-40.)/(65.-40.) )
          ELSE IF ( salt>65. ) THEN
            fSal = zero_
          ENDIF
        ENDIF

        ! Get the light and nutrient limitation.
        ! NITROGEN.
        fNit = 0.0
        IF(data%malgs(mag_i)%simINDynamics /= 0) THEN
         ! IN variable available
         INi = _STATE_VAR_(data%id_in(mag_i))
        ELSE
         ! Assumed constant IN:
         INi = phy*data%malgs(mag_i)%X_ncon
        END IF

        ! Estimate fN limitation from IN or ext N value
        IF(data%malgs(mag_i)%simINDynamics > 1) THEN
         IF (phy > data%malgs(mag_i)%p0) THEN
            fNit = INi / phy
            fNit = phyto_fN(data%malgs,mag_i,IN=fNit)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%malgs(mag_i)%p0) THEN
            fNit = phyto_fN(data%malgs,mag_i,din=no3up+nh4up)
         ENDIF
        ELSE
         fNit = phyto_fN(data%malgs,mag_i,din=no3up+nh4up)
        ENDIF
        IF (data%malgs(mag_i)%simNFixation /= 0) THEN
         ! Nitrogen fixer: apply no N limitation. N Fixation ability
         ! depends on DIN concentration
         a_nfix = (one_ - fNit)
         fNit = one_
        ENDIF

        ! PHOSPHOROUS.
        fPho = zero_
        IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
         ! IP variable available
         IPi = _STATE_VAR_(data%id_ip(mag_i))
        ELSE
         ! Assumed constant IP:
         IPi = phy*data%malgs(mag_i)%X_pcon
        END IF

        ! Estimate fP limitation from IP or ext P value
        IF (data%malgs(mag_i)%simIPDynamics > 1) THEN
         IF (phy > data%malgs(mag_i)%p0) THEN
            fPho = IPi / phy
            fPho = phyto_fP(data%malgs,mag_i,IP=fPho)
         ENDIF
         IF (phy > zero_ .AND. phy <= data%malgs(mag_i)%p0) THEN
            fPho = phyto_fP(data%malgs,mag_i,frp=pup)
         ENDIF
        ELSE
          fPho = phyto_fP(data%malgs,mag_i,frp=pup)
        ENDIF

        ! LIGHT
        extc = _STATE_VAR_(data%id_extc)
        dz = _STATE_VAR_(data%id_dz)        ! water column layer depth
        fI = photosynthesis_irradiance(data%malgs(mag_i)%lightModel, &
               data%malgs(mag_i)%I_K, data%malgs(mag_i)%I_S, par, extc, Io, dz)

        ! METAL AND TOXIC EFFECTS
        fXl = 1.0

        ! Primary production rate
        primprod(mag_i) = data%malgs(mag_i)%R_growth * fT  * fXl * fSal &
                        * findMin(fI,fNit,fPho,fSil)

        ! Adjust primary production rate for nitrogen fixers
        IF (data%malgs(mag_i)%simNFixation /= 0) THEN
         ! Nitrogen fixing species, and the growth rate to  must be reduced
         ! to compensate for the increased metabolic cost of this process
         primprod(mag_i) = primprod(mag_i) * (data%malgs(mag_i)%k_nfix + &
                           (1.0-a_nfix(mag_i))*(1.0-data%malgs(mag_i)%k_nfix))
        ENDIF

        ! Respiration and general metabolic loss
        respiration(mag_i) = bio_respiration(data%malgs(mag_i)%R_resp,data%malgs(mag_i)%theta_resp,temp)

        ! Salinity stress effect on respiration
        !fSal =  phyto_salinity(data%malgs,mag_i,salinity)
        !respiration(mag_i) = respiration(mag_i) * fSal

        ! photo-exudation
        exudation(mag_i) = primprod(mag_i)*data%malgs(mag_i)%f_pr

      ENDIF


      ! Limit respiration if at the min biomass to prevent
      ! leak in the C mass balance
      IF (phy <= data%malgs(mag_i)%p0) THEN
         respiration(mag_i) = zero_
         exudation(mag_i) = zero_
      ENDIF

      ! Carbon uptake and excretion
      cuptake(mag_i)    = -primprod(mag_i) * phy
      cexcretion(mag_i) = (data%malgs(mag_i)%k_fdom*(1.0-data%malgs(mag_i)%k_fres)*respiration(mag_i)+exudation(mag_i)) * phy
      cmortality(mag_i) = ((1.0-data%malgs(mag_i)%k_fdom)*(1.0-data%malgs(mag_i)%k_fres)*respiration(mag_i)) * phy

      ! Nitrogen uptake and excretion
      CALL phyto_internal_nitrogen(data%malgs,mag_i,data%do_N2uptake,phy,INi,primprod(mag_i),&
                             fT,no3up,nh4up,a_nfix(mag_i),respiration(mag_i),exudation(mag_i),PNf,&
                             nuptake(mag_i,:),nexcretion(mag_i),nmortality(mag_i))

      ! Phosphorus uptake and excretion
      CALL phyto_internal_phosphorus(data%malgs,mag_i,data%npup,phy,IPi,primprod(mag_i),&
                             fT,pup,respiration(mag_i),exudation(mag_i),&
                             puptake(mag_i,:),pexcretion(mag_i),pmortality(mag_i))


      ! Diagnostic info
      IF (diag_level>9) THEN
         _DIAG_VAR_(data%id_c2p(mag_i))  =  phy/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)
         _DIAG_VAR_(data%id_c2p(mag_i))  =  phy/MAX(INi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_nmin)
         _DIAG_VAR_(data%id_n2p(mag_i))  =  INi/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)
         _DIAG_VAR_(data%id_fT(mag_i))   =  fT
         _DIAG_VAR_(data%id_fI(mag_i))   =  fI
         _DIAG_VAR_(data%id_fNit(mag_i)) =  fNit
         _DIAG_VAR_(data%id_fPho(mag_i)) =  fPho
         _DIAG_VAR_(data%id_fSal(mag_i)) =  fSal
      ENDIF
   ENDDO


   !-----------------------------------------------------------------
   !-- Check uptake values for availability to prevent -ve numbers

   ! pup   - p available for uptake
   ! no3up - no3 available for uptake
   ! nh4up - nh4 available for uptake
   ! cup   - c available for uptake
   ! rsiup - Si available for uptake

   IF (data%do_Puptake) THEN
      upTot = sum(puptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= pup ) THEN
         DO mag_i=1,data%num_malgae
            puptake(mag_i,1) = (pup*0.99/dtlim) * (puptake(mag_i,1)/upTot)
         ENDDO
      ENDIF
   ENDIF

   IF (data%do_Nuptake) THEN
      upTot = sum(nuptake(:,1))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= no3up ) THEN
         DO mag_i=1,data%num_malgae
            nuptake(mag_i,1) = (no3up*0.99/dtlim) * (nuptake(mag_i,1)/upTot)
         ENDDO
      ENDIF

      upTot = sum(nuptake(:,2))*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= nh4up ) THEN
         DO mag_i=1,data%num_malgae
            nuptake(mag_i,2) = (nh4up*0.99/dtlim) * (nuptake(mag_i,2)/upTot)
         ENDDO
      ENDIF
   ENDIF
   IF (data%do_Cuptake) THEN
      upTot = sum(cuptake)*dtlim
      IF ( abs(upTot) > 1.e-10 .and. upTot >= cup ) THEN
         DO mag_i=1,data%num_malgae
            cuptake(mag_i) = (cup*0.99/dtlim) * (cuptake(mag_i)/upTot)
         ENDDO
      ENDIF
   ENDIF

   !-----------------------------------------------------------------
   ! SET TEMPORAL DERIVATIVES FOR ODE SOLVER
   net_cuptake = zero_
   DO mag_i=1,data%num_malgae

      !# macroalgae PRODUCTION & RESPIRATION
      phy = _STATE_VAR_(data%id_p(mag_i))
      flux = (primprod(mag_i) - respiration(mag_i) - exudation(mag_i)) * phy
      available = MAX(zero_, phy - data%malgs(mag_i)%p0)
      IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
      _FLUX_VAR_(data%id_p(mag_i)) = _FLUX_VAR_(data%id_p(mag_i)) + ( flux)

      !# macroalgae INTERNAL NITROGEN
      IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
         INi = _STATE_VAR_(data%id_in(mag_i))
         flux = (-sum(nuptake(mag_i,:)) - nexcretion(mag_i) - nmortality(mag_i) )
         available = MAX(zero_, INi - data%malgs(mag_i)%X_nmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_in(mag_i)) = _FLUX_VAR_(data%id_in(mag_i)) + ( flux)
      ENDIF

      !# macroalgae INTERNAL PHOSPHORUS
      IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
         IPi = _STATE_VAR_(data%id_ip(mag_i))
         flux = (-sum(puptake(mag_i,:)) - pexcretion(mag_i) - pmortality(mag_i) )
         available = MAX(zero_, IPi - data%malgs(mag_i)%X_pmin*phy)
         IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
         _FLUX_VAR_(data%id_ip(mag_i)) = _FLUX_VAR_(data%id_ip(mag_i)) + ( flux)
      ENDIF

      !# macroalgae CELL DENSITY
      IF ( data%id_rho(mag_i)>0 ) THEN
         ! density increases during carbohydrate creation (daytime)
         flux = zero_
         IF( par>zero_ ) THEN
           flux = data%malgs(mag_i)%c1 * &
              (one_ - EXP(-par/data%malgs(mag_i)%I_K) ) - data%malgs(mag_i)%c3
         ELSE
           ! darkness
           flux = -data%malgs(mag_i)%c3
         ENDIF
         _FLUX_VAR_(data%id_rho(mag_i)) = _FLUX_VAR_(data%id_rho(mag_i)) + flux
         ! check maximum/minimum density are not exceeded
         IF( _STATE_VAR_(data%id_rho(mag_i))>data%max_rho ) THEN
            _FLUX_VAR_(data%id_rho(mag_i)) =zero_
            _STATE_VAR_(data%id_rho(mag_i))=data%max_rho
         ENDIF
         IF( _STATE_VAR_(data%id_rho(mag_i))<data%min_rho ) THEN
             _FLUX_VAR_(data%id_rho(mag_i)) =zero_
             _STATE_VAR_(data%id_rho(mag_i))=data%min_rho
         ENDIF
      ENDIF

      ! BIOGEOCHEMICAL FEEDBACKS
      IF (data%simMalgFeedback) THEN
        ! Now manage uptake of nutrients, CO2 and DO - these cumulative fluxes already limited above loop
        IF (data%do_Puptake) THEN
         DO c = 1,data%npup
            _FLUX_VAR_(data%id_Pupttarget(c)) =    &
                       _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(mag_i,c))
         ENDDO
        ENDIF
        IF (data%do_Nuptake) THEN
         DO c = 1,data%nnup
            _FLUX_VAR_(data%id_Nupttarget(c)) =    &
                       _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(mag_i,c))
         ENDDO
        ENDIF
        IF (data%do_Cuptake) THEN
         _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) +    &
                         (  cuptake(mag_i) - respiration(mag_i)*data%malgs(mag_i)%k_fres*phy )
         net_cuptake = net_cuptake + (cuptake(mag_i) - respiration(mag_i)*data%malgs(mag_i)%k_fres*phy)
        ENDIF
        IF (data%do_DOuptake) THEN
         _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) +    &
                         ( -cuptake(mag_i) + respiration(mag_i)*data%malgs(mag_i)%k_fres*phy )
        ENDIF
        IF (data%do_Siuptake) THEN
         _FLUX_VAR_(data%id_Siupttarget) = _FLUX_VAR_(data%id_Siupttarget) + ( siuptake(mag_i))
        ENDIF
        ! Now manage mortality contributions to POM
        IF (data%do_Pmort) THEN
         _FLUX_VAR_(data%id_Pmor) = _FLUX_VAR_(data%id_Pmor) + (pmortality(mag_i))
        ENDIF
        IF (data%do_Nmort) THEN
         _FLUX_VAR_(data%id_Nmor) = _FLUX_VAR_(data%id_Nmor) + (nmortality(mag_i))
        ENDIF
        IF (data%do_Cmort) THEN
         _FLUX_VAR_(data%id_Cmor) = _FLUX_VAR_(data%id_Cmor) + (cmortality(mag_i))
        ENDIF
        ! Now manage excretion/exudation contributions to DOM
        IF (data%do_Pexc) THEN
         _FLUX_VAR_(data%id_Pexr) = _FLUX_VAR_(data%id_Pexr) + (pexcretion(mag_i))
        ENDIF
        IF (data%do_Nexc) THEN
         _FLUX_VAR_(data%id_Nexr) = _FLUX_VAR_(data%id_Nexr) + (nexcretion(mag_i))
        ENDIF
        IF (data%do_Cexc) THEN
         _FLUX_VAR_(data%id_Cexr) = _FLUX_VAR_(data%id_Cexr) + (cexcretion(mag_i))
        ENDIF
      ENDIF

      !-----------------------------------------------------------------
      ! UPDATE DIAGNOSTIC VARIABLES

      ! total macroalgae carbon
      tphy = tphy + phy
      ! total internal nutrients
      tin = tin + INi
      tip = tip + IPi

      !This overwrites benthic addition (which comes first) - need to sort out order.
      !IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC) = data%malgs(mag_i)%KePHY*phy

   ENDDO

   ! Set diagnostic arrays for combined assemblage properties
   IF (diag_level>0) THEN
     _DIAG_VAR_(data%id_TIN)   =  tin
     _DIAG_VAR_(data%id_TIP)   =  tip
   ENDIF

   ! Summary of productivity
   IF (diag_level>1) THEN
     _DIAG_VAR_(data%id_GPP) = _DIAG_VAR_(data%id_GPP) + sum(cuptake)*secs_per_day
     _DIAG_VAR_(data%id_NUP) = _DIAG_VAR_(data%id_NUP) + sum(nuptake)*secs_per_day
     _DIAG_VAR_(data%id_PUP) = _DIAG_VAR_(data%id_PUP) + sum(puptake)*secs_per_day
     _DIAG_VAR_(data%id_NMP) = _DIAG_VAR_(data%id_NMP) +(sum(primprod)*tphy - sum(respiration)*tphy &
                                                       - sum(exudation)*tphy)*secs_per_day
   ENDIF

END SUBROUTINE aed_calculate_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_calculate_benthic_macroalgae(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate benthic macroalgae (growth, respiration and sloughing)
! Benthic variables are in units per surface area; fluxes are per sec
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mag_i,c
   AED_REAL :: malg, INi, IPi                                        ! State
   AED_REAL :: temp,extc,par,dz,Io,fI,bottom_stress,depth,light,matz
   AED_REAL :: tphy, tin, tip, tchla
   AED_REAL :: salinity, salt

   AED_REAL :: malg_flux,malg_prod,malg_resp                         ! Fluxes
   AED_REAL :: no3up, nh4up, pup, cup, rsiup
   AED_REAL :: primprod(data%num_malgae), exudation(data%num_malgae)
   AED_REAL :: a_nfix(data%num_malgae), respiration(data%num_malgae)
   AED_REAL :: cuptake(data%num_malgae), cexcretion(data%num_malgae), cmortality(data%num_malgae)
   AED_REAL :: nuptake(data%num_malgae,1:4), nexcretion(data%num_malgae), nmortality(data%num_malgae)
   AED_REAL :: puptake(data%num_malgae,1:2), pexcretion(data%num_malgae), pmortality(data%num_malgae)
   AED_REAL :: siuptake(data%num_malgae), siexcretion(data%num_malgae), simortality(data%num_malgae)
   AED_REAL :: fT, fNit, fPho, fSil, fXl, fSal, PNf
   AED_REAL :: upTot, net_cuptake, available, flux
   AED_REAL :: slough_frac = zero_
!
!-------------------------------------------------------------------------------
!BEGIN

  !-- Benthic light fraction and extinction, for diagnostics
  extc  = _STATE_VAR_(data%id_extc)      ! extinction coefficient of bottom cell
  depth = _STATE_VAR_(data%id_dz)        ! water column depth (cell depth if 2D)
  matz  = _STATE_VAR_S_(data%id_sedzone) ! sediment zone / material type

  !-- Initialise the accumulated biomass variables
  IF (diag_level>0) THEN
    _DIAG_VAR_S_(data%id_mag_ben) = zero_
    _DIAG_VAR_S_(data%id_min_ben) = zero_
    _DIAG_VAR_S_(data%id_mip_ben) = zero_
  ENDIF
  IF (diag_level>1) THEN
    _DIAG_VAR_S_(data%id_gpp_ben) = zero_
    _DIAG_VAR_S_(data%id_nmp_ben) = zero_
    _DIAG_VAR_S_(data%id_nup_ben) = zero_
    _DIAG_VAR_S_(data%id_pup_ben) = zero_
    _DIAG_VAR_S_(data%id_slg_ben) = zero_
  ENDIF

  !-- Loop through selected groups, determining if benthic/attached is active
  DO mag_i=1,data%num_malgae

     fI = zero_; fNit = zero_; fPho = zero_; fSil = one_; fSal = one_; fXl = one_

     !-- Process each benthic / attached macroalgae
     IF ( data%malgs(mag_i)%settling == _MOB_ATTACHED_ .AND. &
                                      in_zone_set(matz,data%active_zones) ) THEN

       ! Get local conditions
       salinity = _STATE_VAR_(data%id_sal)            ! local salinity
       temp = _STATE_VAR_(data%id_tem)                ! local temperature
       dz   = _STATE_VAR_(data%id_dz)                 ! cell/layer thickness
       par  = MAX(_STATE_VAR_(data%id_par),zero_)     ! photosynth. active radn
       Io   = MAX(_STATE_VAR_S_(data%id_I_0),zero_)   ! surface shortwave radn
       bottom_stress = _STATE_VAR_S_(data%id_taub)    ! bottom stress

       ! Retrieve current (local) conditions
       pup = 0. ; no3up = 0. ; nh4up = 0. ; cup = 0.
       IF (data%do_Cuptake)  cup   = _STATE_VAR_(data%id_Cupttarget)
       IF (data%do_Puptake)  pup   = _STATE_VAR_(data%id_Pupttarget(1))
       IF (data%do_Nuptake)  no3up = _STATE_VAR_(data%id_Nupttarget(1))
       IF (data%do_Nuptake)  nh4up = _STATE_VAR_(data%id_Nupttarget(2))

       ! Initialise cumualtive rates
       tphy = 0. ; tchla = 0. ; tin = 0. ; tip = 0.  !MOVE?
       INi = 0. ; IPi = 0.

       primprod(mag_i)    = zero_
       exudation(mag_i)   = zero_
       a_nfix(mag_i)      = zero_
       respiration(mag_i) = zero_
       cuptake(mag_i) = zero_ ; cexcretion(mag_i) = zero_ ; cmortality(mag_i) = zero_
       nuptake(mag_i,:) = zero_ ; nexcretion(mag_i) = zero_ ; nmortality(mag_i) = zero_
       puptake(mag_i,:) = zero_ ; pexcretion(mag_i) = zero_ ; pmortality(mag_i) = zero_

       ! Retrieve this macroalgae group
       malg = _STATE_VAR_S_(data%id_pben(mag_i)) ! local malg density



       !-- Use the generic approach or CGM approach
       IF( .NOT. mag_i == data%simCGM ) THEN

         ! Get the temperature limitation function
         fT = fTemp_function(data%malgs(mag_i)%fT_Method,    &
                             data%malgs(mag_i)%T_max,        &
                             data%malgs(mag_i)%T_std,        &
                             data%malgs(mag_i)%theta_growth, &
                             data%malgs(mag_i)%aTn,          &
                             data%malgs(mag_i)%bTn,          &
                             data%malgs(mag_i)%kTn,temp)

         !fSal = fSal_function(salinity, minS, Smin, Smax, maxS )
         !fSal = one_  !fSal_function(salinity, 5., 18., 35., 45. )

         salt = salinity       !COORONG (Ulva) Make me flexible
         IF( salt<=5. ) THEN
           fSal = zero_
         ELSE IF ( salt>5. .AND. salt<=18.  ) THEN
           fSal = 0. + ( (salt-5.)/(18.-5.) )
         ELSE IF ( salt>18. .AND. salt<=40. ) THEN
           fSal = one_
         ELSE IF ( salt>40. .AND. salt<=85. ) THEN
           fSal = 1. - ( (salt-40.)/(85.-40.) )
         ELSE IF ( salt>85. ) THEN
           fSal = zero_
         ENDIF

         ! Get the light and nutrient limitation.
         ! NITROGEN.
         IF(data%malgs(mag_i)%simINDynamics /= 0) THEN
           ! IN variable available
           INi = _STATE_VAR_S_(data%id_inben(mag_i))
         ELSE
           ! Assumed constant IN:
           INi = malg*data%malgs(mag_i)%X_ncon
         END IF

         ! Estimate fN limitation from IN or ext N value
         IF(data%malgs(mag_i)%simINDynamics > 1) THEN
           IF (malg > data%malgs(mag_i)%p0) THEN
             fNit = INi / malg
             fNit = phyto_fN(data%malgs,mag_i,IN=fNit)
           ENDIF
           IF (malg > zero_ .AND. malg <= data%malgs(mag_i)%p0) THEN
             fNit = phyto_fN(data%malgs,mag_i,din=no3up+nh4up)
           ENDIF
         ELSE
           fNit = phyto_fN(data%malgs,mag_i,din=no3up+nh4up)
         ENDIF
         IF (data%malgs(mag_i)%simNFixation /= 0) THEN
           ! Nitrogen fixer: apply no N limitation. N Fixation ability
           ! depends on DIN concentration
           a_nfix = (one_ - fNit)
           fNit = one_
         ENDIF

         ! PHOSPHOROUS.
         IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
           ! IP variable available
           IPi = _STATE_VAR_S_(data%id_ipben(mag_i))
         ELSE
           ! Assumed constant IP:
           IPi = malg*data%malgs(mag_i)%X_pcon
         END IF

         ! Estimate fP limitation from IP or ext P value
         IF (data%malgs(mag_i)%simIPDynamics > 1) THEN
           IF (malg > data%malgs(mag_i)%p0) THEN
             fPho = IPi / malg
             fPho = phyto_fP(data%malgs,mag_i,IP=fPho)
           ENDIF
           IF (malg > zero_ .AND. malg <= data%malgs(mag_i)%p0) THEN
             fPho = phyto_fP(data%malgs,mag_i,frp=pup)
           ENDIF
         ELSE
           fPho = phyto_fP(data%malgs,mag_i,frp=pup)
         ENDIF
         fPho = 1.0

         IF(bottom_stress>0.09) fXl =0.0

         ! Compute P-I
         fI = photosynthesis_irradiance(0, &  ! 0 is vertical integral
              data%malgs(mag_i)%I_K, data%malgs(mag_i)%I_S, par, extc, Io, dz)

         IF (diag_level>9) THEN
           _DIAG_VAR_S_(data%id_fT_ben(mag_i))   =  fT
           _DIAG_VAR_S_(data%id_fI_ben(mag_i))   =  fI
           _DIAG_VAR_S_(data%id_fNit_ben(mag_i)) =  fNit
           _DIAG_VAR_S_(data%id_fPho_ben(mag_i)) =  fPho
           _DIAG_VAR_S_(data%id_fSal_ben(mag_i)) =  fSal
         ENDIF

         ! Primary production rate
         primprod(mag_i) = data%malgs(mag_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl * fSal

         ! Respiration and general metabolic loss
         respiration(mag_i) = bio_respiration(data%malgs(mag_i)%R_resp,data%malgs(mag_i)%theta_resp,temp)

         ! Photo-exudation
         exudation(mag_i) = primprod(mag_i)*data%malgs(mag_i)%f_pr

         ! Limit respiration if at the min biomass to prevent leak in the C mass balance
         IF (malg <= data%malgs(mag_i)%p0) THEN
           respiration(mag_i) = zero_
           exudation(mag_i) = zero_
         ENDIF

         !# Carbon uptake and excretion
         cuptake(mag_i)    = -primprod(mag_i) * malg
         cexcretion(mag_i) = (data%malgs(mag_i)%k_fdom*(1.0-data%malgs(mag_i)%k_fres)*respiration(mag_i)+exudation(mag_i)) * malg
         cmortality(mag_i) = ((1.0-data%malgs(mag_i)%k_fdom)*(1.0-data%malgs(mag_i)%k_fres)*respiration(mag_i)) * malg

         !# Nitrogen uptake and excretion
         CALL phyto_internal_nitrogen(data%malgs,mag_i,data%do_N2uptake,malg,INi,&
                primprod(mag_i),fT,no3up,nh4up,a_nfix(mag_i),respiration(mag_i),&
                exudation(mag_i),PNf,nuptake(mag_i,:),nexcretion(mag_i),nmortality(mag_i))

         !# Phosphorus uptake and excretion
         CALL phyto_internal_phosphorus(data%malgs,mag_i,data%npup,malg,IPi,primprod(mag_i),&
                fT,pup,respiration(mag_i),exudation(mag_i),&
                puptake(mag_i,:),pexcretion(mag_i),pmortality(mag_i))

         malg_flux = ( primprod(mag_i)-respiration(mag_i)-exudation(mag_i) ) * malg

         !# Bottom (nonCGM) macroalgal canopy properties
         macroHgt = 0.1
         macroExt = data%malgs(mag_i)%KePHY*malg

         IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC)  = macroExt
         IF(diag_level>9) _DIAG_VAR_S_(data%id_dHGT) = macroHgt
         IF(diag_level>9) _DIAG_VAR_S_(data%id_dPAR) = par


       ELSE ! group check

         !-- User has selected the special CGM approach (Cladophora Growth Model)
         CALL cgm_calculate_benthic_cladophora(data,column,layer_idx,mag_i,&
             primprod(mag_i),respiration(mag_i),puptake(mag_i,1),nuptake(mag_i,1),INi,IPi)


         exudation(mag_i)  = MAX( zero_,primprod(mag_i)*data%malgs(mag_i)%f_pr )      ! /second

         malg_flux = ( primprod(mag_i)-exudation(mag_i)-respiration(mag_i) ) * malg   ! mmolC/m2/sec

         cuptake(mag_i)    = primprod(mag_i) * malg                                   ! mmolC/m2/sec
         cexcretion(mag_i) = zero_
         cmortality(mag_i) = zero_
         nuptake(mag_i,1)  = -( primprod(mag_i)-exudation(mag_i)-respiration(mag_i) ) * INi  ! Force CGM to adopt a C:N
         nexcretion(mag_i) = zero_
         nmortality(mag_i) = zero_
         puptake(mag_i,1)  = -puptake(mag_i,1) * IPi
         pexcretion(mag_i) = zero_
         pmortality(mag_i) = zero_

         IF (diag_level>9) THEN
          !_DIAG_VAR_S_(data%id_fSal_ben(mag_i)) = one_
         ENDIF

       ENDIF

       !# Record standing biomass of all MAG into diagnostics (mmol/m2)
       IF (diag_level>0) _DIAG_VAR_S_(data%id_mag_ben) = _DIAG_VAR_S_(data%id_mag_ben) + malg
       IF (diag_level>0) _DIAG_VAR_S_(data%id_min_ben) = _DIAG_VAR_S_(data%id_min_ben) + INi
       IF (diag_level>0) _DIAG_VAR_S_(data%id_mip_ben) = _DIAG_VAR_S_(data%id_mip_ben) + IPi
       IF (diag_level>9) _DIAG_VAR_S_(data%id_c2p_ben(mag_i)) = malg/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)
       IF (diag_level>9) _DIAG_VAR_S_(data%id_c2n_ben(mag_i)) = malg/MAX(INi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_nmin)
       IF (diag_level>9) _DIAG_VAR_S_(data%id_n2p_ben(mag_i)) = INi/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)

       !# Record gross productivity and repsiration rates (/day)
       IF (diag_level>1) _DIAG_VAR_S_(data%id_rsp_ben) = -respiration(mag_i) *secs_per_day
       IF (diag_level>1) _DIAG_VAR_S_(data%id_gpp_ben) = &
                                        _DIAG_VAR_S_(data%id_gpp_ben) + cuptake(mag_i) *secs_per_day / MAX(malg,10.)

       !# Record net benthic macroalgae production (mmol/m2/day)
       IF (diag_level>1) _DIAG_VAR_S_(data%id_nmp_ben) = _DIAG_VAR_S_(data%id_nmp_ben) + malg_flux *secs_per_day
       IF (diag_level>1) _DIAG_VAR_S_(data%id_pup_ben) = _DIAG_VAR_S_(data%id_pup_ben) + puptake(mag_i,1) *secs_per_day
       IF (diag_level>1) _DIAG_VAR_S_(data%id_nup_ben) = _DIAG_VAR_S_(data%id_nup_ben) + nuptake(mag_i,1) *secs_per_day

       !# Update the attached (carbon) biomass based on the net growth rate
       _FLUX_VAR_B_(data%id_pben(mag_i)) = _FLUX_VAR_B_(data%id_pben(mag_i)) + malg_flux

       !# Update the INTERNAL NITROGEN biomass of the attached macroalgae
       IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
          INi = _STATE_VAR_S_(data%id_inben(mag_i))
          flux = (sum(-nuptake(mag_i,:)) - nexcretion(mag_i) - nmortality(mag_i) )
          available = MAX(zero_, INi - data%malgs(mag_i)%X_nmin*malg)
          IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
          _FLUX_VAR_B_(data%id_inben(mag_i)) = _FLUX_VAR_B_(data%id_inben(mag_i)) + flux
       ENDIF

       !# Update the INTERNAL PHOSPHORUS biomass of the attached macroalgae
       IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
          !IPi = _STATE_VAR_S_(data%id_ipben(mag_i))
          flux = (sum(-puptake(mag_i,:)) - pexcretion(mag_i) - pmortality(mag_i))
          available = MAX( zero_, IPi-data%malgs(mag_i)%X_pmin*malg )
          IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
          _FLUX_VAR_B_(data%id_ipben(mag_i)) = _FLUX_VAR_B_(data%id_ipben(mag_i)) + flux
       ENDIF

       !# Update the water column vars (O2, CO2, Nuts, OM) based on these fluxes
       IF (data%simMalgFeedback) THEN
         IF (data%do_DOuptake) THEN
           _FLUX_VAR_(data%id_DOupttarget) = _FLUX_VAR_(data%id_DOupttarget) + malg_flux
         ENDIF
         IF (data%do_Cuptake) THEN
           _FLUX_VAR_(data%id_Cupttarget) = _FLUX_VAR_(data%id_Cupttarget) - malg_flux
         ENDIF
         IF (data%do_Puptake) THEN
          DO c = 1,data%npup
             _FLUX_VAR_(data%id_Pupttarget(c)) =         &
                      _FLUX_VAR_(data%id_Pupttarget(c)) + ( puptake(mag_i,c))
          ENDDO
         ENDIF
         IF (data%do_Nuptake) THEN
          DO c = 1,data%nnup
             _FLUX_VAR_(data%id_Nupttarget(c)) =         &
                      _FLUX_VAR_(data%id_Nupttarget(c)) + ( nuptake(mag_i,c))
          ENDDO
         ENDIF

         IF(data%do_Cexc)  _FLUX_VAR_(data%id_Cexr) = &
                                    _FLUX_VAR_(data%id_Cexr) - cexcretion(mag_i) - exudation(mag_i)*malg
         IF(data%do_Nexc)  _FLUX_VAR_(data%id_Nexr) = _FLUX_VAR_(data%id_Nexr) - nexcretion(mag_i)
         IF(data%do_Pexc)  _FLUX_VAR_(data%id_Pexr) = _FLUX_VAR_(data%id_Pexr) - pexcretion(mag_i)
         IF(data%do_Cmort) _FLUX_VAR_(data%id_Cmor) = _FLUX_VAR_(data%id_Cmor) - cmortality(mag_i)
         IF(data%do_Nmort) _FLUX_VAR_(data%id_Nmor) = _FLUX_VAR_(data%id_Nmor) - nmortality(mag_i)
         IF(data%do_Pmort) _FLUX_VAR_(data%id_Pmor) = _FLUX_VAR_(data%id_Pmor) - pmortality(mag_i)
       ENDIF

       !# Redistribute biomass into the water column if sloughing occurs.
       IF( data%simSloughing == 1) THEN
          ! The Coorong Ulva approach
          IF( bottom_stress>data%malgs(mag_i)%tau_0*2. ) THEN
            slough_frac = 0.3
          ELSEIF( bottom_stress>data%malgs(mag_i)%tau_0 .AND. salinity>60.) THEN
            slough_frac = 0.67
          ENDIF
       ELSEIF( data%simSloughing == 2) THEN
          ! The Erie CGM approach
          CALL cgm_slough_cladophora(data,column,layer_idx,mag_i,data%malgs(mag_i)%resuspension,slough_frac)
       ELSE
          ! No sloughing
          slough_frac = zero_
       ENDIF

       _FLUX_VAR_B_(data%id_pben(mag_i)) = _FLUX_VAR_B_(data%id_pben(mag_i)) - slough_frac*malg
       _FLUX_VAR_(data%id_p(mag_i)) = _FLUX_VAR_(data%id_p(mag_i)) + (slough_frac*malg)/depth
       IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
          _FLUX_VAR_B_(data%id_ipben(mag_i)) = _FLUX_VAR_B_(data%id_ipben(mag_i)) - slough_frac*IPi
          _FLUX_VAR_(data%id_ip(mag_i)) = _FLUX_VAR_(data%id_ip(mag_i)) + (slough_frac*IPi)/depth
       ENDIF
       IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
          _FLUX_VAR_B_(data%id_inben(mag_i)) = _FLUX_VAR_B_(data%id_inben(mag_i)) - slough_frac*INi
          _FLUX_VAR_(data%id_in(mag_i)) = _FLUX_VAR_(data%id_in(mag_i)) + (slough_frac*INi)/depth
       ENDIF
       IF (diag_level>1) THEN
         _DIAG_VAR_S_(data%id_slg_ben) = _DIAG_VAR_S_(data%id_slg_ben) - slough_frac*malg*secs_per_day
       ENDIF


       !# TMALG - total macroalgal biomass as gDW/m2
       !  (note this is 1 timestep behind) and only for last group
       IF (diag_level>0) THEN
         _DIAG_VAR_(data%id_TMALG) = (12.*1e-3/data%malgs(mag_i)%Xcc) *        &
                                     (malg + _STATE_VAR_(data%id_p(mag_i))*depth)
       ENDIF

     ENDIF ! end benthic group/zone check

     !# Compute the (growth) HSI for chosen group
     IF( data%simMalgHSI == mag_i) THEN
       _DIAG_VAR_S_(data%id_mhsi) = min( fSal, fI, fNit, fPho ) * fT
     ENDIF

   ENDDO

END SUBROUTINE aed_calculate_benthic_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_mobility_macroalgae(data,column,layer_idx,mobility)
!-------------------------------------------------------------------------------
! Get the vertical movement values for floating macroalgal material
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: mobility(:)
!
!LOCALS
   AED_REAL :: temp, par, rho_p, Io
   AED_REAL :: vvel
   AED_REAL :: pw, pw20, mu, mu20
   AED_REAL :: IN, IC, Q, Qmax
   INTEGER  :: mag_i
!
!-------------------------------------------------------------------------------
!BEGIN

   DO mag_i=1,data%num_malgae
      SELECT CASE (data%malgs(mag_i)%settling)

         CASE ( _MOB_OFF_ )
            ! disable settling by setting vertical velocity to 0
            vvel = zero_

         CASE ( _MOB_CONST_ )
            ! constant settling velocity using user provided value
            vvel = data%malgs(mag_i)%w_p

         CASE ( _MOB_TEMP_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_dens)
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%malgs(mag_i)%w_p*mu20*pw / ( mu*pw20 )

         CASE ( _MOB_STOKES_ )
            ! settling velocity based on Stokes Law calculation and cell density
            pw = _STATE_VAR_(data%id_dens)             ! water density
            temp = _STATE_VAR_(data%id_tem)            ! water temperature
            mu = water_viscosity(temp)                 ! water dynamic viscosity
            IF( data%id_rho(mag_i)>0 ) THEN
              rho_p = _STATE_VAR_(data%id_rho(mag_i))  ! cell density
            ELSE
              rho_p = data%malgs(mag_i)%rho_phy
            ENDIF
            vvel = -9.807*(data%malgs(mag_i)%d_phy**2.)*( rho_p-pw ) / ( 18.*mu )

          CASE ( _MOB_ATTACHED_ )
            ! constant settling velocity @20C corrected for density changes
            pw = _STATE_VAR_(data%id_dens)
            temp = _STATE_VAR_(data%id_tem)
            mu = water_viscosity(temp)
            mu20 = 0.001002  ! N s/m2
            pw20 = 998.2000  ! kg/m3 (assuming freshwater)
            vvel = data%malgs(mag_i)%w_p*mu20*pw / ( mu*pw20 )

         CASE DEFAULT
            ! unknown settling/migration option selection
            vvel =  zero_

      END SELECT
      ! set global mobility array
      mobility(data%id_p(mag_i)) = vvel
      IF(diag_level>9 .AND. data%id_vvel(mag_i)>0) _DIAG_VAR_(data%id_vvel(mag_i)) = vvel
    ENDDO
END SUBROUTINE aed_mobility_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_macroalgae(data,column,layer_idx,extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to macroalgal biomass
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: malg,dz
   INTEGER  :: mag_i
!
!-------------------------------------------------------------------------------
!BEGIN

   DO mag_i=1,data%num_malgae
     ! Retrieve current group floating macroalgal conc (mmol C/m3)
     malg = _STATE_VAR_(data%id_p(mag_i))

     ! Self-shading with explicit contribution from macroalgae concentration
     extinction = extinction + (data%malgs(mag_i)%KePHY*malg)


     ! BENTHIC LIGHT EXTINCTION EFFECT TO BE ADDED (Code below will add to all cells)
!      ! Process benthic / attached macroalgae
!      IF ( data%malgs(mag_i)%settling == _MOB_ATTACHED_ .AND. bottom_cell) THEN
!        dz   = MAX(_STATE_VAR_(data%id_dz),0.05)     ! cell depth
!        malg = _STATE_VAR_S_(data%id_pben(mag_i))/dz ! attached macroalgae
!        extinction = extinction + (data%malgs(mag_i)%KePHY*malg)
!      ENDIF
   ENDDO

END SUBROUTINE aed_light_extinction_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE aed_bio_drag_macroalgae(data,column,layer_idx,drag)
!-------------------------------------------------------------------------------
! Get the effect of macroalgal biomass on benthic drag
!
!  NOTE - THIS IS ADDING THE MALG DRAG EFFECT TO BOTTOM LAYER ONLY
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
!
!LOCALS
   AED_REAL :: dz,malg
   INTEGER  :: mag_i
   AED_REAL, PARAMETER :: KCD_wc  = 0.0  !WARNING
   AED_REAL, PARAMETER :: KCD_ben = 0.0  !WARNING
!
!-------------------------------------------------------------------------------
!BEGIN

 DO mag_i=1,data%num_malgae
   ! Retrieve current group of water macroalgae (mmol C/m3)
   malg = _STATE_VAR_(data%id_p(mag_i))

   ! Drag addition depending on amount of carbon in water volume
   drag = drag + KCD_wc * malg

   ! Process benthic / attached macroalgae
   IF ( data%malgs(mag_i)%settling == _MOB_ATTACHED_ ) THEN
     dz   = MAX(_STATE_VAR_(data%id_dz),0.05)       ! cell depth
     malg = _STATE_VAR_S_(data%id_pben(mag_i))/dz  ! attached macroalgae
     drag = drag + (KCD_ben * (malg/dz) )
   ENDIF
 ENDDO

END SUBROUTINE aed_bio_drag_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE cgm_calculate_benthic_cladophora(data,column,layer_idx,cgm,pf,rf, &
                                                            pu,nu,mag_in,mag_ip)
!-------------------------------------------------------------------------------
! Calculate benthic cladophora productivity, based on the Cladophora Growth
! Model (CGM: Higgens et al 2006). Units per surface area (not volume!)
! and in mmol C/m2, rather than the original g DM/m2.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx,cgm
   AED_REAL, INTENT(inout) :: pf,rf,pu,nu,mag_in,mag_ip
!
!LOCALS
   AED_REAL, PARAMETER :: a01 = -0.7171732e-01 !
   AED_REAL, PARAMETER :: a02 =  0.6610416     !
   AED_REAL, PARAMETER :: a03 =  0.4797417     !
   AED_REAL, PARAMETER :: a04 = -0.2358428e+01 !
   AED_REAL, PARAMETER :: a05 =  0.4662903e+01 !
   AED_REAL, PARAMETER :: a06 = -0.2759110e+01 !
   AED_REAL, PARAMETER :: a07 =  0.3343594e+01 !
   AED_REAL, PARAMETER :: a08 = -0.5609280e+01 !
   AED_REAL, PARAMETER :: a09 = -0.3914500e+01 !
   AED_REAL, PARAMETER :: a10 =  0.3372881e+01 !
   AED_REAL, PARAMETER :: a11 = -0.1742290e+01 !
   AED_REAL, PARAMETER :: a12 =  0.1304549e+01 !
   AED_REAL, PARAMETER :: a13 =  0.2618247e+01 !
   AED_REAL, PARAMETER :: a14 =  0.8705460     !
   AED_REAL, PARAMETER :: a15 = -0.1207573e+01 !
   AED_REAL, PARAMETER :: b01 =  0.8204769E-02 !
   AED_REAL, PARAMETER :: b02 =  0.1622813E+00 !
   AED_REAL, PARAMETER :: b03 =  0.5599809E+00 !
   AED_REAL, PARAMETER :: b04 = -0.9230173E-01 !
   AED_REAL, PARAMETER :: b05 = -0.3701472E+00 !
   AED_REAL, PARAMETER :: b06 = -0.1466722E+01 !
   AED_REAL, PARAMETER :: b07 = -0.1021670E+00 !
   AED_REAL, PARAMETER :: b08 =  0.6366414E-01 !
   AED_REAL, PARAMETER :: b09 =  0.9596720E+00 !
   AED_REAL, PARAMETER :: b10 =  0.1632652E+01 !
   AED_REAL, PARAMETER :: b11 =  0.1431378E+00 !
   AED_REAL, PARAMETER :: b12 =  0.6283502E+00 !
   AED_REAL, PARAMETER :: b13 = -0.1520305E+01 !
   AED_REAL, PARAMETER :: b14 =  0.3120532E+00 !
   AED_REAL, PARAMETER :: b15 = -0.7267257E+00 !

   AED_REAL :: malg
   AED_REAL :: extc, dz, par, Io, temp, salinity
   AED_REAL :: matz, nit, amm, frp
   AED_REAL :: X_maxp, Kq
   AED_REAL :: lght, pplt, prlt, pf_MB
   AED_REAL :: macroPAR_Top, macroPAR_Bot
   AED_REAL :: AvgTemp, AvgLight
   AED_REAL :: sf,tf,lf
!
!-------------------------------------------------------------------------------
!BEGIN

   !----------------------------------------------------------------------------
   !-- Check this cell is in an active zone for cladophoras
   matz = _STATE_VAR_S_(data%id_sedzone)
   IF ( .NOT. in_zone_set(matz, data%active_zones) ) RETURN

   !-- Retrieve current environmental conditions
   salinity = _STATE_VAR_(data%id_sal)        ! local salinity
   temp = _STATE_VAR_(data%id_tem)            ! local temperature
   extc = _STATE_VAR_(data%id_extc)           ! extinction coefficent
   par = MAX(_STATE_VAR_(data%id_par),zero_)  ! photosynthetically active radn
   Io = MAX(_STATE_VAR_S_(data%id_I_0),zero_) ! surface (incident) shortwave radn
   dz = _STATE_VAR_(data%id_dz)               ! dz = cell/layer thickness

   !----------------------------------------------------------------------------
   !-- MACROALGAL BED DEPTH AND EXTINCTION
   !## Resolve the benthic light climate
   !   par = top of cells
   !   macroPAR_Top = top of macroalgal bed
   !   macroPAR_Bot = bottom of macroalgal bed

   ! Retrieve (local) cladophora in gDM/m2 (Xcc=1/0.25 to convert to DM)
   malg = _STATE_VAR_S_(data%id_pben(cgm)) * (12. / 1e3) / data%malgs(cgm)%Xcc

   ! Empirical formulation from Higgin's CGM model
   macroHgt = (1./100.) * 1.467 * ( malg )**0.425
   macroExt = 7.840 * ( malg )**0.240   ! Higgins et al 2006 Fig 6

   macroPAR_Top = MAX( MIN( par*exp(-extc*( MAX(dz-macroHgt,zero_))),Io), zero_)
   macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)
   !macroPAR_Bot = par * exp(-(extc+macroExt)*(dz-macroHgt))

   IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC)  = macroExt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dHGT) = macroHgt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dPAR) = macroPAR_Top

   !----------------------------------------------------------------------------
   !-- MACROALGAL C, N, P
   !## Get the CGM (cladophora) group concentrations
   malg = _STATE_VAR_S_(data%id_pben(cgm))      ! cladophora biomass (mmolC/m2)
   IF(data%malgs(cgm)%simINDynamics/=0) THEN
     mag_in = _STATE_VAR_S_(data%id_inben(cgm)) ! cladophora biomass (mmolN/m2)
   ELSE
     mag_in = data%malgs(cgm)%X_ncon * malg
   ENDIF
   IF(data%malgs(cgm)%simIPDynamics/=0) THEN
     mag_ip = _STATE_VAR_S_(data%id_ipben(cgm)) ! cladophora biomass (mmolP/m2)
   ELSE
     mag_ip = data%malgs(cgm)%X_pcon * malg
   ENDIF

   !----------------------------------------------------------------------------
   !-- Update moving average for daily temp
   !   (averaged over the past 1 day)
   _DIAG_VAR_S_(data%id_tem_avg) = _DIAG_VAR_S_(data%id_tem_avg) &
                                * (1-(DTday/TempAvgTime)) + temp*(DTday/TempAvgTime)
   AvgTemp = _DIAG_VAR_S_(data%id_tem_avg)

   sf = zero_

   !----------------------------------------------------------------------------
   !-- Now check light, and perform specific day or night activities
   IF(Io > 10.0) THEN

      !-------------------------------------------------------------------------
      !-- Calculate the nutrient limitation (phosporus & nitrogen) and
      !   find the most limiting
      IF (malg > zero_) THEN
        sf = mag_ip/malg
      ENDIF
      IF (sf > zero_) THEN
        sf = 1.0 - (data%malgs(cgm)%X_pmin/sf)
      ENDIF
      IF (diag_level>9) _DIAG_VAR_S_(data%id_fPho_ben(cgm)) =  sf

      !-------------------------------------------------------------------------
      !-- Update moving avg for daily light (averaged over the photo-period, PP)

      _DIAG_VAR_S_(data%id_par_avg) = _DIAG_VAR_S_(data%id_par_avg) &
                          * (1-(DTday/LgtAvgTime)) + macroPAR_Top*(DTday/LgtAvgTime)
      AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83   ! AvgLight in uE for CGM

      !-------------------------------------------------------------------------
      !-- Self-shading

      ! Depth (light) based carrying capacity amount, computed empirically with
      ! *Xcc to convert g DM to g C, /12 to get to molC, and 1e3 to get to mmol
      X_maxp = ( 1.18 * AvgLight - 58.7 ) * data%malgs(cgm)%Xcc * 1e3 / 12.

      lf = 1.
      IF(malg>100. .AND. AvgLight>60. ) lf = 1.0 - malg/X_maxp

      IF(lf < zero_) lf = zero_  ;  IF(lf > one_) lf = one_
      IF(diag_level>9) _DIAG_VAR_S_(data%id_fI_ben(cgm)) =  lf

      IF(diag_level>9) _DIAG_VAR_S_(data%id_fSal_ben(cgm)) = macroPAR_Bot/MAX(macroPAR_Top,1.)

      !-------------------------------------------------------------------------
      !-- Now light and temperature function for photosynthesis

      temp = AvgTemp/35.0

      lght = MIN(600.0,AvgLight)/1235.0      ! Capped at 600: Higgins et al 2006

      pplt = a01                             &
           + a02 * temp                      &
           + a03 * lght                      &
           + a04 * temp * temp               &
           + a05 * temp * lght               &
           + a06 * lght * lght               &
           + a07 * temp * temp * temp        &
           + a08 * temp * temp * lght        &
           + a09 * temp * lght * lght        &
           + a10 * lght * lght * lght        &
           + a11 * temp * temp * temp * temp &
           + a12 * temp * temp * temp * lght &
           + a13 * temp * temp * lght * lght &
           + a14 * temp * lght * lght * lght &
           + a15 * lght * lght * lght * lght

      pf = (data%malgs(cgm)%R_growth * pplt * 5.43) * sf * lf

      !-------------------------------------------------------------------------
      !-- Now light and temperature function for daytime respiration

      prlt = b01                             &
           + b02 * temp                      &
           + b03 * lght                      &
           + b04 * temp * temp               &
           + b05 * temp * lght               &
           + b06 * lght * lght               &
           + b07 * temp * temp * temp        &
           + b08 * temp * temp * lght        &
           + b09 * temp * lght * lght        &
           + b10 * lght * lght * lght        &
           + b11 * temp * temp * temp * temp &
           + b12 * temp * temp * temp * lght &
           + b13 * temp * temp * lght * lght &
           + b14 * temp * lght * lght * lght &
           + b15 * lght * lght * lght * lght

     !rf = 0.44 * prlt * 4.52
      rf= zero_
     !IF (diag_level>9) _DIAG_VAR_S_(data%id_fSal_ben(cgm)) = 0.44 * prlt * 4.52

    ELSE
      !-------------------------------------------------------------------------
      !-- Now calculate (night-time) basal respiration....

      rf = data%malgs(cgm)%R_resp * (0.025 * AvgTemp + 0.1)
      pf = zero_

    END IF

    !---------------------------------------------------------------------------
    !-- Get the temperature function for nutrient uptake
    temp = AvgTemp
    IF(temp < 18.0)THEN
      tf = exp((temp-18.0)/39.00)
    ELSE
      tf = exp((18.0-temp)/18.75)
    ENDIF
    IF (diag_level>9) _DIAG_VAR_S_(data%id_fT_ben(cgm)) =  tf

    !---------------------------------------------------------------------------
    !-- Get the INTERNAL PHOSPHORUS stores for the macroalgae groups.
    !   Recall that int. nutrient is in mol and must be converted to molP/molC
    !   by division by macroalgae biomass for the nutrient limitation / uptake

    !-- Compute the internal phosphorus ratio
    IF(malg>zero_) THEN
      pu = mag_ip/malg
    ELSE
      pu = data%malgs(cgm)%X_pmin
    ENDIF

    Kq = 0.0028 * (12e3/31e3) ! 0.07% = 0.0028gP/gC & 12/31 is mol wgt conversion

    frp = _STATE_VAR_(data%id_Pupttarget(1))

    !-- IPmax = Kq; IPmin = Qo; KP = Km; R_puptake = pmax; tau = tf
    pu = data%malgs(cgm)%R_puptake * tf * (frp/(frp + data%malgs(cgm)%K_P)) &
                          * (Kq /(Kq + MAX(pu - data%malgs(cgm)%X_pmin,zero_)))

    !---------------------------------------------------------------------------
    !-- Get the INTERNAL NITROGEN stores for the macroalgae groups.
    !   Recall that internal nutrient is in g and must be converted to g N/g C
    !   by division by macroalgae biomass for the nitrogen limitation

    !-- Get the internal nitrogen ratio.
    IF(malg>zero_) THEN
      nu = mag_in/malg
    ELSE
      nu = data%malgs(cgm)%X_nmin
    ENDIF

    !-- Macroalgae nitrogen uptake
    nu = 16*pu  ! Overwritten above to match C

    IF (diag_level>9) _DIAG_VAR_S_(data%id_fNit_ben(cgm)) =  one_

    !nu = data%malgdata(i)%R_nuptake  * malg * tf * (INmax(bb) - nu) &
    !    / (INmax(bb)-INmin(bb)) * (nit + amm) ) / (nit + amm + KN(bb))

    !-- Update the internal nitrogen store
    !_FLUX_VAR_B_(data%id_mag_in(mag_i)) =  _FLUX_VAR_B_(data%id_mag_in(mag_i))&
    !                                      +  nu * malg
    IF(data%malgs(cgm)%simINDynamics/=0) &
      _STATE_VAR_S_(data%id_inben(cgm)) =  0.9 * data%malgs(cgm)%X_nmax * malg


    !---------------------------------------------------------------------------
    !-- As base of filaments dies, sloughing of live cells into the water occurs
    pf_MB = zero_

    !-- first need to examine light & temp at the bottom of the bed
    temp = AvgTemp / 35.0

    lght = MIN(600.0,macroPAR_Bot) / 1235.0

    pplt = a01                             &
         + a02 * temp                      &
         + a03 * lght                      &
         + a04 * temp * temp               &
         + a05 * temp * lght               &
         + a06 * lght * lght               &
         + a07 * temp * temp * temp        &
         + a08 * temp * temp * lght        &
         + a09 * temp * lght * lght        &
         + a10 * lght * lght * lght        &
         + a11 * temp * temp * temp * temp &
         + a12 * temp * temp * temp * lght &
         + a13 * temp * temp * lght * lght &
         + a14 * temp * lght * lght * lght &
         + a15 * lght * lght * lght * lght

    pf_MB = (data%malgs(cgm)%R_growth * pplt * 5.43) * sf - rf

    IF(malg > data%malgs(cgm)%p0) THEN
      !-- SloughTrigger = SloughTrigger + pf_MB * DTday
      !_FLUX_VAR_B_(data%id_slough_trig) = _FLUX_VAR_B_(data%id_slough_trig) + pf_MB
      _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) + pf_MB*DTday*secs_per_day
    ENDIF

    sf = one_

END SUBROUTINE cgm_calculate_benthic_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE cgm_slough_cladophora(data,column,layer_idx,cgm,slough_rate,hf)
!-------------------------------------------------------------------------------
! Compute the fraction of cladophora biomass that is detatched due to
!  sloughing if the stress is enough
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx,cgm
   AED_REAL,INTENT(in) :: slough_rate
   AED_REAL,INTENT(out) :: hf
!
!LOCALS
   ! Environment
   AED_REAL :: bottom_stress
   ! State
   AED_REAL :: malg
   ! Parameters
   AED_REAL :: slough_trigger, AvgStress, X_maxp, DTsec, AvgLight
!-------------------------------------------------------------------------------
!BEGIN

   hf = zero_
   DTsec = DTday * secs_per_day

   !----------------------------------------------------------------------------
   !-- Retrieve current bottom shear stress condition within the cell
   bottom_stress = MIN( _STATE_VAR_S_(data%id_taub), 100. )

   !-- Update moving average for stress (averaged over the past 2 hrs)
   _DIAG_VAR_S_(data%id_tau_avg) = _DIAG_VAR_S_(data%id_tau_avg) &
                        * (1-(DTday/StrAvgTime)) + bottom_stress *(DTday/StrAvgTime)
   AvgStress = _DIAG_VAR_S_(data%id_tau_avg)

   !-- Retrieve current (local) state variable values.
   malg = _STATE_VAR_S_(data%id_pben(cgm))
   slough_trigger = _DIAG_VAR_S_(data%id_slough_trig)

   !-- Check if growth phase; if so clip to 0
   IF(slough_trigger > zero_) _DIAG_VAR_S_(data%id_slough_trig) = zero_

   !-- Slough off weakened filaments (those with cumulative respiration excess)
   IF(slough_trigger < data%slough_stress) THEN
       hf = 0.95/DTsec
       _DIAG_VAR_S_(data%id_slough_trig) = zero_
       RETURN
   ELSE
       hf = zero_
   ENDIF

   !-- Sloughing of even healthy filaments if the shear stress is high enough
   IF(malg>data%malgs(cgm)%p0 .AND. slough_trigger<-0.001 .AND. AvgStress>data%malgs(cgm)%tau_0)THEN

      !-------------------------------------------------------------------------
      ! Depth (light) based carrying capacity amount, computed empirically with
      ! *0.25 to get from g DM to g C, /12 to get mol C, and 1e3 to get to mmol
      AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83 ! AvgLight in uE for CGM
      X_maxp = ( 1.18 * AvgLight - 58.7 ) * data%malgs(cgm)%Xcc * 1e3 / 12.

      ! Daily slough rate
      ! Lss = Lmax * (T/Tcrit) * (X/Xmax)
      ! hf = 0.242 * exp(-0.3187 * DEPTH) * AvgStress / PCm  * malg / X_maxp
      hf = (AvgStress - data%malgs(cgm)%tau_0) * MIN( malg/MAX(X_maxp,data%malgs(cgm)%p0),one_ )

      hf = slough_rate * hf / secs_per_day  ! daily biomass fraction sloughed per sec

      IF (hf*DTsec>0.95) hf = 0.95/DTsec    ! no more than 95% slough in one interval

      !-- Update the malg and slough variables with following the slough event
      _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) * 0.5
   ENDIF

END SUBROUTINE cgm_slough_cladophora
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_macroalgae
