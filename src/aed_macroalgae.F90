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
!#  Copyright 2017 - 2022 -  The University of Western Australia               #
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
!#   For use in a commercial seeting please contact the authors.               #
!#                                                                             #
!#   -----------------------------------------------------------------------   #
!#                                                                             #
!# Originally created August 2017 by Matthew Hipsey, UWA                       #
!# Follow updates @ https://github.com/AquaticEcoDynamics/libaed-benthic       #
!#                                                                             #
!###############################################################################

#include "aed.h"

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
                             id_fPho_ben(:), id_fSal_ben(:), id_dw_ben(:)
      INTEGER :: id_Pexr, id_Pmor, id_Pupttarget(1:2)
      INTEGER :: id_Nexr, id_Nmor, id_Nupttarget(1:4)
      INTEGER :: id_Cexr, id_Cmor, id_Cupttarget
      INTEGER :: id_Siexctarget,id_Simorttarget,id_Siupttarget
      INTEGER :: id_DOupttarget
      INTEGER :: id_par, id_I_0, id_extc, id_taub, id_sedzone
      INTEGER :: id_tem, id_sal, id_dz, id_dens, id_yearday, id_depth, id_dt
      INTEGER :: id_GPP, id_NMP, id_PPR, id_NPR, id_dPAR, id_bPAR, id_dHGT
      INTEGER :: id_TMALG, id_dEXTC, id_TIN, id_TIP, id_MPB, id_d_MPB, id_d_BPP
      INTEGER :: id_NUP, id_PUP, id_CUP
      INTEGER :: id_mhsi
      INTEGER :: id_mag_ben, id_min_ben, id_mip_ben
      INTEGER :: id_nup_ben, id_pup_ben, id_gpp_ben, id_rsp_ben, id_nmp_ben
      INTEGER :: id_slough_trig, id_slough_tsta, id_slough_days
      INTEGER :: id_tem_avg, id_tau_avg, id_par_avg, id_slg_ben, id_par_bot, id_gpp_bot
      INTEGER :: id_swi_c, id_swi_n, id_swi_p

      !# Model parameters and options
      TYPE(phyto_data_t),DIMENSION(:),ALLOCATABLE :: malgs
      INTEGER  :: num_malgae
      LOGICAL  :: do_Puptake, do_Nuptake, do_Cuptake
      LOGICAL  :: do_Siuptake, do_DOuptake, do_N2uptake
      LOGICAL  :: do_Pmort, do_Nmort, do_Cmort, do_Simort
      LOGICAL  :: do_Pexc, do_Nexc, do_Cexc, do_Siexc
      INTEGER  :: nnup, npup
      INTEGER  :: n_zones

      AED_REAL, ALLOCATABLE :: active_zones(:)
      AED_REAL :: min_rho,max_rho
      AED_REAL :: slough_stress, slough_burial, slough_rate
      LOGICAL  :: simMalgFeedback
      INTEGER  :: simSloughing
      INTEGER  :: simMalgHSI
      INTEGER  :: simCGM

     CONTAINS
         PROCEDURE :: define             => aed_define_macroalgae
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_macroalgae
         PROCEDURE :: calculate          => aed_calculate_macroalgae
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_macroalgae
         PROCEDURE :: mobility           => aed_mobility_macroalgae
         PROCEDURE :: light_extinction   => aed_light_extinction_macroalgae
         PROCEDURE :: bio_drag           => aed_bio_drag_macroalgae
        !PROCEDURE :: delete             => aed_delete_macroalgae

   END TYPE

   ! MODULE GLOBALS
   AED_REAL :: DTday       = (15./60.)/24.  ! 15 minutes
   AED_REAL, PARAMETER :: StrAvgTime  = 2.0/24.0       !  2 hours
   !AED_REAL, PARAMETER :: LgtAvgTime  = 0.5            ! 12 hours
   AED_REAL, PARAMETER :: LgtAvgTime  =  2.0/24.0       !  2 hours
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
INTEGER FUNCTION load_csv(dbase, pd)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(phyto_param_t) :: pd(MAX_PHYTO_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol, dcol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0
   LOGICAL :: meh
   INTEGER :: ret = 0
!
!BEGIN
!-------------------------------------------------------------------------------
   unit = aed_csv_read_header(dbase, csvnames, nccols)
   IF (unit <= 0) THEN
      load_csv = -1
      RETURN !# No file found
   ENDIF

   ALLOCATE(values(nccols))

   DO WHILE ( aed_csv_read_row(unit, values) )
      DO ccol=2,nccols
         dcol = ccol - 1
         pd(dcol)%p_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('p_initial')     ; pd(dcol)%p_initial     = extract_double(values(ccol))
            CASE ('p0')            ; pd(dcol)%p0            = extract_double(values(ccol))
            CASE ('w_p')           ; pd(dcol)%w_p           = extract_double(values(ccol))
            CASE ('Xcc')           ; pd(dcol)%Xcc           = extract_double(values(ccol))
            CASE ('R_growth')      ; pd(dcol)%R_growth      = extract_double(values(ccol))
            CASE ('fT_Method')     ; pd(dcol)%fT_Method     = extract_integer(values(ccol))
            CASE ('theta_growth')  ; pd(dcol)%theta_growth  = extract_double(values(ccol))
            CASE ('T_std')         ; pd(dcol)%T_std         = extract_double(values(ccol))
            CASE ('T_opt')         ; pd(dcol)%T_opt         = extract_double(values(ccol))
            CASE ('T_max')         ; pd(dcol)%T_max         = extract_double(values(ccol))
            CASE ('lightModel')    ; pd(dcol)%lightModel    = extract_integer(values(ccol))
            CASE ('I_K')           ; pd(dcol)%I_K           = extract_double(values(ccol))
            CASE ('I_S')           ; pd(dcol)%I_S           = extract_double(values(ccol))
            CASE ('KePHY')         ; pd(dcol)%KePHY         = extract_double(values(ccol))
            CASE ('f_pr')          ; pd(dcol)%f_pr          = extract_double(values(ccol))
            CASE ('R_resp')        ; pd(dcol)%R_resp        = extract_double(values(ccol))
            CASE ('theta_resp')    ; pd(dcol)%theta_resp    = extract_double(values(ccol))
            CASE ('k_fres')        ; pd(dcol)%k_fres        = extract_double(values(ccol))
            CASE ('k_fdom')        ; pd(dcol)%k_fdom        = extract_double(values(ccol))
            CASE ('salTol')        ; pd(dcol)%salTol        = extract_integer(values(ccol))
            CASE ('S_bep')         ; pd(dcol)%S_bep         = extract_double(values(ccol))
            CASE ('S_maxsp')       ; pd(dcol)%S_maxsp       = extract_double(values(ccol))
            CASE ('S_opt')         ; pd(dcol)%S_opt         = extract_double(values(ccol))
            CASE ('simDINUptake')  ; pd(dcol)%simDINUptake  = extract_integer(values(ccol))
            CASE ('simDONUptake')  ; pd(dcol)%simDONUptake  = extract_integer(values(ccol))
            CASE ('simNFixation')  ; pd(dcol)%simNFixation  = extract_integer(values(ccol))
            CASE ('simINDynamics') ; pd(dcol)%simINDynamics = extract_integer(values(ccol))
            CASE ('N_o')           ; pd(dcol)%N_o           = extract_double(values(ccol))
            CASE ('K_N')           ; pd(dcol)%K_N           = extract_double(values(ccol))
            CASE ('X_ncon')        ; pd(dcol)%X_ncon        = extract_double(values(ccol))
            CASE ('X_nmin')        ; pd(dcol)%X_nmin        = extract_double(values(ccol))
            CASE ('X_nmax')        ; pd(dcol)%X_nmax        = extract_double(values(ccol))
            CASE ('R_nuptake')     ; pd(dcol)%R_nuptake     = extract_double(values(ccol))
            CASE ('k_nfix')        ; pd(dcol)%k_nfix        = extract_double(values(ccol))
            CASE ('R_nfix')        ; pd(dcol)%R_nfix        = extract_double(values(ccol))
            CASE ('simDIPUptake')  ; pd(dcol)%simDIPUptake  = extract_integer(values(ccol))
            CASE ('simIPDynamics') ; pd(dcol)%simIPDynamics = extract_integer(values(ccol))
            CASE ('P_0')           ; pd(dcol)%P_0           = extract_double(values(ccol))
            CASE ('K_P')           ; pd(dcol)%K_P           = extract_double(values(ccol))
            CASE ('X_pcon')        ; pd(dcol)%X_pcon        = extract_double(values(ccol))
            CASE ('X_pmin')        ; pd(dcol)%X_pmin        = extract_double(values(ccol))
            CASE ('X_pmax')        ; pd(dcol)%X_pmax        = extract_double(values(ccol))
            CASE ('R_puptake')     ; pd(dcol)%R_puptake     = extract_double(values(ccol))
            CASE ('simSiUptake')   ; pd(dcol)%simSiUptake   = extract_integer(values(ccol))
            CASE ('Si_0')          ; pd(dcol)%Si_0          = extract_double(values(ccol))
            CASE ('K_Si')          ; pd(dcol)%K_Si          = extract_double(values(ccol))
            CASE ('X_sicon')       ; pd(dcol)%X_sicon       = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit)
   !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_macroalgae_load_params(data, dbase, count, list, settling,      &
                                 growth_form, slough_model, resuspension, tau_0)
!-------------------------------------------------------------------------------
   USE aed_util,ONLY : param_file_type, CSV_TYPE, NML_TYPE
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
   INTEGER,INTENT(in)          :: settling(*)
   INTEGER,INTENT(in)          :: growth_form(*)
   INTEGER,INTENT(in)          :: slough_model(*)
   AED_REAL,INTENT(in)         :: resuspension(*)
   AED_REAL,INTENT(in)         :: tau_0(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(phyto_param_t),ALLOCATABLE :: pd(:)
   NAMELIST /malgae_data/ pd    ! %% phyto_param_t - see aed_bio_utils
!-------------------------------------------------------------------------------
!BEGIN
    ALLOCATE(pd(MAX_PHYTO_TYPES))
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, pd)
       CASE (NML_TYPE)
           print*,"nml format parameter file is deprecated. Please update to CSV format"
           tfil = find_free_lun()
           open(tfil,file=dbase, status='OLD', iostat=status)
           IF (status /= 0) STOP 'Cannot open malgae_data namelist file: ' !,dbase
           read(tfil,nml=malgae_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist malgae_data'

    !---------------------------------------------------------------------------
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
       ALLOCATE(data%id_dw_ben(count)) ; data%id_dw_ben(:) = 0
    ENDIF

    !---------------------------------------------------------------------------
    DO i=1,count
       ! First, assign parameters from database and nml to simulated groups
       data%malgs(i)%p_name       = pd(list(i))%p_name
       data%malgs(i)%p0           = pd(list(i))%p0
       data%malgs(i)%w_p          = pd(list(i))%w_p/secs_per_day
       data%malgs(i)%growth_form  = growth_form(i)
       data%malgs(i)%slough_model = slough_model(i)
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

       !-- Declare variables based on growth form, and sloughing
       !   0 : water only
       !   1 : benthic only
       !   2 : water column
       !   3 : surface

       IF(TRIM(data%malgs(i)%p_name)== 'cgm') THEN
           data%simCGM = i
            ! If sloughing is requested for CGM force to 2 ... for now
           IF(data%simSloughing>0) data%simSloughing = 2
           IF(data%simSloughing>0 .and. (data%malgs(i)%slough_model < 2 .or. data%malgs(i)%slough_model > 4)) data%malgs(i)%slough_model = 2
       ENDIF

       !-- Group requires a water column / pelagic pool
       IF ( growth_form(i)==0 .or. slough_model(i)>0 ) THEN

         ! Register a water column pool for the group as a state variable
         data%id_p(i) = aed_define_variable(                                   &
                              TRIM(data%malgs(i)%p_name),                      &
                              'mmol/m**3',                                     &
                              'macroalgae '//TRIM(data%malgs(i)%p_name),       &
                              pd(list(i))%p0,                                  &
                             !pd(list(i))%p_initial,                           &
                              minimum = pd(list(i))%p0,                        &
                              mobility = data%malgs(i)%w_p)

         ! Register rho (internal density) as a state variable, if required
         IF (data%malgs(i)%settling == _MOB_STOKES_) THEN
           data%id_rho(i) = aed_define_variable(                                &
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

       ENDIF


       !-- Group requires a attached pool, benthic or surface
       IF ( growth_form(i)>0 ) THEN
         ! Register a water column pool for the group as a state variable
        !IF (data%malgs(i)%settling == _MOB_ATTACHED_) THEN
         data%id_pben(i) = aed_define_sheet_variable(                          &
                          TRIM(data%malgs(i)%p_name)//'_ben',                  &
                          'mmolC/m**2',                                        &
                          'macroalgae '//TRIM(data%malgs(i)%p_name),           &
                          pd(list(i))%p_initial,                               &
                          minimum=pd(list(i))%p0,                              &
                          maximum=1e8 )

         IF (data%malgs(i)%simIPDynamics /= 0) THEN
            minNut = data%malgs(i)%p0*data%malgs(i)%X_pmin
            ! Register IP group as a state variable
            data%id_ipben(i) = aed_define_sheet_variable(                      &
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
           data%id_inben(i) = aed_define_sheet_variable(                       &
                          TRIM(data%malgs(i)%p_name)//'_IN_ben',               &
                          'mmol/m**2',                                         &
                          'macroalgae '//TRIM(data%malgs(i)%p_name)//'_IN_ben',&
                          pd(list(i))%p_initial*data%malgs(i)%X_ncon,          &
                          minimum=minNut)
         ENDIF
         ! Group specific diagnostics for the benthic variables
         IF (diag_level>9) THEN
          data%id_fI_ben(i)   = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fI_ben', '-', 'fI (0-1)')
          data%id_fT_ben(i)   = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fT_ben', '-', 'fT (>0)')
          data%id_fNit_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fNit_ben', '-', 'fNit (0-1)')
          data%id_fPho_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fPho_ben', '-', 'fPho (0-1)')
          data%id_fSal_ben(i) = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_fSal_ben', '-', 'fSal (-)')
          data%id_c2p_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2p_ben', '-', 'C:P')
          data%id_c2n_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_c2n_ben', '-', 'C:N')
          data%id_n2p_ben(i)  = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_n2p_ben', '-', 'N:P')
          ! temporary add in so we can see the state var benthic
          data%id_dw_ben(i)   = aed_define_sheet_diag_variable( TRIM(data%malgs(i)%p_name)//'_dw_ben', '-', '-')

         ENDIF

       ENDIF
    ENDDO
    !---------------------------------------------------------------------------
    DEALLOCATE(pd)
END SUBROUTINE aed_macroalgae_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_macroalgae(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macroalgae biogeochemical model
!
!  Here, the aed_macroalgae namelist is read and the variables exported
!  by the model are registered with AED.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macroalgae_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER            :: status, i

!  %% NAMELIST   %%  /aed_macroalgae/
!  %% Last Checked 20/08/2021
   INTEGER            :: num_malgae = 0
   INTEGER            :: the_malgae(MAX_PHYTO_TYPES) = 0
   INTEGER            :: settling(MAX_PHYTO_TYPES) =  _MOB_CONST_
   INTEGER            :: slough_model(MAX_PHYTO_TYPES) = 0
   INTEGER            :: growth_form(MAX_PHYTO_TYPES) = 0
   AED_REAL           :: resuspension(MAX_PHYTO_TYPES) = 0.
   AED_REAL           :: tau_0(MAX_PHYTO_TYPES) = 0.1
   CHARACTER(len=64)  :: p_excretion_target_variable = ''
   CHARACTER(len=64)  :: p_mortality_target_variable = ''
   CHARACTER(len=64)  :: p1_uptake_target_variable = ''
   CHARACTER(len=64)  :: p2_uptake_target_variable = ''
   CHARACTER(len=64)  :: n_excretion_target_variable = ''
   CHARACTER(len=64)  :: n_mortality_target_variable = ''
   CHARACTER(len=64)  :: n1_uptake_target_variable = ''
   CHARACTER(len=64)  :: n2_uptake_target_variable = ''
   CHARACTER(len=64)  :: n3_uptake_target_variable = ''
   CHARACTER(len=64)  :: n4_uptake_target_variable = ''
   CHARACTER(len=64)  :: c_excretion_target_variable = ''
   CHARACTER(len=64)  :: c_mortality_target_variable = ''
   CHARACTER(len=64)  :: c_uptake_target_variable = ''
   CHARACTER(len=64)  :: do_uptake_target_variable = ''
   CHARACTER(len=128) :: dbase = 'aed_malgae_pars.nml'
   AED_REAL           :: zerolimitfudgefactor = 15.*60.
   AED_REAL           :: min_rho = 900.
   AED_REAL           :: max_rho = 1200.
   AED_REAL           :: slough_stress = -1.0
   AED_REAL           :: slough_burial = zero_
   AED_REAL           :: slough_rate = 0.08 ! %/day
   INTEGER            :: simMalgHSI = 0
   INTEGER            :: simSloughing = 0
   INTEGER            :: n_zones = 0
   INTEGER            :: active_zones(1000) = 0
   LOGICAL            :: simMalgFeedback = .true.

! From Module Globals
   LOGICAL  :: extra_debug = .false.  !## Obsolete Use diag_level = 10
!  LOGICAL  :: extra_diag = .false.   !## Obsolete Use diag_level = 10
!  INTEGER  :: diag_level = 10                ! 0 = no diagnostic outputs
                                              ! 1 = basic diagnostic outputs
                                              ! 2 = flux rates, and supporitng
                                              ! 3 = other metrics
                                              !10 = all debug & checking outputs
!  %% END NAMELIST   %%  /aed_macroalgae/

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
                    extra_debug, extra_diag, diag_level, tau_0, dtlim,         &
                    growth_form, slough_model, slough_burial, slough_rate
!-----------------------------------------------------------------------
!BEGIN

   print *,"        aed_macroalgae initialization"

   dtlim = zerolimitfudgefactor

   ! Read the namelist, and set module parameters
   read(namlst,nml=aed_macroalgae,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_macroalgae'
   IF( extra_debug )  extra_diag = .true.       ! legacy use of extra_debug
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
   data%slough_burial = slough_burial / secs_per_day
   data%slough_rate = slough_rate / secs_per_day

   data%simMalgFeedback = simMalgFeedback
   PRINT *,'          NOTE - macroalgae feedbacks to water column properties: ',simMalgFeedback

   ! Store parameter values in a local malgae strcutured type
   ! Note: all rates must be provided in values per day,
   !       but are converted in here to rates per second
   CALL aed_macroalgae_load_params(data,dbase,num_malgae,the_malgae,settling,  &
                                    growth_form,slough_model,resuspension,tau_0)

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
     IF (data%npup>0) &
        data%id_Pupttarget(ifrp) = aed_locate_variable(p1_uptake_target_variable) !; ifrp=1  ! CAB - now constants in bio utils
     IF (data%npup>1) &
        data%id_Pupttarget(idop) = aed_locate_variable(p2_uptake_target_variable) !; idop=2  ! CAB - now constants in bio utils
   ENDIF
   data%nnup = 0
   IF (n1_uptake_target_variable .NE. '') data%nnup = 1
   IF (n2_uptake_target_variable .NE. '') data%nnup = 2
   IF (n3_uptake_target_variable .NE. '') data%nnup = 3
   IF (n4_uptake_target_variable .NE. '') data%nnup = 4
   data%do_Nuptake = .false.
   IF (data%nnup>0) data%do_Nuptake=.true.
   IF (data%do_Nuptake) THEN
     IF (data%nnup>0) data%id_Nupttarget(ino3) = aed_locate_variable( n1_uptake_target_variable)
     IF (data%nnup>1) data%id_Nupttarget(inh4) = aed_locate_variable( n2_uptake_target_variable)
     IF (data%nnup>2) data%id_Nupttarget(idon) = aed_locate_variable( n3_uptake_target_variable)
     IF (data%nnup>3) data%id_Nupttarget(in2)  = aed_locate_variable( n4_uptake_target_variable)
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

   data%id_swi_c = aed_define_sheet_diag_variable('mag_swi_c','mmol C/m2/d', 'MAG C flux across the SWI')
   data%id_swi_n = aed_define_sheet_diag_variable('mag_swi_n','mmol N/m2/d', 'MAG N flux across the SWI')
   data%id_swi_p = aed_define_sheet_diag_variable('mag_swi_p','mmol P/m2/d', 'MAG P flux across the SWI')

   IF ( simMalgHSI>0 ) &
     data%id_mhsi  = aed_define_sheet_diag_variable('HSI','-', 'MAG: macroalgae habitat suitability')

   IF ( data%simCGM >0 ) THEN
     data%id_slough_trig = aed_define_sheet_diag_variable('cgm_sltg','-', 'MAG: cgm slough trigger')
     data%id_slough_tsta = aed_define_sheet_diag_variable('cgm_slst','-', 'MAG: cgm slough start day')
     data%id_slough_days = aed_define_sheet_diag_variable('cgm_sldy','-', 'MAG: cgm slough dark days')
     data%id_tem_avg = aed_define_sheet_diag_variable('cgm_tavg','-', 'MAG: cgm temperature average')
     data%id_par_avg = aed_define_sheet_diag_variable('cgm_lavg','-', 'MAG: cgm light average')
     data%id_par_bot = aed_define_sheet_diag_variable('cgm_parb','-', 'MAG: cgm light bottom')
     data%id_tau_avg = aed_define_sheet_diag_variable('cgm_savg','-', 'MAG: cgm stress average')
     data%id_gpp_bot  = aed_define_sheet_diag_variable('gpp_bot','/day', 'MAG: GPP at canopy base')
   ENDIF

   IF (diag_level>9) THEN
     data%id_dEXTC = aed_define_diag_variable('extc','/m','MAG: extinction due to macroalgae')
     data%id_dPAR  = aed_define_sheet_diag_variable('parc','%', 'MAG: PAR reaching the top of canopy')
     data%id_dHGT  = aed_define_sheet_diag_variable('hgtc','m', 'MAG: height of macroalgal canopy')
   ENDIF

   ! Register the required environmental dependencies
   data%id_dz      = aed_locate_global('layer_ht')
   data%id_tem     = aed_locate_global('temperature')
   data%id_sal     = aed_locate_global('salinity')
   data%id_dens    = aed_locate_global('density')
   data%id_extc    = aed_locate_global('extc_coef')
   data%id_par     = aed_locate_global('par')
   data%id_I_0     = aed_locate_sheet_global('par_sf')
   data%id_taub    = aed_locate_sheet_global('taub')
   data%id_depth   = aed_locate_sheet_global('col_depth') 
   data%id_sedzone = aed_locate_sheet_global('sed_zone')
   data%id_yearday = aed_locate_sheet_global('yearday') 
   data%id_dt      = aed_locate_sheet_global('timestep')  
   
END SUBROUTINE aed_define_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic_macroalgae(data, column, layer_idx)
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
      !IF (data%malgs(mag_i)%settling == _MOB_ATTACHED_) THEN
       IF (data%malgs(mag_i)%growth_form >0) THEN
         _STATE_VAR_S_(data%id_pben(mag_i)) = zero_
         IF (data%malgs(mag_i)%simINDynamics > 0) _STATE_VAR_S_(data%id_inben(mag_i)) = zero_
         IF (data%malgs(mag_i)%simIPDynamics > 0) _STATE_VAR_S_(data%id_ipben(mag_i)) = zero_
       ENDIF
     ENDDO
     RETURN
   ELSE
     DO mag_i=1,data%num_malgae
       IF (data%malgs(mag_i)%growth_form >0) THEN
         IF (data%malgs(mag_i)%simINDynamics > 0) &
          _STATE_VAR_S_(data%id_inben(mag_i)) = _STATE_VAR_S_(data%id_pben(mag_i)) * data%malgs(mag_i)%X_ncon
         IF (data%malgs(mag_i)%simIPDynamics > 0) &
          _STATE_VAR_S_(data%id_ipben(mag_i)) = _STATE_VAR_S_(data%id_pben(mag_i)) * data%malgs(mag_i)%X_pcon
       ENDIF
     ENDDO
     RETURN
   ENDIF
   IF ( data%simCGM >0 ) THEN
     _DIAG_VAR_S_(data%id_tem_avg) = 4.      !_STATE_VAR_S_(data%id_tem)
     _DIAG_VAR_S_(data%id_tau_avg) = 0.001   !_STATE_VAR_S_(data%id_taub)
     _DIAG_VAR_S_(data%id_par_avg) = 0.1     !_STATE_VAR_S_(data%id_par)
     _DIAG_VAR_S_(data%id_slough_days) = 0.0 !
     _DIAG_VAR_S_(data%id_slough_tsta) = 0.0 !
     _DIAG_VAR_S_(data%id_slough_trig) = 0.0 !
    ENDIF

END SUBROUTINE aed_initialize_benthic_macroalgae
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

   !-- Loop through the floating/free macroalgal groups
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

      IF (data%id_p(mag_i)==0) CYCLE

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
         _DIAG_VAR_(data%id_c2n(mag_i))  =  phy/MAX(INi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_nmin)
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

     IF (data%id_p(mag_i)==0) CYCLE

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
        ! Now manage uptake of nutrients, CO2 and DO. Note that these
        ! cumulative fluxes already limited above loop. The uptake arrays
        ! have -ve values in them, so fluxes here will be reducing these pools
        IF (data%do_Puptake) THEN
         DO c = 1,data%npup
            _FLUX_VAR_(data%id_Pupttarget(c)) =    &
                       _FLUX_VAR_(data%id_Pupttarget(c)) + puptake(mag_i,c)
         ENDDO
        ENDIF
        IF (data%do_Nuptake) THEN
         DO c = 1,data%nnup
            _FLUX_VAR_(data%id_Nupttarget(c)) =    &
                       _FLUX_VAR_(data%id_Nupttarget(c)) + nuptake(mag_i,c)
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
         _FLUX_VAR_(data%id_Siupttarget) = _FLUX_VAR_(data%id_Siupttarget) + siuptake(mag_i)
        ENDIF
        ! Now manage mortality contributions to POM
        IF (data%do_Pmort) &
         _FLUX_VAR_(data%id_Pmor) = _FLUX_VAR_(data%id_Pmor) + pmortality(mag_i)
        IF (data%do_Nmort) &
         _FLUX_VAR_(data%id_Nmor) = _FLUX_VAR_(data%id_Nmor) + nmortality(mag_i)
        IF (data%do_Cmort) &
         _FLUX_VAR_(data%id_Cmor) = _FLUX_VAR_(data%id_Cmor) + cmortality(mag_i)
        ! Now manage excretion/exudation contributions to DOM
        IF (data%do_Pexc) &
         _FLUX_VAR_(data%id_Pexr) = _FLUX_VAR_(data%id_Pexr) + pexcretion(mag_i)
        IF (data%do_Nexc) &
         _FLUX_VAR_(data%id_Nexr) = _FLUX_VAR_(data%id_Nexr) + nexcretion(mag_i)
        IF (data%do_Cexc) &
         _FLUX_VAR_(data%id_Cexr) = _FLUX_VAR_(data%id_Cexr) + cexcretion(mag_i)
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
   AED_REAL :: temp,extc,par,dz,Io,bottom_stress,depth,light,matz
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
   AED_REAL :: fI, fT, fNit, fPho, fSil, fXl, fSal, PNf
   AED_REAL :: upTot, net_cuptake, available, flux
   AED_REAL :: slough_frac = zero_
   AED_REAL :: slough_burial, slough, slough_in, slough_ip
   AED_REAL :: aging, aging_rate,R_aging,Kaging,malg_min
   AED_REAL :: DTsec,malgdw,macroPAR_Top,macroPAR_Bot
!
!-------------------------------------------------------------------------------
!BEGIN
   DTday = _STATE_VAR_S_(data%id_dt)/secs_per_day
   DTsec = DTday * secs_per_day

  !-- Benthic light fraction and extinction, for diagnostics
  extc  = _STATE_VAR_  (data%id_extc)    ! extinction coefficient of bottom cell
  depth = _STATE_VAR_  (data%id_dz)      ! water column depth (cell depth if 2D)
  matz  = _STATE_VAR_S_(data%id_sedzone) ! sediment zone / material type

  !-- Initialise the accumulated biomass variables
  IF (diag_level>0) THEN
    _DIAG_VAR_S_(data%id_mag_ben) = zero_
    _DIAG_VAR_S_(data%id_min_ben) = zero_
    _DIAG_VAR_S_(data%id_mip_ben) = zero_
    _DIAG_VAR_(data%id_TMALG) = zero_
   IF (diag_level>1) THEN
     _DIAG_VAR_S_(data%id_gpp_ben) = zero_
     _DIAG_VAR_S_(data%id_nmp_ben) = zero_
     _DIAG_VAR_S_(data%id_nup_ben) = zero_
     _DIAG_VAR_S_(data%id_pup_ben) = zero_
     _DIAG_VAR_S_(data%id_slg_ben) = zero_
    IF(diag_level>9) THEN
      _DIAG_VAR_S_(data%id_dHGT) = zero_
    ENDIF
   ENDIF
  ENDIF
  _DIAG_VAR_(data%id_dEXTC)  = zero_
  _DIAG_VAR_S_(data%id_swi_c) = zero_
  _DIAG_VAR_S_(data%id_swi_n) = zero_
  _DIAG_VAR_S_(data%id_swi_p) = zero_

  !-- Loop through selected groups, determining if benthic/attached is active
  DO mag_i=1,data%num_malgae

     fI = zero_; fNit = zero_; fPho = zero_; fSil = one_; fSal = one_; fXl = one_

     !-- Process each benthic / attached macroalgal group
     IF ( data%malgs(mag_i)%growth_form>0 .AND. &
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

         ! Get the light and nutrient limitation.
         ! NITROGEN.
         IF(data%malgs(mag_i)%simINDynamics /= 0) THEN
           ! Simulated IN variable available
           INi = _STATE_VAR_S_(data%id_inben(mag_i))
         ELSE
           ! Assumed constant IN
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
         !fPho = 1.0

         !IF(bottom_stress>0.09) fXl =0.0  !!!!!!!!!!!

         !----------------------------------------------------------------------------
         !-- MACROALGAL BED DEPTH AND EXTINCTION
         !## Resolve the benthic light climate
         !   par = top of cell/layer
         !   macroPAR_Top = top of macroalgal bed
         !   macroPAR_Bot = bottom of macroalgal bed

         ! Retrieve (local) cladophora in gDM/m2 (Xcc=1/0.25 to convert to DM)
         malgdw = malg * (12./1e3) / data%malgs(mag_i)%Xcc
         IF (diag_level>9)  _DIAG_VAR_S_(data%id_dw_ben(mag_i)) = malgdw

         macroHgt = (1./100.) * 2.4 * ( malgdw )**0.5
         macroExt = data%malgs(mag_i)%KePHY*malg

         IF(data%malgs(mag_i)%growth_form==1) THEN
           macroPAR_Top = MAX( MIN( par*exp(-extc*( MAX(dz-macroHgt,zero_))),Io), zero_)
           macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)
         ELSEIF(data%malgs(mag_i)%growth_form==2) THEN                           ! Assumes contained in bottom cell
           macroPAR_Top = MAX( MIN( par*exp(-extc*( MAX(dz-macroHgt,zero_))),Io), zero_)
           macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)
         ELSEIF(data%malgs(mag_i)%growth_form==3) THEN                           ! Assumes running in 2D (not 3D)
           macroPAR_Top = MAX( Io, zero_)
           macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)
         ENDIF
         _DIAG_VAR_(data%id_dEXTC)  = _DIAG_VAR_(data%id_dEXTC) + macroExt
         IF(diag_level>9) THEN
           _DIAG_VAR_S_(data%id_dHGT) = MAX( _DIAG_VAR_S_(data%id_dHGT), macroHgt )
           IF(data%malgs(mag_i)%growth_form==1) _DIAG_VAR_S_(data%id_dPAR) = macroPAR_Top
         ENDIF

         ! Compute P-I
         fI = photosynthesis_irradiance( data%malgs(mag_i)%lightModel, & ! 0 is vertical integral
                                         data%malgs(mag_i)%I_K,        &
                                         data%malgs(mag_i)%I_S,        &
                                         macroPAR_Top,                 & ! par,
                                         extc+macroExt,                & ! extc,
                                         Io,                           &
                                         macroHgt)                       ! dz)

         ! Primary production rate
         primprod(mag_i) = data%malgs(mag_i)%R_growth * fT * findMin(fI,fNit,fPho,fSil) * fxl !* fSal

         ! Respiration and general metabolic loss
         respiration(mag_i) = bio_respiration(data%malgs(mag_i)%R_resp,data%malgs(mag_i)%theta_resp,temp)

         ! Salinity stress effect on respiration or growth

         !fSal = fSal_function(salinity, minS, Smin, Smax, maxS )
         !fSal = one_  !fSal_function(salinity, 5., 18., 35., 45. )

         !salt = salinity       !COORONG (Ulva) Make me flexible
         !IF( salt<=5. ) THEN
        !   fSal = zero_
        ! ELSE IF ( salt>5. .AND. salt<=18.  ) THEN
        !   fSal = 0. + ( (salt-5.)/(18.-5.) )
         !ELSE IF ( salt>18. .AND. salt<=40. ) THEN
        !   fSal = one_
        ! ELSE IF ( salt>40. .AND. salt<=85. ) THEN
      !     fSal = 1. - ( (salt-40.)/(85.-40.) )
      !   ELSE IF ( salt>85. ) THEN
      !     fSal = zero_
      !   ENDIF
         IF( salt<=5. ) THEN
           fSal = zero_
         ELSE IF ( salt>5. .AND. salt<=15.  ) THEN
           fSal = 0. + ( (salt-5.)/(15.-5.) )
         ELSE IF ( salt>15. .AND. salt<=65. ) THEN
           fSal = one_
         ELSE IF ( salt>65. .AND. salt<=120. ) THEN
           fSal = 1. - ( (salt-65.)/(120.-65.) )
         ELSE IF ( salt>120. ) THEN
           fSal = zero_
         ENDIF


         fSal =  phyto_salinity(data%malgs,mag_i,salinity)
         IF( data%malgs(mag_i)%salTol < 0) THEN
           ! salTol is set for growth supression
           primprod(mag_i) = primprod(mag_i) * fSal
         ELSE
           ! salTol is set for respiration/mortality enhancement
           respiration(mag_i) = respiration(mag_i) * fSal
         ENDIF


         ! Record all environmental limitations
         IF (diag_level>9) THEN
           _DIAG_VAR_S_(data%id_fT_ben(mag_i))   =  fT
           _DIAG_VAR_S_(data%id_fI_ben(mag_i))   =  fI
           _DIAG_VAR_S_(data%id_fNit_ben(mag_i)) =  fNit
           _DIAG_VAR_S_(data%id_fPho_ben(mag_i)) =  fPho
           _DIAG_VAR_S_(data%id_fSal_ben(mag_i)) =  fSal
         ENDIF

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

        ! !# Bottom (nonCGM) macroalgal canopy properties
        ! macroHgt = 0.1
        ! macroExt = data%malgs(mag_i)%KePHY*malg

        ! IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC)  = macroExt
        ! IF(diag_level>9) _DIAG_VAR_S_(data%id_dHGT) = macroHgt
        ! IF(diag_level>9) _DIAG_VAR_S_(data%id_dPAR) = par


       ELSE ! group check for CGM

         !-- User has selected the special CGM approach (Cladophora Growth Model)
         CALL cladophora_calculate_cgm(data,column,layer_idx,mag_i,         &
                                               primprod(mag_i),             &
                                               respiration(mag_i),          &
                                               puptake(mag_i,1),            &
                                               nuptake(mag_i,1),            &
                                               INi,IPi)


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
        ! print *,'sss',malg,matz,primprod(mag_i)

       ENDIF


       !# Record standing biomass of all MAG into diagnostics (mmol/m2)
       IF (diag_level>0) _DIAG_VAR_S_(data%id_mag_ben) = _DIAG_VAR_S_(data%id_mag_ben) + malg
       IF (diag_level>0) _DIAG_VAR_S_(data%id_min_ben) = _DIAG_VAR_S_(data%id_min_ben) + INi
       IF (diag_level>0) _DIAG_VAR_S_(data%id_mip_ben) = _DIAG_VAR_S_(data%id_mip_ben) + IPi
       IF (diag_level>9) _DIAG_VAR_S_(data%id_c2p_ben(mag_i)) = malg/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)
       IF (diag_level>9) _DIAG_VAR_S_(data%id_c2n_ben(mag_i)) = malg/MAX(INi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_nmin)
       IF (diag_level>9) _DIAG_VAR_S_(data%id_n2p_ben(mag_i)) = INi/MAX(IPi,data%malgs(mag_i)%p0*data%malgs(mag_i)%X_pmin)

       !# Record gross productivity and respiration rates (/day)
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
          IPi = _STATE_VAR_S_(data%id_ipben(mag_i))
          flux = (sum(-puptake(mag_i,:)) - pexcretion(mag_i) - pmortality(mag_i))
          available = MAX( zero_, IPi-data%malgs(mag_i)%X_pmin*malg )
          IF ( -flux*dtlim > available  ) flux = -0.99*available/dtlim
          _FLUX_VAR_B_(data%id_ipben(mag_i)) = _FLUX_VAR_B_(data%id_ipben(mag_i)) + flux
       ENDIF

       !# Update the water column vars (O2, CO2, Nuts, OM) based on these fluxes
       !  The uptake arrays have -ve values in them, so fluxes here will be
       !  reducing these pools, but note excretion and mortality are +ve
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
                      _FLUX_VAR_(data%id_Pupttarget(c)) + puptake(mag_i,c)
          ENDDO
         ENDIF
         IF (data%do_Nuptake) THEN
          DO c = 1,data%nnup
             _FLUX_VAR_(data%id_Nupttarget(c)) =         &
                      _FLUX_VAR_(data%id_Nupttarget(c)) + nuptake(mag_i,c)
          ENDDO
         ENDIF

         IF(data%do_Cexc)  _FLUX_VAR_(data%id_Cexr) = _FLUX_VAR_(data%id_Cexr) &
                                     + cexcretion(mag_i) + exudation(mag_i)*malg
         IF(data%do_Nexc)  _FLUX_VAR_(data%id_Nexr) = _FLUX_VAR_(data%id_Nexr) + nexcretion(mag_i)
         IF(data%do_Pexc)  _FLUX_VAR_(data%id_Pexr) = _FLUX_VAR_(data%id_Pexr) + pexcretion(mag_i)
         IF(data%do_Cmort) _FLUX_VAR_(data%id_Cmor) = _FLUX_VAR_(data%id_Cmor) + cmortality(mag_i)
         IF(data%do_Nmort) _FLUX_VAR_(data%id_Nmor) = _FLUX_VAR_(data%id_Nmor) + nmortality(mag_i)
         IF(data%do_Pmort) _FLUX_VAR_(data%id_Pmor) = _FLUX_VAR_(data%id_Pmor) + pmortality(mag_i)
       ENDIF

       !# Redistribute biomass into the water column if sloughing occurs.
       IF( data%simSloughing >0 .and. data%malgs(mag_i)%slough_model >0 ) THEN
         
        IF( data%malgs(mag_i)%slough_model == 1) THEN
          ! The Coorong Ulva approach
          IF( bottom_stress>data%malgs(mag_i)%tau_0*2. ) THEN
            slough_frac = 0.3
          ELSEIF( bottom_stress>data%malgs(mag_i)%tau_0 .AND. salinity>60.) THEN
            slough_frac = 0.67
          ELSE
            slough_frac = 0.0
          ENDIF

         ELSEIF( data%malgs(mag_i)%slough_model == 2) THEN
           ! The Erie CGM approach
           CALL cladophora_slough_cgm(data,column,layer_idx,mag_i,data%slough_rate,slough_frac)
         ELSEIF( data%malgs(mag_i)%slough_model == 3) THEN
            ! The Erie GLCMv3 approach
            CALL cladophora_slough_glcmv3(data,column,layer_idx,mag_i,data%slough_rate,slough_frac)
         ELSEIF( data%malgs(mag_i)%slough_model == 4) THEN
            ! The Erie AED (hybrid) approach
            CALL cladophora_slough_aed(data,column,layer_idx,mag_i,data%slough_rate,slough_frac)
         
         ELSE
           ! No sloughing
           slough_frac = zero_
         ENDIF

         ! Now apply the rate of sloughing to the benthic biomass variable, and add to water pool
         _FLUX_VAR_B_(data%id_pben(mag_i)) = _FLUX_VAR_B_(data%id_pben(mag_i)) - (slough_frac*malg) !/DTsec
         _FLUX_VAR_(data%id_p(mag_i)) = _FLUX_VAR_(data%id_p(mag_i)) + ((slough_frac*malg)/depth) !/DTsec
         IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
          _FLUX_VAR_B_(data%id_ipben(mag_i)) = _FLUX_VAR_B_(data%id_ipben(mag_i)) - (slough_frac*_STATE_VAR_S_(data%id_ipben(mag_i))) !/DTsec
          _FLUX_VAR_(data%id_ip(mag_i)) = _FLUX_VAR_(data%id_ip(mag_i)) + ((slough_frac*IPi)/depth) !/DTsec
         ENDIF
         IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
          _FLUX_VAR_B_(data%id_inben(mag_i)) = _FLUX_VAR_B_(data%id_inben(mag_i)) - (slough_frac*INi) !/DTsec
          _FLUX_VAR_(data%id_in(mag_i)) = _FLUX_VAR_(data%id_in(mag_i)) + ((slough_frac*INi)/depth) !/DTsec
         ENDIF
         IF (diag_level>1) THEN
          !_DIAG_VAR_S_(data%id_slg_ben) = _DIAG_VAR_S_(data%id_slg_ben) - (slough_frac*malg/DTsec)*secs_per_day
!          _DIAG_VAR_S_(data%id_slg_ben) = _DIAG_VAR_S_(data%id_slg_ben) - (slough_frac/DTsec)*secs_per_day
          _DIAG_VAR_S_(data%id_slg_ben) = _DIAG_VAR_S_(data%id_slg_ben) - (slough_frac)*secs_per_day
         ENDIF


         !# Move a fraction of bottom slough biomass into the sediment - slough "burial".
         slough_burial = data%slough_burial ! rate per sec

         slough = _STATE_VAR_(data%id_p(mag_i)) ! local slough density
         _FLUX_VAR_(data%id_p(mag_i)) = _FLUX_VAR_(data%id_p(mag_i)) - (slough_burial*slough)  !/depth

         slough_in = slough * data%malgs(mag_i)%X_ncon
         IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
           slough_in = _STATE_VAR_(data%id_in(mag_i)) ! local slough in
           _FLUX_VAR_(data%id_in(mag_i)) = _FLUX_VAR_(data%id_in(mag_i)) - (slough_burial*slough_in)
         ENDIF

         slough_ip = slough * data%malgs(mag_i)%X_pcon
         IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
           slough_ip = _STATE_VAR_(data%id_ip(mag_i)) ! local slough ip
          _FLUX_VAR_(data%id_ip(mag_i)) = _FLUX_VAR_(data%id_ip(mag_i)) - (slough_burial*slough_ip)
         ENDIF

         ! Increment SWI diagnostics mmol/m2/day (-ve means downward)
         _DIAG_VAR_S_(data%id_swi_c) = _DIAG_VAR_S_(data%id_swi_c) - slough_burial*slough    * secs_per_day
         _DIAG_VAR_S_(data%id_swi_n) = _DIAG_VAR_S_(data%id_swi_n) - slough_burial*slough_in * secs_per_day
         _DIAG_VAR_S_(data%id_swi_p) = _DIAG_VAR_S_(data%id_swi_p) - slough_burial*slough_ip * secs_per_day

       ENDIF




       !# TMALG - total macroalgal biomass as gDW/m2
       !  (note this is 1 timestep behind)
       IF (diag_level>0) THEN
         malg = _STATE_VAR_S_(data%id_pben(mag_i))
         IF( data%id_p(mag_i) >0 ) THEN
          malg = malg + (_STATE_VAR_(data%id_p(mag_i))*depth)
         ENDIF
         _DIAG_VAR_(data%id_TMALG) =  _DIAG_VAR_(data%id_TMALG) &
                                    + malg * (12.*1e-3/data%malgs(mag_i)%Xcc)
       ENDIF

     ENDIF ! end benthic group/zone check

     !# Compute the (growth) HSI for chosen group
     IF( data%simMalgHSI == mag_i) THEN
       _DIAG_VAR_S_(data%id_mhsi) = min( fSal, fI, fNit, fPho ) * fT
     ENDIF


   ENDDO


   !-- Loop through selected groups again, and enact "aging" between groups
   IF(data%num_malgae>1) THEN
     DO mag_i=2,data%num_malgae
       malg = _STATE_VAR_S_(data%id_pben(mag_i-1)) ! local malg density

       IF(mag_i==2)THEN
         R_aging = 0.05/secs_per_day
         Kaging = 1000.
         malg_min = 1000
       ELSEIF(mag_i==3)THEN
         R_aging = 0.08/secs_per_day !/d
         Kaging = 3000.
         malg_min = 3000
       ENDIF

       IF(malg<malg_min)THEN
         aging_rate = zero_
       ELSE
         aging_rate = R_aging * (malg-malg_min)/(Kaging+(malg-malg_min))
       ENDIF

       aging =  aging_rate * malg
       _FLUX_VAR_B_(data%id_pben(mag_i-1)) = _FLUX_VAR_B_(data%id_pben(mag_i-1)) - aging
       _FLUX_VAR_B_(data%id_pben(mag_i)) = _FLUX_VAR_B_(data%id_pben(mag_i)) + aging
       IF (data%malgs(mag_i)%simINDynamics /= 0) THEN
         malg = _STATE_VAR_S_(data%id_inben(mag_i-1))
         aging =  aging_rate * malg
         _FLUX_VAR_B_(data%id_inben(mag_i-1)) = _FLUX_VAR_B_(data%id_inben(mag_i-1)) - aging
         _FLUX_VAR_B_(data%id_inben(mag_i)) = _FLUX_VAR_B_(data%id_inben(mag_i)) + aging
       ENDIF
       IF (data%malgs(mag_i)%simIPDynamics /= 0) THEN
         malg = _STATE_VAR_S_(data%id_ipben(mag_i-1))
         aging =  aging_rate * malg
         _FLUX_VAR_B_(data%id_ipben(mag_i-1)) = _FLUX_VAR_B_(data%id_ipben(mag_i-1)) - aging
         _FLUX_VAR_B_(data%id_ipben(mag_i)) = _FLUX_VAR_B_(data%id_ipben(mag_i)) + aging
       ENDIF
     ENDDO
   ENDIF


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

      IF (data%id_p(mag_i)>0) CYCLE

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

     ! Self-shading with explicit contribution from macroalgae biomass
     extinction = extinction !+ _DIAG_VAR_(data%id_dEXTC)

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

   IF (data%id_p(mag_i)>0) THEN

     ! Retrieve current group of water macroalgae (mmol C/m3)
     malg = _STATE_VAR_(data%id_p(mag_i))

     ! Drag addition depending on amount of carbon in water volume
     drag = drag + KCD_wc * malg
   ENDIF

   ! Process benthic / attached macroalgae
   IF ( data%malgs(mag_i)%growth_form >0 ) THEN
     dz   = MAX(_STATE_VAR_(data%id_dz),0.05)       ! cell depth
     malg = _STATE_VAR_S_(data%id_pben(mag_i))/dz  ! attached macroalgae
     drag = drag + (KCD_ben * (malg/dz) )
   ENDIF

 ENDDO

END SUBROUTINE aed_bio_drag_macroalgae
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE cladophora_calculate_cgm(data,column,layer_idx,cgm,pf,rf, &
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
   AED_REAL :: hour
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

   
   macroHgt = (1./100.) * 1.467 * ( malg )**0.425   ! Formulation from Higgin's CGM model
   !macroHgt = (1./100.) * 1.9415 * ( malg )**0.4138 ! Formulation from S Malkin; GLCMv3 model

   macroExt = 7.840 * ( malg )**0.240   ! Higgins et al 2006 Fig 6

   macroPAR_Top = MAX( MIN( par*exp(-extc*( MAX(dz-macroHgt,zero_))),Io), zero_)
   macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)

   IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC)  = macroExt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dHGT) = macroHgt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dPAR) = macroPAR_Top
   IF(diag_level>9) _DIAG_VAR_S_(data%id_par_bot) = macroPAR_Bot

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
   IF ( data%simCGM >0 ) THEN
     _DIAG_VAR_S_(data%id_tem_avg) = _DIAG_VAR_S_(data%id_tem_avg) &
                                * (1-(DTday/TempAvgTime)) + temp*(DTday/TempAvgTime)
     AvgTemp = _DIAG_VAR_S_(data%id_tem_avg)
   ELSE
     AvgTemp = 0.
   ENDIF

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

      IF ( data%simCGM >0 ) THEN
        _DIAG_VAR_S_(data%id_par_avg) = _DIAG_VAR_S_(data%id_par_avg) &
                          * (1-(DTday/LgtAvgTime)) + macroPAR_Top*(DTday/LgtAvgTime)
        AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83   ! AvgLight in uE for CGM
      ELSE
        AvgLight = 0.
      ENDIF

      !print *,'AvgLight',macroPAR_Top, extc, AvgLight

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

      pplt = pplt*5.43
      
      IF (data%malgs(cgm)%lightModel == 3) pplt = PhotoRate(lght,temp)

      pf = (data%malgs(cgm)%R_growth * pplt) * sf * lf

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

     !print *,'cgm',pplt, lf, AvgLight

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
    pu = data%malgs(cgm)%X_pmin 
    IF( malg>zero_ ) pu = mag_ip/malg

    Kq = 0.0028 * (12e3/31e3) ! 0.07% = 0.0028gP/gC & 12/31 is mol wgt conversion

    frp = _STATE_VAR_(data%id_Pupttarget(1))

    !-- IPmax = Kq; IPmin = Qo; KP = Km; R_puptake = pmax; tau = tf 
    pu = data%malgs(cgm)%R_puptake * tf                                                &
                                   * (frp/(frp + data%malgs(cgm)%K_P))                 &
                                   * (Kq /(Kq + MAX(pu - data%malgs(cgm)%X_pmin,zero_)))

    !---------------------------------------------------------------------------
    !-- Get the INTERNAL NITROGEN stores for the macroalgae groups.
    !   Recall that internal nutrient is in g and must be converted to g N/g C
    !   by division by macroalgae biomass for the nitrogen limitation

    !-- Get the internal nitrogen ratio.
    nu = data%malgs(cgm)%X_nmin
    IF(malg>zero_) nu = mag_in/malg

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

    _DIAG_VAR_S_(data%id_gpp_bot) = pf_MB

!    IF(malg > data%malgs(cgm)%p0) THEN
!      !-- SloughTrigger = SloughTrigger + pf_MB * DTday
!      !_FLUX_VAR_B_(data%id_slough_trig) = _FLUX_VAR_B_(data%id_slough_trig) + pf_MB
!      _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) + pf_MB*DTday*secs_per_day
!    ENDIF

    ! Daily reset of bottom filament checker
    hour = mod(_STATE_VAR_S_(data%id_yearday), 1.0) 

    IF(hour < 0.01) _DIAG_VAR_S_(data%id_slough_trig) = zero_ 

    ! Increment bottom filament checker if pf_MB <0 
    IF(malg > data%malgs(cgm)%p0) THEN
      !-- SloughTrigger = SloughTrigger + pf_MB * DTday
      IF( pf_MB < zero_ ) & 
        _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) + (DTday)
    ENDIF

    !print *,'hr',hour, pf_MB,_DIAG_VAR_S_(data%id_slough_trig)
    sf = one_

END SUBROUTINE cladophora_calculate_cgm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE cladophora_slough_cgm(data,column,layer_idx,cgm,slough_rate,hf)
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
   IF ( data%simCGM >0 ) THEN
      AvgStress = _DIAG_VAR_S_(data%id_tau_avg) &
                * (1-(DTday/StrAvgTime)) + bottom_stress *(DTday/StrAvgTime)
      _DIAG_VAR_S_(data%id_tau_avg) = AvgStress
   ELSE
      AvgStress = 0.
   ENDIF

   !-- Retrieve current (local) state variable values.
   malg = _STATE_VAR_S_(data%id_pben(cgm))

   IF ( data%simCGM >0 ) THEN
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
   ELSE
     slough_trigger = 0.
   ENDIF

   !-- Sloughing of healthy filaments if the shear stress is high enough
   IF(malg>data%malgs(cgm)%p0 .AND. slough_trigger<-0.001 .AND. AvgStress>data%malgs(cgm)%tau_0)THEN

      !-------------------------------------------------------------------------
      ! Depth (light) based carrying capacity amount, computed empirically with
      ! *0.25 to get from g DM to g C, /12 to get mol C, and 1e3 to get to mmol
      IF ( data%simCGM >0 ) THEN
        AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83 ! AvgLight in uE for CGM
      ELSE
        AvgLight = 0.
      ENDIF
      X_maxp = ( 1.18 * AvgLight - 58.7 ) * data%malgs(cgm)%Xcc * 1e3 / 12.

      ! Daily slough rate
      ! Lss = Lmax * (T/Tcrit) * (X/Xmax)
      ! hf = 0.242 * exp(-0.3187 * DEPTH) * AvgStress / PCm  * malg / X_maxp
      hf = (AvgStress - data%malgs(cgm)%tau_0) * MIN( malg/MAX(X_maxp,data%malgs(cgm)%p0),one_ )

      hf = slough_rate * hf   ! daily biomass fraction sloughed per sec

      IF (hf*DTsec>0.95) hf = 0.95/DTsec    ! no more than 95% slough in one interval

      !-- Update the malg and slough variables with following the slough event
      IF ( data%simCGM >0 ) &
        _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) * hf*DTsec/malg ! 0.5
   ENDIF

END SUBROUTINE cladophora_slough_cgm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!###############################################################################
SUBROUTINE cladophora_calculate_glcmv3(data,column,layer_idx,cgm,              & 
                                                           u_canopy,r_canopy,  &
                                                           pu,nu,mag_in,mag_ip )
!-------------------------------------------------------------------------------
! Calculate benthic cladophora productivity, based on the GLCMv3.
!   Units per surface area (not volume!)
!   and switch between mmol C/m2, and the original model which is in g DM/m2.
!-------------------------------------------------------------------------------
!ARGUMENTS
CLASS (aed_macroalgae_data_t),INTENT(in) :: data
TYPE (aed_column_t),INTENT(inout) :: column(:)
INTEGER,INTENT(in) :: layer_idx,cgm
AED_REAL, INTENT(inout) :: u_canopy,r_canopy,pu,nu,mag_in,mag_ip
!
!LOCALS
!
AED_REAL :: malg
AED_REAL :: matz, extc, dz, par, Io, temp, salinity, frp
AED_REAL :: AvgTemp, AvgLight
AED_REAL :: lght, pplt, prlt, pf_MB

AED_REAL :: X_maxp, Kq
AED_REAL :: macroPAR_Top, macroPAR_Bot
AED_REAL :: sf,tf,lf
AED_REAL :: fRIT, Rb, FuIT

AED_REAL :: Q, X, S, Iz, Imat, X_layer, X_left, X_calc, S_calc, Q_layer, fQ, Qmin, rho, U, R, Rmax, pu_canopy, hour, unet_canopy, umax, unet, pf, rf
INTEGER :: layer, total_layers
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

   macroHgt = (1./100.) * 1.9415 * ( malg )**0.4138 ! Formulation from S Malkin; GLCMv3 model
   macroExt = 68. ! 7.840 * ( malg )**0.240   ! Higgins et al 2006 Fig 6
   macroPAR_Top = MAX( MIN( par*exp(-extc*( MAX(dz-macroHgt,zero_))),Io), zero_)
   macroPAR_Bot = MAX( macroPAR_Top * exp(-(extc+macroExt)*macroHgt), zero_)

   IF(diag_level>9) _DIAG_VAR_(data%id_dEXTC)  = macroExt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dHGT) = macroHgt
   IF(diag_level>9) _DIAG_VAR_S_(data%id_dPAR) = macroPAR_Top
   IF(diag_level>9) _DIAG_VAR_S_(data%id_par_bot) = macroPAR_Bot

   !----------------------------------------------------------------------------
   !-- MACROALGAL C, N, P
   !## Get the CGM (cladophora) group concentrations
   malg = _STATE_VAR_S_(data%id_pben(cgm))      ! cladophora biomass (mmolC/m2)
   mag_in = data%malgs(cgm)%X_ncon * malg       ! cladophora biomass (mmolN/m2)
   IF(data%malgs(cgm)%simINDynamics/=0) mag_in = _STATE_VAR_S_(data%id_inben(cgm)) 
   mag_ip = data%malgs(cgm)%X_pcon * malg       ! cladophora biomass (mmolP/m2)
   IF(data%malgs(cgm)%simIPDynamics/=0) mag_ip = _STATE_VAR_S_(data%id_ipben(cgm)) 

   !----------------------------------------------------------------------------
   !-- Update moving average for daily temp
   !   (averaged over the past 1 day)
   IF ( data%simCGM >0 ) THEN
     _DIAG_VAR_S_(data%id_tem_avg) = _DIAG_VAR_S_(data%id_tem_avg) &
                * (1-(DTday/TempAvgTime)) + temp*(DTday/TempAvgTime)
     AvgTemp = _DIAG_VAR_S_(data%id_tem_avg)
   ELSE
     AvgTemp = 0.
   ENDIF

   !----------------------------------------------------------------------------
   ! Calculate current Q and Q limitation factor
   Q = mag_ip/malg * 100;
   X = malg

   ! Light attenuation through water depth - bed height, top of the mat
   Iz = macroPAR_Top

   ! Mat depth / canopy height (in cm) (biomass density of each layer in canopy)
   X_layer = X / MAX( macroHgt*100. ,1.)

   ! Count the number of 1cm layers in the algae mat (canopy height/mat thickness)
   ! Start with 1 layer
   Total_layers =  1
   DO WHILE ( X >= X_layer )
     Total_layers = Total_layers + 1
     X = X - X_layer
   ENDDO
   ! Biomass left after filling each layer to its capacity
   X_left = X

   unet_canopy = zero_
   ! First, fill canopy layers with biomass and P mass
   DO layer = 1,Total_layers
     IF (layer == 1) THEN 
       ! Leftover X goes into the top layer
       X_calc = X_left
     else 
       ! Fill each layer with constant amount
       X_calc = X_layer
     endif

     ! Fill each layer with P mass
     S_calc = Q * X_calc/100.

     ! Calculate Q limitation factor
     Q_layer = S_calc/X_calc * 100
     fQ = 1 - (Qmin/Q_layer)

     ! P uptake: Michaelis-Menten, linear at low Q,
     ! retaining high km = 125 ug/L
     ! P uptake rate (%P/hr)(%P/s)
     rho = (0.012 * (Q_layer**(-2.3))) * (frp / (frp+ 125))   /secs_per_day   ! %P/s

     ! Light extinction through the algae mat
     ! Substract 1 cm for light at top of layer and convert
     ! into meters
     Imat = Iz * exp( -macroExt * (layer - 1) / 100.)

     ! Gross growth rate depending on light,
     ! temperature and store P (1/hr)(1/s)
     fuIT = PhotoRate(Imat, AvgTemp)
     u = umax * fQ * fuIT / secs_per_day;

     ! Total respiration = RIT + Rbasal (1/hr)(1/s)
     fRIT = RespRate(Imat, AvgTemp);
     Rb = BasalResp(AvgTemp);
     r = (Rmax * fRIT + Rb) / secs_per_day;

     ! Interval u and r in the current layer
     u_canopy = u_canopy + u * X_calc 
     r_canopy = r_canopy + r * X_calc 
     pu_canopy = pu_canopy + (rho / 100 * X_calc - r * S_calc) 

enddo 

u_canopy = u_canopy / X
r_canopy = r_canopy / X
pu_canopy = pu_canopy / S


! Daily reset of bottom filament checker
hour = mod(_STATE_VAR_S_(data%id_yearday), 1.0) 

IF(hour < 0.01) _DIAG_VAR_S_(data%id_slough_trig) = zero_ 

! Increment bottom filament checker, if unet <0 
IF(malg > data%malgs(cgm)%p0) THEN
  !-- SloughTrigger = SloughTrigger * DTday
  unet = u_canopy-r_canopy
  IF( (unet) < zero_ ) & 
    _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) + (DTday)
ENDIF



!-------------------------------------------------------------------------
!-- Calculate the nutrient limitation (phosphorus & nitrogen) and
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

IF ( data%simCGM >0 ) THEN
  _DIAG_VAR_S_(data%id_par_avg) = _DIAG_VAR_S_(data%id_par_avg) &
    * (1-(DTday/LgtAvgTime)) + macroPAR_Top*(DTday/LgtAvgTime)
  AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83   ! AvgLight in uE for CGM
ELSE
  AvgLight = 0.
ENDIF

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

temp = AvgTemp
lght = AvgLight  !MIN(600.0,AvgLight)/1235.0      ! Capped at 600: Higgins et al 2006
pplt = PhotoRate(lght,temp)
pf = (data%malgs(cgm)%R_growth * pplt) * sf * lf

!-------------------------------------------------------------------------
!-- Now light and temperature function for daytime respiration

! Total respiration = RIT + Rbasal (1/hr)
fRIT = RespRate(lght,temp);
Rb = BasalResp(temp);
rf = (data%malgs(cgm)%R_Resp * fRIT + Rb) ;


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
pu = data%malgs(cgm)%X_pmin 
IF( malg>zero_ ) pu = mag_ip/malg

Kq = 0.0028 * (12e3/31e3) ! 0.07% = 0.0028gP/gC & 12/31 is mol wgt conversion

frp = _STATE_VAR_(data%id_Pupttarget(1))

!-- IPmax = Kq; IPmin = Qo; KP = Km; R_puptake = pmax; tau = tf 
pu = data%malgs(cgm)%R_puptake * tf                                                &
* (frp/(frp + data%malgs(cgm)%K_P))                 &
* (Kq /(Kq + MAX(pu - data%malgs(cgm)%X_pmin,zero_)))

!---------------------------------------------------------------------------
!-- Get the INTERNAL NITROGEN stores for the macroalgae groups.
!   Recall that internal nutrient is in g and must be converted to g N/g C
!   by division by macroalgae biomass for the nitrogen limitation

!-- Get the internal nitrogen ratio.
nu = data%malgs(cgm)%X_nmin
IF(malg>zero_) nu = mag_in/malg

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


!-------------------------------------------------------------------------
!-- Now light and temperature function for photosynthesis

temp = AvgTemp
lght = macroPAR_Bot  !MIN(600.0,AvgLight)/1235.0      ! Capped at 600: Higgins et al 2006
pplt = PhotoRate(lght,temp)
pf = (data%malgs(cgm)%R_growth * pplt) * sf * lf

!-------------------------------------------------------------------------
!-- Now light and temperature function for daytime respiration

! Total respiration = RIT + Rbasal (1/hr)
fRIT = RespRate(lght,temp);
Rb = BasalResp(temp);
rf = (data%malgs(cgm)%R_resp * fRIT + Rb) ;

pf_MB = pf-rf

_DIAG_VAR_S_(data%id_slough_trig) = pf_MB

IF(malg > data%malgs(cgm)%p0) THEN
  !-- SloughTrigger = SloughTrigger + pf_MB * DTday
  IF( pf_MB < 0 ) & 
    _DIAG_VAR_S_(data%id_slough_trig) = _DIAG_VAR_S_(data%id_slough_trig) + (DTday * secs_per_day)
ENDIF

sf = one_

END SUBROUTINE cladophora_calculate_glcmv3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE cladophora_slough_glcmv3(data,column,layer_idx,cgm,slough_rate,L)
  !-------------------------------------------------------------------------------
  ! Compute the fraction of cladophora biomass that is detatched due to
  !  sloughing if the stress is enough
  !-------------------------------------------------------------------------------
  !ARGUMENTS
     CLASS (aed_macroalgae_data_t),INTENT(in) :: data
     TYPE (aed_column_t),INTENT(inout) :: column(:)
     INTEGER,INTENT(in) :: layer_idx,cgm
     AED_REAL,INTENT(in) :: slough_rate
     AED_REAL,INTENT(out) :: L
  !
  !LOCALS
     ! Environment
     AED_REAL :: yearday, depth
     ! State
     AED_REAL :: tstart, dark_days, slough_trigger
     ! Parameters
     AED_REAL :: Lmax, f1, f2, a, b, tend
     AED_REAL, PARAMETER :: tdur = 30
  !-------------------------------------------------------------------------------
  !BEGIN
  
     yearday= _STATE_VAR_S_(data%id_yearday) 
     depth  = _STATE_VAR_S_(data%id_depth) 

     tstart = _DIAG_VAR_S_(data%id_slough_tsta)
     dark_days = _DIAG_VAR_S_(data%id_slough_days)
     slough_trigger = _DIAG_VAR_S_(data%id_slough_trig) ! Set in calculate_glcmv3/cgm

     L = zero_

     !-- Slough off weakened filaments (those with cumulative respiration excess)
     Lmax = slough_rate  ! 0.08/d 
  
     ! Sloughing from physical factors
     f1 = 0.4635 * exp(-0.3054 * depth)  + 0.5365;
 
     IF( slough_trigger > 0.999 ) THEN
      ! The whole day has been dark, with cumulative respiration deficit 
      dark_days = dark_days + 1.
     ELSE IF ( mod(yearday, 1.0) > (0.999-DTday) ) THEN
      ! Its the end of the day, but slough trigger<1 : the darkness is over!
      dark_days = zero_
      tstart = zero_
    ENDIF 
  
     !-- Check if resp or growth phase; if + clip to 0, else count days
     IF(dark_days > 0.998 .and. dark_days < 1.002 ) THEN
       tstart = yearday;
     ELSEIF(dark_days > 1.) THEN
       tend = tstart + tdur
       a = 12/tdur
       b = 6 * (tstart + tend) / (tstart - tend)
       f2 = 1/ (1 + exp(-(a * yearday + b)))
       L = Lmax * f1 * f2
     ELSEIF(dark_days < 0.99)  THEN
       f2 = zero_
       L = zero_
     ENDIF

     _DIAG_VAR_S_(data%id_slough_days) = dark_days
     _DIAG_VAR_S_(data%id_slough_tsta) = tstart
  
END SUBROUTINE cladophora_slough_glcmv3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!###############################################################################
SUBROUTINE cladophora_slough_aed(data,column,layer_idx,cgm,slough_rate,L)
  !-------------------------------------------------------------------------------
  ! Compute the fraction of cladophora biomass that is detatched due to
  !  sloughing, if the stress is enough and depending on filament health
  !-------------------------------------------------------------------------------
  !ARGUMENTS
     CLASS (aed_macroalgae_data_t),INTENT(in) :: data
     TYPE (aed_column_t),INTENT(inout) :: column(:)
     INTEGER,INTENT(in) :: layer_idx,cgm
     AED_REAL,INTENT(in) :: slough_rate
     AED_REAL,INTENT(out) :: L
  !
  !LOCALS
     ! Environment
     AED_REAL :: yearday, depth
     ! State
     AED_REAL :: malg, tstart, dark_days, slough_trigger
     ! Parameters
     AED_REAL :: Lmax, f1, f2, a, b, tend, tau_crit, X_maxp, bottom_stress, AvgStress
     AED_REAL, PARAMETER :: tdur = 30
     AED_REAL, PARAMETER :: tau_ref = 0.01
  !-------------------------------------------------------------------------------
  !BEGIN
    
    depth  = _STATE_VAR_S_(data%id_depth) 
    tstart = _DIAG_VAR_S_(data%id_slough_tsta)
    yearday = _STATE_VAR_S_(data%id_yearday) 
    dark_days = _DIAG_VAR_S_(data%id_slough_days)
    slough_trigger = _DIAG_VAR_S_(data%id_slough_trig)  ! Set in calculate_glcmv3/cgm

    !----------------------------------------------------------------------------
   !-- Retrieve current bottom shear stress condition within the cell
   bottom_stress = MIN( _STATE_VAR_S_(data%id_taub), 10. )

   !-- Update moving average for stress (averaged over the past 2 hrs)
   IF ( data%simCGM >0 ) THEN
      AvgStress = _DIAG_VAR_S_(data%id_tau_avg) &
                * (1-(DTday/StrAvgTime)) + bottom_stress *(DTday/StrAvgTime)
      _DIAG_VAR_S_(data%id_tau_avg) = AvgStress
   ELSE
      AvgStress = 0.
   ENDIF

     
    !-- Initialise slough rate estimate to zero
    L = zero_

    !-- Define maximum sloughing 
    Lmax = slough_rate  ! 0.08/d 
  
    !-------------------------------------------------------------------------
    !-- Physical controls on sloughing depends on shear and filament health
    tau_crit =  MAX( data%malgs(cgm)%tau_0 * (one_ - MIN(dark_days/tdur,one_)), 0.005)

    !-- Depth/ light based carrying capacity amount, computed empirically with
    ! *0.25 to get from g DM to g C, /12 to get mol C, and 1e3 to get to mmol
    
    ! AvgLight = _DIAG_VAR_S_(data%id_par_avg) * 4.83 ! AvgLight in uE for CGM
    ! X_maxp = ( 1.18 * AvgLight - 58.7 )  * data%malgs(cgm)%Xcc * 1e3 / 12.
   ! Retrieve (local) cladophora in gDM/m2 (Xcc=1/0.25 to convert to DM)

    malg = _STATE_VAR_S_(data%id_pben(cgm)) * (12. / 1e3) / data%malgs(cgm)%Xcc

    X_maxp = 1230 * exp(-0.55 * depth) !* data%malgs(cgm)%Xcc * 1e3 / 12.

    f1 = 0.3
    IF(AvgStress>tau_crit) THEN
      f1 = MIN( f1 + (AvgStress - tau_crit)/tau_ref * MIN( malg/MAX(X_maxp,data%malgs(cgm)%p0),one_ ), 5.)
    ENDIF

    !-------------------------------------------------------------------------
    !-- Physiological controls on filament health and susceptibility to slough
    IF( slough_trigger > 0.999 ) THEN
      ! The whole day has been dark, with cumulative respiration deficit 
      dark_days = dark_days + 1.
     ELSE IF ( mod(yearday, 1.0) > (0.999-DTday)  ) THEN
      ! Its the end of the day, but slough trigger<1 : the darkness is over!
      dark_days = zero_
      tstart = zero_
    ENDIF 
  
     !-- Check if resp or growth phase; if + clip to 0, else count days
     IF(dark_days > 0.998 .and. dark_days < 1.002 ) THEN
       tstart = yearday;
     ELSEIF(dark_days > 1.) THEN
       tend = tstart + tdur
       a = 12/tdur
       b = 6 * (tstart + tend) / (tstart - tend)
       f2 = 1/ (1 + exp(-(a * yearday + b)))
       L = Lmax * f1 * f2
     ELSEIF(dark_days < 0.99)  THEN
       f2 = zero_
       L = zero_
     ENDIF

     _DIAG_VAR_S_(data%id_slough_days) = dark_days
     _DIAG_VAR_S_(data%id_slough_tsta) = tstart
  
END SUBROUTINE cladophora_slough_aed
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION BasalResp(temp) RESULT(Rb)
  !-------------------------------------------------------------------------------
  !% Kuzynski et al. (2020)
  !-------------------------------------------------------------------------------
  !ARGUMENTS
  AED_REAL,INTENT(in) :: temp
  AED_REAL :: Rb

  AED_REAL, PARAMETER :: Rb_a = 0.07
  AED_REAL, PARAMETER :: Rb_Theta = 1.04

  Rb = Rb_a * Rb_Theta ** (temp - 20.)

END FUNCTION BasalResp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
FUNCTION RespRate(lgt,temp) RESULT(fRIT)
  !-------------------------------------------------------------------------------
  !% Kuzynski et al. (2020)
  !-------------------------------------------------------------------------------
  !ARGUMENTS
  AED_REAL,INTENT(in) :: lgt,temp
  AED_REAL :: fRIT

  AED_REAL :: Rmax, R_alpha
  AED_REAL, PARAMETER :: R_a = 0.1
  AED_REAL, PARAMETER :: R_Theta = 1.03
  AED_REAL, PARAMETER :: R_T = 20
  AED_REAL, PARAMETER :: R_alpha_a = 0.00168
  AED_REAL, PARAMETER :: R_alpha_b = 2.5

  AED_REAL, PARAMETER :: Rmax_modeled = 0.187

  Rmax = R_a * R_Theta ** (temp - R_T);
  R_alpha = R_alpha_a * (temp / (temp + R_alpha_b));

  fRIT = (Rmax * (1 - exp(-R_alpha * lgt / Rmax))) / Rmax_modeled;

END FUNCTION RespRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
FUNCTION PhotoRate(lgt,temp) RESULT(fuIT)
!-------------------------------------------------------------------------------
!% Kuzynski et al. (2020)
!-------------------------------------------------------------------------------
!ARGUMENTS
  AED_REAL,INTENT(in) :: lgt,temp
  AED_REAL :: fuIT
!
!LOCALS
  AED_REAL :: Pmax,u_alpha,u_beta

  AED_REAL, PARAMETER :: Pmax_a    = 0.34;
  AED_REAL, PARAMETER :: Pmax_b    = 3.1;
  AED_REAL, PARAMETER :: u_alpha_a = 0.55;
  AED_REAL, PARAMETER :: u_alpha_b = 0.001;
  AED_REAL, PARAMETER :: u_alpha_c = 0.048;
  AED_REAL, PARAMETER :: u_beta_a  = 1E-17;
  AED_REAL, PARAMETER :: u_beta_b  = 9.3;

  AED_REAL, PARAMETER :: umax_modeled = 0.28778;

!-------------------------------------------------------------------------------
!BEGIN

  Pmax = Pmax_a * (temp / (temp + Pmax_b));
  u_alpha = u_alpha_a * (1 - exp(-u_alpha_b * temp / u_alpha_a))  &
          * exp(-u_alpha_c * temp / u_alpha_a);
  u_beta = u_beta_a * temp ** u_beta_b;

  fuIT = (Pmax * (1 - exp(-u_alpha * lgt / Pmax)) * exp(-u_beta * lgt / Pmax)) / umax_modeled;

END FUNCTION PhotoRate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE aed_macroalgae
