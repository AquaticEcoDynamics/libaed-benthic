!###############################################################################
!#                                                                             #
!# aed_macrophyte.F90                                                          #
!#                                                                             #
!#  Developed by :                                                             #
!#      AquaticEcoDynamics (AED) Group                                         #
!#      School of Agriculture and Environment                                  #
!#      The University of Western Australia                                    #
!#                                                                             #
!#      http://aquatic.science.uwa.edu.au/                                     #
!#                                                                             #
!#  Copyright 2015 - 2025 -  The University of Western Australia               #
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
!# Created July 2024                                                           #
!#                                                                             #
!###############################################################################

#include "aed.h"


MODULE aed_macrophyte
!===============================================================================
!  aed_macrophyte --- multi-group macrophyte (/seagrass) model, with options
!===============================================================================
   USE aed_core
   USE aed_util
  !USE aed_maths
   USE aed_bio_utils
  !USE kobe_sav_posture

   IMPLICIT NONE

   PRIVATE   ! By default make everything private
!
   PUBLIC aed_macrophyte_data_t

!-------------------------------------------------------------------------------
!  MODULE SPECIES PARAMETER TYPE DEFINITION

!  %% NAMELIST    %% type :  macrophyte_params_t
   TYPE :: macrophyte_params_t

       CHARACTER(64) :: m_name
       INTEGER       :: growth_form
       INTEGER       :: zone_lock
       AED_REAL      :: m_initial
       AED_REAL      :: m_initial_minb
       AED_REAL      :: m_initial_mind
       AED_REAL      :: m0

       ! Canopy/physical characteristics
       INTEGER       :: drag_model
       AED_REAL      :: K_CD
       AED_REAL      :: shoot_diameter   ! shoot/leaf diameter or blade width
       AED_REAL      :: shoot_length     ! shoot/leaf typical length
       AED_REAL      :: tissue_density   ! shoot/leaf density (kg/m3 of tissue)
       AED_REAL      :: elastic_modulus  ! shoot/leaf elasticity
       INTEGER       :: shading_model
       AED_REAL      :: KeMAC
       AED_REAL      :: kA               ! shelf shading effect
       AED_REAL      :: k_omega          ! Omega_MAC?
       AED_REAL      :: sine_blade       ! sine blade shape
       AED_REAL      :: X_ldw            ! leaf length to biomass
       AED_REAL      :: K_ldw            ! half-biomass of leaf to biomass curve
       AED_REAL      :: X_sdw            ! shoot to biomass ratio

       ! Growth rate options/parameters
       AED_REAL      :: R_growth
       INTEGER       :: temp_method
       AED_REAL      :: theta_growth
       AED_REAL      :: T_std
       AED_REAL      :: T_opt
       AED_REAL      :: T_max
       AED_REAL      :: kTn, aTn, bTn
       INTEGER       :: light_model         ! Light_model
       AED_REAL      :: I_K
       AED_REAL      :: I_S
       AED_REAL      :: f_pr
       AED_REAL      :: K_epi

       ! Loss rate options/parameters
       AED_REAL      :: R_resp             ! maximum AG respiration rate, /day
       AED_REAL      :: theta_resp
       AED_REAL      :: E_comp             ! compensation scalar PAR irradiance, mol photon/m2/d;
       AED_REAL      :: R_mort             ! maximum AG mortality rate, /day
       INTEGER       :: sal_method
       AED_REAL      :: S_bep
       AED_REAL      :: S_maxsp
       AED_REAL      :: S_opt

       ! Belowground options/parameters
       AED_REAL      :: m_initial_bg
       AED_REAL      :: f_bg               ! f_below
       AED_REAL      :: tau_tran           ! AG/BG translocation rate, /day
       AED_REAL      :: R_resp_bg          ! maximum BG respiration rate, /day
       AED_REAL      :: R_mort_bg          ! maximum BG mortality rate, /day

       ! Stoichiometry and limitation parameters
       AED_REAL      :: X_cchl
       AED_REAL      :: X_cdw
       AED_REAL      :: K_N
       AED_REAL      :: X_ncon
       AED_REAL      :: K_P
       AED_REAL      :: X_pcon
       AED_REAL      :: K_S
       AED_REAL      :: K_NPP

       ! Fruiting/recruitment options/parameters
       INTEGER       :: fruit_model         ! Fruiting
       AED_REAL      :: m_initial_f
       AED_REAL      :: tau_tran_fruit  ! translocation rate from leaves to fruit, /day
       AED_REAL      :: r_release         ! fruit releasing rate, /day (assume all fruit released in 10 days)
       AED_REAL      :: f_seed          ! seed biomass/(AG+seed biomass) fraction
       AED_REAL      :: t_start_g       ! calenda day to grow, days of the year
       AED_REAL      :: t_dur_g          ! duration to reach maximum growth, days
       AED_REAL      :: t_start_r       ! calenda day to release fruits, days of the year
       AED_REAL      :: t_dur_r           ! duration to reach maximum release, days

   END TYPE
!  %% END NAMELIST    %% type :  macrophyte_params_t
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  MODULE DATA STRUCTURE DEFINITION

   TYPE,extends(aed_model_data_t) :: aed_macrophyte_data_t
      !# Variable identifiers (local and linked)
      INTEGER,ALLOCATABLE :: id_mphya(:),id_mphyb(:),id_mphyf(:)
      INTEGER :: id_epi, id_epib, id_epig, id_epir
      INTEGER :: id_par, id_tem, id_sal, id_extc, id_I_0, id_yearday, id_vel
      INTEGER :: id_depth, id_dz
      INTEGER :: id_atem,id_theta
      INTEGER :: id_mac_ag, id_mac_bg, id_root_d, id_root_o, id_mac_tr
      INTEGER :: id_d_par, id_gpp, id_p2r, id_mac, id_sed_zone
      INTEGER :: id_shdiam(MAX_PHYTO_TYPES)
      INTEGER :: id_shdens(MAX_PHYTO_TYPES)
      INTEGER :: id_shhght(MAX_PHYTO_TYPES)
      INTEGER :: id_shabve(MAX_PHYTO_TYPES)
      INTEGER :: id_shbvol(MAX_PHYTO_TYPES)
      INTEGER :: id_canopy_stem_diam, id_canopy_stem_dens, id_canopy_height, id_canopy_biovolume
      INTEGER :: id_canopy_sh_dens, id_canopy_sh_diam, id_canopy_frarea
      INTEGER :: id_canopy_velocity, id_canopy_drag, id_canopy_blockage, id_canopy_lai
      INTEGER :: id_parc, id_aeff, id_kemac, id_mac_fr, id_mac_ft, id_mac_gt
      INTEGER :: id_oxy, id_dic, id_nox, id_nh4, id_po4

      !# Macrophyte simulation settings
      INTEGER  :: num_mphy
      INTEGER  :: n_zones
      LOGICAL  :: simMacFeedback, simStaticBiomass
      LOGICAL  :: simEpiphytes, simFruiting, simMacDrag
      INTEGER  :: drag_model, mac_pattern, mac_mode, mac_initial
      AED_REAL :: water_nutrient_frac
      AED_REAL,ALLOCATABLE :: active_zones(:)

      !# Macrophyte module species parameters
      TYPE(macrophyte_params_t),DIMENSION(:),ALLOCATABLE :: mpars

      !# Epiphyte sub-module option & parameters
      INTEGER  :: epi_model
      AED_REAL :: R_epig, R_epir, R_epib, I_Kepi, epi_max
      AED_REAL :: theta_epi_growth, theta_epi_resp
      AED_REAL :: epi_initial, epi_Xpc, epi_Xnc, epi_K_N, epi_K_P
      LOGICAL  :: do_Cuptake, do_Nuptake, do_Puptake


     CONTAINS      ! Selected AED methods to activate:
         PROCEDURE :: define             => aed_define_macrophyte
         PROCEDURE :: initialize_benthic => aed_initialize_benthic_macrophyte
         PROCEDURE :: calculate_column   => aed_calculate_column_macrophyte
         PROCEDURE :: calculate_benthic  => aed_calculate_benthic_macrophyte
         PROCEDURE :: light_extinction   => aed_light_extinction_macrophyte
       ! PROCEDURE :: bio_drag           => aed_bio_drag_macrophyte

   END TYPE
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  MODULE WORKING VARIABLES

   INTEGER :: i
   AED_REAL, parameter :: h = 6.626e-34        ! Planck constant, Js
   AED_REAL, parameter :: c = 2.998e8          ! light speed, m/s
   AED_REAL, parameter :: Av= 6.02e23          ! Avagadro number, /mol, convert light capture to units of mol photon/m2/s
   AED_REAL, parameter :: WL_range(501) = (/ (I, I = 300, 800) /) ![300 800] ! wavelength range
   AED_REAL, parameter :: WL_mean = 550        ! mean wavelength, nm
   AED_REAL, ALLOCATABLE :: WLint(:), landaiint(:), ALlint(:)
   AED_REAL :: landaiSUM, ALl

   ! clear-sky wave lengths and irradiance at a particular wave-band W/m2/nm, from EMS model
   AED_REAL, parameter :: wavei(150) = (/ &
       140.00, 150.00, 160.00, 170.00, &
       180.00, 190.00, 200.00, 205.00, 210.00, 215.00, 220.00, &
       225.00, 230.00, 235.00, 240.00, 245.00, 250.00, 255.00, &
       260.00, 265.00, 270.00, 275.00, 280.00, 285.00, 290.00, &
       295.00, 300.00, 305.00, 310.00, 315.00, 320.00, 325.00, &
       330.00, 335.00, 340.00, 345.00, 350.00, 355.00, 360.00, &
       365.00, 370.00, 375.00, 380.00, 385.00, 390.00, 395.00, &
       400.00, 405.00, 410.00, 415.00, 420.00, 425.00, 430.00, &
       435.00, 440.00, 445.00, 450.00, 455.00, 460.00, 465.00, &
       470.00, 475.00, 480.00, 485.00, 490.00, 495.00, 500.00, &
       505.00, 510.00, 515.00, 520.00, 525.00, 530.00, 535.00, &
       540.00, 545.00, 550.00, 555.00, 560.00, 565.00, 570.00, &
       575.00, 580.00, 585.00, 590.00, 595.00, 600.00, 610.00, &
       620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, &
       690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, &
       800.00, 850.00, 900.00, 950.00, 1000.00, 1100.00, 1200.00, &
       1300.00, 1400.00, 1500.00, 1600.00, 1700.00, 1800.00, 1900.00, &
       2000.00, 2100.00, 2200.00, 2300.00, 2400.00, 2500.00, 2600.00, &
       2700.00, 2800.00, 2900.00, 3000.00, 3100.00, 3200.00, 3300.00, &
       3400.00, 3500.00, 3600.00, 3700.00, 3800.00, 3900.00, 4000.00, &
       4100.00, 4200.00, 4300.00, 4400.00, 4500.00, 4600.00, 4700.00, &
       4800.00, 4900.00, 5000.00, 6000.00, 7000.00, 8000.00   /)

   AED_REAL, parameter :: landai(150) = (/ &
       0.0000, 0.0001, 0.0000, 0.0004, &
       0.0009,0.0017,0.0030,0.0050,0.0100,0.0180,0.0300,0.0420,0.0520, &
       0.0540,0.0580,0.0640,0.0640,0.1000,0.1300,0.2000,0.2500,0.2200, &
       0.2400,0.3400,0.5200,0.6300,0.6100,0.6700,0.7600,0.8200,0.8500, &
       1.0200,1.1500,1.1100,1.1100,1.1700,1.1800,1.1600,1.1600,1.2900, &
       1.3300,1.3200,1.2300,1.1500,1.1200,1.2000,1.5400,1.8800,1.9400, &
       1.9200,1.9200,1.8900,1.7800,1.8200,2.0300,2.1500,2.2000,2.1900, &
       2.1600,2.1500,2.1700,2.2000,2.1600,2.0300,1.9900,2.0400,1.9800, &
       1.9700,1.9600,1.8900,1.8700,1.9200,1.9500,1.9700,1.9800,1.9800, &
       1.9500,1.9200,1.9000,1.8900,1.8700,1.8700,1.8700,1.8500,1.8400, &
       1.8300,1.8100,1.7700,1.7400,1.7000,1.6600,1.6200,1.5900,1.5500, &
       1.5100,1.4800,1.4400,1.4100,1.3700,1.3400,1.3000,1.2700,1.1270, &
       1.0030,0.8950,0.8030,0.7250,0.6060,0.5010,0.4060,0.3280,0.2670, &
       0.2200,0.1820,0.1520,0.1274,0.1079,0.0917,0.0785,0.0676,0.0585, &
       0.0509,0.0445,0.0390,0.0343,0.0303,0.0268,0.0230,0.0214,0.0191, &
       0.0171,0.0153,0.0139,0.0125,0.0114,0.0103,0.0095,0.0087,0.0080, &
       0.0073,0.0067,0.0061,0.0056,0.0051,0.0048,0.0044,0.0042,0.0021, &
       0.0012,0.0006 /)

   INTEGER, PARAMETER :: SUBMERGED = 1
   INTEGER, PARAMETER :: EMERGENT  = 2
   INTEGER, PARAMETER :: FLOATING  = 3
   INTEGER            :: diag_level = 10     ! 0 = no diagnostic outputs
                                             ! 1 = basic diagnostic outputs
                                             ! 2 = flux rates, and supporitng
                                             ! 3 = other metrics
                                             !10 = all debug & checking outputs
!-------------------------------------------------------------------------------



CONTAINS
!===============================================================================




!###############################################################################
INTEGER FUNCTION load_csv(dbase, md)
!-------------------------------------------------------------------------------
   USE aed_csv_reader
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*),INTENT(in) :: dbase
   TYPE(macrophyte_params_t) :: md(MAX_PHYTO_TYPES)
!
!LOCALS
   INTEGER :: unit, nccols, ccol, dcol
   CHARACTER(len=32),POINTER,DIMENSION(:) :: csvnames
   CHARACTER(len=32) :: name
   TYPE(AED_SYMBOL),DIMENSION(:),ALLOCATABLE :: values
   INTEGER :: idx_col = 0, ret = 0
   LOGICAL :: meh
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
         dcol = ccol-1
         md(dcol)%m_name = csvnames(ccol)

         CALL copy_name(values(1), name)
         SELECT CASE (name)
            CASE ('growth_form')    ; md(dcol)%growth_form  = extract_integer(values(ccol))
            CASE ('zone_lock')      ; md(dcol)%zone_lock    = extract_integer(values(ccol))
            CASE ('m_initial')      ; md(dcol)%m_initial    = extract_double(values(ccol))
            CASE ('m_initial_minb') ; md(dcol)%m_initial_minb = extract_double(values(ccol))
            CASE ('m_initial_mind') ; md(dcol)%m_initial_mind = extract_double(values(ccol))
            CASE ('m0')             ; md(dcol)%m0           = extract_double(values(ccol))

            CASE ('drag_model')     ; md(dcol)%drag_model   = extract_integer(values(ccol))
            CASE ('K_CD')           ; md(dcol)%K_CD         = extract_double(values(ccol))
            CASE ('shoot_diameter') ; md(dcol)%shoot_diameter = extract_double(values(ccol))
            CASE ('shoot_length')   ; md(dcol)%shoot_length = extract_double(values(ccol))
            CASE ('tissue_density') ; md(dcol)%tissue_density = extract_double(values(ccol))
            CASE ('elastic_modulus'); md(dcol)%elastic_modulus = extract_double(values(ccol))

            CASE ('shading_model')  ; md(dcol)%shading_model= extract_integer(values(ccol))
            CASE ('KeMAC')          ; md(dcol)%KeMAC        = extract_double(values(ccol))
            CASE ('kA')             ; md(dcol)%kA           = extract_double(values(ccol))
            CASE ('k_omega')        ; md(dcol)%k_omega      = extract_double(values(ccol))
            CASE ('sine_blade')     ; md(dcol)%sine_blade   = extract_double(values(ccol))
            CASE ('X_ldw')          ; md(dcol)%X_ldw        = extract_double(values(ccol))
            CASE ('K_ldw')          ; md(dcol)%K_ldw        = extract_double(values(ccol))
            CASE ('X_sdw')          ; md(dcol)%X_sdw        = extract_double(values(ccol))

            CASE ('R_growth')       ; md(dcol)%R_growth     = extract_double(values(ccol))
            CASE ('temp_method')    ; md(dcol)%temp_method  = extract_integer(values(ccol))
            CASE ('theta_growth')   ; md(dcol)%theta_growth = extract_double(values(ccol))
            CASE ('T_std')          ; md(dcol)%T_std        = extract_double(values(ccol))
            CASE ('T_opt')          ; md(dcol)%T_opt        = extract_double(values(ccol))
            CASE ('T_max')          ; md(dcol)%T_max        = extract_double(values(ccol))
            CASE ('light_model')    ; md(dcol)%light_model  = extract_integer(values(ccol))
            CASE ('I_K')            ; md(dcol)%I_K          = extract_double(values(ccol))
            CASE ('I_S')            ; md(dcol)%I_S          = extract_double(values(ccol))
            CASE ('f_pr')           ; md(dcol)%f_pr         = extract_double(values(ccol))
            CASE ('K_epi')          ; md(dcol)%K_epi        = extract_double(values(ccol))

            CASE ('R_resp')         ; md(dcol)%R_resp       = extract_double(values(ccol))
            CASE ('theta_resp')     ; md(dcol)%theta_resp   = extract_double(values(ccol))
            CASE ('E_comp')         ; md(dcol)%E_comp       = extract_double(values(ccol))
            CASE ('R_mort')         ; md(dcol)%R_mort       = extract_double(values(ccol))
            CASE ('sal_method')     ; md(dcol)%sal_method   = extract_integer(values(ccol))
            CASE ('S_bep')          ; md(dcol)%S_bep        = extract_double(values(ccol))
            CASE ('S_maxsp')        ; md(dcol)%S_maxsp      = extract_double(values(ccol))
            CASE ('S_opt')          ; md(dcol)%S_opt        = extract_double(values(ccol))

            CASE ('m_initial_bg')   ; md(dcol)%m_initial_bg = extract_double(values(ccol))
            CASE ('f_bg')           ; md(dcol)%f_bg         = extract_double(values(ccol))
            CASE ('R_resp_bg')      ; md(dcol)%R_resp_bg    = extract_double(values(ccol))
            CASE ('R_mort_bg')      ; md(dcol)%R_mort_bg    = extract_double(values(ccol))
            CASE ('tau_tran')       ; md(dcol)%tau_tran     = extract_double(values(ccol))
            CASE ('X_cchl')         ; md(dcol)%X_cchl       = extract_double(values(ccol))
            CASE ('X_cdw')          ; md(dcol)%X_cdw        = extract_double(values(ccol))
            CASE ('K_N')            ; md(dcol)%K_N          = extract_double(values(ccol))
            CASE ('X_ncon')         ; md(dcol)%X_ncon       = extract_double(values(ccol))
            CASE ('K_P')            ; md(dcol)%K_P          = extract_double(values(ccol))
            CASE ('X_pcon')         ; md(dcol)%X_pcon       = extract_double(values(ccol))
            CASE ('K_S')            ; md(dcol)%K_S          = extract_double(values(ccol))
            CASE ('K_NPP')          ; md(dcol)%K_NPP        = extract_double(values(ccol))

            CASE ('fruit_model')    ; md(dcol)%fruit_model  = extract_integer(values(ccol))
            CASE ('m_initial_f')    ; md(dcol)%m_initial_f  = extract_double(values(ccol))
            CASE ('tau_tran_fruit') ; md(dcol)%tau_tran_fruit = extract_double(values(ccol))
            CASE ('r_release')      ; md(dcol)%r_release    = extract_double(values(ccol))
            CASE ('f_seed')         ; md(dcol)%f_seed       = extract_double(values(ccol))
            CASE ('t_start_g')      ; md(dcol)%t_start_g    = extract_double(values(ccol))
            CASE ('t_dur_g')        ; md(dcol)%t_dur_g      = extract_double(values(ccol))
            CASE ('t_start_r')      ; md(dcol)%t_start_r    = extract_double(values(ccol))
            CASE ('t_dur_r')        ; md(dcol)%t_dur_r      = extract_double(values(ccol))

            CASE DEFAULT ; print *, 'Unknown row "', TRIM(name), '"'
         END SELECT
      ENDDO
   ENDDO

   meh = aed_csv_close(unit) !# don't care if close fails

   IF (ASSOCIATED(csvnames)) DEALLOCATE(csvnames)
   IF (ALLOCATED(values))    DEALLOCATE(values)

   load_csv = ret
END FUNCTION load_csv
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_macrophyte_load_params(data, dbase, count, list)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(inout) :: data
   CHARACTER(len=*),INTENT(in) :: dbase
   INTEGER,INTENT(in)          :: count
   INTEGER,INTENT(in)          :: list(*)
!
!LOCALS
   INTEGER  :: status
   INTEGER  :: i,tfil
   AED_REAL :: minNut

   TYPE(macrophyte_params_t),ALLOCATABLE :: md(:)
   NAMELIST /macrophyte_data/ md  ! %% type : macrophyte_params_t - see above
!-------------------------------------------------------------------------------
!BEGIN
    ALLOCATE(md(MAX_PHYTO_TYPES))
    SELECT CASE (param_file_type(dbase))
       CASE (CSV_TYPE)
           status = load_csv(dbase, md)
       CASE (NML_TYPE)
           open(NEWUNIT=tfil,file=dbase, status='OLD', iostat=status)
           IF (status /= 0) STOP 'Cannot open macrophyte_data namelist file for macrophytes'
           read(tfil,nml=macrophyte_data,iostat=status)
           close(tfil)
       CASE DEFAULT
           print *,'Unknown file type "',TRIM(dbase),'"'; status=1
    END SELECT
    IF (status /= 0) STOP 'Error reading namelist macrophyte_data for macrophytes'

    data%num_mphy = count
    ALLOCATE(data%mpars(count))

    ALLOCATE(data%id_mphya(count))
    ALLOCATE(data%id_mphyb(count))
    ALLOCATE(data%id_mphyf(count))

    DO i=1,count
       ! Assign parameters from database to simulated groups
       data%mpars(i)%growth_form    = md(list(i))%growth_form
       data%mpars(i)%zone_lock      = md(list(i))%zone_lock
       data%mpars(i)%m_name         = md(list(i))%m_name
       data%mpars(i)%m_initial      = md(list(i))%m_initial
       data%mpars(i)%m_initial_minb = md(list(i))%m_initial_minb
       data%mpars(i)%m_initial_mind = md(list(i))%m_initial_mind
       data%mpars(i)%m0             = md(list(i))%m0

       data%mpars(i)%drag_model     = md(list(i))%drag_model
       data%mpars(i)%K_CD           = md(list(i))%K_CD
       data%mpars(i)%shoot_diameter = md(list(i))%shoot_diameter
       data%mpars(i)%shoot_length   = md(list(i))%shoot_length
       data%mpars(i)%tissue_density = md(list(i))%tissue_density
       data%mpars(i)%elastic_modulus= md(list(i))%elastic_modulus

       data%mpars(i)%shading_model  = md(list(i))%shading_model
       data%mpars(i)%KeMAC          = md(list(i))%KeMAC
       data%mpars(i)%kA             = md(list(i))%kA
       data%mpars(i)%k_omega        = md(list(i))%k_omega
       data%mpars(i)%sine_blade     = md(list(i))%sine_blade
       data%mpars(i)%X_ldw          = md(list(i))%X_ldw
       data%mpars(i)%K_ldw          = md(list(i))%K_ldw
       data%mpars(i)%X_sdw          = md(list(i))%X_sdw

       data%mpars(i)%R_growth       = md(list(i))%R_growth/secs_per_day
       data%mpars(i)%temp_method    = md(list(i))%temp_method
       data%mpars(i)%theta_growth   = md(list(i))%theta_growth
       data%mpars(i)%T_std          = md(list(i))%T_std
       data%mpars(i)%T_opt          = md(list(i))%T_opt
       data%mpars(i)%T_max          = md(list(i))%T_max
       data%mpars(i)%light_model    = md(list(i))%light_model
       data%mpars(i)%I_K            = md(list(i))%I_K
       data%mpars(i)%I_S            = md(list(i))%I_S
       data%mpars(i)%f_pr           = md(list(i))%f_pr
       data%mpars(i)%K_epi          = md(list(i))%K_epi

       data%mpars(i)%R_resp         = md(list(i))%R_resp/secs_per_day
       data%mpars(i)%theta_resp     = md(list(i))%theta_resp
       data%mpars(i)%E_comp         = md(list(i))%E_comp
       data%mpars(i)%R_mort         = md(list(i))%R_mort/secs_per_day
       data%mpars(i)%sal_method     = md(list(i))%sal_method
       data%mpars(i)%S_bep          = md(list(i))%S_bep
       data%mpars(i)%S_maxsp        = md(list(i))%S_maxsp
       data%mpars(i)%S_opt          = md(list(i))%S_opt

       data%mpars(i)%m_initial_bg   = md(list(i))%m_initial_bg
       data%mpars(i)%f_bg           = md(list(i))%f_bg
       data%mpars(i)%R_resp_bg      = md(list(i))%R_resp_bg/secs_per_day
       data%mpars(i)%R_mort_bg      = md(list(i))%R_mort_bg/secs_per_day
       data%mpars(i)%tau_tran       = md(list(i))%tau_tran/secs_per_day
       data%mpars(i)%X_cchl         = md(list(i))%X_cchl
       data%mpars(i)%X_cdw          = md(list(i))%X_cdw
       data%mpars(i)%K_N            = md(list(i))%K_N
       data%mpars(i)%X_ncon         = md(list(i))%X_ncon
       data%mpars(i)%K_P            = md(list(i))%K_P
       data%mpars(i)%X_pcon         = md(list(i))%X_pcon
       data%mpars(i)%K_S            = md(list(i))%K_S
       data%mpars(i)%K_NPP          = md(list(i))%K_NPP

       data%mpars(i)%fruit_model    = md(list(i))%fruit_model
       data%mpars(i)%m_initial_f    = md(list(i))%m_initial_f
       data%mpars(i)%tau_tran_fruit = md(list(i))%tau_tran_fruit/secs_per_day
       data%mpars(i)%r_release      = md(list(i))%r_release/secs_per_day
       data%mpars(i)%f_seed         = md(list(i))%f_seed
       data%mpars(i)%t_start_g      = md(list(i))%t_start_g
       data%mpars(i)%t_dur_g        = md(list(i))%t_dur_g
       data%mpars(i)%t_start_r      = md(list(i))%t_start_r
       data%mpars(i)%t_dur_r        = md(list(i))%t_dur_r
       IF (.NOT. data%simFruiting)  data%mpars(i)%fruit_model = 0


       ! Register group's above ground biomass as a state variable
       data%id_mphya(i) = aed_define_sheet_variable(                         &
                              TRIM(md(list(i))%m_name)//'_ag', 'mmol C/m2',  &
                               'macrophyte group above-ground biomass',      &
                              md(list(i))%m_initial,                         &
                              minimum=md(list(i))%m0)
       ! Register group's below ground biomass as a state variable
       data%id_mphyb(i) = aed_define_sheet_variable(                         &
                              TRIM(md(list(i))%m_name)//'_bg', 'mmol C/m2',  &
                              'macrophyte group below-ground biomass',       &
                              md(list(i))%m_initial_bg,                      &
                              minimum=md(list(i))%m0)
       ! Register group's fruit biomass as a state variable
       IF( data%mpars(i)%fruit_model >0 ) THEN
          data%id_mphyf(i) = aed_define_sheet_variable(                      &
                              TRIM(md(list(i))%m_name)//'_fr', 'mmol C/m2',  &
                               'macrophyte group seed/fruit biomass',        &
                              md(list(i))%m_initial_f,                       &
                              minimum=md(list(i))%m0/2.)
       ENDIF
!       ! Register group's shoot number as diagnostic
!       IF( data%mpars(i)%drag_model >0 ) THEN
!         data%id_shdens(i) = aed_define_sheet_diag_variable(                 &
!                             TRIM(md(list(i))%m_name)//'_shdens', '# sh/m2', &
!                              'macrophyte group shoot density')
!         data%id_shhght(i) = aed_define_sheet_diag_variable(                 &
!                             TRIM(md(list(i))%m_name)//'_shhght', '# m', &
!                              'macrophyte group shoot/leaf length')
!       ENDIF

   ENDDO
    DEALLOCATE(md)
END SUBROUTINE aed_macrophyte_load_params
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_define_macrophyte(data, namlst)
!-------------------------------------------------------------------------------
! Initialise the macrophyte/seagrass  model
!
!  Here, the aed_ namelist is read and the variables exported
!  by the model are registered
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(inout) :: data
   INTEGER,INTENT(in) :: namlst
!
!LOCALS
   INTEGER  :: status, i, ww
   AED_REAL :: landaimid

   AED_REAL :: Wli(8), ALli(8)

!  %% NAMELIST   %%  /aed_macrophyte/
!  %% Last Checked 20/08/2021
   INTEGER            :: num_mphy = 0
   INTEGER            :: the_mphy(MAX_PHYTO_TYPES) = 0
   CHARACTER(len=128) :: dbase = 'aed_macrophyte_pars.nml'

   INTEGER            :: n_zones = 0
   INTEGER            :: active_zones(MAX_ZONES)

   LOGICAL            :: simStaticBiomass = .FALSE.
   LOGICAL            :: simMacFeedback = .FALSE.
   LOGICAL            :: simEpiphytes = .FALSE.
   LOGICAL            :: simFruiting = .FALSE.
   LOGICAL            :: simMacDrag = .FALSE.

   INTEGER            :: mac_initial = 0

   AED_REAL           :: water_nutrient_frac = zero_

   !AED_REAL           :: bg_gpp_frac = 0.1
   !AED_REAL           :: coef_bm_hgt = 1e-5

   INTEGER            :: epi_model = 1
   AED_REAL           :: epi_initial = 0.1
   AED_REAL           :: R_epig = zero_
   AED_REAL           :: R_epib = zero_ 
   AED_REAL           :: R_epir = zero_ 
   AED_REAL           :: I_Kepi = 100.
   AED_REAL           :: epi_Xnc = 16./106.
   AED_REAL           :: epi_Xpc = 1./106.
   AED_REAL           :: epi_max = 10000.
   AED_REAL           :: theta_epi_growth = 1.06
   AED_REAL           :: theta_epi_resp = 1.06
   AED_REAL           :: epi_K_N = 10.
   AED_REAL           :: epi_K_P = 1.0

!  %% END NAMELIST   %%  /aed_macrophyte/

   NAMELIST /aed_macrophyte/ num_mphy, the_mphy, dbase,               &
                             n_zones, active_zones,                   &
                             simMacFeedback, simStaticBiomass,        &
                             simEpiphytes, simFruiting, simMacDrag,   &
                             water_nutrient_frac,                     &
                             epi_model,R_epig,R_epib,R_epir,I_Kepi,epi_max,  &
                             epi_Xnc,epi_Xpc,epi_K_N,epi_K_P,         &
                             theta_epi_growth,theta_epi_resp,         &
                             mac_initial, epi_initial, diag_level

!-----------------------------------------------------------------------
!BEGIN
   print *,"        aed_macrophyte configuration"

   ! Read the namelist
   read(namlst,nml=aed_macrophyte,iostat=status)
   IF (status /= 0) STOP 'Error reading namelist aed_macrophyte'

   data%simStaticBiomass = simStaticBiomass
   data%simMacFeedback = simMacFeedback
   data%simEpiphytes = simEpiphytes
   data%simFruiting = simFruiting
   data%simMacDrag = simMacDrag
   data%water_nutrient_frac = water_nutrient_frac

   data%mac_initial = mac_initial
   data%epi_initial = epi_initial

   IF( .NOT. simEpiphytes) epi_model = 0; data%epi_model = epi_model
   data%theta_epi_growth = theta_epi_growth ; data%theta_epi_resp = theta_epi_resp
   data%epi_Xpc = epi_Xpc ; data%epi_K_P = epi_K_P
   data%epi_K_N = epi_K_N ; data%epi_Xnc = epi_Xnc
   data%I_Kepi = I_Kepi ; data%epi_max = epi_max
   data%R_epig = R_epig /secs_per_day
   data%R_epib = R_epib /secs_per_day
   data%R_epir = R_epir /secs_per_day

  ! !remove
  ! data%bg_gpp_frac = bg_gpp_frac  ! make this species specific
  ! data%coef_bm_hgt = coef_bm_hgt  ! make this species specific

   data%n_zones = n_zones
   IF (n_zones > 0) THEN
      ALLOCATE(data%active_zones(n_zones))
      DO i=1,n_zones
         data%active_zones(i) = active_zones(i)
      ENDDO
   ENDIF

   ! Store parameter values in modules own data structure (macrophyte_params_t)
   !   NB: all rates are provided in values per day,
   !       and are converted here to values per second for internal AED use.
   !  The requested macrophyte state variables are registered in this function
   CALL aed_macrophyte_load_params(data, dbase, num_mphy, the_mphy)

   CALL aed_bio_temp_function(  data%num_mphy,               &
                                data%mpars%theta_growth,     &
                                data%mpars%T_std,            &
                                data%mpars%T_opt,            &
                                data%mpars%T_max,            &
                                data%mpars%aTn,              &
                                data%mpars%bTn,              &
                                data%mpars%kTn,              &
                                data%mpars%m_name)

   ! Variables related to epiphytes
   IF ( simEpiphytes ) THEN
      data%id_epi = aed_define_sheet_variable( 'epiphyte', 'mmol C/m2',       &
                                               'epiphyte biomass',            &
                                                epi_initial, 0.001            )
     !data%id_epib = aed_define_sheet_diag_variable('epi_ben','mmol C/m2','total epiphyte biomass')
      IF( diag_level > 0) THEN
        data%id_epib = aed_define_sheet_diag_variable('epi_ben','mmol C/m2','total epiphyte biomass')
        IF( diag_level > 1) THEN
         data%id_epig = aed_define_sheet_diag_variable('epi_gpp','mmol C/m2/d','benthic gross productivity')
         data%id_epir = aed_define_sheet_diag_variable('epi_rsp','mmol C/m2/d','benthic phyto respiration')
        ENDIF
       ENDIF
   ENDIF


   ! Register diagnostic variables, for macrophyte the community canopy
   data%id_d_par  = aed_define_sheet_diag_variable('par','W/m2','light intensity at canopy top')
!  data%id_p2r    = aed_define_sheet_diag_variable('npp','/d','macrophyte net productivity')
   data%id_gpp    = aed_define_sheet_diag_variable('gpp','mmol C/m2/d','benthic plant productivity')
   data%id_mac    = aed_define_sheet_diag_variable('mac_ben','mmol C/m2','total macrophyte biomass')
   data%id_canopy_lai = aed_define_sheet_diag_variable('mac_lai','m2/m2','macrophyte leaf area density')
   data%id_mac_ag = aed_define_sheet_diag_variable('mac_ag','mmol C/m2','total above ground macrophyte biomass')
   data%id_mac_bg = aed_define_sheet_diag_variable('mac_bg','mmol C/m2','total below ground macrophyte biomass')
   data%id_root_d = aed_define_sheet_diag_variable('mac_root_depth','m','mean depth of roots below the sediment surface')
   data%id_root_o = aed_define_sheet_diag_variable('mac_root_o2','mmol O2/m2','mean O2 injection rate of roots into sediment')
   data%id_mac_tr = aed_define_sheet_diag_variable('trn_bg','mmol C/m2/day', 'macrophyte biomass translocation')

   ! Variables related to fruiting
   IF ( ANY(data%mpars(:)%fruit_model >0) ) THEN
      data%id_mac_fr = aed_define_sheet_diag_variable('mac_fr','mmol C/m2', 'macrophyte community fruit biomass')
      data%id_mac_ft = aed_define_sheet_diag_variable('trn_fr','mmol C/m2/day', 'macrophyte fruit biomass translocation')
      data%id_mac_gt = aed_define_sheet_diag_variable('tri_fr','mmol C/m2/day', 'macrophyte fruit growth trigger')
   ENDIF

   ! Variables for macrophyte canopy light attenuation
   data%id_kemac = aed_define_diag_variable('ke_mac','/m', 'macrophyte light attenuation')
   IF (ANY(data%mpars(:)%light_model == 2) .or. ANY(data%mpars(:)%light_model == 3)) THEN
      ! For spectral light capture
      data%id_parc = aed_define_sheet_diag_variable('mac_ag_parc','', 'macrophyte leaf photon capture')
      data%id_aeff = aed_define_sheet_diag_variable('mac_ag_aeff','', 'macrophyte leaf effective area')
   ENDIF

   ! Variables for drag feedback, computed based on macrophyte community canopy properties
   data%id_canopy_height  = aed_define_sheet_diag_variable('canopy_height','m', '')
   data%id_canopy_sh_dens = aed_define_sheet_diag_variable('canopy_shoot_dens','shoot/m2', '')
   data%id_canopy_sh_diam = aed_define_sheet_diag_variable('canopy_shoot_diam','m/shoot', '')
 !  IF (data%drag_model > 1) THEN
      data%id_canopy_drag      = aed_define_sheet_diag_variable('canopy_drag','-', '')       !2D becoming 3D
      data%id_canopy_velocity= aed_define_sheet_diag_variable('canopy_velocity','m/s', '')
      data%id_canopy_biovolume = aed_define_sheet_diag_variable('canopy_biovolume','m3 leaf/m2', '')
      data%id_canopy_blockage  = aed_define_diag_variable('canopy_blockage','m3/m2', '')
      data%id_canopy_frarea    = aed_define_diag_variable('canopy_frarea','m2/m2?', '')
 !  ENDIF

   ! Link to variables from other modules
   IF ( simMacFeedback ) THEN
      data%id_oxy = aed_locate_variable('OXY_oxy')
      data%id_nox = aed_locate_variable('NIT_nit')
      data%id_nh4 = aed_locate_variable('NIT_amm')
      data%id_po4 = aed_locate_variable('PHS_frp')
      data%id_dic = 0 !aed_locate_sheet_variable('CAR_dic')
   ENDIF
        
   ! Register environmental dependenciess
   data%id_tem      = aed_locate_global('temperature')
   data%id_extc     = aed_locate_global('extc_coef')
   data%id_sal      = aed_locate_global('salinity')
   data%id_dz       = aed_locate_global('layer_ht')
   data%id_par      = aed_locate_global('par')
   data%id_sed_zone = aed_locate_sheet_global('sed_zone')
   data%id_atem     = aed_locate_sheet_global('air_temp')
   data%id_yearday  = aed_locate_sheet_global('yearday')
   data%id_I_0      = aed_locate_sheet_global('par_sf')
   data%id_depth    = aed_locate_sheet_global('col_depth')


   ! Light spectra pre-processing tasks
   ! calculate the integrated irradiance over the 300-800 nm wave lengths
   WLint=wavei(27:103)             ! 300-800 nm
   landaiint=landai(27:103)        ! corresponding irradiance at 300-800 nm
   ! total irradiance integrate over the 300-800 nm
   landaiSUM=0
   DO ww=1,SIZE(WLint)-1
     landaimid = (landaiint(ww)+landaiint(ww+1))/2
     !WLmid = (WLint(ww)+WLint(ww+1))/2 !NOT USED
     landaiSUM=landaiSUM+landaimid*(WLint(ww+1)-WLint(ww))
   ENDDO
   ! leaf light absorbance at selected wavelengths; from Baird et al (2016)
   Wli = (/ 300,  390,  500,  530,  640,  680,  705,  800 /)
   ALli= (/ 0.72, 0.72, 0.68, 0.38, 0.38, 0.60, 0.04, 0.0 /)
   ! interpolate into fine-resolution wave lengths for further processing
   ALLOCATE (ALlint(501))
   CALL interp_0d(SIZE(Wli),Wli,ALli,SIZE(WLint),WLint,ALlint) ! ALlint=interp1(Wli,ALli,WLint)
   ALl=0.5412  !SUM(ALlint)/SIZE(ALlint)  ! overall (mean) leaf absorbance, fixed to 0.5412 as SUM(ALlint)/SIZE(ALlint) gives 0.083


END SUBROUTINE aed_define_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_initialize_benthic_macrophyte(data, column, layer_idx)
!-------------------------------------------------------------------------------
! Routine to initialize bottom diagnostics, and other checks
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mi
   AED_REAL :: matz, mphy

   AED_REAL :: lai_to_biomass, depth
   AED_REAL :: n_stems, v_diameter, v_height
   AED_REAL :: bio_volume, wet_biovolume, dry_biovolume
   AED_REAL :: n_shoots,sh_diameter, sh_height
   AED_REAL :: minBiomass, maxBiomass, minDepth, maxDepth, scale, offset, a ,biomass

   AED_REAL, PARAMETER :: coef_bm_hgt = 0.5

   AED_REAL :: MAC_A, MAC_B
!-------------------------------------------------------------------------------
!BEGIN

   !---------------------------------------------------------------------------
   ! Check to ensure this cell/zone is colonisable
   matz = _STATE_VAR_S_(data%id_sed_zone)

   IF ( .NOT. in_zone_set(matz, data%active_zones) ) THEN
     _DIAG_VAR_S_(data%id_mac) = zero_
     _DIAG_VAR_S_(data%id_mac_ag) = zero_
     _DIAG_VAR_S_(data%id_mac_bg) = zero_
     _DIAG_VAR_S_(data%id_canopy_lai) = zero_
     _DIAG_VAR_S_(data%id_gpp) = zero_
     _DIAG_VAR_S_(data%id_root_o) = zero_
     _DIAG_VAR_S_(data%id_root_d) = zero_
     _DIAG_VAR_S_(data%id_canopy_sh_dens) = zero_

     DO mi=1,data%num_mphy
      _STATE_VAR_S_(data%id_mphya(mi))  =  zero_
      _STATE_VAR_S_(data%id_mphyb(mi))  =  zero_
      IF( data%mpars(mi)%fruit_model >0 ) THEN
       _STATE_VAR_S_(data%id_mphyf(mi))  =  zero_
      ENDIF
     ENDDO
     RETURN
   ENDIF


   !---------------------------------------------------------------------------
   ! Set initial diagnostics
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_

   _DIAG_VAR_S_(data%id_gpp) = zero_
   _DIAG_VAR_S_(data%id_root_o) = zero_
   _DIAG_VAR_S_(data%id_root_d) = MAX( MIN(_DIAG_VAR_S_(data%id_mac_ag) * coef_bm_hgt,0.25),0.01)

   _DIAG_VAR_S_(data%id_canopy_biovolume) = zero_

   IF ( ANY(data%mpars(:)%fruit_model >0) ) _DIAG_VAR_S_(data%id_mac_gt) = one_

   !---------------------------------------------------------------------------
   ! (Re)set local bottom cell macrophyte details with data read in from benthic maps


   SELECT CASE (data%mac_initial)

      CASE ( 0 )
         ! Re-inforce "m_initial" etc, as read in from parameters
         DO mi=1,data%num_mphy
           _STATE_VAR_S_(data%id_mphya(mi))  =  data%mpars(mi)%m_initial         ! leaf biomass
           _STATE_VAR_S_(data%id_mphyb(mi))  =  data%mpars(mi)%m_initial_bg      ! below-ground biomass
           IF( data%mpars(mi)%fruit_model >0 ) THEN
            _STATE_VAR_S_(data%id_mphyf(mi))  =  data%mpars(mi)%m_initial_f      ! fruit biomass
           ENDIF
         ENDDO

      CASE ( 1 )
         ! compute canopy variables from input biomass and depth (Cockburn Sound)
         depth = _STATE_VAR_S_(data%id_depth)      ! incoming water depth above cell

         DO mi=1,data%num_mphy
           ! Check if this species/group is locked to grow in a specific zone
           IF( data%mpars(mi)%zone_lock >0 .and. data%mpars(mi)%zone_lock .ne. INT(matz) ) CYCLE

            !minBiomass = 1.                                    ! min biomass (log10 g DW/m2)
            minBiomass = data%mpars(mi)%m_initial_minb          ! min biomass from pars file is mmol C/m2
            minBiomass = LOG10(minBiomass/data%mpars(mi)%X_cdw) ! min biomass (log10 g DW/m2)
            maxBiomass = data%mpars(mi)%m_initial               ! max biomass from pars file is mmol C/m2
            maxBiomass = LOG10(maxBiomass/data%mpars(mi)%X_cdw) ! max biomass (log10 g DW/m2)
            minDepth   = 1.2                                    ! min depth below which is max biomass;
            maxDepth   = data%mpars(mi)%m_initial_mind          ! max depth above which is 0 biomass; was 13.8

            scale      = 12. / (maxDepth-minDepth)
            offset     = -6./scale - minDepth
            a          = (depth+offset)*scale
            biomass    = minBiomass + exp(-a)/(1.+exp(-a))*(maxBiomass-minBiomass)
            biomass    = 10.**(biomass)                         ! g DW/m2
            n_shoots   = biomass * data%mpars(mi)%X_sdw         ! shoots/m2

            _STATE_VAR_S_(data%id_mphya(mi)) = biomass*data%mpars(mi)%X_cdw        ! mmol C/m2  !0.5/12.*1000.
            _STATE_VAR_S_(data%id_mphyb(mi)) =  data%mpars(mi)%f_bg * _STATE_VAR_S_(data%id_mphya(mi)) ! below-ground biomass
            IF( data%mpars(mi)%fruit_model >0 ) THEN
               _STATE_VAR_S_(data%id_mphyf(mi)) = data%mpars(mi)%f_seed * _STATE_VAR_S_(data%id_mphya(mi)) ! fruit biomass
            ENDIF
            _DIAG_VAR_S_(data%id_canopy_sh_dens) = _DIAG_VAR_S_(data%id_canopy_sh_dens) + n_shoots
         ENDDO

      CASE ( 2 )
         ! compute macrophyte group biomass from a single LAI and unit conversion
         lai_to_biomass = 1e6 ! mmol C/m2 / unit LAI
         DO mi=1,data%num_mphy
            _STATE_VAR_S_(data%id_mphya(mi)) =  &
                           lai_to_biomass * _DIAG_VAR_S_(data%id_canopy_lai) / data%num_mphy
         END DO

      CASE ( 3 )
         ! compute macrophyte group biomass from input canopy variables

         DO mi=1,data%num_mphy
            n_shoots    = INT(_DIAG_VAR_S_(data%id_shdens(mi)))
            sh_diameter = _DIAG_VAR_S_(data%id_shdiam(mi))
            sh_height   = _DIAG_VAR_S_(data%id_shhght(mi))

            ! actual volume of vegetation
            bio_volume = n_shoots * 3.1418 * sh_diameter**2 / 4 * sh_height

            ! volume above the water
            wet_biovolume = bio_volume * MIN( depth/sh_height, one_ )
            dry_biovolume = bio_volume - wet_biovolume

            ! cumulated volume of all macrophyte groups in this bottom cell
            _DIAG_VAR_S_(data%id_canopy_biovolume) = _DIAG_VAR_S_(data%id_canopy_biovolume) + bio_volume

           ! amount of carbon biomass (within bottom cell)
           _STATE_VAR_S_(data%id_mphya(mi)) =  bio_volume * data%mpars(mi)%X_cchl  ! check stoich
        END DO


      CASE ( 4 )


      CASE DEFAULT
         ! no re-initialisation required


   END SELECT

   DO mi=1,data%num_mphy
      mphy = _STATE_VAR_S_(data%id_mphya(mi))
      _DIAG_VAR_S_(data%id_mac) = _DIAG_VAR_S_(data%id_mac) + mphy
      _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + mphy*(one_-data%mpars(mi)%f_bg)
      _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + mphy*(data%mpars(mi)%f_bg)
      _DIAG_VAR_S_(data%id_canopy_lai) = _DIAG_VAR_S_(data%id_canopy_lai) +       &
                    (one_ - exp(-data%mpars(mi)%k_omega * mphy*(one_-data%mpars(mi)%f_bg)))
  ENDDO

  _STATE_VAR_S_(data%id_epi)  =  data%epi_initial

END SUBROUTINE aed_initialize_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_column_macrophyte(data,column,layer_map)
!-------------------------------------------------------------------------------
! Vertical column loop, to compute epiphytes and macrophyte biovolume
! Note that variables are set here before computations are completed in benthic
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_map(:)
!
!LOCALS
   INTEGER  :: layer_idx, layer, mi, done
   AED_REAL :: yearday, dz, matz
   AED_REAL :: hgt_layer_bot, hgt_canopy, hgt_layer_top, layer_frac, canopy_leaf_area
   AED_REAL :: epi, epi_prod, epi_resp, epi_flux, fI, par, extc, temp, Io, epi_max
   AED_REAL :: canopy_par, canopy_frac, ratio

   !-------------------------------------------------------------------------------
   !BEGIN

   ! Check this column is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   IF ( .NOT. in_zone_set(matz, data%active_zones) ) RETURN

   layer_idx = layer_map(1)                    ! Set layer to the top of the column
   layer_frac = zero_

   ! Set column environmental conditions (from host)
   yearday  = _STATE_VAR_S_(data%id_yearday)
   extc = _STATE_VAR_(data%id_extc)            ! extinction
   par = _STATE_VAR_(data%id_par)              ! photosynthetically active radiation
   dz = _STATE_VAR_(data%id_dz)                ! cell thickness
   Io = _STATE_VAR_S_(data%id_I_0)             ! surface short wave radiation

   ! Refer to total canopy height, and if exceeds the benthic/bottom layer,
   !  then distribute relevant properties vertically
   mi = 1
   done = 0


   ! If submerged plant group, loop up through the water column
   IF(data%mpars(mi)%growth_form == SUBMERGED) THEN

         ! Initialise all layers canopy properties to 0 (these are on the 3D grid)
         DO layer = SIZE(layer_map),1,-1
            layer_idx = layer_map(layer)
            _DIAG_VAR_(data%id_canopy_blockage) = zero_
            _DIAG_VAR_(data%id_canopy_frarea) = zero_
            _DIAG_VAR_(data%id_kemac) = zero_
         END DO

      _DIAG_VAR_S_(data%id_epib) = zero_

      ! Canopy height in this column
      hgt_canopy = _DIAG_VAR_S_(data%id_canopy_height)
      IF (hgt_canopy<0.005) THEN
         _DIAG_VAR_S_(data%id_d_par) = MAX( MIN( par * exp(-extc*( MAX(dz-0.005,zero_))),Io), zero_)
         print *,'plants too small: hgt_canopy= ',hgt_canopy
         RETURN ! plants too small, nothing to do here
      ENDIF

      ! Work up through the layers
      hgt_layer_top = zero_
      DO layer = SIZE(layer_map),1,-1

       layer_idx = layer_map(layer)
       dz  = _STATE_VAR_(data%id_dz)
       hgt_layer_top = hgt_layer_top + dz

       canopy_par = par
       canopy_frac = zero_

       ! If canopy height < layer thickness, then all is in bottom layer
       IF (hgt_canopy<=hgt_layer_top)THEN
         ! canopy ends in this layer
         hgt_layer_bot = hgt_layer_top - dz
         layer_frac = (hgt_canopy-hgt_layer_bot)/dz   ! vertical fraction canopy occupies in this layer
         canopy_frac=  (hgt_canopy-hgt_layer_bot)/hgt_canopy ! fraction of total canopy in this layer
         canopy_leaf_area = layer_frac * _DIAG_VAR_S_(data%id_canopy_lai) / dz ! m2 leaf / m3 water
         canopy_par = MAX( MIN( par * exp(-extc*( MAX(dz-(hgt_canopy-hgt_layer_bot),zero_))),Io), zero_)
         _DIAG_VAR_S_(data%id_d_par) = canopy_par
          done = 1
       ELSE
         ! canopy fully spans through layer
         layer_frac = one_
         canopy_leaf_area = layer_frac * _DIAG_VAR_S_(data%id_canopy_lai) / dz ! m2 leaf / m3 water
         canopy_par = par
         canopy_frac = dz/hgt_canopy
       ENDIF


       ! Set layer's macrophyte biovolume/blockage, and frontal area (e.g. for layer specific drag)
       _DIAG_VAR_(data%id_canopy_blockage) = _DIAG_VAR_S_(data%id_canopy_biovolume) * canopy_frac
       _DIAG_VAR_(data%id_canopy_frarea) = _DIAG_VAR_S_(data%id_canopy_sh_dens) * _DIAG_VAR_S_(data%id_canopy_sh_diam) * (dz*layer_frac)
       _DIAG_VAR_(data%id_kemac) = _DIAG_VAR_(data%id_canopy_blockage) * data%mpars(mi)%KeMAC ! or maybe use Aeff?

          ! Allow epiphyte growth in this layer
          IF( data%epi_model >2 .AND. _DIAG_VAR_(data%id_canopy_blockage) > zero_ ) THEN

         epi_prod = zero_ ; epi_resp = zero_

         ! Compute photosynthesis and respiration
         epi =  _STATE_VAR_S_(data%id_epi) * canopy_frac ! biomass of epiphytes in this layer
         fI = photosynthesis_irradiance(10,data%I_Kepi,data%I_Kepi,canopy_par,extc,Io,dz)
         epi_prod = data%R_epig*fI*(data%theta_epi_growth**(temp-20.))*(1.-(MIN(epi,data%epi_max*canopy_frac)/data%epi_max*canopy_frac))
         epi_resp = (data%R_epir*(data%theta_epi_resp**(temp-20.)))
         epi_flux = (epi_prod-epi_resp)*epi * (canopy_leaf_area/dz)

         ! Increment this layers epi productivity into the "bulk" epi (benthic) pool
         _FLUX_VAR_B_(data%id_epi) = _FLUX_VAR_B_(data%id_epi) + epi_flux

         ! Update epi diagnostics
         IF (diag_level>0) _DIAG_VAR_S_(data%id_epib) = _DIAG_VAR_S_(data%id_epib) + epi
         IF (diag_level>1) _DIAG_VAR_S_(data%id_epig) = _DIAG_VAR_S_(data%id_epig) + epi_prod *epi*(canopy_leaf_area/dz)
         IF (diag_level>1) _DIAG_VAR_S_(data%id_epir) = _DIAG_VAR_S_(data%id_epir) + epi_resp *epi*(canopy_leaf_area/dz)

         ! Add/remove any metabolism products to this water layer
         IF( data%simMacFeedback ) THEN
           ! Update flux terms for O2 (mmol O2/m2/s)
           _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + epi_flux
           IF (data%do_Cuptake ) THEN
             ! Update flux terms for CO2 (mmol C/m2/s)
             _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - epi_flux
           ENDIF
           IF (data%do_Nuptake ) THEN
             ! Update flux terms for nitrogen (mmol N/m2/s)
             ratio = data%epi_Xnc
             _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) - epi_flux * ratio * 0.5
             _FLUX_VAR_(data%id_nh4) = _FLUX_VAR_(data%id_nh4) - epi_flux * ratio * 0.5
             ! log this uptake into the bulk community N uptake diagnostic (mmol N/m3/d)
             !   _DIAG_VAR_(data%id_NUP1)= _DIAG_VAR_(data%id_NUP1)- (epi_flux) * ratio *0.5 * secs_per_day
             !   _DIAG_VAR_(data%id_NUP2)= _DIAG_VAR_(data%id_NUP2)- (epi_flux) * ratio *0.5 * secs_per_day
           ENDIF
           IF (data%do_Puptake ) THEN
             ! Update flux terms for phosphate (mmol P/m2/s)
             ratio = data%epi_Xpc
             _FLUX_VAR_(data%id_po4) = _FLUX_VAR_(data%id_po4) - epi_flux * ratio
             ! log this uptake into the bulk community P uptake diagnostic (mmol P/m3/d)
             !   _DIAG_VAR_(data%id_PUP) = _DIAG_VAR_(data%id_PUP) - mpb_flux * ratio * secs_per_day
           ENDIF

           ENDIF
        ENDIF

        IF(done==1) EXIT ! top of canopy was found so stop looping up
     END DO
     ! End vertical/column loop

     print *,'canopy_par',_DIAG_VAR_S_(data%id_d_par)
     ! UPDATE CANOPY POSTURE & DRAG
     ! Based on cell velocity and canopy geometry, compute in-canopy velocity
!     CALL sav_canopy_velocity( U , canopy_height , dz , u1 , u2 , drag )
     ! Use in canopy velocity to update canopy posture, and drag
!     CALL sav_posture_model( u1 , canopy_height , drag)

   ENDIF
END SUBROUTINE aed_calculate_column_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_calculate_benthic_macrophyte(data,column,layer_idx)
!-------------------------------------------------------------------------------
! Calculate submerged macrophyte production and respiration etc
! Everything in units per surface area (not volume!) per time.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
!
!LOCALS
   INTEGER  :: mi
   AED_REAL :: mphy, MAC_A, MAC_B, MAC_F       ! State
   AED_REAL :: mphy_flux, mphy_flux_a, mphy_flux_b
   AED_REAL :: fT, fI, fSal, fDO, fEpi
   AED_REAL :: extc, dz, par, Io, temp, salinity, day, matz

   AED_REAL :: canopy_hgt, canopy_extc, par_canopy, par_swi

   AED_REAL :: sh_diameter, sh_height, n_shoot, biovolume

   AED_REAL :: primprod(data%num_mphy), respiration(data%num_mphy)

   AED_REAL :: landaimid, WLmid
   AED_REAL :: gpp, npp
   AED_REAL :: resp, respiration_A, respiration_B, R_resp_A, R_resp_B, theta_resp, k_resp
   AED_REAL :: mort_A, mort_B, R_mort_A, R_mort_B
   AED_REAL :: f_tran, tau_tran, f_tran_fruit, tau_tran_fruit, f1, f2, f_seed, f_below
   AED_REAL :: Omega_MAC, sine_blade, E_comp,  R_growth, I_K, kA  , A_eff
   AED_REAL :: t_start_g, t_dur_g,t_max_g, t_max_r, t_start_r, t_dur_r, light_int, lterm1, factor, kI, term1, term2, factor2, npp0, x, tmp1, tmp2, tmp12, tmp22
   AED_REAL :: MAC_F_release, r_release, f_release, trigger_fruit_growth
   AED_REAL :: dw, rho_v
   INTEGER ::  ll
   AED_REAL :: epi, epi_prod, epi_resp, epi_flux, epi_max, ratio, din, dip, fN, fP, K_epi, leaf_area
!
!-------------------------------------------------------------------------------
!BEGIN

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)

   IF ( .NOT. in_zone_set(matz, data%active_zones) ) return

   A_eff = 0. ; f_tran = 0. ; mphy_flux_a = 0. ; mphy_flux_b = 0.

   !--- BENTHIC CELL SETUP TASKS
   ! Retrieve current (local) environmental conditions in benthic cell
   salinity = _STATE_VAR_(data%id_sal)         ! salinity
   temp  = _STATE_VAR_(data%id_tem)            ! temperature
   extc = _STATE_VAR_(data%id_extc)            ! extinction
   par = _STATE_VAR_(data%id_par)              ! photosynthetically active radiation
   dz = _STATE_VAR_(data%id_dz)                ! cell thickness
   Io = _STATE_VAR_S_(data%id_I_0)             ! surface short wave radiation

   ! Adjust PAR to top of canopy
   !canopy_extc = MAX( 0.01,_DIAG_VAR_(data%id_kemac))       !0.1
   !canopy_hgt = _DIAG_VAR_S_(data%id_canopy_height)         !0.1 ;
   !par_canopy = _DIAG_VAR_S_(data%id_d_par) ! MAX( MIN( par * exp(-extc*( MAX(dz-canopy_hgt,zero_))),Io), zero_)
   canopy_extc = 0.1
   canopy_hgt = 0.1
   par_canopy = MAX( MIN( par * exp(-extc*( MAX(dz-canopy_hgt,zero_))),Io), zero_)
   par_swi = MAX( par_canopy * exp(-(extc+canopy_extc)*canopy_hgt), zero_)

   ! Initialise epiphyte density (for macrophyte fI limitation)
   epi = zero_
   IF ( data%simEpiphytes ) epi = _STATE_VAR_S_(data%id_epi)      
   leaf_area = (3.142*(_DIAG_VAR_S_(data%id_canopy_sh_diam))*_DIAG_VAR_S_(data%id_canopy_height))*_DIAG_VAR_S_(data%id_canopy_sh_dens)
   epi = epi / MAX(leaf_area, 1e-3)  ! ( mmolC epiphytes/m2 leaf = mmolC epiphytes /m2 benthos / m2leaf/m2 benthos )

   ! Reset cumulative biomass/canopy diagnostics (used to sum over all groups)
   _DIAG_VAR_S_(data%id_canopy_biovolume) = zero_
   _DIAG_VAR_S_(data%id_canopy_sh_dens) = zero_
   _DIAG_VAR_S_(data%id_canopy_sh_diam) = zero_
   _DIAG_VAR_S_(data%id_canopy_height) = zero_
   _DIAG_VAR_S_(data%id_canopy_lai) = zero_
   _DIAG_VAR_S_(data%id_mac_ag) = zero_
   _DIAG_VAR_S_(data%id_mac_bg) = zero_
   _DIAG_VAR_S_(data%id_mac_tr) = zero_
   _DIAG_VAR_S_(data%id_mac) = zero_
   _DIAG_VAR_S_(data%id_gpp) = zero_

   !IF (data%mpars(mi)%fruit_model == 1) THEN
   IF( data%simFruiting ) THEN
      _DIAG_VAR_S_(data%id_mac_fr) = zero_
      _DIAG_VAR_S_(data%id_mac_ft) = zero_
   ENDIF

   ! Initialise fruit growth control variables
   trigger_fruit_growth = _DIAG_VAR_S_(data%id_mac_gt)    ! =0 when start releasing
   MAC_F = zero_                                          ! default fruit mass


   !--- LOOP THROUGH EACH GROUP/SPECIES
   ! Compute productivity, respiration, translocation, fruiting
   DO mi=1,data%num_mphy

      ! Check if this species/group is locked to grow in a specific zone
      IF( data%mpars(mi)%zone_lock >0 .and. data%mpars(mi)%zone_lock .ne. INT(matz) ) CYCLE

      ! Set local pars for species
      R_growth  = data%mpars(mi)%R_growth
      I_K       = data%mpars(mi)%I_K
      Omega_MAC = data%mpars(mi)%k_omega
      theta_resp= data%mpars(mi)%theta_resp
      f_below   = data%mpars(mi)%f_bg
      tau_tran  = data%mpars(mi)%tau_tran
      sine_blade= data%mpars(mi)%sine_blade
      E_comp    = data%mpars(mi)%E_comp
      kA        = data%mpars(mi)%kA
      R_mort_A  = data%mpars(mi)%R_mort
      R_mort_B  = data%mpars(mi)%R_mort_bg
      R_resp_A  = data%mpars(mi)%R_resp
      R_resp_B  = data%mpars(mi)%R_resp_bg
      K_epi     = data%mpars(mi)%K_epi

      !  Retrieve current (local cell) biomass for AG, BG pools
      MAC_A = _STATE_VAR_S_(data%id_mphya(mi))            ! leaf biomass
      MAC_B = _STATE_VAR_S_(data%id_mphyb(mi))            ! below-ground biomass

      !--- GROWTH
      IF ( data%mpars(mi)%light_model == 1 ) THEN
         ! use Baird approach with total/bulk light approach

         ! gross primary production
         fT           = 1.0                               ! temperature limitation
         fSal         = 1.0                               ! salinity limitation
         fEpi         = 1.0                               ! epiphyte limitation
         IF ( data%simEpiphytes ) fEpi = epi/(epi+K_epi)  ! epiphyte shading
         x            = (fEpi*par_canopy)/I_K             ! light limitation
         A_eff        = 1 - exp(-Omega_MAC*MAC_A)         ! effective area
         fI           = x/(1 + x) * (kA/(kA+A_eff))       ! light limitation+shelf shading
         primprod(mi) = R_growth*(min(fI,one_)) * MAC_A   ! GPP rate (mmol C/m2/s)

         ! respiration; note the definition of respiration in total light
         ! model is different to the one in the spectral light model. Here
         ! the respiration includes pure respiratory fraction, mortality and excretion
         respiration_A = R_resp_A * theta_resp**(temp-20.0) * MAC_A
         respiration_B = R_resp_B * theta_resp**(temp-20.0) * MAC_B

         respiration(mi) = respiration_A + respiration_B  ! Resp rate (mmol C/m2/s)

         ! translocation
         f_tran = (f_below-(MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! net flux rates (mmol C/m2/s)
         mphy_flux_a = (primprod(mi)  - respiration_A - f_tran)
         mphy_flux_b = (- respiration_B + f_tran)
         npp = primprod(mi)  - respiration_A


      ELSEIF ( data%mpars(mi)%light_model == 2 ) THEN
         ! use Baird approach with assumed (fixed) light spectra

         ! Integrate over 300-800nm to calculate photons being captured (Eq 8 in Baird et al 2016)
         fEpi         = 1.0                               ! epiphyte limitation
         IF ( data%simEpiphytes ) fEpi = 1.0              ! epiphyte shading
         light_int = zero_
         DO ll = 1,SIZE(WLint)-1
             lterm1    = 1. - exp(-ALlint(ll)*Omega_MAC*MAC_A*sine_blade)
             landaimid = (landaiint(ll)+landaiint(ll+1))/2.
             WLmid     = (WLint(ll)+WLint(ll+1))/2.
             light_int = light_int + landaimid*WLmid*(WLint(ll+1)-WLint(ll))*lterm1
         ENDDO
         light_int = light_int*(par_canopy*fEpi)/landaiSUM ! proportion of incoming par to clear-sky irradiance
         factor    = 1./(h*c*Av*1e9)                      ! conversion constant of photons from W/m2 to photon/m2/s
         kI        = factor*light_int                     ! rate of photon capture, mol photon/m2/s;
         primprod(mi) = kI  !?

         ! respiration
         ! compensation light (sum of respiration and mortality), mol photon/m2/day,
         term1 = E_comp*ALl*Omega_MAC*sine_blade

         term2 = 5500./550./1000. * R_mort_A * 86400.       ! mortality, converted to mol photon/m2/day;
         k_resp = 2.*(term1 - term2) * MAC_A /86400.      ! respiration rate in mol photon/m2/s, Eq(9) of Baird et al. 2016

         ! net production
         factor2 = 550./5500.*1000.                       ! factor to convert photon to carbon, mmol C/m2/s
         resp = k_resp * factor2                          ! respiration rate in mmol C/m2/s
         respiration(mi) = resp
         npp0 = max(zero_,(kI*factor2-resp))              ! net production rate
         npp = min(R_growth*MAC_A,npp0)                   ! cross-check of NPP npp0; %

         ! mortality
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A
         mort_B = R_mort_B * theta_resp**(temp-20.0) * MAC_B

         ! translocation between AG/BG
         f_tran = (f_below - (MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! update AG and BG biomass rates, and effective projected area fraction, A_eff
         mphy_flux_a = (npp - mort_A - f_tran)
         mphy_flux_b = (- mort_B + f_tran)

         A_eff = 1 - exp(-Omega_MAC * MAC_A)
         _DIAG_VAR_S_(data%id_aeff) = A_eff

         !PRINT *,'LIGHT2',npp, term1, term2, kI, k_resp !, ALl, E_comp, Omega_MAC, sine_blade


      ELSEIF ( data%mpars(mi)%light_model == 3 ) THEN
         ! use Baird approach with AED spectrally-resolved light

         light_int = _DIAG_VAR_S_(data%id_parc)           ! light capture by plant leaves (returned from OASIM)
         factor    = 1./(h*c*Av*1e9)                      ! conversion constant of photons from W/m2 to photon/m2/s
         kI        = factor*light_int                     ! rate of photon capture, mol photon/m2/s;
         primprod(mi) = kI  !?

         ! respiration
         ! compensation light (sum of respiration and mortality), mol photon/m2/day, divided by 86400 to mol photon/m2/s
         term1 = E_comp*ALl*Omega_MAC*sine_blade/86400.

         term2 = 5500./550./1000. * R_mort_A              ! mortality, converted to mol photon/m2/s;
         k_resp = 2.*(term1 - term2) * MAC_A              ! respiration rate in mol photon/m2/s, Eq(9) of Baird et al. 2016

         ! net production
         factor2 = 550./5500.*1000.                       ! factor to convert photon to carbon, mmol C/m2/s
         resp = k_resp * factor2                          ! respiration rate in mmol C/m2/s
         respiration(mi) = resp
         npp0 = max(zero_,(kI*factor2-resp))              ! net production rate
         npp = min(R_growth*MAC_A,npp0)                   ! cross-check of NPP npp0; %

         ! mortality
         mort_A = R_mort_A * theta_resp**(temp-20.0) * MAC_A
         mort_B = R_mort_B * theta_resp**(temp-20.0) * MAC_B

         ! translocation between AG/BG
         f_tran = (f_below - (MAC_B)/(MAC_A+MAC_B)) * (MAC_A+MAC_B) * tau_tran

         ! update AG and BG biomass rates, and effective projected area fraction, A_eff
         mphy_flux_a = (npp - mort_A - f_tran)
         mphy_flux_b = (- mort_B + f_tran)

         A_eff = 1 - exp(-Omega_MAC * MAC_A)
         _DIAG_VAR_S_(data%id_aeff) = A_eff


     ELSEIF ( data%mpars(mi)%light_model == 10 ) THEN
         ! use legacy AED approach with bulk light

         fT   = 1.0   ! fTemp_function(temp,     )
         fSal = 1.0   ! fSal_function(salinity,10.,72.,123.,230.)
         fI   = photosynthesis_irradiance(data%mpars(mi)%light_model, &
                                 I_K, data%mpars(mi)%I_S, par_canopy, extc, Io, dz)

         ! Primary productivity
         primprod(mi) = R_growth * (fI * fT * fSal) *  MAC_A

         ! Respiration and general metabolic loss
         respiration(mi) = bio_respiration(R_resp_A, theta_resp, temp)

         ! Salinity stress effect on respiration
         fSal = 1.0 ! phyto_salinity(data%mpars,mi,salinity)
         fDO  = 1.0 ! phyto_oxygen(data%mpars,mi,oxygen)

         respiration(mi) = respiration(mi) * (fSal * fDO) *  MAC_A

         mphy_flux_a = (primprod(mi) - respiration(mi))
         mphy_flux_b = zero_
         npp = mphy_flux_a

     ENDIF


     !--- SET GROWTH/RESPIRATION RATES
     IF( .NOT. data%simStaticBiomass ) THEN
          ! Set bottom fluxes for the benthic pools (mmol C/m2/s)
          _FLUX_VAR_B_(data%id_mphya(mi)) = _FLUX_VAR_B_(data%id_mphya(mi)) + mphy_flux_a
          _FLUX_VAR_B_(data%id_mphyb(mi)) = _FLUX_VAR_B_(data%id_mphyb(mi)) + mphy_flux_b
     ENDIF


     !--- FRUITING : seagrass fruit/seed development, & release trigger
     f_tran_fruit = zero_
     IF (data%mpars(mi)%fruit_model == 1) THEN

        f_seed         = data%mpars(mi)%f_seed
        tau_tran_fruit = data%mpars(mi)%tau_tran_fruit
        r_release      = data%mpars(mi)%r_release
        t_start_g      = data%mpars(mi)%t_start_g
        t_dur_g        = data%mpars(mi)%t_dur_g
        t_start_r      = data%mpars(mi)%t_start_r
        t_dur_r        = data%mpars(mi)%t_dur_r

          ! 1. get the day of the year, and prior fruit biomass
         MAC_F = max(zero_, _STATE_VAR_S_(data%id_mphyf(mi)))     ! fruits biomass
         day = FLOOR(_STATE_VAR_S_(data%id_yearday))  ! day number
         ! reset trigger for fruit to grow if the day<t_start_g

         IF (day<t_start_g) trigger_fruit_growth=one_

         ! 2. translocation from AG to fruit/seeds
         t_max_g = t_start_g+t_dur_g
         tmp1    = 12./t_dur_g*day+6.*(t_start_g+t_max_g)/(t_start_g-t_max_g)
         tmp2    = exp(-tmp1)
         f1      = 1./(1.+tmp2)

         f_tran_fruit = max(zero_,(f_seed-MAC_F/(MAC_F+MAC_A)) * (MAC_F+MAC_A)*tau_tran_fruit*f1)

         ! ...assume that if fruit ratio reach the 90% of maximum (mature) value,
         ! ...then stop fruit growth & start releasing at constant rate
         IF ( (MAC_F)/(MAC_F+MAC_A) > f_seed*0.9 ) THEN
            trigger_fruit_growth=zero_

         ENDIF

         ! 3. releasing
         t_max_r = t_start_r+t_dur_r
         tmp12   = 12. /t_dur_r*day + 6.*(t_start_r+t_max_r)/(t_start_r-t_max_r)
         tmp22   = exp(-tmp12)
         f2      = 1./(1.+tmp22)
         IF (trigger_fruit_growth>0.) THEN
           _FLUX_VAR_B_(data%id_mphyf(mi)) = _FLUX_VAR_B_(data%id_mphyf(mi)) + f_tran_fruit
           _FLUX_VAR_B_(data%id_mphya(mi)) = _FLUX_VAR_B_(data%id_mphya(mi)) - f_tran_fruit
         ELSE
            f_tran_fruit = zero_
            MAC_F_release = MAC_F*r_release
            f_release = MAC_F_release*f2
           ! _STATE_VAR_S_(data%id_mphyf(mi)) = max(zero_,MAC_F-f_release)

            _FLUX_VAR_B_(data%id_mphyf(mi)) = _FLUX_VAR_B_(data%id_mphyf(mi)) - f_release
           !_STATE_VAR_S_(data%id_mphyf(mi)) = max(zero_,_STATE_VAR_S_(data%id_mphyf(mi)))

            !MAC_F=max(0,MAC_F-f_release)
            !MAC_A(dd)=MAC_A(dd)   ! PH: is this missing something?
         ENDIF
         _DIAG_VAR_S_(data%id_mac_gt) = trigger_fruit_growth
         _DIAG_VAR_S_(data%id_mac_fr) = _DIAG_VAR_S_(data%id_mac_fr) + MAC_F
         _DIAG_VAR_S_(data%id_mac_ft) = _DIAG_VAR_S_(data%id_mac_ft) - f_tran_fruit

     ENDIF

     !--- FEEDBACK TO WATER COLUMN (NUTRIENTS, OXYGEN)
     IF( data%simMacFeedback .and. dz>0.05 ) THEN
        IF(data%id_oxy>0)&
          _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + mphy_flux_a/dz
        IF(data%id_dic>0)&
          _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - mphy_flux_a/dz

        ! ASSUMED NUTRIENT STOICHIOMETRY FOR NOW - NEED TO ADD SPECIES SPECIFIC VALUES FOR N and P
        IF(data%id_nox>0)&
          _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) - (mphy_flux_a/dz) * (16./106.) * 0.5 * data%water_nutrient_frac
        IF(data%id_nh4>0)&
          _FLUX_VAR_(data%id_nh4) = _FLUX_VAR_(data%id_nh4) - (mphy_flux_a/dz) * (16./106.) * 0.5 * data%water_nutrient_frac
        IF(data%id_po4>0)&
          _FLUX_VAR_(data%id_po4) = _FLUX_VAR_(data%id_po4) - (mphy_flux_a/dz) * (1./106.) * data%water_nutrient_frac
     ENDIF

     !--- CANOPY CHARACTERISTICS
     rho_v= 150e3  !!CHECK, unit in g DW/m3
     dw = MAC_A / data%mpars(mi)%X_cdw
     n_shoot = dw * data%mpars(mi)%X_sdw
     sh_height =  max(0.05, data%mpars(mi)%X_ldw * (dw/(dw+120.)))   ! set min height 5cm to avoid unrealistic large diameter
     sh_diameter = 2.* sqrt( (dw/rho_v/(3.142*sh_height*n_shoot)) )  
     biovolume = n_shoot* sh_height * 3.142*(sh_diameter/2.)**2.

     !--- ADD GROUP INFO TO TOTAL MAC COMMUNITY / CANOPY
     ! Increment updates to the plant community level diagnistics for this group
     _DIAG_VAR_S_(data%id_mac_ag) = _DIAG_VAR_S_(data%id_mac_ag) + MAC_A
     _DIAG_VAR_S_(data%id_mac_bg) = _DIAG_VAR_S_(data%id_mac_bg) + MAC_B
     _DIAG_VAR_S_(data%id_mac_tr) = _DIAG_VAR_S_(data%id_mac_tr) - f_tran
     _DIAG_VAR_S_(data%id_mac)    = _DIAG_VAR_S_(data%id_mac) + (MAC_A+MAC_B+MAC_F)
     _DIAG_VAR_S_(data%id_gpp)    = _DIAG_VAR_S_(data%id_gpp) + primprod(mi)*secs_per_day ! (mmol C/m2/day)

     ! Increment canopy characteristics
     _DIAG_VAR_S_(data%id_canopy_sh_dens) = _DIAG_VAR_S_(data%id_canopy_sh_dens) + n_shoot
     _DIAG_VAR_S_(data%id_canopy_sh_diam) = _DIAG_VAR_S_(data%id_canopy_sh_diam) + sh_diameter
     _DIAG_VAR_S_(data%id_canopy_biovolume) = _DIAG_VAR_S_(data%id_canopy_biovolume) + biovolume
     _DIAG_VAR_S_(data%id_canopy_height) = _DIAG_VAR_S_(data%id_canopy_height) + sh_height
     _DIAG_VAR_S_(data%id_canopy_lai) = _DIAG_VAR_S_(data%id_canopy_lai) + A_eff  ! TBC+       &
                     ! (one_ - exp(-data%mpars(mi)%k_omega * mphy*(one_-data%mpars(mi)%f_bg)))

   ENDDO ! Finish looping through groups

   !--- FINALISE CANOPY & COMMUNITY CALCULATIONS

   ! Finalise canopy averaging
   _DIAG_VAR_S_(data%id_canopy_sh_diam) = _DIAG_VAR_S_(data%id_canopy_sh_diam) / REAL(data%num_mphy)
   _DIAG_VAR_S_(data%id_canopy_height) = _DIAG_VAR_S_(data%id_canopy_height) / REAL(data%num_mphy)

   ! Approximate root depth and oxygen input (eg., for sediment bgc model)
   _DIAG_VAR_S_(data%id_root_o) = _DIAG_VAR_S_(data%id_gpp) * 0.2 ! data%bg_gpp_frac     ! mmolO2/m2/d
   _DIAG_VAR_S_(data%id_root_d) = 0.05 !MAX( MIN(_DIAG_VAR_S_(data%id_mac_ag) * data%coef_bm_hgt,0.25),0.01) !m

   ! Export additional diagnostic variables
   _DIAG_VAR_S_(data%id_d_par)= par_canopy






    ! BULK epiphyte model

    IF( data%epi_model == 1 ) THEN

      epi_prod = zero_ ; epi_resp = zero_

      ! Compute maximum epiphyte capacity, based on leaf area
      leaf_area = (3.142*(_DIAG_VAR_S_(data%id_canopy_sh_diam))*_DIAG_VAR_S_(data%id_canopy_height))*_DIAG_VAR_S_(data%id_canopy_sh_dens) 
      epi_max = data%epi_max * leaf_area  ! ( mmolC epiphytes /m2 benthos = mmolC epiphytes/m2 leaf * m2leaf/m2 benthos )

      ! Compute nutrient limitation 
      din  = _STATE_VAR_(data%id_nox)+_STATE_VAR_(data%id_nh4)
      dip  = _STATE_VAR_(data%id_po4)
      fN = MAX( MIN(  (din)/(din+data%epi_K_N) ,one_), zero_ )
      fP = MAX( MIN(  (dip)/(dip+data%epi_K_P) ,one_), zero_ )

      ! Compute photosynthesis and respiration
      epi =  _STATE_VAR_S_(data%id_epi)  ! biomass of epiphytes in the benthos
      fI = photosynthesis_irradiance(10,data%I_Kepi,data%I_Kepi,par_canopy,extc,Io,dz)
      epi_prod = data%R_epig*fI*(data%theta_epi_growth**(temp-20.)*MIN(fN,fP))
      epi_prod = epi_prod * (1.-(MIN(epi,epi_max)/epi_max))
      epi_resp = (data%R_epir*(data%theta_epi_resp**(temp-20.)))
      epi_flux = (epi_prod-epi_resp)*epi

      ! Increment this layers epi productivity into the "bulk" epi (benthic) pool
      _FLUX_VAR_B_(data%id_epi) = _FLUX_VAR_B_(data%id_epi) + epi_flux

      ! Update epi diagnostics
      IF (diag_level>0) _DIAG_VAR_S_(data%id_epib) = epi
      IF (diag_level>1) _DIAG_VAR_S_(data%id_epig) = _DIAG_VAR_S_(data%id_epig) + epi_prod *epi
      IF (diag_level>1) _DIAG_VAR_S_(data%id_epir) = _DIAG_VAR_S_(data%id_epir) + epi_resp *epi

      ! Add/remove any metabolism products to this water layer
      IF( data%simMacFeedback ) THEN
        ! Update flux terms for O2 (mmol O2/m2/s)

         _FLUX_VAR_(data%id_oxy) = _FLUX_VAR_(data%id_oxy) + epi_flux/dz
        IF (data%id_dic>0 ) THEN
          ! Update flux terms for CO2 (mmol C/m2/s)
          _FLUX_VAR_(data%id_dic) = _FLUX_VAR_(data%id_dic) - epi_flux/dz
        ENDIF
        IF (data%id_nox>0 ) THEN
          ! Update flux terms for nitrogen (mmol N/m2/s)
          ratio = data%epi_Xnc
          _FLUX_VAR_(data%id_nox) = _FLUX_VAR_(data%id_nox) - epi_flux * ratio * 0.5 /dz
          _FLUX_VAR_(data%id_nh4) = _FLUX_VAR_(data%id_nh4) - epi_flux * ratio * 0.5 /dz
          ! log this uptake into the bulk community N uptake diagnostic (mmol N/m3/d)
          !   _DIAG_VAR_(data%id_NUP1)= _DIAG_VAR_(data%id_NUP1)- (epi_flux) * ratio *0.5 * secs_per_day
          !   _DIAG_VAR_(data%id_NUP2)= _DIAG_VAR_(data%id_NUP2)- (epi_flux) * ratio *0.5 * secs_per_day
        ENDIF
        IF (data%id_po4>0 ) THEN
          ! Update flux terms for phosphate (mmol P/m2/s)
          ratio = data%epi_Xpc
          _FLUX_VAR_(data%id_po4) = _FLUX_VAR_(data%id_po4) - epi_flux * ratio /dz
          ! log this uptake into the bulk community P uptake diagnostic (mmol P/m3/d)
          !   _DIAG_VAR_(data%id_PUP) = _DIAG_VAR_(data%id_PUP) - mpb_flux * ratio * secs_per_day
        ENDIF

        ENDIF
     ENDIF


END SUBROUTINE aed_calculate_benthic_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_light_extinction_macrophyte(data, column, layer_idx, extinction)
!-------------------------------------------------------------------------------
! Get the light extinction coefficient due to macrophyte biomass
!
!  WARNING - THIS IS ADDDING MACROPHYTE EFFECT TO ALL CELLS IN WATER COLUMN
!
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: extinction
!
!LOCALS
   AED_REAL :: dz, mphy = 0., matz
   INTEGER  :: mi
   AED_REAL, PARAMETER :: max_extc = 2.0
!
!-------------------------------------------------------------------------------
!BEGIN
   extinction = zero_

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! Recall extinction value (for this cell) recorded during column loop calculations
   extinction = MIN( _DIAG_VAR_(data%id_kemac) , max_extc)

END SUBROUTINE aed_light_extinction_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed_bio_drag_macrophyte(data, column, layer_idx, drag)
!-------------------------------------------------------------------------------
! Get the effect of macrophyte biomass on benthic drag
!
!
!-------------------------------------------------------------------------------
   CLASS (aed_macrophyte_data_t),INTENT(in) :: data
   TYPE (aed_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: layer_idx
   AED_REAL,INTENT(inout) :: drag
  ! AED_REAL,INTENT(inout) :: blockage
!
!LOCALS
   INTEGER  :: mi, n_shoot
   AED_REAL :: sh_height, sh_diameter, m_biomass, blockage
   AED_REAL :: dz, matz, vel, depth, mphy
   AED_REAL :: Sx = one_, Sy = one_
   AED_REAL :: K_CD, Uc
!-------------------------------------------------------------------------------
!BEGIN

   drag = zero_

   ! Check this cell is in an active zone for macrophytes
   matz = _STATE_VAR_S_(data%id_sed_zone)
   if ( .NOT. in_zone_set(matz, data%active_zones) ) return

   ! get the necessary environment information - cell thickness, water velocity
   vel  = _STATE_VAR_(data%id_vel)
   dz   = _STATE_VAR_(data%id_dz)
   depth= _STATE_VAR_S_(data%id_depth)

   IF( data%drag_model==1 ) THEN

     DO mi=1,data%num_mphy
      !# above ground density of macrophyte group i
      mphy = _STATE_VAR_S_(data%id_mphya(mi))

      !# additional drag due to biomass
      drag = drag + (data%mpars(mi)%K_CD * (mphy/dz) )
     ENDDO

   ELSEIF( data%drag_model==2 ) THEN
     ! set canopy average conditions for drag computation
     n_shoot = INT(_DIAG_VAR_S_(data%id_canopy_stem_dens))
     sh_diameter = _DIAG_VAR_S_(data%id_canopy_stem_diam)
     sh_height = _DIAG_VAR_S_(data%id_canopy_height)

     ! now do computation of velocity and drag within canopy
     Uc = vegetation_drag(n_shoot, data%mac_pattern, &
                          sh_diameter, Sx, Sy, vel, sh_height, K_CD)
     _DIAG_VAR_S_(data%id_canopy_velocity) = Uc

     !# additional drag due to canopy
     drag = K_CD

   ELSEIF( data%drag_model==3 ) THEN
     ! Work out leaf effective length etc and buoyancy
   ENDIF

  _DIAG_VAR_S_(data%id_canopy_drag) = drag

END SUBROUTINE aed_bio_drag_macrophyte
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION vegetation_drag(n,vp,d,Sx,Sy,Ub,hv,Cd) RESULT(Uc)
   !-------------------------------------------------------------------------------
   ! Return the drag on the water flow given vegetation canopy characteristics
   !-------------------------------------------------------------------------------
   !#ARGUMENTS
      INTEGER, INTENT(in) :: n     !# Number of stems per unit area (1/m^2)
      INTEGER, INTENT(in) :: vp    !# Vegetation pattern/arrangement
      AED_REAL,INTENT(in) :: d     !# Mean stem diameter (m)
      AED_REAL,INTENT(in) :: Sx    !# Distance between elemetns in x-direction (streamwise direction) (m)
      AED_REAL,INTENT(in) :: Sy    !# Distance between elemetns in y-direction (lateral direction) (m)
      AED_REAL,INTENT(in) :: Ub    !# Bulk velocity that is the average cross-sectional velocity outside the canopy (m/s)
      AED_REAL,INTENT(in) :: hv    !# Mean height of vegetation elements
      AED_REAL,INTENT(inout) :: Cd !# Drag coefficient (computed here)
   !
   !#RETURN VAL
      AED_REAL :: Uc               !# Returned value of within-canopy velocity
   !
   !#LOCALS
      AED_REAL :: nu, Lambdap, Up
      AED_REAL :: nSqrt, Beta, Rec
   !
   !-------------------------------------------------------------------------------
   !Environmental parameters; Units: metric;

      nu = 10**(-6)                ! Kinematic viscosity of water
      Lambdap = n * 3.1418  * d**2 / 4  ! Planar area (dimensionless)
      Up = Ub / (1-Lambdap);       ! Pore velocity (m/s)

      IF (vp == 1) THEN
        !(a) If linear or staggered arrangement, use a constricted velocity
        Beta  = Sx / Sy
        Uc    = (1-Lambdap) * Up/(1-sqrt(4*Lambdap/(3.1418 *Beta)))
      ELSE
        !(b) If random array of stems, use the following relationship for Uc :
        nSqrt = sqrt(REAL(n))
        Uc    = (1-Lambdap) * Up/(1-d*nSqrt)
      ENDIF

      Rec = Uc * d / nu            ! Stem Reynolds number
      Cd  = 1 + Rec**(-2/3)        ! Drag coefficient

END FUNCTION vegetation_drag
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE interp_0d(nsource,x,y,ntarget,targetx,targety)
   !-------------------------------------------------------------------------------
   ! 0D interpolation, extrapolates beyond boundaries
   !-------------------------------------------------------------------------------
   !ARGUMENTS
      INTEGER,INTENT(in)                      :: nsource, ntarget
      AED_REAL,DIMENSION(nsource),INTENT(in)  :: x, y
      AED_REAL,DIMENSION(ntarget),INTENT(in)  :: targetx
      AED_REAL,DIMENSION(ntarget),INTENT(out) :: targety
   !
   !LOCALS
      INTEGER  :: i, j
      AED_REAL :: frac

   !-------------------------------------------------------------------------------
   !BEGIN
      i = 1
      DO j = 1,ntarget
         DO while (i+1<nsource)
            IF (x(i+1)>=targetx(j)) exit
            i = i+1
         ENDDO
         frac = (targetx(j)-x(i))/(x(i+1)-x(i))
         targety(j) = y(i) + frac*(y(i+1)-y(i))
      ENDDO
   END SUBROUTINE interp_0d
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE aed_macrophyte


#if 0



%% configuration for general parameters

% initial condition
MAC_A_INI = 25000;  % above-ground (AG, leaf without seeds) biomass, mmol C/m2
MAC_B_INI = 20000;  % below-ground (BG) biomass, mmol C/m2

% options
Light_model = 2;   % 1: total light model; 2: spectral-resolved model
Fruiting    = 1;   % 0: no fruiting; 1: include fruiting

% BGC parameters
R_growth = 0.4;    % maximum AG growth rate, /day
theta_resp=1.04;   % theta coefficient for respiration
T_standard=20;     % standard temperature for respiration
f_below=0.5;       % BG/(AG+BG) fraction
Omega_MAC = 1e-4;  %0.0011; % carbon-specific area of seagrass (mmol C/m2)^-1, converted from (240 mgC/m2)^-1
tau_tran = 0.033;  % AG/BG translocation rate, /day
sine_blade = 0.5;  % sine blade shape

% fruiting parameters
MAC_F_INI = 0;    % initial fruit biomass, mmol C/m2
tau_tran_fruit=0.01;  % translocation rate from leaves to fruit, /day
r_release=0.1;    % fruit releasing rate, /day (assume all fruit released in 10 days)

f_seed =0.4;       % seed biomass/(AG+seed biomass) fraction
t_start_g = 120;   % calenda day to grow, days of the year
t_dur_g = 30;      % duration to reach maximum growth, days

t_start_r = 210;   % calenda day to release fruits, days of the year
t_dur_r = 15;       % duration to reach maximum release, days




%% total light model parameters

I_K=80; % Half saturation constant for light limitation of growth
%I_S=120;  % Saturating light intensity for optimum photosynthesis
kA = 1;   % shelf shading effect

R_resp_A = 0.06;   % maximum AG respiration rate, /day
R_resp_B = 0.004;  % maximum BG respiration rate, /day

% fT = 1; % no temperature limiation
% fS = 1; % no salinity limiation
% fN = 1; % no sediment nutrient limiation

%% spectral-resolved parameters

WL_range = [300 800]; % wavelength range
WL_mean  = 550;       % mean wavelength, nm

h = 6.626e-34;  % Planck constant, Js;
c = 2.998e8;    % light speed, m/s;
Av= 6.02e23;    % Avagadro number, /mol, convert light capture to units of mol photon/m2/s;

E_comp = 15; %4.5*3.5; % compensation scalar PAR irradiance, mol photon/m2/d;
R_mort_A = 0.03;   % maximum AG mortality rate, /day
R_mort_B = 0.004;  % maximum BG mortality rate, /day


%% clear-sky wave lengths and irradiance at a particular wave-band W/m2/nm, from EMS model
% wave lengths
wavei = [140.00, 150.00, 160.00, 170.00, ...
    180.00, 190.00, 200.00, 205.00, 210.00, 215.00, 220.00, ...
    225.00, 230.00, 235.00, 240.00, 245.00, 250.00, 255.00, ...
    260.00, 265.00, 270.00, 275.00, 280.00, 285.00, 290.00, ...
    295.00, 300.00, 305.00, 310.00, 315.00, 320.00, 325.00, ...
    330.00, 335.00, 340.00, 345.00, 350.00, 355.00, 360.00, ...
    365.00, 370.00, 375.00, 380.00, 385.00, 390.00, 395.00, ...
    400.00, 405.00, 410.00, 415.00, 420.00, 425.00, 430.00, ...
    435.00, 440.00, 445.00, 450.00, 455.00, 460.00, 465.00, ...
    470.00, 475.00, 480.00, 485.00, 490.00, 495.00, 500.00, ...
    505.00, 510.00, 515.00, 520.00, 525.00, 530.00, 535.00, ...
    540.00, 545.00, 550.00, 555.00, 560.00, 565.00, 570.00, ...
    575.00, 580.00, 585.00, 590.00, 595.00, 600.00, 610.00, ...
    620.00, 630.00, 640.00, 650.00, 660.00, 670.00, 680.00, ...
    690.00, 700.00, 710.00, 720.00, 730.00, 740.00, 750.00, ...
    800.00, 850.00, 900.00, 950.00, 1000.00, 1100.00, 1200.00, ...
    1300.00, 1400.00, 1500.00, 1600.00, 1700.00, 1800.00, 1900.00, ...
    2000.00, 2100.00, 2200.00, 2300.00, 2400.00, 2500.00, 2600.00, ...
    2700.00, 2800.00, 2900.00, 3000.00, 3100.00, 3200.00, 3300.00, ...
    3400.00, 3500.00, 3600.00, 3700.00, 3800.00, 3900.00, 4000.00, ...
    4100.00, 4200.00, 4300.00, 4400.00, 4500.00, 4600.00, 4700.00, ...
    4800.00, 4900.00, 5000.00, 6000.00, 7000.00, 8000.00 ];

% clear-sky irradiance at a particular wave-band W/m2/nm, from EMS model
landai = [ 0.0000, 0.0001, 0.0000, 0.0004, ...
    0.0009,0.0017,0.0030,0.0050,0.0100,0.0180,0.0300,0.0420,0.0520, ...
    0.0540,0.0580,0.0640,0.0640,0.1000,0.1300,0.2000,0.2500,0.2200, ...
    0.2400,0.3400,0.5200,0.6300,0.6100,0.6700,0.7600,0.8200,0.8500, ...
    1.0200,1.1500,1.1100,1.1100,1.1700,1.1800,1.1600,1.1600,1.2900, ...
    1.3300,1.3200,1.2300,1.1500,1.1200,1.2000,1.5400,1.8800,1.9400, ...
    1.9200,1.9200,1.8900,1.7800,1.8200,2.0300,2.1500,2.2000,2.1900, ...
    2.1600,2.1500,2.1700,2.2000,2.1600,2.0300,1.9900,2.0400,1.9800, ...
    1.9700,1.9600,1.8900,1.8700,1.9200,1.9500,1.9700,1.9800,1.9800, ...
    1.9500,1.9200,1.9000,1.8900,1.8700,1.8700,1.8700,1.8500,1.8400, ...
    1.8300,1.8100,1.7700,1.7400,1.7000,1.6600,1.6200,1.5900,1.5500, ...
    1.5100,1.4800,1.4400,1.4100,1.3700,1.3400,1.3000,1.2700,1.1270, ...
    1.0030,0.8950,0.8030,0.7250,0.6060,0.5010,0.4060,0.3280,0.2670, ...
    0.2200,0.1820,0.1520,0.1274,0.1079,0.0917,0.0785,0.0676,0.0585, ...
    0.0509,0.0445,0.0390,0.0343,0.0303,0.0268,0.0230,0.0214,0.0191, ...
    0.0171,0.0153,0.0139,0.0125,0.0114,0.0103,0.0095,0.0087,0.0080, ...
    0.0073,0.0067,0.0061,0.0056,0.0051,0.0048,0.0044,0.0042,0.0021, ...
    0.0012,0.0006];

% calculate the integrated irradiance over the 300-800 nm wave lengths
WLint=wavei(27:103);              % 300-800 nm
landaiint=landai(27:103);         % corresponding irradiance at 300-800 nm

% total irradiance integrate over the 300-800 nm
landaiSUM=0;
for ww=1:length(WLint)-1
    landaimid = (landaiint(ww)+landaiint(ww+1))/2;
    WLmid = (WLint(ww)+WLint(ww+1))/2;
    landaiSUM=landaiSUM+landaimid*(WLint(ww+1)-WLint(ww));
end

%% leaf light absorbance at wave lengths, from Baird et al., 2016
Wli =[300  390  500  530  640  680  705  800];
ALli=[0.72 0.72 0.68 0.38 0.38 0.60 0.04 0.0];

% interpolate into fine-resolution wave lengths for further processing
ALlint=interp1(Wli,ALli,WLint);
ALl=mean(ALlint);  % overall leaf absorbance





%% This is a Matlab script for the seagrass dynamic model, only consider
%  the light limitation, e.g. no temperature, salinity, and sediment
%  nutrient limitation is included.
%
%  configurations: in Seagrass_model_config.m;
%  forcing data: surface PAR data at Cockburn Sound centre, exported
%                CSIEM model, then bottom PAR is calculated based on
%                pre-set water depth and extinction coefficient

clear; close all;

% load in configuration file
run('./Seagrass_model_config.m');

% load in light data
load('.\environmental_condition\light_CScenter.mat');
PAR_time = output_site.PAR_time;   % model times
timestep = 4/24;                   % timestpes
extc = 0.1;                        % extinction coefficient, /m
d = 5;                             % depth
par=output_site.PAR_surf.*exp(-extc*d);   % calculated bottom PAR


% allocate biomass for AG, BG, and Fruits
MAC_A = zeros(size(PAR_time));  % leaf biomass
MAC_B = zeros(size(PAR_time));  % below-ground biomass
MAC_F = zeros(size(PAR_time));  % fruits biomass

% allocate parameters for spectral-resoved model
npp = zeros(size(PAR_time));    % net production
resp = zeros(size(PAR_time));   % respiration of AG due to compensation light
mort_A = zeros(size(PAR_time)); % mortality of leaf
mort_B = zeros(size(PAR_time)); % mortality of roots
f_tran = zeros(size(PAR_time)); % translocation rates between AG/BG
A_eff = zeros(size(PAR_time));  % effective area

% allocate parameters for total light model
gpp = zeros(size(PAR_time));           % gross production
respiration_A = zeros(size(PAR_time)); % respiraton of leaf
respiration_B = zeros(size(PAR_time)); % respiration of roots

% allocate parameters for fruiting processes
f_tran_fruit = zeros(size(PAR_time));  % translocation rate of leaf to fruits
f1 = zeros(size(PAR_time));            % f1 function for growth
f2 = zeros(size(PAR_time));            % f2 function for release

% allocate parameters for spectral-resoved model

% initialize the biomass for AG, BG, and Fruits
MAC_A(1)=MAC_A_INI;
MAC_B(1)=MAC_B_INI;
MAC_F(1)=MAC_F_INI;

% a control on fruit growth, becomes 0 when start releasing
trigger_fruit_growth = 1;

% loop through the timesteps
for dd=2:length(PAR_time)

    % growth
    if Light_model == 2 % spectral-resolved model

        % integrated over 300-800nm to calculate photons being captured
        % Eq(8) of Baird et al. 2016
        light_int=0;
        for ll=1:length(WLint)-1
            lterm1 = 1 - exp(-ALlint(ll)*Omega_MAC.*MAC_A(dd-1)*sine_blade);
            landaimid = (landaiint(ll)+landaiint(ll+1))/2;
            WLmid = (WLint(ll)+WLint(ll+1))/2;

            light_int=light_int+landaimid*WLmid*(WLint(ll+1)-WLint(ll))*lterm1;
        end
        light_int=light_int*par(dd)/landaiSUM; %proportion of incoming irradiation par(dd) to clear-sky irradiance

        factor    = 1/(h*c*Av*1e9);    % conversion constant of photons from W/m2 to photon/m2/s
        kI = factor*light_int;         % rate of photon capture, mol photon/m2/s;

        % respiration

        term1 = E_comp*ALl*Omega_MAC.*sine_blade; % compensation light
        term2 = 5500/550/1000*R_mort_A;           % respiration, converted to mol photon/m2/s;
        k_resp = 2*(term1 - term2)*MAC_A(dd-1);   % respiration rate in photon, Eq(9) of Baird et al. 2016


        % net production
        factor2 = 550/5500*1000*86400; % factor to converting photon to carbon, mmol C/m2/day

        resp(dd) = k_resp*factor2/86400;        % respiration rate in mmol C/m2/day
        npp0 = max(0,(kI*factor2-resp(dd)));    % net production rate
        npp(dd) = min(R_growth*MAC_A(dd-1),npp0); % cross-check of NPP npp0; %

        % mortality
        mort_A(dd) = MAC_A(dd-1)*R_mort_A * theta_resp^(T_standard-20.0);
        mort_B(dd) = MAC_B(dd-1)*R_mort_B * theta_resp^(T_standard-20.0);

        % translocation between AG/BG
        f_tran(dd)=(f_below - (MAC_B(dd-1))/(MAC_A(dd-1)+MAC_B(dd-1)))*(MAC_A(dd-1)+MAC_B(dd-1))*tau_tran;

        % update AG and BG biomass, and effective projected area fraction
        % A_eff
        MAC_A(dd)= MAC_A(dd-1) + (npp(dd) - mort_A(dd) - f_tran(dd))*timestep;
        MAC_B(dd)= MAC_B(dd-1) + (- mort_B(dd) + f_tran(dd))*timestep;
        A_eff(dd) = 1 - exp(-Omega_MAC.*MAC_A(dd));

    elseif Light_model == 1 % total light model
        A_eff(dd) = 1 - exp(-Omega_MAC.*MAC_A(dd-1)); % effective area
        x = par(dd)/I_K;                             % light limitation
        fI = x ./ (1 + x)*(kA/(kA+A_eff(dd))); % light limitation+shelf shading
        gpp(dd) = MAC_A(dd-1)*R_growth*(min(fI,1)); % GPP rate per day

        % respiration; note the definition of respiration in total light
        % model is different to the one in the spectral light model. Here
        % the respiration includes pure respiratory fraction, mortality and excretion
        respiration_A(dd) = MAC_A(dd-1)*R_resp_A * theta_resp^(T_standard-20.0);
        respiration_B(dd) = MAC_B(dd-1)*R_resp_B * theta_resp^(T_standard-20.0);

        f_tran(dd)=(f_below - (MAC_B(dd-1))/(MAC_A(dd-1)+MAC_B(dd-1)))*(MAC_A(dd-1)+MAC_B(dd-1))*tau_tran;

        MAC_A(dd)= MAC_A(dd-1) + (gpp(dd) - respiration_A(dd) - f_tran(dd))*timestep;
        MAC_B(dd)= MAC_B(dd-1) + (- respiration_B(dd) + f_tran(dd))*timestep;

    else

        error('Light model option can be recognized');
    end

    % seagrass fruiting growth and release
    if Fruiting == 1

        % 1. get the day of the year
        t_vec = datevec(PAR_time(1));
        t = PAR_time(dd)-datenum(t_vec(1),1,1);

        % reset trigger for fruit to grow if t<t_start_g
        if t<t_start_g
            trigger_fruit_growth=1;
        end

        % 1. translocation from AG to seeds
        t_max_g=t_start_g+t_dur_g;
        tmp1=12/t_dur_g*t+6*(t_start_g+t_max_g)/(t_start_g-t_max_g);
        tmp2=exp(-tmp1);
        f1(dd)=1/(1+tmp2);

        f_tran_fruit(dd)=(f_seed-MAC_F(dd-1)/(MAC_F(dd-1)+MAC_A(dd-1)))*...
            (MAC_F(dd-1)+MAC_A(dd-1))*tau_tran_fruit*f1(dd);

        % assume if fruit ratio reach the 90% of maximum value, then
        % stop fruit growing and start releasing at constant rate
        if MAC_F(dd-1)/(MAC_F(dd-1)+MAC_A(dd-1))>f_seed*0.9
            trigger_fruit_growth=0;
            MAC_F_release=MAC_F(dd-1)*r_release;
        end

        % 2. releasing

        t_max_r=t_start_r+t_dur_r;
        tmp12=12/t_dur_r*t+6*(t_start_r+t_max_r)/(t_start_r-t_max_r);
        tmp22=exp(-tmp12);
        f2(dd)=1/(1+tmp22);
        if trigger_fruit_growth
            MAC_F(dd)=MAC_F(dd-1)+f_tran_fruit(dd);
            MAC_A(dd)=MAC_A(dd)-f_tran_fruit(dd);
        else
            f_tran_fruit(dd)=0;
            f_release=MAC_F_release*f2(dd);
            MAC_F(dd)=max(0,MAC_F(dd-1)-f_release);
            MAC_A(dd)=MAC_A(dd);
        end
    end



end

#endif
