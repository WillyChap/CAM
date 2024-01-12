module damlining

!=====================================================================
!
! Purpose: Implement DAMLining of the model state of U,V,T,Q, and/or PS
!          toward specified values from random analyses increments.
!
! Author: Wiliiam Chapman -- wchapman@ucar.edu
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    DAMLin the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the DAMLining
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the DAMLining module is set up to accomodate DAMLining of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The DAMLining of the model toward the analyses data is controlled by
!    the 'DAMLining_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which DAMLining is applied, the strength of the DAMLining
!    tendencies, and its spatial distribution.
!
!    FORCING:
!    --------
!    DAMLining tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'DAMLin_Force_Opt'
!    option, the timescale of the forcing is determined from the given
!    'DAMLin_TimeScale_Opt', and the DAMLining strength Alpha=[0.,1.] for each
!    variable is specified by the 'DAMLin_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_DAMLin = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied DAMLining can be limited using Horizontal/Vertical
!    window functions that are constructed using a parameterization of the
!    Heaviside step function.
!
!    The Heaviside window function is the product of separate horizonal and vertical
!    windows that are controled via 12 parameters:
!
!        DAMLin_Hwin_lat0:     Specify the horizontal center of the window in degrees.
!        DAMLin_Hwin_lon0:     The longitude must be in the range [0,360] and the
!                             latitude should be [-90,+90].
!        DAMLin_Hwin_latWidth: Specify the lat and lon widths of the window as positive
!        DAMLin_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999)
!                             renders the window a constant in that direction.
!        DAMLin_Hwin_latDelta: Controls the sharpness of the window transition with a
!        DAMLin_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step
!                             function while a large value yeilds a smoother transition.
!        DAMLin_Hwin_Invert  : A logical flag used to invert the horizontal window function
!                             to get its compliment.(e.g. to DAMLin outside a given window).
!
!        DAMLin_Vwin_Lindex:   In the vertical, the window is specified in terms of model
!        DAMLin_Vwin_Ldelta:   level indcies. The High and Low transition levels should
!        DAMLin_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also
!        DAMLin_Vwin_Hdelta:   specified in terms of model indices. For a window function
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition
!                             lengths should be set to 0.001
!        DAMLin_Vwin_Invert  : A logical flag used to invert the vertical window function
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent
!                 of the vertical (30 levels):
!                        DAMLin_Hwin_lat0     = 0.         DAMLin_Vwin_Lindex = 0.
!                        DAMLin_Hwin_latWidth = 30.        DAMLin_Vwin_Ldelta = 0.001
!                        DAMLin_Hwin_latDelta = 5.0        DAMLin_Vwin_Hindex = 31.
!                        DAMLin_Hwin_lon0     = 180.       DAMLin_Vwin_Hdelta = 0.001
!                        DAMLin_Hwin_lonWidth = 999.       DAMLin_Vwin_Invert = .false.
!                        DAMLin_Hwin_lonDelta = 1.0
!                        DAMLin_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply DAMLining at the poles and
!                 not at the equator, the settings would be similar but with:
!                        DAMLin_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before
!    running the model. Lookat_DAMLinWindow.ncl is a script avalable in the tools directory
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'DAMLin_U','DAMLin_V','DAMLin_T',and 'DAMLin_Q'}
!    The target values that the model state is DAMLind toward are available for history
!    file output via the variables:  {'Target_daml_U','Target_daml_V','Target_daml_T',and 'Target_daml_Q'}
!
!    &DAMLining_nl
!      DAMLin_Model         - LOGICAL toggle to activate DAMLining.
!                              TRUE  -> DAMLining is on.
!                              FALSE -> DAMLining is off.                            [DEFAULT]
!
!      DAMLin_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/DAMLining/ERAI-Data/')
!
!      DAMLin_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      DAMLin_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Model_Times_Per_Day_DAMLin - INT Number of times to update the model state (used for DAMLining)
!                                each day. The value is restricted to be longer than the
!                                current model timestep and shorter than the analyses
!                                timestep. As this number is increased, the DAMLining
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      DAMLin_Beg_Year      - INT DAMLining begining year.  [1979- ]
!      DAMLin_Beg_Month     - INT DAMLining begining month. [1-12]
!      DAMLin_Beg_Day       - INT DAMLining begining day.   [1-31]
!      DAMLin_End_Year      - INT DAMLining ending year.    [1979-]
!      DAMLin_End_Month     - INT DAMLining ending month.   [1-12]
!      DAMLin_End_Day       - INT DAMLining ending day.     [1-31]
!
!      DAMLin_Force_Opt     - INT Index to select the DAMLining Target for a relaxation
!                                forcing of the form:
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      DAMLin_TimeScale_Opt - INT Index to select the timescale for DAMLining.
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      DAMLin_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      DAMLin_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      DAMLin_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      DAMLin_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      DAMLin_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No DAMLining of this variable)
!                                         1 == CONSTANT (Spatially Uniform DAMLining)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      DAMLin_Ucoef         - REAL fractional DAMLining coeffcient for U.
!      DAMLin_Vcoef         - REAL fractional DAMLining coeffcient for V.
!      DAMLin_Tcoef         - REAL fractional DAMLining coeffcient for T.
!      DAMLin_Qcoef         - REAL fractional DAMLining coeffcient for Q.
!      DAMLin_PScoef        - REAL fractional DAMLining coeffcient for PS.
!
!                                 The strength of the DAMLining is specified as a fractional
!                                 coeffcient between [0,1].
!
!      DAMLin_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      DAMLin_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      DAMLin_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      DAMLin_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      DAMLin_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      DAMLin_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      DAMLin_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      DAMLin_Vwin_Lindex   - REAL LO model index of transition
!      DAMLin_Vwin_Hindex   - REAL HI model index of transition
!      DAMLin_Vwin_Ldelta   - REAL LO transition length
!      DAMLin_Vwin_Hdelta   - REAL HI transition length
!      DAMLin_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps DAMLining????
!
!=====================================================================
! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc, mstrid=>masterprocid, mpicom
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: DAMLin_Model,DAMLin_Model_Regress,DAMLin_ON
  public:: damlining_readnl
  public:: damlining_init
  public:: regress_daml_timestep_tend
  public:: regress_diurnal_daml_timestep_tend
  public:: damlining_timestep_init
  public:: damlining_timestep_tend
  public:: damlining_diurnal_timestep_tend
  
  private::DAMLining_readin_weights_all
  private::DAMLining_readin_weights
  private::DAMLining_readin_weights_diurnal
  private::DAMLining_readin_normdict
  private::DAMLining_set_PSprofile
  private::DAMLining_set_profile
  
  !private::calc_DryStaticEnergy

  ! DAMLining Parameters
  !--------------------
  logical          :: DAMLin_Model       =.false.
  logical          :: DAMLin_Model_Regress =.false.
  logical          :: DAMLin_ON          =.false.
  logical          :: DAMLin_Initialized =.false.
  character(len=cl):: DAMLin_Path
  character(len=cs):: DAMLin_File,DAMLin_File_Template
  character(len=cs):: DAMLin_File_Diurnal,DAMLin_File_Anomaly,DAMLin_File_Norm
  character(len=cs):: Ran_File
  integer          :: DAMLin_Force_Opt
  integer          :: DAMLin_TimeScale_Opt
  integer          :: DAMLin_TSmode
  integer          :: DAMLin_Times_Per_Day
  integer          :: Model_Times_Per_Day_DAMLin
  real(r8)         :: DAMLin_Ucoef,DAMLin_Vcoef
  integer          :: DAMLin_Uprof,DAMLin_Vprof
  real(r8)         :: DAMLin_Qcoef,DAMLin_Tcoef
  integer          :: DAMLin_Qprof,DAMLin_Tprof
  real(r8)         :: DAMLin_PScoef
  integer          :: DAMLin_PSprof
  integer          :: DAMLin_Beg_Year ,DAMLin_Beg_Month
  integer          :: DAMLin_Beg_Day  ,DAMLin_Beg_Sec
  integer          :: DAMLin_End_Year ,DAMLin_End_Month
  integer          :: DAMLin_End_Day  ,DAMLin_End_Sec
  integer          :: DAMLin_Curr_Year,DAMLin_Curr_Month
  integer          :: DAMLin_Curr_Day ,DAMLin_Curr_Sec
  integer          :: DAMLin_Next_Year,DAMLin_Next_Month
  integer          :: DAMLin_Next_Day ,DAMLin_Next_Sec
  integer          :: Random_Next_Year,Random_Next_Month !++WEC
  integer          :: Random_Next_Day ,Random_Next_Sec !++WEC
  integer          :: Memory_Rand !++WEC
  integer          :: DAMLin_Step
  integer          :: DAMLin_Climo_Year  = 1000
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: DAMLin_Hwin_lat0
  real(r8)         :: DAMLin_Hwin_latWidth
  real(r8)         :: DAMLin_Hwin_latDelta
  real(r8)         :: DAMLin_Hwin_lon0
  real(r8)         :: DAMLin_Hwin_lonWidth
  real(r8)         :: DAMLin_Hwin_lonDelta
  logical          :: DAMLin_Hwin_Invert = .false.
  real(r8)         :: DAMLin_Hwin_lo
  real(r8)         :: DAMLin_Hwin_hi
  real(r8)         :: DAMLin_Vwin_Hindex
  real(r8)         :: DAMLin_Vwin_Hdelta
  real(r8)         :: DAMLin_Vwin_Lindex
  real(r8)         :: DAMLin_Vwin_Ldelta
  logical          :: DAMLin_Vwin_Invert =.false.
  real(r8)         :: DAMLin_Vwin_lo
  real(r8)         :: DAMLin_Vwin_hi
  real(r8)         :: DAMLin_Hwin_latWidthH
  real(r8)         :: DAMLin_Hwin_lonWidthH
  real(r8)         :: DAMLin_Hwin_max
  real(r8)         :: DAMLin_Hwin_min
  integer          :: DAMLin_Do !adding control [WEC]
  integer          :: DAMLin_Pred !adding control [WEC]
  
  integer :: &
    upwp_clubb_idx, &          ! Cloud fraction
    taux_idx,       &
    pblh_idx,       &
    ustarwec_idx,   &
    flntc_idx,      &
    fsntoa_idx

  ! DAMLining State Arrays
  !-----------------------
  integer DAMLin_nlon,DAMLin_nlat,DAMLin_ncol,DAMLin_nlev,DAMLin_npredvar
  
  real(r8),allocatable:: DAMLin_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_PStau (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_PSstep(:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Ustep_avg (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Vstep_avg (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Sstep_avg (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_Qstep_avg (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: DAMLin_PSstep_avg(:,:)    !(pcols,begchunk:endchunk)

  ! DAMLining Observation Arrays and ML Weight Arrays
  !-----------------------------
  integer               DAMLin_NumObs
  integer,allocatable:: DAMLin_ObsInd(:)
  logical ,allocatable::DAMLin_File_Present(:)
  real(r8),allocatable::Weights_U (:,:,:,:,:) !(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_V (:,:,:,:,:) !(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_T (:,:,:,:,:) !(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_Q (:,:,:,:,:) !(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_PS(:,:,:)   !(pcols,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_U (:,:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_V (:,:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_T (:,:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_Q (:,:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_PS(:,:,:)   !(pcols,begchunk:endchunk,DAMLin_NumObs)gchunk:endchunk,:DAMLin_NumObs)=0._r8
  real(r8),allocatable::Weights_U_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_V_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_T_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_Q_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::Weights_PS_int(:,:,:)   !(pcols,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_U_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_V_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_T_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_Q_int (:,:,:,:) !(pcols,pver,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::WeightsDiurnal_PS_int(:,:,:)   !(pcols,begchunk:endchunk,DAMLin_NumObs)gchunk:endchunk,:DAMLin_NumObs)=0._r8
  real(r8),allocatable::mean_preds(:,:,:,:) !(Rpvar,pcols,begchunk:endchunk,DAMLin_NumObs)
  real(r8),allocatable::std_preds(:,:,:,:) !(Rpvar,pcols,begchunk:endchunk,DAMLin_NumObs)
  
contains



  subroutine damlining_readnl(nlfile)
   !
   ! DAMLinING_READNL: Initialize default values controlling the DAMLining
   !                 process. Then read namelist values to override
   !                 them.
   !===============================================================
   use ppgrid        ,only: pver
   use namelist_utils,only:find_group_name
   use units         ,only:getunit,freeunit
   
   !
   ! Arguments
   !-------------
   character(len=*),intent(in)::nlfile
   !
   ! Local Values
   !---------------
   integer ierr,unitn

   namelist /DAMLining_nl/ DAMLin_Model,DAMLin_Model_Regress,DAMLin_Path,  &
                         DAMLin_File_Template,DAMLin_File_Anomaly,         &
                         DAMLin_File_Diurnal,DAMLin_File_Norm, DAMLin_Force_Opt,&
                         DAMLin_TimeScale_Opt,                          &
                         DAMLin_Times_Per_Day,Model_Times_Per_Day_DAMLin,      &
                         DAMLin_Ucoef ,DAMLin_Uprof,                     &
                         DAMLin_Vcoef ,DAMLin_Vprof,                     &
                         DAMLin_Qcoef ,DAMLin_Qprof,                     &
                         DAMLin_Tcoef ,DAMLin_Tprof,                     &
                         DAMLin_PScoef,DAMLin_PSprof,                    &
                         DAMLin_Beg_Year,DAMLin_Beg_Month,DAMLin_Beg_Day, &
                         DAMLin_End_Year,DAMLin_End_Month,DAMLin_End_Day, &
                         DAMLin_Hwin_lat0,DAMLin_Hwin_lon0,              &
                         DAMLin_Hwin_latWidth,DAMLin_Hwin_lonWidth,      &
                         DAMLin_Hwin_latDelta,DAMLin_Hwin_lonDelta,      &
                         DAMLin_Hwin_Invert,                            &
                         DAMLin_Vwin_Lindex,DAMLin_Vwin_Hindex,          &
                         DAMLin_Vwin_Ldelta,DAMLin_Vwin_Hdelta,          &
                         DAMLin_Vwin_Invert

   ! DAMLining is NOT initialized yet, For now
   ! DAMLining will always begin/end at midnight.
   !--------------------------------------------
   DAMLin_Initialized =.false.
   DAMLin_ON          =.false.
   DAMLin_Beg_Sec=0
   DAMLin_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   DAMLin_Model         = .false.
   DAMLin_Model_Regress = .false.
   DAMLin_Path          = '/glade/scratch/wchapman/DA_ML/CESML_AI/Data/LINREG_weights/'
   DAMLin_File_Template = 'YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   DAMLin_File_Anomaly  = 'linreg_DiAn_anomaly_regression_weights_UVT.nc'
   DAMLin_File_Diurnal  = 'linreg_DiAn_diurnal_regression_weights_UVT.nc'
   DAMLin_File_Norm     = 'linreg_DiAn_diurnal_normalization_dictionary_UVT.nc'
   DAMLin_Force_Opt     = 0
   DAMLin_TimeScale_Opt = 0
   DAMLin_TSmode        = 0
   DAMLin_Times_Per_Day = 24
   Model_Times_Per_Day_DAMLin = 48
   DAMLin_Ucoef         = 1._r8
   DAMLin_Vcoef         = 1._r8
   DAMLin_Qcoef         = 1._r8
   DAMLin_Tcoef         = 1._r8
   DAMLin_PScoef        = 1._r8
   DAMLin_Uprof         = 1
   DAMLin_Vprof         = 1
   DAMLin_Qprof         = 1
   DAMLin_Tprof         = 1
   DAMLin_PSprof        = 1
   DAMLin_Beg_Year      = 1969
   DAMLin_Beg_Month     = 1
   DAMLin_Beg_Day       = 1
   DAMLin_End_Year      = 2008
   DAMLin_End_Month     = 9
   DAMLin_End_Day       = 1
   DAMLin_Hwin_lat0     = 0._r8
   DAMLin_Hwin_latWidth = 9999._r8
   DAMLin_Hwin_latDelta = 1.0_r8
   DAMLin_Hwin_lon0     = 180._r8
   DAMLin_Hwin_lonWidth = 9999._r8
   DAMLin_Hwin_lonDelta = 1.0_r8
   DAMLin_Hwin_Invert   = .false.
   DAMLin_Hwin_lo       = 0.0_r8
   DAMLin_Hwin_hi       = 1.0_r8
   DAMLin_Vwin_Hindex   = float(pver+1)
   DAMLin_Vwin_Hdelta   = 0.001_r8
   DAMLin_Vwin_Lindex   = 0.0_r8
   DAMLin_Vwin_Ldelta   = 0.001_r8
   DAMLin_Vwin_Invert   = .false.


   DAMLin_Vwin_lo       = 0.0_r8
   DAMLin_Vwin_hi       = 1.0_r8
   DAMLin_Do            = 1 !adding control to start ML [WEC]
   DAMLin_Pred          = 1 !adding control TODO WEC -- whether to include predictor field.
   Memory_Rand          = 5 !++WEC

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     write(iulog,*) 'Namelist file ',trim(nlfile)
     call find_group_name(unitn,'damlining_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,damlining_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('damlining_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(DAMLin_Hwin_Invert) then
     DAMLin_Hwin_lo = 1.0_r8
     DAMLin_Hwin_hi = 0.0_r8
   else
     DAMLin_Hwin_lo = 0.0_r8
     DAMLin_Hwin_hi = 1.0_r8
   endif

   if(DAMLin_Vwin_Invert) then
     DAMLin_Vwin_lo = 1.0_r8
     DAMLin_Vwin_hi = 0.0_r8
   else
     DAMLin_Vwin_lo = 0.0_r8
     DAMLin_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values
   !----------------------------------
   if((DAMLin_Hwin_lat0.lt.-90._r8).or.(DAMLin_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'DAMLinING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_lat0=',DAMLin_Hwin_lat0
     call endrun('damlining_readnl:: ERROR in namelist')
   endif

   if((DAMLin_Hwin_lon0.lt.0._r8).or.(DAMLin_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'DAMLinING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_lon0=',DAMLin_Hwin_lon0
     call endrun('damlining_readnl:: ERROR in namelist')
   endif

   if((DAMLin_Vwin_Lindex.gt.DAMLin_Vwin_Hindex)                         .or. &
      (DAMLin_Vwin_Hindex.gt.float(pver+1)).or.(DAMLin_Vwin_Hindex.lt.0._r8).or. &
      (DAMLin_Vwin_Lindex.gt.float(pver+1)).or.(DAMLin_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'DAMLinING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'DAMLinING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'DAMLinING: Lindex must be LE than Hindex'
     write(iulog,*) 'DAMLinING:  DAMLin_Vwin_Lindex=',DAMLin_Vwin_Lindex
     write(iulog,*) 'DAMLinING:  DAMLin_Vwin_Hindex=',DAMLin_Vwin_Hindex
     call endrun('damlining_readnl:: ERROR in namelist')
   endif

   if((DAMLin_Hwin_latDelta.le.0._r8).or.(DAMLin_Hwin_lonDelta.le.0._r8).or. &
      (DAMLin_Vwin_Hdelta  .le.0._r8).or.(DAMLin_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'DAMLinING: Window Deltas must be positive'
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_latDelta=',DAMLin_Hwin_latDelta
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_lonDelta=',DAMLin_Hwin_lonDelta
     write(iulog,*) 'DAMLinING:  DAMLin_Vwin_Hdelta=',DAMLin_Vwin_Hdelta
     write(iulog,*) 'DAMLinING:  DAMLin_Vwin_Ldelta=',DAMLin_Vwin_Ldelta
     call endrun('damlining_readnl:: ERROR in namelist')

   endif

   if((DAMLin_Hwin_latWidth.le.0._r8).or.(DAMLin_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'DAMLinING: Window widths must be positive'
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_latWidth=',DAMLin_Hwin_latWidth
     write(iulog,*) 'DAMLinING:  DAMLin_Hwin_lonWidth=',DAMLin_Hwin_lonWidth
     call endrun('damlining_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(DAMLin_Path         ,len(DAMLin_Path)         ,mpichar,0,mpicom)
   call mpibcast(DAMLin_File_Template,len(DAMLin_File_Template),mpichar,0,mpicom)
   call mpibcast(DAMLin_File_Diurnal,len(DAMLin_File_Diurnal),mpichar,0,mpicom)
   call mpibcast(DAMLin_File_Anomaly,len(DAMLin_File_Anomaly),mpichar,0,mpicom)
   call mpibcast(DAMLin_File_Norm,len(DAMLin_File_Norm),mpichar,0,mpicom)
   call mpibcast(DAMLin_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day_DAMLin, 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Vwin_Invert,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! damlining_readnl
  !================================================================

  !================================================================
  subroutine damlining_init
    !
    ! DAMLinING_INIT: Allocate space and initialize DAMLining values
    !===============================================================
    use ppgrid        ,only: pver,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec

    ! Local values
    !----------------
    integer  Year,Month,Day,Sec,Climo_Year
    integer  YMD1,YMD
    logical  After_Beg,Before_End
    integer  istat,lchnk,ncol,icol,ilev
    integer  hdim1_d,hdim2_d
    integer  dtime
    real(r8) rlat,rlon
    real(r8) Wprof(pver)
    real(r8) lonp,lon0,lonn,latp,lat0,latn
    real(r8) Val1_p,Val2_p,Val3_p,Val4_p
    real(r8) Val1_0,Val2_0,Val3_0,Val4_0
    real(r8) Val1_n,Val2_n,Val3_n,Val4_n
    integer               nn
    ! WEC TODO: make these into namelist variables: 
    integer, parameter :: pvar = 182
    integer, parameter :: Rpvar = 177
    integer, parameter :: Dpvar = 12

    ! Get the time step size
    !------------------------
    dtime = get_step_size()
    
    ! Allocate Space for DAMLining data arrays
    !-----------------------------------------
    
    !write(iulog,*) 'DAML 1'
    ! Allocate Space for spatial dependence of
    ! DAMLining Coefs and DAMLining Forcing.
    !-------------------------------------------
    allocate(DAMLin_Utau(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Utau',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Vtau',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Stau(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Stau',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Qtau',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_PStau(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_PStau',pcols*((endchunk-begchunk)+1))

    allocate(DAMLin_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Ustep',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Vstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Sstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Qstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_PSstep(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_PSstep',pcols*((endchunk-begchunk)+1))
    
    allocate(DAMLin_Ustep_avg(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Ustep_avg',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Vstep_avg(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Vstep_avg',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Sstep_avg(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Sstep_avg',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_Qstep_avg(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_Qstep_avg',pcols*pver*((endchunk-begchunk)+1))
    allocate(DAMLin_PSstep_avg(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'DAMLining_init','DAMLin_PSstep_avg',pcols*((endchunk-begchunk)+1))
    !write(iulog,*) 'DAML 2'

    ! Register output fields with the cam history module
     !-----------------------------------------------------
    call addfld( 'DAMLin_U',(/ 'lev' /),'A','m/s/s'  ,'U DAMLining Tendency')
    call addfld( 'DAMLin_V',(/ 'lev' /),'A','m/s/s'  ,'V DAMLining Tendency')
    call addfld( 'DAMLin_T',(/ 'lev' /),'A','K/s'    ,'T DAMLining Tendency')
    call addfld( 'DAMLin_Q',(/ 'lev' /),'A','kg/kg/s','Q DAMLining Tendency')
    
    call addfld( 'DAMLin_U_avg',(/ 'lev' /),'A','m/s/s'  ,'U DAMLining Tendency')
    call addfld( 'DAMLin_V_avg',(/ 'lev' /),'A','m/s/s'  ,'V DAMLining Tendency')
    call addfld( 'DAMLin_T_avg',(/ 'lev' /),'A','K/s'    ,'T DAMLining Tendency')
    call addfld( 'DAMLin_Q_avg',(/ 'lev' /),'A','kg/kg/s','Q DAMLining Tendency')

    !if(DAMLin_Do==1) then
    !   call addfld('Target_daml_U',(/ 'lev' /),'A','m/s/s'    ,'U DAMLining Target'  ) !need to edit units if we DAMLin ML 'm/s/s' [WEC]
    !   call addfld('Target_daml_V',(/ 'lev' /),'A','m/s/s'    ,'V DAMLining Target'  ) !need to edit units if we DAMLin ML 'm/s/s' [WEC]
    !   call addfld('Target_daml_T',(/ 'lev' /),'A','K/s'      ,'T DAMLining Target'  ) !need to edit units if we DAMLin ML 'K/s' [WEC
    !   call addfld('Target_daml_Q',(/ 'lev' /),'A','kg/kg/s'  ,'Q DAMLining Target  ') !need to edit units if we DAMLin ML 'kg/kg/s' [WEC]
    !else
    !   call addfld('Target_daml_U',(/ 'lev' /),'A','m/s'    ,'U DAMLining Target'  )
    !   call addfld('Target_daml_V',(/ 'lev' /),'A','m/s'    ,'V DAMLining Target'  )
    !   call addfld('Target_daml_T',(/ 'lev' /),'A','K'      ,'T DAMLining Target'  )
    !   call addfld('Target_daml_Q',(/ 'lev' /),'A','kg/kg'  ,'Q DAMLining Target  ')
    !endif
    !-----------------------------------------
    ! Values initialized only by masterproc
    !-----------------------------------------

    if(masterproc) then
      ! Set the Stepping intervals for Model and DAMLining values
      ! Ensure that the Model_Step is not smaller then one timestep
      !  and not larger then the DAMLin_Step.
      !--------------------------------------------------------
      Model_Step=86400/Model_Times_Per_Day_DAMLin
      DAMLin_Step=86400/DAMLin_Times_Per_Day
      if(Model_Step.lt.dtime) then
        write(iulog,*) ' '
        write(iulog,*) 'DAMLinING: Model_Step cannot be less than a model timestep'
        write(iulog,*) 'DAMLinING:  Setting Model_Step=dtime , dtime=',dtime
        write(iulog,*) ' '
        Model_Step=dtime
      endif
      if(Model_Step.gt.DAMLin_Step) then
        write(iulog,*) ' '
        write(iulog,*) 'DAMLinING: Model_Step cannot be more than DAMLin_Step'
        write(iulog,*) 'DAMLinING:  Setting Model_Step=DAMLin_Step, DAMLin_Step=',DAMLin_Step
        write(iulog,*) ' '
        Model_Step=DAMLin_Step
      endif
      ! Initialize column and level dimensions
      !--------------------------------------------------------
      call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
      DAMLin_nlon=hdim1_d
      DAMLin_nlat=hdim2_d
      DAMLin_ncol=hdim1_d*hdim2_d
      DAMLin_nlev=pver
      DAMLin_npredvar=pvar
      
      ! Check the time relative to the DAMLining window
      !------------------------------------------------
      call get_curr_date(Year,Month,Day,Sec)
      YMD=(Year*10000) + (Month*100) + Day
      YMD1=(DAMLin_Beg_Year*10000) + (DAMLin_Beg_Month*100) + DAMLin_Beg_Day
      call timemgr_time_ge(YMD1,DAMLin_Beg_Sec,         &
                            YMD ,Sec          ,After_Beg)
      YMD1=(DAMLin_End_Year*10000) + (DAMLin_End_Month*100) + DAMLin_End_Day
      call timemgr_time_ge(YMD ,Sec          ,          &
                            YMD1,DAMLin_End_Sec,Before_End)
      if((After_Beg).and.(Before_End)) then
        ! Set Time indicies so that the next call to
        ! timestep_init will initialize the data arrays.
        !--------------------------------------------
        Model_Next_Year =Year
        Model_Next_Month=Month
        Model_Next_Day  =Day
        Model_Next_Sec  =(Sec/Model_Step)*Model_Step
        DAMLin_Next_Year =Year
        DAMLin_Next_Month=Month
        DAMLin_Next_Day  =Day
        DAMLin_Next_Sec  =(Sec/DAMLin_Step)*DAMLin_Step

      elseif(.not.After_Beg) then
        ! Set Time indicies to DAMLining start,
        ! timestep_init will initialize the data arrays.
        !--------------------------------------------
        Model_Next_Year =DAMLin_Beg_Year
        Model_Next_Month=DAMLin_Beg_Month
        Model_Next_Day  =DAMLin_Beg_Day
        Model_Next_Sec  =DAMLin_Beg_Sec
        DAMLin_Next_Year =DAMLin_Beg_Year
        DAMLin_Next_Month=DAMLin_Beg_Month
        DAMLin_Next_Day  =DAMLin_Beg_Day
        DAMLin_Next_Sec  =DAMLin_Beg_Sec

      elseif(.not.Before_End) then
        ! DAMLining will never occur, so switch it off
        !--------------------------------------------
        DAMLin_Model=.false.
        DAMLin_Model_Regress=.false.
        
        DAMLin_ON   =.false.
        write(iulog,*) ' '
        write(iulog,*) 'DAMLinING: WARNING - DAMLining has been requested but it will'
        write(iulog,*) 'DAMLinING:           never occur for the given time values'
        write(iulog,*) ' '
      endif
       

      ! Initialize values for window function
      !----------------------------------------
      lonp= 180._r8
      lon0=   0._r8
      lonn=-180._r8
      latp=  90._r8-DAMLin_Hwin_lat0
      lat0=   0._r8
      latn= -90._r8-DAMLin_Hwin_lat0

      DAMLin_Hwin_lonWidthH=DAMLin_Hwin_lonWidth/2._r8
      DAMLin_Hwin_latWidthH=DAMLin_Hwin_latWidth/2._r8

      Val1_p=(1._r8+tanh((DAMLin_Hwin_lonWidthH+lonp)/DAMLin_Hwin_lonDelta))/2._r8
      Val2_p=(1._r8+tanh((DAMLin_Hwin_lonWidthH-lonp)/DAMLin_Hwin_lonDelta))/2._r8
      Val3_p=(1._r8+tanh((DAMLin_Hwin_latWidthH+latp)/DAMLin_Hwin_latDelta))/2._r8
      Val4_p=(1._r8+tanh((DAMLin_Hwin_latWidthH-latp)/DAMLin_Hwin_latDelta))/2_r8
      Val1_0=(1._r8+tanh((DAMLin_Hwin_lonWidthH+lon0)/DAMLin_Hwin_lonDelta))/2._r8
      Val2_0=(1._r8+tanh((DAMLin_Hwin_lonWidthH-lon0)/DAMLin_Hwin_lonDelta))/2._r8
      Val3_0=(1._r8+tanh((DAMLin_Hwin_latWidthH+lat0)/DAMLin_Hwin_latDelta))/2._r8
      Val4_0=(1._r8+tanh((DAMLin_Hwin_latWidthH-lat0)/DAMLin_Hwin_latDelta))/2._r8

      Val1_n=(1._r8+tanh((DAMLin_Hwin_lonWidthH+lonn)/DAMLin_Hwin_lonDelta))/2._r8
      Val2_n=(1._r8+tanh((DAMLin_Hwin_lonWidthH-lonn)/DAMLin_Hwin_lonDelta))/2._r8
      Val3_n=(1._r8+tanh((DAMLin_Hwin_latWidthH+latn)/DAMLin_Hwin_latDelta))/2._r8
      Val4_n=(1._r8+tanh((DAMLin_Hwin_latWidthH-latn)/DAMLin_Hwin_latDelta))/2._r8
       
      DAMLin_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
      DAMLin_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                          (Val1_p*Val2_p*Val3_p*Val4_p), &
                          (Val1_n*Val2_n*Val3_n*Val4_n), &
                          (Val1_n*Val2_n*Val3_p*Val4_p))

      ! Initialize number of DAMLining observation values to keep track of.
      ! Allocate and initialize observation indices
      !-----------------------------------------------------------------
      if((DAMLin_Force_Opt.ge.0).and.(DAMLin_Force_Opt.le.1)) then
        DAMLin_NumObs=2
      else
        ! Additional Options may need OBS values at more times.
        !------------------------------------------------------
        DAMLin_NumObs=2
        write(iulog,*) 'DAMLinING: Setting DAMLin_NumObs=2'
        write(iulog,*) 'DAMLinING: WARNING: Unknown DAMLin_Force_Opt=',DAMLin_Force_Opt
        call endrun('DAMLinING: Unknown Forcing Option')
      endif
      allocate(DAMLin_ObsInd(DAMLin_NumObs),stat=istat)
      call alloc_err(istat,'DAMLining_init','DAMLin_ObsInd',DAMLin_NumObs)
      allocate(DAMLin_File_Present(DAMLin_NumObs),stat=istat)
      call alloc_err(istat,'DAMLining_init','DAMLin_File_Present',DAMLin_NumObs)
      do nn=1,DAMLin_NumObs
        DAMLin_ObsInd(nn) = DAMLin_NumObs+1-nn
      end do
      DAMLin_File_Present(:)=.false.
       
      !write(iulog,*) 'DAML 8'

      ! Initialization is done,
      !--------------------------
      DAMLin_Initialized=.true.

      ! Check that this is a valid DYCORE model
      !------------------------------------------
      if((.not.dycore_is('UNSTRUCTURED')).and. &
         (.not.dycore_is('EUL')         ).and. &
         (.not.dycore_is('LR')          )      ) then
        call endrun('DAMLinING IS CURRENTLY ONLY CONFIGURED FOR FV')
      endif
       
      !write(iulog,*) 'DAML 9'
        
      ! Informational Output
      !---------------------------
      write(iulog,*) ' '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) '  MODEL DAMLinING INITIALIZED WITH THE FOLLOWING SETTINGS: '
      write(iulog,*) '---------------------------------------------------------'
      write(iulog,*) 'DAMLinING: DAMLin_Model=',DAMLin_Model
      write(iulog,*) 'DAMLinING: DAMLin_Model_Regress=',DAMLin_Model_Regress
      write(iulog,*) 'DAMLinING: DAMLin_Path=',DAMLin_Path
      write(iulog,*) 'DAMLinING: DAMLin_File_Template =',DAMLin_File_Template
      write(iulog,*) 'DAMLinING: Diurnal Weights File =',DAMLin_File_Diurnal !WEC
      write(iulog,*) 'DAMLinING: Anomaly Weights File =',DAMLin_File_Anomaly !WEC
      write(iulog,*) 'DAMLinING: Normalization Dictionary =',DAMLin_File_Norm !WEC
      write(iulog,*) 'DAMLinING: DAMLin_Force_Opt=',DAMLin_Force_Opt
      write(iulog,*) 'DAMLinING: DAMLin_TimeScale_Opt=',DAMLin_TimeScale_Opt
      write(iulog,*) 'DAMLinING: DAMLin_TSmode=',DAMLin_TSmode
      write(iulog,*) 'DAMLinING: DAMLin_Times_Per_Day=',DAMLin_Times_Per_Day
      write(iulog,*) 'DAMLinING: Model_Times_Per_Day_DAMLin=',Model_Times_Per_Day_DAMLin
      write(iulog,*) 'DAMLinING: DAMLin_Step=',DAMLin_Step
      write(iulog,*) 'DAMLinING: Model_Step=',Model_Step
      write(iulog,*) 'DAMLinING: DAMLin_Ucoef  =',DAMLin_Ucoef
      write(iulog,*) 'DAMLinING: DAMLin_Vcoef  =',DAMLin_Vcoef
      write(iulog,*) 'DAMLinING: DAMLin_Qcoef  =',DAMLin_Qcoef
      write(iulog,*) 'DAMLinING: DAMLin_Tcoef  =',DAMLin_Tcoef
      write(iulog,*) 'DAMLinING: DAMLin_PScoef =',DAMLin_PScoef
      write(iulog,*) 'DAMLinING: DAMLin_Uprof  =',DAMLin_Uprof
      write(iulog,*) 'DAMLinING: DAMLin_Vprof  =',DAMLin_Vprof
      write(iulog,*) 'DAMLinING: DAMLin_Qprof  =',DAMLin_Qprof
      write(iulog,*) 'DAMLinING: DAMLin_Tprof  =',DAMLin_Tprof
      write(iulog,*) 'DAMLinING: DAMLin_PSprof =',DAMLin_PSprof
      write(iulog,*) 'DAMLinING: DAMLin_Beg_Year =',DAMLin_Beg_Year
      write(iulog,*) 'DAMLinING: DAMLin_Beg_Month=',DAMLin_Beg_Month
      write(iulog,*) 'DAMLinING: DAMLin_Beg_Day  =',DAMLin_Beg_Day
      write(iulog,*) 'DAMLinING: DAMLin_End_Year =',DAMLin_End_Year
      write(iulog,*) 'DAMLinING: DAMLin_End_Month=',DAMLin_End_Month
      write(iulog,*) 'DAMLinING: DAMLin_End_Day  =',DAMLin_End_Day
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lat0     =',DAMLin_Hwin_lat0
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_latWidth =',DAMLin_Hwin_latWidth
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_latDelta =',DAMLin_Hwin_latDelta
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lon0     =',DAMLin_Hwin_lon0
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lonWidth =',DAMLin_Hwin_lonWidth
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lonDelta =',DAMLin_Hwin_lonDelta
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_Invert   =',DAMLin_Hwin_Invert
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lo       =',DAMLin_Hwin_lo
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_hi       =',DAMLin_Hwin_hi
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_Hindex   =',DAMLin_Vwin_Hindex
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_Hdelta   =',DAMLin_Vwin_Hdelta
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_Lindex   =',DAMLin_Vwin_Lindex
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_Ldelta   =',DAMLin_Vwin_Ldelta
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_Invert   =',DAMLin_Vwin_Invert
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_lo       =',DAMLin_Vwin_lo
      write(iulog,*) 'DAMLinING: DAMLin_Vwin_hi       =',DAMLin_Vwin_hi
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_latWidthH=',DAMLin_Hwin_latWidthH
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_lonWidthH=',DAMLin_Hwin_lonWidthH
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_max      =',DAMLin_Hwin_max
      write(iulog,*) 'DAMLinING: DAMLin_Hwin_min      =',DAMLin_Hwin_min
      write(iulog,*) 'DAMLinING: DAMLin_Initialized   =',DAMLin_Initialized
      write(iulog,*) ' '
      write(iulog,*) 'DAMLinING: DAMLin_NumObs=',DAMLin_NumObs
      write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Random_Next_Sec     ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Year    ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Month   ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Day     ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Memory_Rand         ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(DAMLin_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_Model_Regress ,            1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(DAMLin_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_npredvar       ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(DAMLin_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(DAMLin_NumObs        ,            1, mpiint, 0, mpicom)
#endif

    ! All non-masterproc processes also need to allocate space
    ! before the broadcast of DAMLin_NumObs dependent data.
    !------------------------------------------------------------
    if(.not.masterproc) then
      allocate(DAMLin_ObsInd(DAMLin_NumObs),stat=istat)
      call alloc_err(istat,'DAMLining_init','DAMLin_ObsInd',DAMLin_NumObs)
      allocate(DAMLin_File_Present(DAMLin_NumObs),stat=istat)
      call alloc_err(istat,'DAMLining_init','DAMLin_File_Present',DAMLin_NumObs)
    endif
#ifdef SPMD
   call mpibcast(DAMLin_ObsInd        , DAMLin_NumObs, mpiint, 0, mpicom)
   call mpibcast(DAMLin_File_Present  , DAMLin_NumObs, mpilog, 0, mpicom)
#endif

  ! Allocate Space for DAMLining observation arrays, initialize with 0's
  !---------------------------------------------------------------------
  allocate(Weights_U(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','Weights_U',pcols*pvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(Weights_V(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','Weights_V',pcols*pvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(Weights_T(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','Weights_T',pcols*pvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(Weights_Q(pvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','Weights_Q',pcols*pvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  !allocate(Weights_PS(pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  !call alloc_err(istat,'DAMLining_init','Weights_PS',pcols*((endchunk-begchunk)+1)*DAMLin_NumObs)

  Weights_U(:pvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  Weights_V(:pvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  Weights_T(:pvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  Weights_Q(:pvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  !Weights_PS(:pcols     ,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   
   
  allocate(mean_preds(Rpvar,pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','norm_preds',pcols*Rpvar*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(std_preds(Rpvar,pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','std_preds',pcols*Rpvar*((endchunk-begchunk)+1)*DAMLin_NumObs)
   
  mean_preds(:Rpvar,:pcols,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  std_preds(:Rpvar,:pcols,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  ! Allocate Space for DAMLining observation arrays, initialize with 0's
  !---------------------------------------------------------------------
  allocate(WeightsDiurnal_U(Dpvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','WeightsDiurnal_U',pcols*Dpvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(WeightsDiurnal_V(Dpvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','WeightsDiurnal_V',pcols*Dpvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(WeightsDiurnal_T(Dpvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','WeightsDiurnal_T',pcols*Dpvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  allocate(WeightsDiurnal_Q(Dpvar,pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  call alloc_err(istat,'DAMLining_init','WeightsDiurnal_Q',pcols*Dpvar*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
  !allocate(WeightsDiurnal_PS(pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
  !call alloc_err(istat,'DAMLining_init','WeightsDiurnal_PS',pcols*((endchunk-begchunk)+1)*DAMLin_NumObs)

  WeightsDiurnal_U(:Dpvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  WeightsDiurnal_V(:Dpvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  WeightsDiurnal_T(:Dpvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  WeightsDiurnal_Q(:Dpvar,:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
  !WeightsDiurnal_PS(:pcols    ,begchunk:endchunk,:DAMLin_NumObs)=0._r8
     
     
   ! Allocate Space for DAMLining observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   allocate(Weights_U_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','Weights_U_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(Weights_V_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','Weights_V_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(Weights_T_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','Weights_T_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(Weights_Q_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','Weights_Q_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(Weights_PS_int(pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','_int_PS_int',pcols*((endchunk-begchunk)+1)*DAMLin_NumObs)

   Weights_U_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   Weights_V_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   Weights_T_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   Weights_Q_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   Weights_PS_int(:pcols     ,begchunk:endchunk,:DAMLin_NumObs)=0._r8


   ! Allocate Space for DAMLining observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   allocate(WeightsDiurnal_U_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','WeightsDiurnal_U_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(WeightsDiurnal_V_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','WeightsDiurnal_V_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(WeightsDiurnal_T_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','WeightsDiurnal_T_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(WeightsDiurnal_Q_int(pcols,pver,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','WeightsDiurnal_Q_int',pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
   allocate(WeightsDiurnal_PS_int(pcols,begchunk:endchunk,DAMLin_NumObs),stat=istat)
   call alloc_err(istat,'DAMLining_init','WeightsDiurnal_PS_int',pcols*((endchunk-begchunk)+1)*DAMLin_NumObs)

   WeightsDiurnal_U_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   WeightsDiurnal_V_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   WeightsDiurnal_T_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   WeightsDiurnal_Q_int(:pcols,:pver,begchunk:endchunk,:DAMLin_NumObs)=0._r8
   WeightsDiurnal_PS_int(:pcols    ,begchunk:endchunk,:DAMLin_NumObs)=0._r8

!!DIAG
   if(masterproc) then
     write(iulog,*) 'DAMLinING: DAMLining_init() OBS arrays allocated and initialized'
     write(iulog,*) 'DAMLinING: DAMLining_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)
     write(iulog,*) 'DAMLinING: DAMLining_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*DAMLin_NumObs)/(1024._r8*1024._r8)
     write(iulog,*) 'DAMLinING: DAMLining_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'DAMLinING: DAMLining_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'DAMLinING: DAMLining_init() chunk:',(endchunk-begchunk+1),' DAMLin_NumObs=',DAMLin_NumObs
     write(iulog,*) 'DAMLinING: DAMLining_init() DAMLin_ObsInd=',DAMLin_ObsInd
     write(iulog,*) 'DAMLinING: DAMLining_init() DAMLin_File_Present=',DAMLin_File_Present
   endif
!!DIAG

    
   if(masterproc) then
    write(iulog,*) 'DAMLinING: Reading anomaly weights:',trim(DAMLin_Path)//trim(DAMLin_File_Anomaly)
    write(iulog,*) 'DAMLinING: Reading diurnal weights:',trim(DAMLin_Path)//trim(DAMLin_File_Diurnal)
    write(iulog,*) 'DAMLinING: Reading norm weights:',trim(DAMLin_Path)//trim(DAMLin_File_Norm)
   endif

   ! Rotate DAMLin_ObsInd() indices for new data, then update
   ! the DAMLin observation arrays with analysis data at the
   ! NEXT==DAMLin_ObsInd(1) time.
   !----------------------------------------------------------
   !WEC
   if(dycore_is('UNSTRUCTURED')) then
     
     call DAMLining_readin_weights_all(trim(DAMLin_Path)//trim(DAMLin_File_Anomaly),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Diurnal),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Norm))
    
   elseif(dycore_is('EUL')) then
    
     call DAMLining_readin_weights_all(trim(DAMLin_Path)//trim(DAMLin_File_Anomaly),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Diurnal),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Norm))
     
   else !if(dycore_is('LR')) then
     
     
     call DAMLining_readin_weights_all(trim(DAMLin_Path)//trim(DAMLin_File_Anomaly),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Diurnal),&
                                       trim(DAMLin_Path)//trim(DAMLin_File_Norm))
     
   endif
    
    
   !write(iulog,*) 'Wabba4:::',Weights_U(10,10,4:5,begchunk:endchunk,1) 
   ! Initialize DAMLining Coeffcient profiles in local arrays
   ! Load zeros into DAMLining arrays
   !------------------------------------------------------
   !write(iulog,*) 'done with that stuff.....'
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call DAMLining_set_profile(rlat,rlon,DAMLin_Uprof,Wprof,pver)
       DAMLin_Utau(icol,:,lchnk)=Wprof(:)
       call DAMLining_set_profile(rlat,rlon,DAMLin_Vprof,Wprof,pver)
       DAMLin_Vtau(icol,:,lchnk)=Wprof(:)
       call DAMLining_set_profile(rlat,rlon,DAMLin_Tprof,Wprof,pver)
       DAMLin_Stau(icol,:,lchnk)=Wprof(:)
       call DAMLining_set_profile(rlat,rlon,DAMLin_Qprof,Wprof,pver)
       DAMLin_Qtau(icol,:,lchnk)=Wprof(:)
            
       !+++ WEC set   
       DAMLin_PStau(icol,lchnk)=DAMLining_set_PSprofile(rlat,rlon,DAMLin_PSprof)
     end do

     if(DAMLin_Do==1) then
       DAMLin_Utau(:ncol,:pver,lchnk) =                             &
       DAMLin_Utau(:ncol,:pver,lchnk) * DAMLin_Ucoef
       DAMLin_Vtau(:ncol,:pver,lchnk) =                             &
       DAMLin_Vtau(:ncol,:pver,lchnk) * DAMLin_Vcoef
       DAMLin_Stau(:ncol,:pver,lchnk) =                             &
       DAMLin_Stau(:ncol,:pver,lchnk) * DAMLin_Tcoef
       DAMLin_Qtau(:ncol,:pver,lchnk) =                             &
       DAMLin_Qtau(:ncol,:pver,lchnk) * DAMLin_Qcoef
       DAMLin_PStau(:ncol,lchnk)=                             &
       DAMLin_PStau(:ncol,lchnk)* DAMLin_PScoef
    else
      DAMLin_Utau(:ncol,:pver,lchnk) =                             &
      DAMLin_Utau(:ncol,:pver,lchnk) * DAMLin_Ucoef/float(DAMLin_Step)
      DAMLin_Vtau(:ncol,:pver,lchnk) =                             &
      DAMLin_Vtau(:ncol,:pver,lchnk) * DAMLin_Vcoef/float(DAMLin_Step)
      DAMLin_Stau(:ncol,:pver,lchnk) =                             &
      DAMLin_Stau(:ncol,:pver,lchnk) * DAMLin_Tcoef/float(DAMLin_Step)
      DAMLin_Qtau(:ncol,:pver,lchnk) =                             &
      DAMLin_Qtau(:ncol,:pver,lchnk) * DAMLin_Qcoef/float(DAMLin_Step)
      DAMLin_PStau(:ncol,lchnk)=                             &
      DAMLin_PStau(:ncol,lchnk)* DAMLin_PScoef/float(DAMLin_Step)
    endif
 
     DAMLin_Ustep(:pcols,:pver,lchnk)=0._r8
     DAMLin_Vstep(:pcols,:pver,lchnk)=0._r8
     DAMLin_Sstep(:pcols,:pver,lchnk)=0._r8
     DAMLin_Qstep(:pcols,:pver,lchnk)=0._r8
     DAMLin_PSstep(:pcols,lchnk)=0._r8
     
     DAMLin_Ustep_avg(:pcols,:pver,lchnk)=0._r8
     DAMLin_Vstep_avg(:pcols,:pver,lchnk)=0._r8
     DAMLin_Sstep_avg(:pcols,:pver,lchnk)=0._r8
     DAMLin_Qstep_avg(:pcols,:pver,lchnk)=0._r8
     DAMLin_PSstep_avg(:pcols,lchnk)=0._r8

   end do
   
   !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine ! damlining_init
  !================================================================

  subroutine DAMLining_readin_weights_all(anal_file,anal_file2,anal_file3)
   !
   ! NUDGING_UPDATE_ANALYSES:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid          ,only: pcols,pver,begchunk,endchunk
   use cam_pio_utils   ,only: cam_pio_openfile
   use pio             ,only: PIO_BCAST_ERROR,PIO_INTERNAL_ERROR,pio_inq_dimlen,pio_inq_dimid
   use pio             ,only: pio_inq_att, pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, pio_inq_vardimid, pio_inq_varid,pio_get_var
   use pio             ,only: pio_closefile,pio_seterrorhandling,file_desc_t
   use ncdio_atm       ,only: infld
   use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
   use ioFileMod       ,only: getfil
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file
   character(len=*),intent(in):: anal_file2
   character(len=*),intent(in):: anal_file3

   ! Local values
   !-------------
   type(file_desc_t)  :: fileID,fileID2,fileID3
   integer            :: nn,Nindex,pvar,c
   integer            :: dimid,dimid2,dimid3
   logical            :: VARflag
   integer            :: grid_id
   integer            :: ierr
   integer            :: err_handling,err_handling2,err_handling3
   
   
   character(len=256) :: locfn,locfn2,locfn3
   character(len=256) :: filespec,filespec2,filespec3
   integer nlon,nlat,nlev,npredvar,npredvar2,npredvar3

   real(r8),allocatable:: Tmp4D(:,:,:,:)
   real(r8),allocatable:: Tmp3D(:,:,:)
   
   real(r8),allocatable:: Tmp4D_2(:,:,:,:)
   real(r8),allocatable:: Tmp3D_2(:,:,:)
   
   real(r8),allocatable:: Tmp3D_3(:,:,:)

   character(len=*), parameter :: prefix = 'nudging_update_analyses: '
   character(len=*), parameter :: sub = 'DAMLdiurnal_weights_init'

   ! Rotate DAMLin_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=DAMLin_ObsInd(DAMLin_NumObs)
     do nn=DAMLin_NumObs,2,-1
       DAMLin_ObsInd(nn)=DAMLin_ObsInd(nn-1)
     end do
     DAMLin_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=DAMLin_File_Present(DAMLin_ObsInd(1)))
     write(iulog,*)'NUDGING: DAMLin_ObsInd=',DAMLin_ObsInd
     write(iulog,*)'NUDGING: DAMLin_File_Present=',DAMLin_File_Present
   endif

   call MPI_bcast(DAMLin_File_Present, DAMLin_NumObs, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_File_Present')
   call MPI_bcast(DAMLin_ObsInd      , DAMLin_NumObs, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_ObsInd')

   if(.not. DAMLin_File_Present(DAMLin_ObsInd(1))) then
      return
   end if

   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec = trim( anal_file )
   call getfil( filespec, locfn, 0 )
   call cam_pio_openfile(fileID,trim(locfn), 0)
   
   call pio_seterrorhandling(fileID,PIO_BCAST_ERROR,oldmethod=err_handling)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file)
   
   ierr = pio_inq_dimid( fileID, 'predvar', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, npredvar )
   
   ierr = pio_inq_dimid( fileID, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlev )
   
   allocate(Tmp4D(npredvar,pcols,pver,begchunk:endchunk))
   allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'npredvar:',npredvar,'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('U',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: U not found in anomaly file')
             
   Weights_U(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U::',Weights_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('U_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: U int not found in anomaly file')
              
   Weights_U_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _U_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U_int::',WeightsDiurnal_U_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('V',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: V not found in anomaly file')
              
   Weights_V(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U::',Weights_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('V_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: V int not found in anomaly file')
              
   Weights_V_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _V_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_V_int::',WeightsDiurnal_V_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   
   call infld('T',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: T not found in anomaly file')
              
   Weights_T(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_T::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_T::',Weights_T(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('T_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: T int not found in anomaly file')
              
   Weights_T_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _T_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_T_int::',WeightsDiurnal_T_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('Q',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q not found in anomaly file')
              
   Weights_Q(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_Q::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_Q::',Weights_Q(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   
   call infld('Q_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q int not found in anomaly file')
              
   Weights_Q_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _Q_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_Q_int::',WeightsDiurnal_Q_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   !write(iulog,*) 'Wabba3:::',Weights_U(10,10,4:5,begchunk:endchunk,1)
   
   !##################################################################################################
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec2 = trim( anal_file2 )
   call getfil( filespec2, locfn2, 0 )
   call cam_pio_openfile(fileID2,trim(locfn2), 0)
   
   call pio_seterrorhandling(fileID2,PIO_BCAST_ERROR,oldmethod=err_handling2)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file2)
   
   ierr = pio_inq_dimid( fileID2, 'predvar', dimid2 )
   ierr = pio_inq_dimlen( fileID2, dimid2, npredvar2 )
   
   ierr = pio_inq_dimid( fileID2, 'lat', dimid2 )
   ierr = pio_inq_dimlen( fileID2, dimid2, nlat )
   
   ierr = pio_inq_dimid( fileID2, 'lon', dimid2 )
   ierr = pio_inq_dimlen( fileID2, dimid2, nlon )
   
   ierr = pio_inq_dimid( fileID2, 'lev', dimid2 )
   ierr = pio_inq_dimlen( fileID2, dimid2, nlev )
   
   allocate(Tmp4D_2(npredvar2,pcols,pver,begchunk:endchunk))
   allocate(Tmp3D_2(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'npredvar:',npredvar2,'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !------------------------------------

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('dU',fileID2,'predvar','lat','lev','lon', &
              1,npredvar2,1,pcols,1,pver,begchunk,endchunk,Tmp4D_2, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: dU not found in anomaly file')
   
   WeightsDiurnal_U(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D_2(:,:,:,begchunk:endchunk)
   
   call infld('dU_int',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D_2, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: U int not found in anomaly file')
              
   WeightsDiurnal_U_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_2(:,:,begchunk:endchunk)
   
   
   call infld('dV',fileID2,'predvar','lat','lev','lon', &
              1,npredvar2,1,pcols,1,pver,begchunk,endchunk,Tmp4D_2, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: V not found in anomaly file')
              
   WeightsDiurnal_V(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D_2(:,:,:,begchunk:endchunk)

   call infld('dV_int',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D_2, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: V int not found in anomaly file')
              
   WeightsDiurnal_V_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_2(:,:,begchunk:endchunk)
   
   call infld('dT',fileID2,'predvar','lat','lev','lon', &
              1,npredvar2,1,pcols,1,pver,begchunk,endchunk,Tmp4D_2, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: T not found in anomaly file')
              
   WeightsDiurnal_T(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D_2(:,:,:,begchunk:endchunk)
   
   call infld('dT_int',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D_2, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: T int not found in anomaly file')
              
   WeightsDiurnal_T_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_2(:,:,begchunk:endchunk)
   
   call infld('dQ',fileID2,'predvar','lat','lev','lon', &
              1,npredvar2,1,pcols,1,pver,begchunk,endchunk,Tmp4D_2, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q not found in anomaly file')
              
   WeightsDiurnal_Q(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D_2(:,:,:,begchunk:endchunk)
   
   
   call infld('dQ_int',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D_2, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q int not found in anomaly file')
              
   WeightsDiurnal_Q_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_2(:,:,begchunk:endchunk)
   !############################################################################################
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec3 = trim( anal_file3 )
   call getfil( filespec3, locfn3, 0 )
   call cam_pio_openfile(fileID3,trim(locfn3), 0)
   
   call pio_seterrorhandling(fileID3,PIO_BCAST_ERROR,oldmethod=err_handling3)
   
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file)
   
   ierr = pio_inq_dimid( fileID3, 'predvar', dimid3 )
   ierr = pio_inq_dimlen( fileID3, dimid3, npredvar3 )
   
   ierr = pio_inq_dimid( fileID3, 'lat', dimid3 )
   ierr = pio_inq_dimlen( fileID3, dimid3, nlat )
   
   ierr = pio_inq_dimid( fileID3, 'lon', dimid3 )
   ierr = pio_inq_dimlen( fileID3, dimid3, nlon )
   
   allocate(Tmp3D_3(npredvar3,pcols,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'npredvar:',npredvar3,'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   call infld('mean_vars',fileID3,'predvar','lat','lon', &
              1,npredvar3,1,pcols,begchunk,endchunk,Tmp3D_3, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: mean_vars not found in anomaly file')        
   mean_preds(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_3(:,:,begchunk:endchunk)
   
   call infld('std_vars',fileID3,'predvar','lat','lon', &
              1,npredvar3,1,pcols,begchunk,endchunk,Tmp3D_3, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: std_vars int not found in anomaly file')
   std_preds(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_3(:,:,begchunk:endchunk)
   
   if(masterproc) write(iulog,*)'weights written'   
   
  
   ! Restore old error handling
   !----------------------------
   call pio_seterrorhandling(fileID,err_handling)
   call pio_seterrorhandling(fileID2,err_handling2)
   call pio_seterrorhandling(fileID3,err_handling3)

   ! Close the analyses file
   !-----------------------
   deallocate(Tmp3D)
   deallocate(Tmp4D)
   deallocate(Tmp3D_2)
   deallocate(Tmp4D_2)
   deallocate(Tmp3D_3)
   call pio_closefile(fileID)
   call pio_closefile(fileID2)
   call pio_closefile(fileID3)

   ! End Routine
   !------------

  end subroutine DAMLining_readin_weights_all


    !================================================================
  subroutine DAMLining_readin_weights(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid          ,only: pcols,pver,begchunk,endchunk
   use cam_pio_utils   ,only: cam_pio_openfile
   use pio             ,only: PIO_BCAST_ERROR,PIO_INTERNAL_ERROR,pio_inq_dimlen,pio_inq_dimid
   use pio             ,only: pio_inq_att, pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, pio_inq_vardimid, pio_inq_varid,pio_get_var
   use pio             ,only: pio_closefile,pio_seterrorhandling,file_desc_t
   use ncdio_atm       ,only: infld
   use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
   use ioFileMod       ,only: getfil
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   type(file_desc_t)  :: fileID
   integer            :: nn,Nindex,pvar,c
   integer            :: dimid,vid
   logical            :: VARflag
   integer            :: grid_id
   integer            :: ierr
   character(len=DLEN):: dim1name,dim2name
   integer            :: err_handling
   character(len=256) :: locfn
   character(len=256) :: filespec
   integer nlon,nlat,nlev,npredvar

   real(r8),allocatable:: Tmp4D(:,:,:,:)
   real(r8),allocatable:: Tmp3D(:,:,:)

   character(len=*), parameter :: prefix = 'nudging_update_analyses: '
   character(len=*), parameter :: sub = 'DAMLdiurnal_weights_init'

   ! Rotate DAMLin_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=DAMLin_ObsInd(DAMLin_NumObs)
     do nn=DAMLin_NumObs,2,-1
       DAMLin_ObsInd(nn)=DAMLin_ObsInd(nn-1)
     end do
     DAMLin_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=DAMLin_File_Present(DAMLin_ObsInd(1)))
     write(iulog,*)'NUDGING: DAMLin_ObsInd=',DAMLin_ObsInd
     write(iulog,*)'NUDGING: DAMLin_File_Present=',DAMLin_File_Present
   endif

   call MPI_bcast(DAMLin_File_Present, DAMLin_NumObs, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_File_Present')
   call MPI_bcast(DAMLin_ObsInd      , DAMLin_NumObs, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_ObsInd')

   if(.not. DAMLin_File_Present(DAMLin_ObsInd(1))) then
      return
   end if

   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec = trim( anal_file )
   call getfil( filespec, locfn, 0 )
   call cam_pio_openfile(fileID,trim(locfn), 0)
   
   call pio_seterrorhandling(fileID,PIO_BCAST_ERROR,oldmethod=err_handling)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file)
   
   ierr = pio_inq_dimid( fileID, 'predvar', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, npredvar )
   
   ierr = pio_inq_dimid( fileID, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlev )
   
   allocate(Tmp4D(npredvar,pcols,pver,begchunk:endchunk))
   allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'npredvar:',npredvar,'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('U',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: U not found in anomaly file')
             
   Weights_U(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U::',Weights_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('U_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: U int not found in anomaly file')
              
   Weights_U_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _U_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U_int::',WeightsDiurnal_U_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('V',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: V not found in anomaly file')
              
   Weights_V(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U::',Weights_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('V_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: V int not found in anomaly file')
              
   Weights_V_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _V_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_V_int::',WeightsDiurnal_V_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   
   call infld('T',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: T not found in anomaly file')
              
   Weights_T(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_T::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_T::',Weights_T(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('T_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: T int not found in anomaly file')
              
   Weights_T_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _T_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_T_int::',WeightsDiurnal_T_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('Q',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q not found in anomaly file')
              
   Weights_Q(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP4D_Q::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_Q::',Weights_Q(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   
   call infld('Q_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q int not found in anomaly file')
              
   Weights_Q_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Pawg TEMP3D _Q_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_Q_int::',WeightsDiurnal_Q_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   !write(iulog,*) 'Wabba3:::',Weights_U(10,10,4:5,begchunk:endchunk,1)

   
   if(masterproc) write(iulog,*)'weights written'  

   ! Restore old error handling
   !----------------------------
   call pio_seterrorhandling(fileID,err_handling)

   ! Close the analyses file
   !-----------------------
   deallocate(Tmp3D)
   deallocate(Tmp4D)
   call pio_closefile(fileID)

   ! End Routine
   !------------

  end subroutine DAMLining_readin_weights
  !================================================================
  !================================================================
 
  
  
  
  
  
  
  
  
  
  
  subroutine DAMLining_readin_weights_diurnal(anal_file)
   !
   ! NUDGING_UPDATE_ANALYSES:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid          ,only: pcols,pver,begchunk,endchunk
   use cam_pio_utils   ,only: cam_pio_openfile
   use pio             ,only: PIO_BCAST_ERROR,PIO_INTERNAL_ERROR,pio_inq_dimlen,pio_inq_dimid
   use pio             ,only: pio_inq_att, pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, pio_inq_vardimid, pio_inq_varid,pio_get_var
   use pio             ,only: pio_closefile,pio_seterrorhandling,file_desc_t
   use ncdio_atm       ,only: infld
   use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
   use ioFileMod       ,only: getfil
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   type(file_desc_t)  :: fileID
   integer            :: nn,Nindex,pvar,c
   integer            :: dimid,vid
   logical            :: VARflag
   integer            :: grid_id
   integer            :: ierr
   character(len=DLEN):: dim1name,dim2name
   integer            :: err_handling
   character(len=256) :: locfn
   character(len=256) :: filespec
   integer nlon,nlat,nlev,npredvar

   real(r8),allocatable:: Tmp4D(:,:,:,:)
   real(r8),allocatable:: Tmp3D(:,:,:)

   character(len=*), parameter :: prefix = 'nudging_update_analyses: '
   character(len=*), parameter :: sub = 'DAMLdiurnal_weights_init'

   ! Rotate DAMLin_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=DAMLin_ObsInd(DAMLin_NumObs)
     do nn=DAMLin_NumObs,2,-1
       DAMLin_ObsInd(nn)=DAMLin_ObsInd(nn-1)
     end do
     DAMLin_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=DAMLin_File_Present(DAMLin_ObsInd(1)))
     write(iulog,*)'NUDGING: DAMLin_ObsInd=',DAMLin_ObsInd
     write(iulog,*)'NUDGING: DAMLin_File_Present=',DAMLin_File_Present
   endif

   call MPI_bcast(DAMLin_File_Present, DAMLin_NumObs, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_File_Present')
   call MPI_bcast(DAMLin_ObsInd      , DAMLin_NumObs, mpi_integer, mstrid, mpicom, ierr)
   if (ierr /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_ObsInd')

   if(.not. DAMLin_File_Present(DAMLin_ObsInd(1))) then
      return
   end if

   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec = trim( anal_file )
   call getfil( filespec, locfn, 0 )
   call cam_pio_openfile(fileID,trim(locfn), 0)
   
   call pio_seterrorhandling(fileID,PIO_BCAST_ERROR,oldmethod=err_handling)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file)
   
   ierr = pio_inq_dimid( fileID, 'predvar', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, npredvar )
   
   ierr = pio_inq_dimid( fileID, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID, dimid, nlev )
   
   allocate(Tmp4D(npredvar,pcols,pver,begchunk:endchunk))
   allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('dU',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: dU not found in anomaly file')
   
   WeightsDiurnal_U(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP Diurnal_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_U::',WeightsDiurnal_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('dU_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: U int not found in anomaly file')
              
   WeightsDiurnal_U_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP3D Diurnal_U_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_U_int::',WeightsDiurnal_U_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('dV',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: V not found in anomaly file')
              
   WeightsDiurnal_V(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP Diurnal_V::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_V::',WeightsDiurnal_V(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('dV_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: V int not found in anomaly file')
              
   WeightsDiurnal_V_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP3D Diurnal_V_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_V_int::',WeightsDiurnal_V_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   
   call infld('dT',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: T not found in anomaly file')
              
   WeightsDiurnal_T(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP Diurnal_T::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_T::',WeightsDiurnal_T(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('dT_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: T int not found in anomaly file')
              
   WeightsDiurnal_T_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP3D Diurnal_T_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_T_int::',WeightsDiurnal_T_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('dQ',fileID,'predvar','lat','lev','lon', &
              1,npredvar,1,pcols,1,pver,begchunk,endchunk,Tmp4D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q not found in anomaly file')
              
   WeightsDiurnal_Q(:,:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp4D(:,:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP Diurnal_Q::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_Q::',WeightsDiurnal_Q(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   
   call infld('dQ_int',fileID,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1)
              
   if (.not. VARflag) call endrun(sub//': ERROR: Q int not found in anomaly file')
              
   WeightsDiurnal_Q_int(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D(:,:,begchunk:endchunk)
   !write(iulog,*),'Dawg TEMP3D Diurnal_Q_int::',Tmp3D(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Dawg WeightsDiurnal_Q_int::',WeightsDiurnal_Q_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   if(masterproc) write(iulog,*)'weights written'   
   
   

   ! Restore old error handling
   !----------------------------
   call pio_seterrorhandling(fileID,err_handling)

   ! Close the analyses file
   !-----------------------
   deallocate(Tmp3D)
   deallocate(Tmp4D)
   call pio_closefile(fileID)

   ! End Routine
   !------------

  end subroutine DAMLining_readin_weights_diurnal

  
  subroutine DAMLining_readin_normdict(anal_nd_file)
   ! Toadstool
   ! NUDGING_UPDATE_ANALYSES:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid          ,only: pcols,pver,begchunk,endchunk
   use cam_pio_utils   ,only: cam_pio_openfile
   use pio             ,only: PIO_BCAST_ERROR,PIO_INTERNAL_ERROR,pio_inq_dimlen,pio_inq_dimid
   use pio             ,only: pio_inq_att, pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, pio_inq_vardimid, pio_inq_varid,pio_get_var
   use pio             ,only: pio_closefile,pio_seterrorhandling,file_desc_t
   use ncdio_atm       ,only: infld
   use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
   use ioFileMod       ,only: getfil
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_nd_file

   ! Local values
   !-------------
   type(file_desc_t)  :: fileID_nd
   integer            :: nn_nd,Nindex_nd
   integer            :: dimid_nd
   logical            :: VARflag_nd
   integer            :: ierr_nd
   integer            :: err_handling_nd
   character(len=256) :: locfn_nd
   character(len=256) :: filespec_nd
   integer nlon,nlat,nlev,npredvar_nd

   real(r8),allocatable:: Tmp3D_nd(:,:,:)

   character(len=*), parameter :: prefix = 'nudging_update_analyses: '
   character(len=*), parameter :: sub = 'DAMLnorm_weights_init'

   ! Rotate DAMLin_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex_nd=DAMLin_ObsInd(DAMLin_NumObs)
     do nn_nd=DAMLin_NumObs,2,-1
       DAMLin_ObsInd(nn_nd)=DAMLin_ObsInd(nn_nd-1)
     end do
     DAMLin_ObsInd(1)=Nindex_nd
     inquire(FILE=trim(anal_nd_file),EXIST=DAMLin_File_Present(DAMLin_ObsInd(1)))
     write(iulog,*)'NUDGING: DAMLin_ObsInd=',DAMLin_ObsInd
     write(iulog,*)'NUDGING: DAMLin_File_Present=',DAMLin_File_Present
   endif

   call MPI_bcast(DAMLin_File_Present, DAMLin_NumObs, mpi_logical, mstrid, mpicom, ierr_nd)
   if (ierr_nd /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_File_Present')
   call MPI_bcast(DAMLin_ObsInd      , DAMLin_NumObs, mpi_integer, mstrid, mpicom, ierr_nd)
   if (ierr_nd /= mpi_success) call endrun(prefix//'FATAL: mpi_bcast: DAMLin_ObsInd')

   if(.not. DAMLin_File_Present(DAMLin_ObsInd(1))) then
      return
   end if

   ! Open the file and get the fileID_nd, and define the predvar length
   !-------------------------------------
   filespec_nd = trim( anal_nd_file )
   call getfil( filespec_nd, locfn_nd, 0 )
   call cam_pio_openfile(fileID_nd,trim(locfn_nd), 0)
   
   call pio_seterrorhandling(fileID_nd,PIO_BCAST_ERROR,oldmethod=err_handling_nd)
   
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_nd_file)
   
   ierr_nd = pio_inq_dimid( fileID_nd, 'predvar', dimid_nd )
   ierr_nd = pio_inq_dimlen( fileID_nd, dimid_nd, npredvar_nd )
   
   ierr_nd = pio_inq_dimid( fileID_nd, 'lat', dimid_nd )
   ierr_nd = pio_inq_dimlen( fileID_nd, dimid_nd, nlat )
   
   ierr_nd = pio_inq_dimid( fileID_nd, 'lon', dimid_nd )
   ierr_nd = pio_inq_dimlen( fileID_nd, dimid_nd, nlon )
   
   allocate(Tmp3D_nd(npredvar_nd,pcols,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'npredvar_nd:',npredvar_nd,'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('mean_vars',fileID_nd,'predvar','lat','lon', &
              1,npredvar_nd,1,pcols,begchunk,endchunk,Tmp3D_nd, &
              VARflag_nd,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag_nd) call endrun(sub//': ERROR: mean_vars not found in anomaly file')        
   mean_preds(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_nd(:,:,begchunk:endchunk)
   
   !write(iulog,*),'Pawg TEMP4D_U::',Tmp4D(10,10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U::',Weights_U(10,10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   
   call infld('std_vars',fileID_nd,'predvar','lat','lon', &
              1,npredvar_nd,1,pcols,begchunk,endchunk,Tmp3D_nd, &
              VARflag_nd,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag_nd) call endrun(sub//': ERROR: std_vars int not found in anomaly file')
   std_preds(:,:,begchunk:endchunk,DAMLin_ObsInd(1)) = Tmp3D_nd(:,:,begchunk:endchunk)
   
   !write(iulog,*),'Pawg TEMP3D _U_int::',Tmp3D_nd(10,4:5,begchunk:endchunk)
   !write(iulog,*),'Pawg Weights_U_int::',WeightsDiurnal_U_int(10,4:5,begchunk:endchunk,DAMLin_ObsInd(1))
   !write(iulog,*) 'Wabba3:::',Weights_U(10,10,4:5,begchunk:endchunk,1)
   
   if(masterproc) write(iulog,*)'weights written'  
   
   ! Restore old error handling
   !----------------------------
   call pio_seterrorhandling(fileID_nd,err_handling_nd)
   
   ! Close the analyses file
   !-----------------------
   deallocate(Tmp3D_nd)
   call pio_closefile(fileID_nd)
   ! End Routine
   !------------

  end subroutine DAMLining_readin_normdict
  !================================================================
  !================================================================
  
  

  subroutine regress_diurnal_daml_timestep_tend(phys_state,phys_tend,pbuf,cam_in,cam_out)
  !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the DAMLin
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pverp,pcols,begchunk,endchunk
   use cam_history  ,only: outfld
   use physconst,    only: zvir, gravit, cpair, rair
   use phys_grid    ,only: get_rlat_p, get_rlon_p
   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             physics_buffer_desc
                             
   use camsrfexch,     only: cam_in_t,cam_out_t


   
   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(cam_in_t),      intent(in)    :: cam_in
   type(cam_out_t),     intent(in)    :: cam_out
   type(physics_ptend), intent(out):: phys_tend
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state) :: state1                ! Local copy of state variable
   
   

   ! Local values
   !--------------------
   integer :: lchnk,ncol,indw,icol,ipver,ipp,ipverp,tango,itim_old
   integer  Year,Month,Day,Sec,Hour,calday
   logical lq(pcnst)
   real(r8) hilat,lolat,hilon,lolon,PI
   real(r8) lattmp, lontmp
   real(r8) DOY1,DOY2,DOY3,DOY4,DOY5,DOY6,DOY7,DOY8
   real(r8) HOD1,HOD2,HOD3,HOD4
   real(r8) tend_u,tend_v,tend_q,tend_t
   real(r8) tend_ur,tend_vr,tend_qr,tend_tr
   
   real(r8), pointer, dimension(:,:) :: upwp      ! upwp                               [no idea]
   real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height                [m]
   real(r8), pointer, dimension(:) :: taux     ! Surface zonal momentum flux                [m]
   real(r8), pointer, dimension(:) :: tauy     ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: ustarwec ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: fsntoa ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: flntc  ! Surface meridional momentum flux
   
   hilat = 65.0_r8
   lolat = 63.0_r8
   hilon = 20.0_r8
   lolon = 18.0_r8
   PI    = 4.D0*DATAN(1.D0)
   
   call physics_state_copy(phys_state, state1)
   
   !TODO WEC - make a namelist selection for just regressing out Seasonal + Diurnal Cycle. 
   
   ncol =phys_state%ncol
   lchnk=phys_state%lchnk

   !call cnst_get_ind('Q',indw)
   !lq(:)   =.false.
   !lq(indw)=.true.
   !call physics_ptend_init(phys_tend,phys_state%psetcols,'daml',lu=.true.,lv=.true.,ls=.true.,lq=lq)
   
  !Call hour to regress T.O.D. and D.O.Y. 
  call get_curr_date(Year,Month,Day,Sec)
  Hour = Sec/3600.0_r8
  
  calday = get_curr_calday()
  
  !write(iulog,*)'YMDH:',Year,'|',Month,'|',Day,'|',Hour,'|'
  !Set DOY/TOD vars: 
  DOY1= SIN(2.0_r8*PI*calday/365.0_r8)
  DOY2= COS(2.0_r8*PI*calday/365.0_r8)
  DOY3= SIN(4.0_r8*PI*calday/365.0_r8)
  DOY4= COS(4.0_r8*PI*calday/365.0_r8)
  DOY5= SIN(6.0_r8*PI*calday/365.0_r8)
  DOY6= COS(6.0_r8*PI*calday/365.0_r8)
  DOY7= SIN(8.0_r8*PI*calday/365.0_r8)
  DOY8= COS(8.0_r8*PI*calday/365.0_r8)
  HOD1= SIN(2.0_r8*PI*Hour/24.0_r8)
  HOD2= COS(2.0_r8*PI*Hour/24.0_r8)
  HOD3= SIN(4.0_r8*PI*Hour/24.0_r8)
  HOD4= COS(4.0_r8*PI*Hour/24.0_r8)
  
  
  !write(iulog,*)'YMDH:',Year,'|',Month,'|',Day,'|',Hour,'|'
  !write(iulog,*)'Harmonic SINE/COS DOY:',DOY1,'|',DOY2,'|',DOY3,'|',DOY4,'|',DOY5,'|',DOY6,'|',DOY7,'|',DOY8,'|'
  !write(iulog,*)'Harmonic SINE/COS HOD:',HOD1,'|',HOD2,'|',HOD3,'|',HOD4,'|'
  
  
  !check the weights: 
  
  !do icol=1,ncol
  !    lattmp = get_rlat_p(lchnk,icol)*57.2957795131_r8
  !    lontmp = get_rlon_p(lchnk,icol)*57.2957795131_r8
  !    if (lattmp < 65.0_r8 .and. lattmp > 63.0_r8 .and. lontmp < 20.0_r8 .and. lontmp > 18.0_r8) then
  !        write(iulog,*),'WDU_ plane ','lat:',lattmp,'lon:',lontmp,WeightsDiurnal_U(1,icol,:,lchnk,DAMLin_ObsInd(1))
  !        write(iulog,*),'W_U train ','lat:',lattmp,'lon:',lontmp,Weights_U(1,icol,:,lchnk,DAMLin_ObsInd(1))    
  !    end if 
  !    
  !    if (lattmp < 25.0_r8 .and. lattmp > 23.0_r8 .and. lontmp < 105.0_r8 .and. lontmp > 103.0_r8) then
  !        write(iulog,*),'WDU_ plane ','lat:',lattmp,'lon:',lontmp,WeightsDiurnal_U(1,icol,:,lchnk,DAMLin_ObsInd(1))
  !        write(iulog,*),'W_U train ','lat:',lattmp,'lon:',lontmp,Weights_U(1,icol,:,lchnk,DAMLin_ObsInd(1))    
  !    end if 
  !end do
  

  do icol=1,ncol
    lattmp = get_rlat_p(lchnk,icol)*57.2957795131_r8
    lontmp = get_rlon_p(lchnk,icol)*57.2957795131_r8
    do ipver=1,pver
        tend_u = WeightsDiurnal_U(1,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD1 + &
                        WeightsDiurnal_U(2,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD2 + &
                        WeightsDiurnal_U(3,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD3 + &
                        WeightsDiurnal_U(4,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD4 + &
                        WeightsDiurnal_U(5,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY1 + &
                        WeightsDiurnal_U(6,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY2 + &
                        WeightsDiurnal_U(7,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY3 + &
                        WeightsDiurnal_U(8,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY4 + &
                        WeightsDiurnal_U(9,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY5 + &
                        WeightsDiurnal_U(10,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY6 + &
                        WeightsDiurnal_U(11,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY7 + &
                        WeightsDiurnal_U(12,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY8 + WeightsDiurnal_U_int(icol,ipver,lchnk,DAMLin_ObsInd(1))
                        
        !if(masterproc) write(iulog,*)'chicken tend_U=',WeightsDiurnal_U(3,icol,ipver,lchnk,DAMLin_ObsInd(1))
                        
        DAMLin_Ustep_avg(icol,ipver,lchnk)=tend_u
        
        tend_v = WeightsDiurnal_V(1,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD1 + &
                        WeightsDiurnal_V(2,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD2 + &
                        WeightsDiurnal_V(3,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD3 + &
                        WeightsDiurnal_V(4,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD4 + &
                        WeightsDiurnal_V(5,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY1 + &
                        WeightsDiurnal_V(6,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY2 + &
                        WeightsDiurnal_V(7,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY3 + &
                        WeightsDiurnal_V(8,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY4 + &
                        WeightsDiurnal_V(9,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY5 + &
                        WeightsDiurnal_V(10,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY6 + &
                        WeightsDiurnal_V(11,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY7 + &
                        WeightsDiurnal_V(12,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY8 + WeightsDiurnal_V_int(icol,ipver,lchnk,DAMLin_ObsInd(1))
                        
        DAMLin_Vstep_avg(icol,ipver,lchnk)=tend_v
        
        tend_t = WeightsDiurnal_T(1,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD1 + &
                        WeightsDiurnal_T(2,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD2 + &
                        WeightsDiurnal_T(3,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD3 + &
                        WeightsDiurnal_T(4,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD4 + &
                        WeightsDiurnal_T(5,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY1 + &
                        WeightsDiurnal_T(6,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY2 + &
                        WeightsDiurnal_T(7,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY3 + &
                        WeightsDiurnal_T(8,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY4 + &
                        WeightsDiurnal_T(9,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY5 + &
                        WeightsDiurnal_T(10,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY6 + &
                        WeightsDiurnal_T(11,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY7 + &
                        WeightsDiurnal_T(12,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY8 + WeightsDiurnal_T_int(icol,ipver,lchnk,DAMLin_ObsInd(1))
                        
        DAMLin_Sstep_avg(icol,ipver,lchnk)=tend_t*cpair
        
        tend_q = WeightsDiurnal_Q(1,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD1 + &
                        WeightsDiurnal_Q(2,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD2 + &
                        WeightsDiurnal_Q(3,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD3 + &
                        WeightsDiurnal_Q(4,icol,ipver,lchnk,DAMLin_ObsInd(1))*HOD4 + &
                        WeightsDiurnal_Q(5,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY1 + &
                        WeightsDiurnal_Q(6,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY2 + &
                        WeightsDiurnal_Q(7,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY3 + &
                        WeightsDiurnal_Q(8,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY4 + &
                        WeightsDiurnal_Q(9,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY5 + &
                        WeightsDiurnal_Q(10,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY6 + &
                        WeightsDiurnal_Q(11,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY7 + &
                        WeightsDiurnal_Q(12,icol,ipver,lchnk,DAMLin_ObsInd(1))*DOY8 + WeightsDiurnal_Q_int(icol,ipver,lchnk,DAMLin_ObsInd(1))
                        
        DAMLin_Qstep_avg(icol,ipver,lchnk)=tend_q
    end do
  end do
  
  do lchnk=begchunk,endchunk !Remove:
       DAMLin_Ustep_avg(:ncol,:pver,lchnk)=(DAMLin_Ustep_avg(:ncol,:pver,lchnk))*DAMLin_Utau(:ncol,:pver,lchnk) !These are all tendencies already [WEC]
       DAMLin_Vstep_avg(:ncol,:pver,lchnk)=(DAMLin_Vstep_avg(:ncol,:pver,lchnk))*DAMLin_Vtau(:ncol,:pver,lchnk)
       DAMLin_Sstep_avg(:ncol,:pver,lchnk)=(DAMLin_Sstep_avg(:ncol,:pver,lchnk))*DAMLin_Stau(:ncol,:pver,lchnk) !cant do Temp [WEC]
       DAMLin_Qstep_avg(:ncol,:pver,lchnk)=(DAMLin_Qstep_avg(:ncol,:pver,lchnk))*DAMLin_Qtau(:ncol,:pver,lchnk)
       !DAMLin_PSstep(:ncol,     lchnk)=(Target_stoch_PS(:ncol,lchnk))*Stochai_PStau(:ncol,lchnk) !cant do PS [WEC]
  end do
  
  
  end subroutine regress_diurnal_daml_timestep_tend
  
  
  
  subroutine regress_daml_timestep_tend(phys_state,phys_tend,pbuf,cam_in,cam_out)
  !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the DAMLin
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pverp,pcols,begchunk,endchunk
   use cam_history  ,only: outfld
   use physconst,    only: zvir, gravit, cpair, rair
   use phys_grid    ,only: get_rlat_p, get_rlon_p
   use physics_buffer, only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                             physics_buffer_desc                          
   use camsrfexch,     only: cam_in_t,cam_out_t

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(cam_in_t),      intent(in)    :: cam_in
   type(cam_out_t),     intent(in)    :: cam_out
   type(physics_ptend), intent(out):: phys_tend
   type(physics_buffer_desc), pointer :: pbuf(:)
   type(physics_state) :: state1                ! Local copy of state variable
   
   ! Local values
   !--------------------
   integer :: lchnk,ncol,indw,icol,ipver,ipp,ipverp,tango,itim_old
   integer  Year,Month,Day,Sec,Hour,DP,calday
   logical lq(pcnst)
   real(r8) hilat,lolat,hilon,lolon,PI
   real(r8) lattmp, lontmp
   real(r8) DOY1,DOY2,DOY3,DOY4,DOY5,DOY6,DOY7,DOY8
   real(r8) HOD1,HOD2,HOD3,HOD4
   real(r8) tend_u,tend_v,tend_q,tend_t
   real(r8) threshy
   real(r8) tend_ur_u,tend_vr_u,tend_qr_u,tend_tr_u,tend_or_u,tend_auxr_u
   real(r8) tend_ur_v,tend_vr_v,tend_qr_v,tend_tr_v,tend_or_v,tend_auxr_v
   real(r8) tend_ur_t,tend_vr_t,tend_qr_t,tend_tr_t,tend_or_t,tend_auxr_t,ttend_val
   real(r8) tend_ur_q,tend_vr_q,tend_qr_q,tend_tr_q,tend_or_q,tend_auxr_q
   
   real(r8), pointer, dimension(:,:) :: upwp      ! upwp                               [no idea]
   real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height       [m]
   real(r8), pointer, dimension(:) :: taux     ! Surface zonal momentum flux           [m]
   real(r8), pointer, dimension(:) :: tauy     ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: ustarwec ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: fsntoa ! Surface meridional momentum flux
   real(r8), pointer, dimension(:) :: flntc  ! Surface meridional momentum flux
   
   hilat = 65.0_r8
   lolat = 63.0_r8
   hilon = 20.0_r8
   lolon = 18.0_r8
   threshy = 0.1
   PI    = 4.D0*DATAN(1.D0)
   
   call physics_state_copy(phys_state, state1)
   
   !TODO WEC - make a namelist selection for just regressing out Seasonal + Diurnal Cycle. 
   
   ncol =phys_state%ncol
   lchnk=phys_state%lchnk

   !call cnst_get_ind('Q',indw)
   !lq(:)   =.false.
   !lq(indw)=.true.
   !call physics_ptend_init(phys_tend,phys_state%psetcols,'daml',lu=.true.,lv=.true.,ls=.true.,lq=lq)
   
  !Call hour to regress T.O.D. and D.O.Y. 
  call get_curr_date(Year,Month,Day,Sec)
  calday = get_curr_calday()
  Hour = Sec/3600.0_r8
  
  DP = DAMLin_ObsInd(1)
  
  !write(iulog,*)'YMDH:',Year,'|',Month,'|',Day,'|',Hour,'|'
  !Set DOY/TOD vars: 
  DOY1= SIN(2.0_r8*PI*calday/365.0_r8)
  DOY2= COS(2.0_r8*PI*calday/365.0_r8)
  DOY3= SIN(4.0_r8*PI*calday/365.0_r8)
  DOY4= COS(4.0_r8*PI*calday/365.0_r8)
  DOY5= SIN(6.0_r8*PI*calday/365.0_r8)
  DOY6= COS(6.0_r8*PI*calday/365.0_r8)
  DOY7= SIN(8.0_r8*PI*calday/365.0_r8)
  DOY8= COS(8.0_r8*PI*calday/365.0_r8)
  HOD1= SIN(2.0_r8*PI*Hour/24.0_r8)
  HOD2= COS(2.0_r8*PI*Hour/24.0_r8)
  HOD3= SIN(4.0_r8*PI*Hour/24.0_r8)
  HOD4= COS(4.0_r8*PI*Hour/24.0_r8)
 
  !write(iulog,*)'YMDH:',Year,'|',Month,'|',Day,'|',Hour,'|'
  !write(iulog,*)'Harmonic SINE/COS DOY:',DOY1,'|',DOY2,'|',DOY3,'|',DOY4,'|',DOY5,'|',DOY6,'|',DOY7,'|',DOY8,'|'
  !write(iulog,*)'Harmonic SINE/COS HOD:',HOD1,'|',HOD2,'|',HOD3,'|',HOD4,'|'
  !check the weights: 
  
  !do icol=1,ncol
  !    lattmp = get_rlat_p(lchnk,icol)*57.2957795131_r8
  !    lontmp = get_rlon_p(lchnk,icol)*57.2957795131_r8
  !    if (lattmp < 65.0_r8 .and. lattmp > 63.0_r8 .and. lontmp < 20.0_r8 .and. lontmp > 18.0_r8) then
  !        write(iulog,*),'WDU_ plane ','lat:',lattmp,'lon:',lontmp,WeightsDiurnal_U(1,icol,:,lchnk,DAMLin_ObsInd(1))
  !        write(iulog,*),'W_U train ','lat:',lattmp,'lon:',lontmp,Weights_U(1,icol,:,lchnk,DAMLin_ObsInd(1))    
  !    end if 
  !    
  !    if (lattmp < 25.0_r8 .and. lattmp > 23.0_r8 .and. lontmp < 105.0_r8 .and. lontmp > 103.0_r8) then
  !        write(iulog,*),'WDU_ plane ','lat:',lattmp,'lon:',lontmp,WeightsDiurnal_U(1,icol,:,lchnk,DAMLin_ObsInd(1))
  !        write(iulog,*),'W_U train ','lat:',lattmp,'lon:',lontmp,Weights_U(1,icol,:,lchnk,DAMLin_ObsInd(1))    
  !    end if 
  !end do
  
  !ISSUE2 this is like repeadting or some shit.
  tend_ur_u = 0.0_r8
  tend_vr_u = 0.0_r8
  tend_tr_u = 0.0_r8
  tend_qr_u = 0.0_r8
  tend_or_u = 0.0_r8
  tend_auxr_u = 0.0_r8
  
  tend_ur_v = 0.0_r8
  tend_vr_v = 0.0_r8
  tend_tr_v = 0.0_r8
  tend_qr_v = 0.0_r8
  tend_or_v = 0.0_r8
  tend_auxr_v = 0.0_r8
  
  tend_ur_t = 0.0_r8
  tend_vr_t = 0.0_r8
  tend_tr_t = 0.0_r8
  tend_qr_t = 0.0_r8
  tend_or_t = 0.0_r8
  tend_auxr_t = 0.0_r8
  ttend_val = 0.0_r8
  
  tend_ur_q = 0.0_r8
  tend_vr_q = 0.0_r8
  tend_tr_q = 0.0_r8
  tend_qr_q = 0.0_r8
  tend_or_q = 0.0_r8
  tend_auxr_q = 0.0_r8
  
  itim_old = pbuf_old_tim_idx() 
      
  upwp_clubb_idx     = pbuf_get_index('UPWP')         ! UPWP
  call pbuf_get_field(pbuf, upwp_clubb_idx,    upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))
      
  ustarwec_idx     = pbuf_get_index('ustarwec')         ! Ustar !16|33
  call pbuf_get_field(pbuf, ustarwec_idx,    ustarwec)
      
  pblh_idx     = pbuf_get_index('pblh')         ! PBL height  
  call pbuf_get_field(pbuf, pblh_idx,    pblh) 
      
  fsntoa_idx     = pbuf_get_index('FSNTOA')         !16
  call pbuf_get_field(pbuf, fsntoa_idx,    fsntoa)
      
  fsntoa_idx     = pbuf_get_index('FSNTOA')         !16
  call pbuf_get_field(pbuf, fsntoa_idx,    fsntoa)
      
  flntc_idx     = pbuf_get_index('FLNTC')         !16
  call pbuf_get_field(pbuf, flntc_idx,    flntc)

  ncol =phys_state%ncol
  lchnk=phys_state%lchnk
  !!regression loop:
      
  do icol=1,ncol
  do ipver=1,pver
  
  
     tend_ur_u = Weights_U(1,icol,ipver,lchnk,DP)*((state1%u(icol,1)-mean_preds(1,icol,lchnk,DP))/std_preds(1,icol,lchnk,DP))+ &
                                     Weights_U(2,icol,ipver,lchnk,DP)*((state1%u(icol,2)-mean_preds(2,icol,lchnk,DP))/std_preds(2,icol,lchnk,DP))+ &
                                     Weights_U(3,icol,ipver,lchnk,DP)*((state1%u(icol,3)-mean_preds(3,icol,lchnk,DP))/std_preds(3,icol,lchnk,DP))+ &
                                     Weights_U(4,icol,ipver,lchnk,DP)*((state1%u(icol,4)-mean_preds(4,icol,lchnk,DP))/std_preds(4,icol,lchnk,DP))+ &
                                     Weights_U(5,icol,ipver,lchnk,DP)*((state1%u(icol,5)-mean_preds(5,icol,lchnk,DP))/std_preds(5,icol,lchnk,DP))+ &
                                     Weights_U(6,icol,ipver,lchnk,DP)*((state1%u(icol,6)-mean_preds(6,icol,lchnk,DP))/std_preds(6,icol,lchnk,DP))+ &
                                     Weights_U(7,icol,ipver,lchnk,DP)*((state1%u(icol,7)-mean_preds(7,icol,lchnk,DP))/std_preds(7,icol,lchnk,DP))+ &
                                     Weights_U(8,icol,ipver,lchnk,DP)*((state1%u(icol,8)-mean_preds(8,icol,lchnk,DP))/std_preds(8,icol,lchnk,DP))+ &
                                     Weights_U(9,icol,ipver,lchnk,DP)*((state1%u(icol,9)-mean_preds(9,icol,lchnk,DP))/std_preds(9,icol,lchnk,DP))+ &
                                     Weights_U(10,icol,ipver,lchnk,DP)*((state1%u(icol,10)-mean_preds(10,icol,lchnk,DP))/std_preds(10,icol,lchnk,DP))+ &
                                     Weights_U(11,icol,ipver,lchnk,DP)*((state1%u(icol,11)-mean_preds(11,icol,lchnk,DP))/std_preds(11,icol,lchnk,DP))+ &
                                     Weights_U(12,icol,ipver,lchnk,DP)*((state1%u(icol,12)-mean_preds(12,icol,lchnk,DP))/std_preds(12,icol,lchnk,DP))+ &
                                     Weights_U(13,icol,ipver,lchnk,DP)*((state1%u(icol,13)-mean_preds(13,icol,lchnk,DP))/std_preds(13,icol,lchnk,DP))+ &
                                     Weights_U(14,icol,ipver,lchnk,DP)*((state1%u(icol,14)-mean_preds(14,icol,lchnk,DP))/std_preds(14,icol,lchnk,DP))+ &
                                     Weights_U(15,icol,ipver,lchnk,DP)*((state1%u(icol,15)-mean_preds(15,icol,lchnk,DP))/std_preds(15,icol,lchnk,DP))+ &
                                     Weights_U(16,icol,ipver,lchnk,DP)*((state1%u(icol,16)-mean_preds(16,icol,lchnk,DP))/std_preds(16,icol,lchnk,DP))+ &
                                     Weights_U(17,icol,ipver,lchnk,DP)*((state1%u(icol,17)-mean_preds(17,icol,lchnk,DP))/std_preds(17,icol,lchnk,DP))+ &
                                     Weights_U(18,icol,ipver,lchnk,DP)*((state1%u(icol,18)-mean_preds(18,icol,lchnk,DP))/std_preds(18,icol,lchnk,DP))+ &
                                     Weights_U(19,icol,ipver,lchnk,DP)*((state1%u(icol,19)-mean_preds(19,icol,lchnk,DP))/std_preds(19,icol,lchnk,DP))+ &
                                     Weights_U(20,icol,ipver,lchnk,DP)*((state1%u(icol,20)-mean_preds(20,icol,lchnk,DP))/std_preds(20,icol,lchnk,DP))+ &
                                     Weights_U(21,icol,ipver,lchnk,DP)*((state1%u(icol,21)-mean_preds(21,icol,lchnk,DP))/std_preds(21,icol,lchnk,DP))+ &
                                     Weights_U(22,icol,ipver,lchnk,DP)*((state1%u(icol,22)-mean_preds(22,icol,lchnk,DP))/std_preds(22,icol,lchnk,DP))+ &
                                     Weights_U(23,icol,ipver,lchnk,DP)*((state1%u(icol,23)-mean_preds(23,icol,lchnk,DP))/std_preds(23,icol,lchnk,DP))+ &
                                     Weights_U(24,icol,ipver,lchnk,DP)*((state1%u(icol,24)-mean_preds(24,icol,lchnk,DP))/std_preds(24,icol,lchnk,DP))+ &
                                     Weights_U(25,icol,ipver,lchnk,DP)*((state1%u(icol,25)-mean_preds(25,icol,lchnk,DP))/std_preds(25,icol,lchnk,DP))+ &
                                     Weights_U(26,icol,ipver,lchnk,DP)*((state1%u(icol,26)-mean_preds(26,icol,lchnk,DP))/std_preds(26,icol,lchnk,DP))+ &
                                     Weights_U(27,icol,ipver,lchnk,DP)*((state1%u(icol,27)-mean_preds(27,icol,lchnk,DP))/std_preds(27,icol,lchnk,DP))+ &
                                     Weights_U(28,icol,ipver,lchnk,DP)*((state1%u(icol,28)-mean_preds(28,icol,lchnk,DP))/std_preds(28,icol,lchnk,DP))+ &
                                     Weights_U(29,icol,ipver,lchnk,DP)*((state1%u(icol,29)-mean_preds(29,icol,lchnk,DP))/std_preds(29,icol,lchnk,DP))+ &
                                     Weights_U(30,icol,ipver,lchnk,DP)*((state1%u(icol,30)-mean_preds(30,icol,lchnk,DP))/std_preds(30,icol,lchnk,DP))+ &
                                     Weights_U(31,icol,ipver,lchnk,DP)*((state1%u(icol,31)-mean_preds(31,icol,lchnk,DP))/std_preds(31,icol,lchnk,DP))+ &
                                     Weights_U(32,icol,ipver,lchnk,DP)*((state1%u(icol,32)-mean_preds(32,icol,lchnk,DP))/std_preds(32,icol,lchnk,DP)) 
                                     
                                     
     tend_vr_u = Weights_U(33,icol,ipver,lchnk,DP)*((state1%v(icol,1)-mean_preds(33,icol,lchnk,DP))/std_preds(33,icol,lchnk,DP))+ &
                                     Weights_U(34,icol,ipver,lchnk,DP)*((state1%v(icol,2)-mean_preds(34,icol,lchnk,DP))/std_preds(34,icol,lchnk,DP))+ &
                                     Weights_U(35,icol,ipver,lchnk,DP)*((state1%v(icol,3)-mean_preds(35,icol,lchnk,DP))/std_preds(35,icol,lchnk,DP))+ &
                                     Weights_U(36,icol,ipver,lchnk,DP)*((state1%v(icol,4)-mean_preds(36,icol,lchnk,DP))/std_preds(36,icol,lchnk,DP))+ &
                                     Weights_U(37,icol,ipver,lchnk,DP)*((state1%v(icol,5)-mean_preds(37,icol,lchnk,DP))/std_preds(37,icol,lchnk,DP))+ &
                                     Weights_U(38,icol,ipver,lchnk,DP)*((state1%v(icol,6)-mean_preds(38,icol,lchnk,DP))/std_preds(38,icol,lchnk,DP))+ &
                                     Weights_U(39,icol,ipver,lchnk,DP)*((state1%v(icol,7)-mean_preds(39,icol,lchnk,DP))/std_preds(39,icol,lchnk,DP))+ &
                                     Weights_U(40,icol,ipver,lchnk,DP)*((state1%v(icol,8)-mean_preds(40,icol,lchnk,DP))/std_preds(40,icol,lchnk,DP))+ &
                                     Weights_U(41,icol,ipver,lchnk,DP)*((state1%v(icol,9)-mean_preds(41,icol,lchnk,DP))/std_preds(41,icol,lchnk,DP))+ &
                                     Weights_U(42,icol,ipver,lchnk,DP)*((state1%v(icol,10)-mean_preds(42,icol,lchnk,DP))/std_preds(42,icol,lchnk,DP))+ &
                                     Weights_U(43,icol,ipver,lchnk,DP)*((state1%v(icol,11)-mean_preds(43,icol,lchnk,DP))/std_preds(43,icol,lchnk,DP))+ &
                                     Weights_U(44,icol,ipver,lchnk,DP)*((state1%v(icol,12)-mean_preds(44,icol,lchnk,DP))/std_preds(44,icol,lchnk,DP))+ &
                                     Weights_U(45,icol,ipver,lchnk,DP)*((state1%v(icol,13)-mean_preds(45,icol,lchnk,DP))/std_preds(45,icol,lchnk,DP))+ &
                                     Weights_U(46,icol,ipver,lchnk,DP)*((state1%v(icol,14)-mean_preds(46,icol,lchnk,DP))/std_preds(46,icol,lchnk,DP))+ &
                                     Weights_U(47,icol,ipver,lchnk,DP)*((state1%v(icol,15)-mean_preds(47,icol,lchnk,DP))/std_preds(47,icol,lchnk,DP))+ &
                                     Weights_U(48,icol,ipver,lchnk,DP)*((state1%v(icol,16)-mean_preds(48,icol,lchnk,DP))/std_preds(48,icol,lchnk,DP))+ &
                                     Weights_U(49,icol,ipver,lchnk,DP)*((state1%v(icol,17)-mean_preds(49,icol,lchnk,DP))/std_preds(49,icol,lchnk,DP))+ &
                                     Weights_U(50,icol,ipver,lchnk,DP)*((state1%v(icol,18)-mean_preds(50,icol,lchnk,DP))/std_preds(50,icol,lchnk,DP))+ &
                                     Weights_U(51,icol,ipver,lchnk,DP)*((state1%v(icol,19)-mean_preds(51,icol,lchnk,DP))/std_preds(51,icol,lchnk,DP))+ &
                                     Weights_U(52,icol,ipver,lchnk,DP)*((state1%v(icol,20)-mean_preds(52,icol,lchnk,DP))/std_preds(52,icol,lchnk,DP))+ &
                                     Weights_U(53,icol,ipver,lchnk,DP)*((state1%v(icol,21)-mean_preds(53,icol,lchnk,DP))/std_preds(53,icol,lchnk,DP))+ &
                                     Weights_U(54,icol,ipver,lchnk,DP)*((state1%v(icol,22)-mean_preds(54,icol,lchnk,DP))/std_preds(54,icol,lchnk,DP))+ &
                                     Weights_U(55,icol,ipver,lchnk,DP)*((state1%v(icol,23)-mean_preds(55,icol,lchnk,DP))/std_preds(55,icol,lchnk,DP))+ &
                                     Weights_U(56,icol,ipver,lchnk,DP)*((state1%v(icol,24)-mean_preds(56,icol,lchnk,DP))/std_preds(56,icol,lchnk,DP))+ &
                                     Weights_U(57,icol,ipver,lchnk,DP)*((state1%v(icol,25)-mean_preds(57,icol,lchnk,DP))/std_preds(57,icol,lchnk,DP))+ &
                                     Weights_U(58,icol,ipver,lchnk,DP)*((state1%v(icol,26)-mean_preds(58,icol,lchnk,DP))/std_preds(58,icol,lchnk,DP))+ &
                                     Weights_U(59,icol,ipver,lchnk,DP)*((state1%v(icol,27)-mean_preds(59,icol,lchnk,DP))/std_preds(59,icol,lchnk,DP))+ &
                                     Weights_U(60,icol,ipver,lchnk,DP)*((state1%v(icol,28)-mean_preds(60,icol,lchnk,DP))/std_preds(60,icol,lchnk,DP))+ &
                                     Weights_U(61,icol,ipver,lchnk,DP)*((state1%v(icol,29)-mean_preds(61,icol,lchnk,DP))/std_preds(61,icol,lchnk,DP))+ &
                                     Weights_U(62,icol,ipver,lchnk,DP)*((state1%v(icol,30)-mean_preds(62,icol,lchnk,DP))/std_preds(62,icol,lchnk,DP))+ &
                                     Weights_U(63,icol,ipver,lchnk,DP)*((state1%v(icol,31)-mean_preds(63,icol,lchnk,DP))/std_preds(63,icol,lchnk,DP))+ &
                                     Weights_U(64,icol,ipver,lchnk,DP)*((state1%v(icol,32)-mean_preds(64,icol,lchnk,DP))/std_preds(64,icol,lchnk,DP))
                                                                  
                                
      tend_tr_u = Weights_U(65,icol,ipver,lchnk,DP)*((state1%t(icol,1)-mean_preds(65,icol,lchnk,DP))/std_preds(65,icol,lchnk,DP))+ &
                                     Weights_U(66,icol,ipver,lchnk,DP)*((state1%t(icol,2)-mean_preds(66,icol,lchnk,DP))/std_preds(66,icol,lchnk,DP))+ &
                                     Weights_U(67,icol,ipver,lchnk,DP)*((state1%t(icol,3)-mean_preds(67,icol,lchnk,DP))/std_preds(67,icol,lchnk,DP))+ &
                                     Weights_U(68,icol,ipver,lchnk,DP)*((state1%t(icol,4)-mean_preds(68,icol,lchnk,DP))/std_preds(68,icol,lchnk,DP))+ &
                                     Weights_U(69,icol,ipver,lchnk,DP)*((state1%t(icol,5)-mean_preds(69,icol,lchnk,DP))/std_preds(69,icol,lchnk,DP))+ &
                                     Weights_U(70,icol,ipver,lchnk,DP)*((state1%t(icol,6)-mean_preds(70,icol,lchnk,DP))/std_preds(70,icol,lchnk,DP))+ &
                                     Weights_U(71,icol,ipver,lchnk,DP)*((state1%t(icol,7)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_U(71,icol,ipver,lchnk,DP)*((state1%t(icol,8)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_U(73,icol,ipver,lchnk,DP)*((state1%t(icol,9)-mean_preds(73,icol,lchnk,DP))/std_preds(73,icol,lchnk,DP))+ &
                                     Weights_U(74,icol,ipver,lchnk,DP)*((state1%t(icol,10)-mean_preds(74,icol,lchnk,DP))/std_preds(74,icol,lchnk,DP))+ &
                                     Weights_U(75,icol,ipver,lchnk,DP)*((state1%t(icol,11)-mean_preds(75,icol,lchnk,DP))/std_preds(75,icol,lchnk,DP))+ &
                                     Weights_U(76,icol,ipver,lchnk,DP)*((state1%t(icol,12)-mean_preds(76,icol,lchnk,DP))/std_preds(76,icol,lchnk,DP))+ &
                                     Weights_U(77,icol,ipver,lchnk,DP)*((state1%t(icol,13)-mean_preds(77,icol,lchnk,DP))/std_preds(77,icol,lchnk,DP))+ &
                                     Weights_U(78,icol,ipver,lchnk,DP)*((state1%t(icol,14)-mean_preds(78,icol,lchnk,DP))/std_preds(78,icol,lchnk,DP))+ &
                                     Weights_U(79,icol,ipver,lchnk,DP)*((state1%t(icol,15)-mean_preds(79,icol,lchnk,DP))/std_preds(79,icol,lchnk,DP))+ &
                                     Weights_U(80,icol,ipver,lchnk,DP)*((state1%t(icol,16)-mean_preds(80,icol,lchnk,DP))/std_preds(80,icol,lchnk,DP))+ &
                                     Weights_U(81,icol,ipver,lchnk,DP)*((state1%t(icol,17)-mean_preds(81,icol,lchnk,DP))/std_preds(81,icol,lchnk,DP))+ &
                                     Weights_U(82,icol,ipver,lchnk,DP)*((state1%t(icol,18)-mean_preds(82,icol,lchnk,DP))/std_preds(82,icol,lchnk,DP))+ &
                                     Weights_U(83,icol,ipver,lchnk,DP)*((state1%t(icol,19)-mean_preds(83,icol,lchnk,DP))/std_preds(83,icol,lchnk,DP))+ &
                                     Weights_U(84,icol,ipver,lchnk,DP)*((state1%t(icol,20)-mean_preds(84,icol,lchnk,DP))/std_preds(84,icol,lchnk,DP))+ &
                                     Weights_U(85,icol,ipver,lchnk,DP)*((state1%t(icol,21)-mean_preds(85,icol,lchnk,DP))/std_preds(85,icol,lchnk,DP))+ &
                                     Weights_U(86,icol,ipver,lchnk,DP)*((state1%t(icol,22)-mean_preds(86,icol,lchnk,DP))/std_preds(86,icol,lchnk,DP))+ &
                                     Weights_U(87,icol,ipver,lchnk,DP)*((state1%t(icol,23)-mean_preds(87,icol,lchnk,DP))/std_preds(87,icol,lchnk,DP))+ &
                                     Weights_U(88,icol,ipver,lchnk,DP)*((state1%t(icol,24)-mean_preds(88,icol,lchnk,DP))/std_preds(88,icol,lchnk,DP))+ &
                                     Weights_U(89,icol,ipver,lchnk,DP)*((state1%t(icol,25)-mean_preds(89,icol,lchnk,DP))/std_preds(89,icol,lchnk,DP))+ &
                                     Weights_U(90,icol,ipver,lchnk,DP)*((state1%t(icol,26)-mean_preds(90,icol,lchnk,DP))/std_preds(90,icol,lchnk,DP))+ &
                                     Weights_U(91,icol,ipver,lchnk,DP)*((state1%t(icol,27)-mean_preds(91,icol,lchnk,DP))/std_preds(91,icol,lchnk,DP))+ &
                                     Weights_U(92,icol,ipver,lchnk,DP)*((state1%t(icol,28)-mean_preds(92,icol,lchnk,DP))/std_preds(92,icol,lchnk,DP))+ &
                                     Weights_U(93,icol,ipver,lchnk,DP)*((state1%t(icol,29)-mean_preds(93,icol,lchnk,DP))/std_preds(93,icol,lchnk,DP))+ &
                                     Weights_U(94,icol,ipver,lchnk,DP)*((state1%t(icol,30)-mean_preds(94,icol,lchnk,DP))/std_preds(94,icol,lchnk,DP))+ &
                                     Weights_U(95,icol,ipver,lchnk,DP)*((state1%t(icol,31)-mean_preds(95,icol,lchnk,DP))/std_preds(95,icol,lchnk,DP))+ &
                                     Weights_U(96,icol,ipver,lchnk,DP)*((state1%t(icol,32)-mean_preds(96,icol,lchnk,DP))/std_preds(96,icol,lchnk,DP))
                                     
      tend_qr_u = Weights_U(97,icol,ipver,lchnk,DP)*((state1%q(icol,1,1)-mean_preds(97,icol,lchnk,DP))/std_preds(97,icol,lchnk,DP))+ &
                                     Weights_U(98,icol,ipver,lchnk,DP)*((state1%q(icol,2,1)-mean_preds(98,icol,lchnk,DP))/std_preds(98,icol,lchnk,DP))+ &
                                     Weights_U(99,icol,ipver,lchnk,DP)*((state1%q(icol,3,1)-mean_preds(99,icol,lchnk,DP))/std_preds(99,icol,lchnk,DP))+ &
                                     Weights_U(100,icol,ipver,lchnk,DP)*((state1%q(icol,4,1)-mean_preds(100,icol,lchnk,DP))/std_preds(100,icol,lchnk,DP))+ &
                                     Weights_U(101,icol,ipver,lchnk,DP)*((state1%q(icol,5,1)-mean_preds(101,icol,lchnk,DP))/std_preds(101,icol,lchnk,DP))+ &
                                     Weights_U(102,icol,ipver,lchnk,DP)*((state1%q(icol,6,1)-mean_preds(102,icol,lchnk,DP))/std_preds(102,icol,lchnk,DP))+ &
                                     Weights_U(103,icol,ipver,lchnk,DP)*((state1%q(icol,7,1)-mean_preds(103,icol,lchnk,DP))/std_preds(103,icol,lchnk,DP))+ &
                                     Weights_U(104,icol,ipver,lchnk,DP)*((state1%q(icol,8,1)-mean_preds(104,icol,lchnk,DP))/std_preds(104,icol,lchnk,DP))+ &
                                     Weights_U(105,icol,ipver,lchnk,DP)*((state1%q(icol,9,1)-mean_preds(105,icol,lchnk,DP))/std_preds(105,icol,lchnk,DP))+ &
                                     Weights_U(106,icol,ipver,lchnk,DP)*((state1%q(icol,10,1)-mean_preds(106,icol,lchnk,DP))/std_preds(106,icol,lchnk,DP))+ &
                                     Weights_U(107,icol,ipver,lchnk,DP)*((state1%q(icol,11,1)-mean_preds(107,icol,lchnk,DP))/std_preds(107,icol,lchnk,DP))+ &
                                     Weights_U(108,icol,ipver,lchnk,DP)*((state1%q(icol,12,1)-mean_preds(108,icol,lchnk,DP))/std_preds(108,icol,lchnk,DP))+ &
                                     Weights_U(109,icol,ipver,lchnk,DP)*((state1%q(icol,13,1)-mean_preds(109,icol,lchnk,DP))/std_preds(109,icol,lchnk,DP))+ &
                                     Weights_U(110,icol,ipver,lchnk,DP)*((state1%q(icol,14,1)-mean_preds(110,icol,lchnk,DP))/std_preds(110,icol,lchnk,DP))+ &
                                     Weights_U(111,icol,ipver,lchnk,DP)*((state1%q(icol,15,1)-mean_preds(111,icol,lchnk,DP))/std_preds(111,icol,lchnk,DP))+ &
                                     Weights_U(112,icol,ipver,lchnk,DP)*((state1%q(icol,16,1)-mean_preds(112,icol,lchnk,DP))/std_preds(112,icol,lchnk,DP))+ &
                                     Weights_U(113,icol,ipver,lchnk,DP)*((state1%q(icol,17,1)-mean_preds(113,icol,lchnk,DP))/std_preds(113,icol,lchnk,DP))+ &
                                     Weights_U(114,icol,ipver,lchnk,DP)*((state1%q(icol,18,1)-mean_preds(114,icol,lchnk,DP))/std_preds(114,icol,lchnk,DP))+ &
                                     Weights_U(115,icol,ipver,lchnk,DP)*((state1%q(icol,19,1)-mean_preds(115,icol,lchnk,DP))/std_preds(115,icol,lchnk,DP))+ &
                                     Weights_U(116,icol,ipver,lchnk,DP)*((state1%q(icol,20,1)-mean_preds(116,icol,lchnk,DP))/std_preds(116,icol,lchnk,DP))+ &
                                     Weights_U(117,icol,ipver,lchnk,DP)*((state1%q(icol,21,1)-mean_preds(117,icol,lchnk,DP))/std_preds(117,icol,lchnk,DP))+ &
                                     Weights_U(118,icol,ipver,lchnk,DP)*((state1%q(icol,22,1)-mean_preds(118,icol,lchnk,DP))/std_preds(118,icol,lchnk,DP))+ &
                                     Weights_U(119,icol,ipver,lchnk,DP)*((state1%q(icol,23,1)-mean_preds(119,icol,lchnk,DP))/std_preds(119,icol,lchnk,DP))+ &
                                     Weights_U(120,icol,ipver,lchnk,DP)*((state1%q(icol,24,1)-mean_preds(120,icol,lchnk,DP))/std_preds(120,icol,lchnk,DP))+ &
                                     Weights_U(121,icol,ipver,lchnk,DP)*((state1%q(icol,25,1)-mean_preds(121,icol,lchnk,DP))/std_preds(121,icol,lchnk,DP))+ &
                                     Weights_U(122,icol,ipver,lchnk,DP)*((state1%q(icol,26,1)-mean_preds(122,icol,lchnk,DP))/std_preds(122,icol,lchnk,DP))+ &
                                     Weights_U(123,icol,ipver,lchnk,DP)*((state1%q(icol,27,1)-mean_preds(123,icol,lchnk,DP))/std_preds(123,icol,lchnk,DP))+ &
                                     Weights_U(124,icol,ipver,lchnk,DP)*((state1%q(icol,28,1)-mean_preds(124,icol,lchnk,DP))/std_preds(124,icol,lchnk,DP))+ &
                                     Weights_U(125,icol,ipver,lchnk,DP)*((state1%q(icol,29,1)-mean_preds(125,icol,lchnk,DP))/std_preds(125,icol,lchnk,DP))+ &
                                     Weights_U(126,icol,ipver,lchnk,DP)*((state1%q(icol,30,1)-mean_preds(126,icol,lchnk,DP))/std_preds(126,icol,lchnk,DP))+ &
                                     Weights_U(127,icol,ipver,lchnk,DP)*((state1%q(icol,31,1)-mean_preds(127,icol,lchnk,DP))/std_preds(127,icol,lchnk,DP))+ &
                                     Weights_U(128,icol,ipver,lchnk,DP)*((state1%q(icol,32,1)-mean_preds(128,icol,lchnk,DP))/std_preds(128,icol,lchnk,DP))
                                     
      tend_or_u = Weights_U(129,icol,ipver,lchnk,DP)*((state1%omega(icol,1)-mean_preds(129,icol,lchnk,DP))/std_preds(129,icol,lchnk,DP))+ &
                                     Weights_U(130,icol,ipver,lchnk,DP)*((state1%omega(icol,2)-mean_preds(130,icol,lchnk,DP))/std_preds(130,icol,lchnk,DP))+ &
                                     Weights_U(131,icol,ipver,lchnk,DP)*((state1%omega(icol,3)-mean_preds(131,icol,lchnk,DP))/std_preds(131,icol,lchnk,DP))+ &
                                     Weights_U(132,icol,ipver,lchnk,DP)*((state1%omega(icol,4)-mean_preds(132,icol,lchnk,DP))/std_preds(132,icol,lchnk,DP))+ &
                                     Weights_U(133,icol,ipver,lchnk,DP)*((state1%omega(icol,5)-mean_preds(133,icol,lchnk,DP))/std_preds(133,icol,lchnk,DP))+ &
                                     Weights_U(134,icol,ipver,lchnk,DP)*((state1%omega(icol,6)-mean_preds(134,icol,lchnk,DP))/std_preds(134,icol,lchnk,DP))+ &
                                     Weights_U(135,icol,ipver,lchnk,DP)*((state1%omega(icol,7)-mean_preds(135,icol,lchnk,DP))/std_preds(135,icol,lchnk,DP))+ &
                                     Weights_U(136,icol,ipver,lchnk,DP)*((state1%omega(icol,8)-mean_preds(136,icol,lchnk,DP))/std_preds(136,icol,lchnk,DP))+ &
                                     Weights_U(137,icol,ipver,lchnk,DP)*((state1%omega(icol,9)-mean_preds(137,icol,lchnk,DP))/std_preds(137,icol,lchnk,DP))+ &
                                     Weights_U(138,icol,ipver,lchnk,DP)*((state1%omega(icol,10)-mean_preds(138,icol,lchnk,DP))/std_preds(138,icol,lchnk,DP))+ &
                                     Weights_U(139,icol,ipver,lchnk,DP)*((state1%omega(icol,11)-mean_preds(139,icol,lchnk,DP))/std_preds(139,icol,lchnk,DP))+ &
                                     Weights_U(140,icol,ipver,lchnk,DP)*((state1%omega(icol,12)-mean_preds(140,icol,lchnk,DP))/std_preds(140,icol,lchnk,DP))+ &
                                     Weights_U(141,icol,ipver,lchnk,DP)*((state1%omega(icol,13)-mean_preds(141,icol,lchnk,DP))/std_preds(141,icol,lchnk,DP))+ &
                                     Weights_U(142,icol,ipver,lchnk,DP)*((state1%omega(icol,14)-mean_preds(142,icol,lchnk,DP))/std_preds(142,icol,lchnk,DP))+ &
                                     Weights_U(143,icol,ipver,lchnk,DP)*((state1%omega(icol,15)-mean_preds(143,icol,lchnk,DP))/std_preds(143,icol,lchnk,DP))+ &
                                     Weights_U(144,icol,ipver,lchnk,DP)*((state1%omega(icol,16)-mean_preds(144,icol,lchnk,DP))/std_preds(144,icol,lchnk,DP))+ &
                                     Weights_U(145,icol,ipver,lchnk,DP)*((state1%omega(icol,17)-mean_preds(145,icol,lchnk,DP))/std_preds(145,icol,lchnk,DP))+ &
                                     Weights_U(146,icol,ipver,lchnk,DP)*((state1%omega(icol,18)-mean_preds(146,icol,lchnk,DP))/std_preds(146,icol,lchnk,DP))+ &
                                     Weights_U(147,icol,ipver,lchnk,DP)*((state1%omega(icol,19)-mean_preds(147,icol,lchnk,DP))/std_preds(147,icol,lchnk,DP))+ &
                                     Weights_U(148,icol,ipver,lchnk,DP)*((state1%omega(icol,20)-mean_preds(148,icol,lchnk,DP))/std_preds(148,icol,lchnk,DP))+ &
                                     Weights_U(149,icol,ipver,lchnk,DP)*((state1%omega(icol,21)-mean_preds(149,icol,lchnk,DP))/std_preds(149,icol,lchnk,DP))+ &
                                     Weights_U(150,icol,ipver,lchnk,DP)*((state1%omega(icol,22)-mean_preds(150,icol,lchnk,DP))/std_preds(150,icol,lchnk,DP))+ &
                                     Weights_U(151,icol,ipver,lchnk,DP)*((state1%omega(icol,23)-mean_preds(151,icol,lchnk,DP))/std_preds(151,icol,lchnk,DP))+ &
                                     Weights_U(152,icol,ipver,lchnk,DP)*((state1%omega(icol,24)-mean_preds(152,icol,lchnk,DP))/std_preds(152,icol,lchnk,DP))+ &
                                     Weights_U(153,icol,ipver,lchnk,DP)*((state1%omega(icol,25)-mean_preds(153,icol,lchnk,DP))/std_preds(153,icol,lchnk,DP))+ &
                                     Weights_U(154,icol,ipver,lchnk,DP)*((state1%omega(icol,26)-mean_preds(154,icol,lchnk,DP))/std_preds(154,icol,lchnk,DP))+ &
                                     Weights_U(155,icol,ipver,lchnk,DP)*((state1%omega(icol,27)-mean_preds(155,icol,lchnk,DP))/std_preds(155,icol,lchnk,DP))+ &
                                     Weights_U(156,icol,ipver,lchnk,DP)*((state1%omega(icol,28)-mean_preds(156,icol,lchnk,DP))/std_preds(156,icol,lchnk,DP))+ &
                                     Weights_U(157,icol,ipver,lchnk,DP)*((state1%omega(icol,29)-mean_preds(157,icol,lchnk,DP))/std_preds(157,icol,lchnk,DP))+ &
                                     Weights_U(158,icol,ipver,lchnk,DP)*((state1%omega(icol,30)-mean_preds(158,icol,lchnk,DP))/std_preds(158,icol,lchnk,DP))+ &
                                     Weights_U(159,icol,ipver,lchnk,DP)*((state1%omega(icol,31)-mean_preds(159,icol,lchnk,DP))/std_preds(159,icol,lchnk,DP))+ &
                                     Weights_U(160,icol,ipver,lchnk,DP)*((state1%omega(icol,32)-mean_preds(160,icol,lchnk,DP))/std_preds(160,icol,lchnk,DP))
                                     
                                     
      tend_auxr_u = Weights_U(161,icol,ipver,lchnk,DP)*((upwp(icol,24)-mean_preds(161,icol,lchnk,DP))/std_preds(161,icol,lchnk,DP))+ &
                                     Weights_U(162,icol,ipver,lchnk,DP)*((upwp(icol,25)-mean_preds(162,icol,lchnk,DP))/std_preds(162,icol,lchnk,DP))+ &
                                     Weights_U(163,icol,ipver,lchnk,DP)*((upwp(icol,26)-mean_preds(163,icol,lchnk,DP))/std_preds(163,icol,lchnk,DP))+ &
                                     Weights_U(164,icol,ipver,lchnk,DP)*((upwp(icol,27)-mean_preds(164,icol,lchnk,DP))/std_preds(164,icol,lchnk,DP))+ &
                                     Weights_U(165,icol,ipver,lchnk,DP)*((upwp(icol,28)-mean_preds(165,icol,lchnk,DP))/std_preds(165,icol,lchnk,DP))+ &
                                     Weights_U(166,icol,ipver,lchnk,DP)*((upwp(icol,29)-mean_preds(166,icol,lchnk,DP))/std_preds(166,icol,lchnk,DP))+ &
                                     Weights_U(167,icol,ipver,lchnk,DP)*((upwp(icol,30)-mean_preds(167,icol,lchnk,DP))/std_preds(167,icol,lchnk,DP))+ &
                                     Weights_U(168,icol,ipver,lchnk,DP)*((upwp(icol,31)-mean_preds(168,icol,lchnk,DP))/std_preds(168,icol,lchnk,DP))+ &
                                     Weights_U(169,icol,ipver,lchnk,DP)*((upwp(icol,32)-mean_preds(169,icol,lchnk,DP))/std_preds(169,icol,lchnk,DP))+ &
                                     Weights_U(170,icol,ipver,lchnk,DP)*((upwp(icol,33)-mean_preds(170,icol,lchnk,DP))/std_preds(170,icol,lchnk,DP))+ &
                                     Weights_U(171,icol,ipver,lchnk,DP)*((cam_in%wsx(icol)-mean_preds(171,icol,lchnk,DP))/std_preds(171,icol,lchnk,DP))+ &
                                     Weights_U(172,icol,ipver,lchnk,DP)*((cam_in%wsy(icol)-mean_preds(172,icol,lchnk,DP))/std_preds(172,icol,lchnk,DP))+ &
                                     Weights_U(173,icol,ipver,lchnk,DP)*((flntc(icol)-mean_preds(173,icol,lchnk,DP))/std_preds(173,icol,lchnk,DP))+ &
                                     Weights_U(174,icol,ipver,lchnk,DP)*((fsntoa(icol)-mean_preds(174,icol,lchnk,DP))/std_preds(174,icol,lchnk,DP))+ &
                                     Weights_U(175,icol,ipver,lchnk,DP)*((pblh(icol)-mean_preds(175,icol,lchnk,DP))/std_preds(175,icol,lchnk,DP))+ &
                                     Weights_U(176,icol,ipver,lchnk,DP)*((cam_out%tbot(icol)-mean_preds(176,icol,lchnk,DP))/std_preds(176,icol,lchnk,DP))+ &
                                     Weights_U(177,icol,ipver,lchnk,DP)*((ustarwec(icol)-mean_preds(177,icol,lchnk,DP))/std_preds(177,icol,lchnk,DP))+ &
                                     Weights_U(178,icol,ipver,lchnk,DP)*1.5_r8*DOY2+ &
                                     Weights_U(179,icol,ipver,lchnk,DP)*1.5_r8*DOY1+ &
                                     Weights_U(180,icol,ipver,lchnk,DP)*1.5_r8*HOD2+ &
                                     Weights_U(181,icol,ipver,lchnk,DP)*1.5_r8*HOD1+ &
                                     Weights_U(182,icol,ipver,lchnk,DP)
      
      
      if (abs(tend_ur_u + tend_vr_u + tend_tr_u + tend_qr_u + tend_or_u + tend_auxr_u) > threshy) then
          DAMLin_Ustep(icol,ipver,lchnk) = sign(threshy, tend_ur_u + tend_vr_u + tend_tr_u + tend_qr_u + tend_or_u + tend_auxr_u)
      else 
          DAMLin_Ustep(icol,ipver,lchnk) = tend_ur_u + tend_vr_u + tend_tr_u + tend_qr_u + tend_or_u + tend_auxr_u
      end if       
      
      tend_ur_v = Weights_V(1,icol,ipver,lchnk,DP)*((state1%u(icol,1)-mean_preds(1,icol,lchnk,DP))/std_preds(1,icol,lchnk,DP))+ &
                                     Weights_V(2,icol,ipver,lchnk,DP)*((state1%u(icol,2)-mean_preds(2,icol,lchnk,DP))/std_preds(2,icol,lchnk,DP))+ &
                                     Weights_V(3,icol,ipver,lchnk,DP)*((state1%u(icol,3)-mean_preds(3,icol,lchnk,DP))/std_preds(3,icol,lchnk,DP))+ &
                                     Weights_V(4,icol,ipver,lchnk,DP)*((state1%u(icol,4)-mean_preds(4,icol,lchnk,DP))/std_preds(4,icol,lchnk,DP))+ &
                                     Weights_V(5,icol,ipver,lchnk,DP)*((state1%u(icol,5)-mean_preds(5,icol,lchnk,DP))/std_preds(5,icol,lchnk,DP))+ &
                                     Weights_V(6,icol,ipver,lchnk,DP)*((state1%u(icol,6)-mean_preds(6,icol,lchnk,DP))/std_preds(6,icol,lchnk,DP))+ &
                                     Weights_V(7,icol,ipver,lchnk,DP)*((state1%u(icol,7)-mean_preds(7,icol,lchnk,DP))/std_preds(7,icol,lchnk,DP))+ &
                                     Weights_V(8,icol,ipver,lchnk,DP)*((state1%u(icol,8)-mean_preds(8,icol,lchnk,DP))/std_preds(8,icol,lchnk,DP))+ &
                                     Weights_V(9,icol,ipver,lchnk,DP)*((state1%u(icol,9)-mean_preds(9,icol,lchnk,DP))/std_preds(9,icol,lchnk,DP))+ &
                                     Weights_V(10,icol,ipver,lchnk,DP)*((state1%u(icol,10)-mean_preds(10,icol,lchnk,DP))/std_preds(10,icol,lchnk,DP))+ &
                                     Weights_V(11,icol,ipver,lchnk,DP)*((state1%u(icol,11)-mean_preds(11,icol,lchnk,DP))/std_preds(11,icol,lchnk,DP))+ &
                                     Weights_V(12,icol,ipver,lchnk,DP)*((state1%u(icol,12)-mean_preds(12,icol,lchnk,DP))/std_preds(12,icol,lchnk,DP))+ &
                                     Weights_V(13,icol,ipver,lchnk,DP)*((state1%u(icol,13)-mean_preds(13,icol,lchnk,DP))/std_preds(13,icol,lchnk,DP))+ &
                                     Weights_V(14,icol,ipver,lchnk,DP)*((state1%u(icol,14)-mean_preds(14,icol,lchnk,DP))/std_preds(14,icol,lchnk,DP))+ &
                                     Weights_V(15,icol,ipver,lchnk,DP)*((state1%u(icol,15)-mean_preds(15,icol,lchnk,DP))/std_preds(15,icol,lchnk,DP))+ &
                                     Weights_V(16,icol,ipver,lchnk,DP)*((state1%u(icol,16)-mean_preds(16,icol,lchnk,DP))/std_preds(16,icol,lchnk,DP))+ &
                                     Weights_V(17,icol,ipver,lchnk,DP)*((state1%u(icol,17)-mean_preds(17,icol,lchnk,DP))/std_preds(17,icol,lchnk,DP))+ &
                                     Weights_V(18,icol,ipver,lchnk,DP)*((state1%u(icol,18)-mean_preds(18,icol,lchnk,DP))/std_preds(18,icol,lchnk,DP))+ &
                                     Weights_V(19,icol,ipver,lchnk,DP)*((state1%u(icol,19)-mean_preds(19,icol,lchnk,DP))/std_preds(19,icol,lchnk,DP))+ &
                                     Weights_V(20,icol,ipver,lchnk,DP)*((state1%u(icol,20)-mean_preds(20,icol,lchnk,DP))/std_preds(20,icol,lchnk,DP))+ &
                                     Weights_V(21,icol,ipver,lchnk,DP)*((state1%u(icol,21)-mean_preds(21,icol,lchnk,DP))/std_preds(21,icol,lchnk,DP))+ &
                                     Weights_V(22,icol,ipver,lchnk,DP)*((state1%u(icol,22)-mean_preds(22,icol,lchnk,DP))/std_preds(22,icol,lchnk,DP))+ &
                                     Weights_V(23,icol,ipver,lchnk,DP)*((state1%u(icol,23)-mean_preds(23,icol,lchnk,DP))/std_preds(23,icol,lchnk,DP))+ &
                                     Weights_V(24,icol,ipver,lchnk,DP)*((state1%u(icol,24)-mean_preds(24,icol,lchnk,DP))/std_preds(24,icol,lchnk,DP))+ &
                                     Weights_V(25,icol,ipver,lchnk,DP)*((state1%u(icol,25)-mean_preds(25,icol,lchnk,DP))/std_preds(25,icol,lchnk,DP))+ &
                                     Weights_V(26,icol,ipver,lchnk,DP)*((state1%u(icol,26)-mean_preds(26,icol,lchnk,DP))/std_preds(26,icol,lchnk,DP))+ &
                                     Weights_V(27,icol,ipver,lchnk,DP)*((state1%u(icol,27)-mean_preds(27,icol,lchnk,DP))/std_preds(27,icol,lchnk,DP))+ &
                                     Weights_V(28,icol,ipver,lchnk,DP)*((state1%u(icol,28)-mean_preds(28,icol,lchnk,DP))/std_preds(28,icol,lchnk,DP))+ &
                                     Weights_V(29,icol,ipver,lchnk,DP)*((state1%u(icol,29)-mean_preds(29,icol,lchnk,DP))/std_preds(29,icol,lchnk,DP))+ &
                                     Weights_V(30,icol,ipver,lchnk,DP)*((state1%u(icol,30)-mean_preds(30,icol,lchnk,DP))/std_preds(30,icol,lchnk,DP))+ &
                                     Weights_V(31,icol,ipver,lchnk,DP)*((state1%u(icol,31)-mean_preds(31,icol,lchnk,DP))/std_preds(31,icol,lchnk,DP))+ &
                                     Weights_V(32,icol,ipver,lchnk,DP)*((state1%u(icol,32)-mean_preds(32,icol,lchnk,DP))/std_preds(32,icol,lchnk,DP)) 
                                     
                                     
     tend_vr_v = Weights_V(33,icol,ipver,lchnk,DP)*((state1%v(icol,1)-mean_preds(33,icol,lchnk,DP))/std_preds(33,icol,lchnk,DP))+ &
                                     Weights_V(34,icol,ipver,lchnk,DP)*((state1%v(icol,2)-mean_preds(34,icol,lchnk,DP))/std_preds(34,icol,lchnk,DP))+ &
                                     Weights_V(35,icol,ipver,lchnk,DP)*((state1%v(icol,3)-mean_preds(35,icol,lchnk,DP))/std_preds(35,icol,lchnk,DP))+ &
                                     Weights_V(36,icol,ipver,lchnk,DP)*((state1%v(icol,4)-mean_preds(36,icol,lchnk,DP))/std_preds(36,icol,lchnk,DP))+ &
                                     Weights_V(37,icol,ipver,lchnk,DP)*((state1%v(icol,5)-mean_preds(37,icol,lchnk,DP))/std_preds(37,icol,lchnk,DP))+ &
                                     Weights_V(38,icol,ipver,lchnk,DP)*((state1%v(icol,6)-mean_preds(38,icol,lchnk,DP))/std_preds(38,icol,lchnk,DP))+ &
                                     Weights_V(39,icol,ipver,lchnk,DP)*((state1%v(icol,7)-mean_preds(39,icol,lchnk,DP))/std_preds(39,icol,lchnk,DP))+ &
                                     Weights_V(40,icol,ipver,lchnk,DP)*((state1%v(icol,8)-mean_preds(40,icol,lchnk,DP))/std_preds(40,icol,lchnk,DP))+ &
                                     Weights_V(41,icol,ipver,lchnk,DP)*((state1%v(icol,9)-mean_preds(41,icol,lchnk,DP))/std_preds(41,icol,lchnk,DP))+ &
                                     Weights_V(42,icol,ipver,lchnk,DP)*((state1%v(icol,10)-mean_preds(42,icol,lchnk,DP))/std_preds(42,icol,lchnk,DP))+ &
                                     Weights_V(43,icol,ipver,lchnk,DP)*((state1%v(icol,11)-mean_preds(43,icol,lchnk,DP))/std_preds(43,icol,lchnk,DP))+ &
                                     Weights_V(44,icol,ipver,lchnk,DP)*((state1%v(icol,12)-mean_preds(44,icol,lchnk,DP))/std_preds(44,icol,lchnk,DP))+ &
                                     Weights_V(45,icol,ipver,lchnk,DP)*((state1%v(icol,13)-mean_preds(45,icol,lchnk,DP))/std_preds(45,icol,lchnk,DP))+ &
                                     Weights_V(46,icol,ipver,lchnk,DP)*((state1%v(icol,14)-mean_preds(46,icol,lchnk,DP))/std_preds(46,icol,lchnk,DP))+ &
                                     Weights_V(47,icol,ipver,lchnk,DP)*((state1%v(icol,15)-mean_preds(47,icol,lchnk,DP))/std_preds(47,icol,lchnk,DP))+ &
                                     Weights_V(48,icol,ipver,lchnk,DP)*((state1%v(icol,16)-mean_preds(48,icol,lchnk,DP))/std_preds(48,icol,lchnk,DP))+ &
                                     Weights_V(49,icol,ipver,lchnk,DP)*((state1%v(icol,17)-mean_preds(49,icol,lchnk,DP))/std_preds(49,icol,lchnk,DP))+ &
                                     Weights_V(50,icol,ipver,lchnk,DP)*((state1%v(icol,18)-mean_preds(50,icol,lchnk,DP))/std_preds(50,icol,lchnk,DP))+ &
                                     Weights_V(51,icol,ipver,lchnk,DP)*((state1%v(icol,19)-mean_preds(51,icol,lchnk,DP))/std_preds(51,icol,lchnk,DP))+ &
                                     Weights_V(52,icol,ipver,lchnk,DP)*((state1%v(icol,20)-mean_preds(52,icol,lchnk,DP))/std_preds(52,icol,lchnk,DP))+ &
                                     Weights_V(53,icol,ipver,lchnk,DP)*((state1%v(icol,21)-mean_preds(53,icol,lchnk,DP))/std_preds(53,icol,lchnk,DP))+ &
                                     Weights_V(54,icol,ipver,lchnk,DP)*((state1%v(icol,22)-mean_preds(54,icol,lchnk,DP))/std_preds(54,icol,lchnk,DP))+ &
                                     Weights_V(55,icol,ipver,lchnk,DP)*((state1%v(icol,23)-mean_preds(55,icol,lchnk,DP))/std_preds(55,icol,lchnk,DP))+ &
                                     Weights_V(56,icol,ipver,lchnk,DP)*((state1%v(icol,24)-mean_preds(56,icol,lchnk,DP))/std_preds(56,icol,lchnk,DP))+ &
                                     Weights_V(57,icol,ipver,lchnk,DP)*((state1%v(icol,25)-mean_preds(57,icol,lchnk,DP))/std_preds(57,icol,lchnk,DP))+ &
                                     Weights_V(58,icol,ipver,lchnk,DP)*((state1%v(icol,26)-mean_preds(58,icol,lchnk,DP))/std_preds(58,icol,lchnk,DP))+ &
                                     Weights_V(59,icol,ipver,lchnk,DP)*((state1%v(icol,27)-mean_preds(59,icol,lchnk,DP))/std_preds(59,icol,lchnk,DP))+ &
                                     Weights_V(60,icol,ipver,lchnk,DP)*((state1%v(icol,28)-mean_preds(60,icol,lchnk,DP))/std_preds(60,icol,lchnk,DP))+ &
                                     Weights_V(61,icol,ipver,lchnk,DP)*((state1%v(icol,29)-mean_preds(61,icol,lchnk,DP))/std_preds(61,icol,lchnk,DP))+ &
                                     Weights_V(62,icol,ipver,lchnk,DP)*((state1%v(icol,30)-mean_preds(62,icol,lchnk,DP))/std_preds(62,icol,lchnk,DP))+ &
                                     Weights_V(63,icol,ipver,lchnk,DP)*((state1%v(icol,31)-mean_preds(63,icol,lchnk,DP))/std_preds(63,icol,lchnk,DP))+ &
                                     Weights_V(64,icol,ipver,lchnk,DP)*((state1%v(icol,32)-mean_preds(64,icol,lchnk,DP))/std_preds(64,icol,lchnk,DP))
                                                                  
                                
      tend_tr_v = Weights_V(65,icol,ipver,lchnk,DP)*((state1%t(icol,1)-mean_preds(65,icol,lchnk,DP))/std_preds(65,icol,lchnk,DP))+ &
                                     Weights_V(66,icol,ipver,lchnk,DP)*((state1%t(icol,2)-mean_preds(66,icol,lchnk,DP))/std_preds(66,icol,lchnk,DP))+ &
                                     Weights_V(67,icol,ipver,lchnk,DP)*((state1%t(icol,3)-mean_preds(67,icol,lchnk,DP))/std_preds(67,icol,lchnk,DP))+ &
                                     Weights_V(68,icol,ipver,lchnk,DP)*((state1%t(icol,4)-mean_preds(68,icol,lchnk,DP))/std_preds(68,icol,lchnk,DP))+ &
                                     Weights_V(69,icol,ipver,lchnk,DP)*((state1%t(icol,5)-mean_preds(69,icol,lchnk,DP))/std_preds(69,icol,lchnk,DP))+ &
                                     Weights_V(70,icol,ipver,lchnk,DP)*((state1%t(icol,6)-mean_preds(70,icol,lchnk,DP))/std_preds(70,icol,lchnk,DP))+ &
                                     Weights_V(71,icol,ipver,lchnk,DP)*((state1%t(icol,7)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_V(71,icol,ipver,lchnk,DP)*((state1%t(icol,8)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_V(73,icol,ipver,lchnk,DP)*((state1%t(icol,9)-mean_preds(73,icol,lchnk,DP))/std_preds(73,icol,lchnk,DP))+ &
                                     Weights_V(74,icol,ipver,lchnk,DP)*((state1%t(icol,10)-mean_preds(74,icol,lchnk,DP))/std_preds(74,icol,lchnk,DP))+ &
                                     Weights_V(75,icol,ipver,lchnk,DP)*((state1%t(icol,11)-mean_preds(75,icol,lchnk,DP))/std_preds(75,icol,lchnk,DP))+ &
                                     Weights_V(76,icol,ipver,lchnk,DP)*((state1%t(icol,12)-mean_preds(76,icol,lchnk,DP))/std_preds(76,icol,lchnk,DP))+ &
                                     Weights_V(77,icol,ipver,lchnk,DP)*((state1%t(icol,13)-mean_preds(77,icol,lchnk,DP))/std_preds(77,icol,lchnk,DP))+ &
                                     Weights_V(78,icol,ipver,lchnk,DP)*((state1%t(icol,14)-mean_preds(78,icol,lchnk,DP))/std_preds(78,icol,lchnk,DP))+ &
                                     Weights_V(79,icol,ipver,lchnk,DP)*((state1%t(icol,15)-mean_preds(79,icol,lchnk,DP))/std_preds(79,icol,lchnk,DP))+ &
                                     Weights_V(80,icol,ipver,lchnk,DP)*((state1%t(icol,16)-mean_preds(80,icol,lchnk,DP))/std_preds(80,icol,lchnk,DP))+ &
                                     Weights_V(81,icol,ipver,lchnk,DP)*((state1%t(icol,17)-mean_preds(81,icol,lchnk,DP))/std_preds(81,icol,lchnk,DP))+ &
                                     Weights_V(82,icol,ipver,lchnk,DP)*((state1%t(icol,18)-mean_preds(82,icol,lchnk,DP))/std_preds(82,icol,lchnk,DP))+ &
                                     Weights_V(83,icol,ipver,lchnk,DP)*((state1%t(icol,19)-mean_preds(83,icol,lchnk,DP))/std_preds(83,icol,lchnk,DP))+ &
                                     Weights_V(84,icol,ipver,lchnk,DP)*((state1%t(icol,20)-mean_preds(84,icol,lchnk,DP))/std_preds(84,icol,lchnk,DP))+ &
                                     Weights_V(85,icol,ipver,lchnk,DP)*((state1%t(icol,21)-mean_preds(85,icol,lchnk,DP))/std_preds(85,icol,lchnk,DP))+ &
                                     Weights_V(86,icol,ipver,lchnk,DP)*((state1%t(icol,22)-mean_preds(86,icol,lchnk,DP))/std_preds(86,icol,lchnk,DP))+ &
                                     Weights_V(87,icol,ipver,lchnk,DP)*((state1%t(icol,23)-mean_preds(87,icol,lchnk,DP))/std_preds(87,icol,lchnk,DP))+ &
                                     Weights_V(88,icol,ipver,lchnk,DP)*((state1%t(icol,24)-mean_preds(88,icol,lchnk,DP))/std_preds(88,icol,lchnk,DP))+ &
                                     Weights_V(89,icol,ipver,lchnk,DP)*((state1%t(icol,25)-mean_preds(89,icol,lchnk,DP))/std_preds(89,icol,lchnk,DP))+ &
                                     Weights_V(90,icol,ipver,lchnk,DP)*((state1%t(icol,26)-mean_preds(90,icol,lchnk,DP))/std_preds(90,icol,lchnk,DP))+ &
                                     Weights_V(91,icol,ipver,lchnk,DP)*((state1%t(icol,27)-mean_preds(91,icol,lchnk,DP))/std_preds(91,icol,lchnk,DP))+ &
                                     Weights_V(92,icol,ipver,lchnk,DP)*((state1%t(icol,28)-mean_preds(92,icol,lchnk,DP))/std_preds(92,icol,lchnk,DP))+ &
                                     Weights_V(93,icol,ipver,lchnk,DP)*((state1%t(icol,29)-mean_preds(93,icol,lchnk,DP))/std_preds(93,icol,lchnk,DP))+ &
                                     Weights_V(94,icol,ipver,lchnk,DP)*((state1%t(icol,30)-mean_preds(94,icol,lchnk,DP))/std_preds(94,icol,lchnk,DP))+ &
                                     Weights_V(95,icol,ipver,lchnk,DP)*((state1%t(icol,31)-mean_preds(95,icol,lchnk,DP))/std_preds(95,icol,lchnk,DP))+ &
                                     Weights_V(96,icol,ipver,lchnk,DP)*((state1%t(icol,32)-mean_preds(96,icol,lchnk,DP))/std_preds(96,icol,lchnk,DP))
                                     
      tend_qr_v = Weights_V(97,icol,ipver,lchnk,DP)*((state1%q(icol,1,1)-mean_preds(97,icol,lchnk,DP))/std_preds(97,icol,lchnk,DP))+ &
                                     Weights_V(98,icol,ipver,lchnk,DP)*((state1%q(icol,2,1)-mean_preds(98,icol,lchnk,DP))/std_preds(98,icol,lchnk,DP))+ &
                                     Weights_V(99,icol,ipver,lchnk,DP)*((state1%q(icol,3,1)-mean_preds(99,icol,lchnk,DP))/std_preds(99,icol,lchnk,DP))+ &
                                     Weights_V(100,icol,ipver,lchnk,DP)*((state1%q(icol,4,1)-mean_preds(100,icol,lchnk,DP))/std_preds(100,icol,lchnk,DP))+ &
                                     Weights_V(101,icol,ipver,lchnk,DP)*((state1%q(icol,5,1)-mean_preds(101,icol,lchnk,DP))/std_preds(101,icol,lchnk,DP))+ &
                                     Weights_V(102,icol,ipver,lchnk,DP)*((state1%q(icol,6,1)-mean_preds(102,icol,lchnk,DP))/std_preds(102,icol,lchnk,DP))+ &
                                     Weights_V(103,icol,ipver,lchnk,DP)*((state1%q(icol,7,1)-mean_preds(103,icol,lchnk,DP))/std_preds(103,icol,lchnk,DP))+ &
                                     Weights_V(104,icol,ipver,lchnk,DP)*((state1%q(icol,8,1)-mean_preds(104,icol,lchnk,DP))/std_preds(104,icol,lchnk,DP))+ &
                                     Weights_V(105,icol,ipver,lchnk,DP)*((state1%q(icol,9,1)-mean_preds(105,icol,lchnk,DP))/std_preds(105,icol,lchnk,DP))+ &
                                     Weights_V(106,icol,ipver,lchnk,DP)*((state1%q(icol,10,1)-mean_preds(106,icol,lchnk,DP))/std_preds(106,icol,lchnk,DP))+ &
                                     Weights_V(107,icol,ipver,lchnk,DP)*((state1%q(icol,11,1)-mean_preds(107,icol,lchnk,DP))/std_preds(107,icol,lchnk,DP))+ &
                                     Weights_V(108,icol,ipver,lchnk,DP)*((state1%q(icol,12,1)-mean_preds(108,icol,lchnk,DP))/std_preds(108,icol,lchnk,DP))+ &
                                     Weights_V(109,icol,ipver,lchnk,DP)*((state1%q(icol,13,1)-mean_preds(109,icol,lchnk,DP))/std_preds(109,icol,lchnk,DP))+ &
                                     Weights_V(110,icol,ipver,lchnk,DP)*((state1%q(icol,14,1)-mean_preds(110,icol,lchnk,DP))/std_preds(110,icol,lchnk,DP))+ &
                                     Weights_V(111,icol,ipver,lchnk,DP)*((state1%q(icol,15,1)-mean_preds(111,icol,lchnk,DP))/std_preds(111,icol,lchnk,DP))+ &
                                     Weights_V(112,icol,ipver,lchnk,DP)*((state1%q(icol,16,1)-mean_preds(112,icol,lchnk,DP))/std_preds(112,icol,lchnk,DP))+ &
                                     Weights_V(113,icol,ipver,lchnk,DP)*((state1%q(icol,17,1)-mean_preds(113,icol,lchnk,DP))/std_preds(113,icol,lchnk,DP))+ &
                                     Weights_V(114,icol,ipver,lchnk,DP)*((state1%q(icol,18,1)-mean_preds(114,icol,lchnk,DP))/std_preds(114,icol,lchnk,DP))+ &
                                     Weights_V(115,icol,ipver,lchnk,DP)*((state1%q(icol,19,1)-mean_preds(115,icol,lchnk,DP))/std_preds(115,icol,lchnk,DP))+ &
                                     Weights_V(116,icol,ipver,lchnk,DP)*((state1%q(icol,20,1)-mean_preds(116,icol,lchnk,DP))/std_preds(116,icol,lchnk,DP))+ &
                                     Weights_V(117,icol,ipver,lchnk,DP)*((state1%q(icol,21,1)-mean_preds(117,icol,lchnk,DP))/std_preds(117,icol,lchnk,DP))+ &
                                     Weights_V(118,icol,ipver,lchnk,DP)*((state1%q(icol,22,1)-mean_preds(118,icol,lchnk,DP))/std_preds(118,icol,lchnk,DP))+ &
                                     Weights_V(119,icol,ipver,lchnk,DP)*((state1%q(icol,23,1)-mean_preds(119,icol,lchnk,DP))/std_preds(119,icol,lchnk,DP))+ &
                                     Weights_V(120,icol,ipver,lchnk,DP)*((state1%q(icol,24,1)-mean_preds(120,icol,lchnk,DP))/std_preds(120,icol,lchnk,DP))+ &
                                     Weights_V(121,icol,ipver,lchnk,DP)*((state1%q(icol,25,1)-mean_preds(121,icol,lchnk,DP))/std_preds(121,icol,lchnk,DP))+ &
                                     Weights_V(122,icol,ipver,lchnk,DP)*((state1%q(icol,26,1)-mean_preds(122,icol,lchnk,DP))/std_preds(122,icol,lchnk,DP))+ &
                                     Weights_V(123,icol,ipver,lchnk,DP)*((state1%q(icol,27,1)-mean_preds(123,icol,lchnk,DP))/std_preds(123,icol,lchnk,DP))+ &
                                     Weights_V(124,icol,ipver,lchnk,DP)*((state1%q(icol,28,1)-mean_preds(124,icol,lchnk,DP))/std_preds(124,icol,lchnk,DP))+ &
                                     Weights_V(125,icol,ipver,lchnk,DP)*((state1%q(icol,29,1)-mean_preds(125,icol,lchnk,DP))/std_preds(125,icol,lchnk,DP))+ &
                                     Weights_V(126,icol,ipver,lchnk,DP)*((state1%q(icol,30,1)-mean_preds(126,icol,lchnk,DP))/std_preds(126,icol,lchnk,DP))+ &
                                     Weights_V(127,icol,ipver,lchnk,DP)*((state1%q(icol,31,1)-mean_preds(127,icol,lchnk,DP))/std_preds(127,icol,lchnk,DP))+ &
                                     Weights_V(128,icol,ipver,lchnk,DP)*((state1%q(icol,32,1)-mean_preds(128,icol,lchnk,DP))/std_preds(128,icol,lchnk,DP))
                                     
      tend_or_v = Weights_V(129,icol,ipver,lchnk,DP)*((state1%omega(icol,1)-mean_preds(129,icol,lchnk,DP))/std_preds(129,icol,lchnk,DP))+ &
                                     Weights_V(130,icol,ipver,lchnk,DP)*((state1%omega(icol,2)-mean_preds(130,icol,lchnk,DP))/std_preds(130,icol,lchnk,DP))+ &
                                     Weights_V(131,icol,ipver,lchnk,DP)*((state1%omega(icol,3)-mean_preds(131,icol,lchnk,DP))/std_preds(131,icol,lchnk,DP))+ &
                                     Weights_V(132,icol,ipver,lchnk,DP)*((state1%omega(icol,4)-mean_preds(132,icol,lchnk,DP))/std_preds(132,icol,lchnk,DP))+ &
                                     Weights_V(133,icol,ipver,lchnk,DP)*((state1%omega(icol,5)-mean_preds(133,icol,lchnk,DP))/std_preds(133,icol,lchnk,DP))+ &
                                     Weights_V(134,icol,ipver,lchnk,DP)*((state1%omega(icol,6)-mean_preds(134,icol,lchnk,DP))/std_preds(134,icol,lchnk,DP))+ &
                                     Weights_V(135,icol,ipver,lchnk,DP)*((state1%omega(icol,7)-mean_preds(135,icol,lchnk,DP))/std_preds(135,icol,lchnk,DP))+ &
                                     Weights_V(136,icol,ipver,lchnk,DP)*((state1%omega(icol,8)-mean_preds(136,icol,lchnk,DP))/std_preds(136,icol,lchnk,DP))+ &
                                     Weights_V(137,icol,ipver,lchnk,DP)*((state1%omega(icol,9)-mean_preds(137,icol,lchnk,DP))/std_preds(137,icol,lchnk,DP))+ &
                                     Weights_V(138,icol,ipver,lchnk,DP)*((state1%omega(icol,10)-mean_preds(138,icol,lchnk,DP))/std_preds(138,icol,lchnk,DP))+ &
                                     Weights_V(139,icol,ipver,lchnk,DP)*((state1%omega(icol,11)-mean_preds(139,icol,lchnk,DP))/std_preds(139,icol,lchnk,DP))+ &
                                     Weights_V(140,icol,ipver,lchnk,DP)*((state1%omega(icol,12)-mean_preds(140,icol,lchnk,DP))/std_preds(140,icol,lchnk,DP))+ &
                                     Weights_V(141,icol,ipver,lchnk,DP)*((state1%omega(icol,13)-mean_preds(141,icol,lchnk,DP))/std_preds(141,icol,lchnk,DP))+ &
                                     Weights_V(142,icol,ipver,lchnk,DP)*((state1%omega(icol,14)-mean_preds(142,icol,lchnk,DP))/std_preds(142,icol,lchnk,DP))+ &
                                     Weights_V(143,icol,ipver,lchnk,DP)*((state1%omega(icol,15)-mean_preds(143,icol,lchnk,DP))/std_preds(143,icol,lchnk,DP))+ &
                                     Weights_V(144,icol,ipver,lchnk,DP)*((state1%omega(icol,16)-mean_preds(144,icol,lchnk,DP))/std_preds(144,icol,lchnk,DP))+ &
                                     Weights_V(145,icol,ipver,lchnk,DP)*((state1%omega(icol,17)-mean_preds(145,icol,lchnk,DP))/std_preds(145,icol,lchnk,DP))+ &
                                     Weights_V(146,icol,ipver,lchnk,DP)*((state1%omega(icol,18)-mean_preds(146,icol,lchnk,DP))/std_preds(146,icol,lchnk,DP))+ &
                                     Weights_V(147,icol,ipver,lchnk,DP)*((state1%omega(icol,19)-mean_preds(147,icol,lchnk,DP))/std_preds(147,icol,lchnk,DP))+ &
                                     Weights_V(148,icol,ipver,lchnk,DP)*((state1%omega(icol,20)-mean_preds(148,icol,lchnk,DP))/std_preds(148,icol,lchnk,DP))+ &
                                     Weights_V(149,icol,ipver,lchnk,DP)*((state1%omega(icol,21)-mean_preds(149,icol,lchnk,DP))/std_preds(149,icol,lchnk,DP))+ &
                                     Weights_V(150,icol,ipver,lchnk,DP)*((state1%omega(icol,22)-mean_preds(150,icol,lchnk,DP))/std_preds(150,icol,lchnk,DP))+ &
                                     Weights_V(151,icol,ipver,lchnk,DP)*((state1%omega(icol,23)-mean_preds(151,icol,lchnk,DP))/std_preds(151,icol,lchnk,DP))+ &
                                     Weights_V(152,icol,ipver,lchnk,DP)*((state1%omega(icol,24)-mean_preds(152,icol,lchnk,DP))/std_preds(152,icol,lchnk,DP))+ &
                                     Weights_V(153,icol,ipver,lchnk,DP)*((state1%omega(icol,25)-mean_preds(153,icol,lchnk,DP))/std_preds(153,icol,lchnk,DP))+ &
                                     Weights_V(154,icol,ipver,lchnk,DP)*((state1%omega(icol,26)-mean_preds(154,icol,lchnk,DP))/std_preds(154,icol,lchnk,DP))+ &
                                     Weights_V(155,icol,ipver,lchnk,DP)*((state1%omega(icol,27)-mean_preds(155,icol,lchnk,DP))/std_preds(155,icol,lchnk,DP))+ &
                                     Weights_V(156,icol,ipver,lchnk,DP)*((state1%omega(icol,28)-mean_preds(156,icol,lchnk,DP))/std_preds(156,icol,lchnk,DP))+ &
                                     Weights_V(157,icol,ipver,lchnk,DP)*((state1%omega(icol,29)-mean_preds(157,icol,lchnk,DP))/std_preds(157,icol,lchnk,DP))+ &
                                     Weights_V(158,icol,ipver,lchnk,DP)*((state1%omega(icol,30)-mean_preds(158,icol,lchnk,DP))/std_preds(158,icol,lchnk,DP))+ &
                                     Weights_V(159,icol,ipver,lchnk,DP)*((state1%omega(icol,31)-mean_preds(159,icol,lchnk,DP))/std_preds(159,icol,lchnk,DP))+ &
                                     Weights_V(160,icol,ipver,lchnk,DP)*((state1%omega(icol,32)-mean_preds(160,icol,lchnk,DP))/std_preds(160,icol,lchnk,DP))
                                     
                                     
      tend_auxr_v = Weights_V(161,icol,ipver,lchnk,DP)*((upwp(icol,24)-mean_preds(161,icol,lchnk,DP))/std_preds(161,icol,lchnk,DP))+ &
                                     Weights_V(162,icol,ipver,lchnk,DP)*((upwp(icol,25)-mean_preds(162,icol,lchnk,DP))/std_preds(162,icol,lchnk,DP))+ &
                                     Weights_V(163,icol,ipver,lchnk,DP)*((upwp(icol,26)-mean_preds(163,icol,lchnk,DP))/std_preds(163,icol,lchnk,DP))+ &
                                     Weights_V(164,icol,ipver,lchnk,DP)*((upwp(icol,27)-mean_preds(164,icol,lchnk,DP))/std_preds(164,icol,lchnk,DP))+ &
                                     Weights_V(165,icol,ipver,lchnk,DP)*((upwp(icol,28)-mean_preds(165,icol,lchnk,DP))/std_preds(165,icol,lchnk,DP))+ &
                                     Weights_V(166,icol,ipver,lchnk,DP)*((upwp(icol,29)-mean_preds(166,icol,lchnk,DP))/std_preds(166,icol,lchnk,DP))+ &
                                     Weights_V(167,icol,ipver,lchnk,DP)*((upwp(icol,30)-mean_preds(167,icol,lchnk,DP))/std_preds(167,icol,lchnk,DP))+ &
                                     Weights_V(168,icol,ipver,lchnk,DP)*((upwp(icol,31)-mean_preds(168,icol,lchnk,DP))/std_preds(168,icol,lchnk,DP))+ &
                                     Weights_V(169,icol,ipver,lchnk,DP)*((upwp(icol,32)-mean_preds(169,icol,lchnk,DP))/std_preds(169,icol,lchnk,DP))+ &
                                     Weights_V(170,icol,ipver,lchnk,DP)*((upwp(icol,33)-mean_preds(170,icol,lchnk,DP))/std_preds(170,icol,lchnk,DP))+ &
                                     Weights_V(171,icol,ipver,lchnk,DP)*((cam_in%wsx(icol)-mean_preds(171,icol,lchnk,DP))/std_preds(171,icol,lchnk,DP))+ &
                                     Weights_V(172,icol,ipver,lchnk,DP)*((cam_in%wsy(icol)-mean_preds(172,icol,lchnk,DP))/std_preds(172,icol,lchnk,DP))+ &
                                     Weights_V(173,icol,ipver,lchnk,DP)*((flntc(icol)-mean_preds(173,icol,lchnk,DP))/std_preds(173,icol,lchnk,DP))+ &
                                     Weights_V(174,icol,ipver,lchnk,DP)*((fsntoa(icol)-mean_preds(174,icol,lchnk,DP))/std_preds(174,icol,lchnk,DP))+ &
                                     Weights_V(175,icol,ipver,lchnk,DP)*((pblh(icol)-mean_preds(175,icol,lchnk,DP))/std_preds(175,icol,lchnk,DP))+ &
                                     Weights_V(176,icol,ipver,lchnk,DP)*((cam_out%tbot(icol)-mean_preds(176,icol,lchnk,DP))/std_preds(176,icol,lchnk,DP))+ &
                                     Weights_V(177,icol,ipver,lchnk,DP)*((ustarwec(icol)-mean_preds(177,icol,lchnk,DP))/std_preds(177,icol,lchnk,DP))+ &
                                     Weights_V(178,icol,ipver,lchnk,DP)*1.5_r8*DOY2+ &
                                     Weights_V(179,icol,ipver,lchnk,DP)*1.5_r8*DOY1+ &
                                     Weights_V(180,icol,ipver,lchnk,DP)*1.5_r8*HOD2+ &
                                     Weights_V(181,icol,ipver,lchnk,DP)*1.5_r8*HOD1+ &
                                     Weights_V(182,icol,ipver,lchnk,DP)
                       
      
      if (abs(tend_ur_v + tend_vr_v + tend_tr_v + tend_qr_v + tend_or_v + tend_auxr_v) > threshy) then
          DAMLin_Vstep(icol,ipver,lchnk) = sign(threshy, tend_ur_v + tend_vr_v + tend_tr_v + tend_qr_v + tend_or_v + tend_auxr_v)
      else 
          DAMLin_Vstep(icol,ipver,lchnk) = tend_ur_v + tend_vr_v + tend_tr_v + tend_qr_v + tend_or_v + tend_auxr_v
      end if  
      
      tend_ur_t = Weights_T(1,icol,ipver,lchnk,DP)*((state1%u(icol,1)-mean_preds(1,icol,lchnk,DP))/std_preds(1,icol,lchnk,DP))+ &
                                     Weights_T(2,icol,ipver,lchnk,DP)*((state1%u(icol,2)-mean_preds(2,icol,lchnk,DP))/std_preds(2,icol,lchnk,DP))+ &
                                     Weights_T(3,icol,ipver,lchnk,DP)*((state1%u(icol,3)-mean_preds(3,icol,lchnk,DP))/std_preds(3,icol,lchnk,DP))+ &
                                     Weights_T(4,icol,ipver,lchnk,DP)*((state1%u(icol,4)-mean_preds(4,icol,lchnk,DP))/std_preds(4,icol,lchnk,DP))+ &
                                     Weights_T(5,icol,ipver,lchnk,DP)*((state1%u(icol,5)-mean_preds(5,icol,lchnk,DP))/std_preds(5,icol,lchnk,DP))+ &
                                     Weights_T(6,icol,ipver,lchnk,DP)*((state1%u(icol,6)-mean_preds(6,icol,lchnk,DP))/std_preds(6,icol,lchnk,DP))+ &
                                     Weights_T(7,icol,ipver,lchnk,DP)*((state1%u(icol,7)-mean_preds(7,icol,lchnk,DP))/std_preds(7,icol,lchnk,DP))+ &
                                     Weights_T(8,icol,ipver,lchnk,DP)*((state1%u(icol,8)-mean_preds(8,icol,lchnk,DP))/std_preds(8,icol,lchnk,DP))+ &
                                     Weights_T(9,icol,ipver,lchnk,DP)*((state1%u(icol,9)-mean_preds(9,icol,lchnk,DP))/std_preds(9,icol,lchnk,DP))+ &
                                     Weights_T(10,icol,ipver,lchnk,DP)*((state1%u(icol,10)-mean_preds(10,icol,lchnk,DP))/std_preds(10,icol,lchnk,DP))+ &
                                     Weights_T(11,icol,ipver,lchnk,DP)*((state1%u(icol,11)-mean_preds(11,icol,lchnk,DP))/std_preds(11,icol,lchnk,DP))+ &
                                     Weights_T(12,icol,ipver,lchnk,DP)*((state1%u(icol,12)-mean_preds(12,icol,lchnk,DP))/std_preds(12,icol,lchnk,DP))+ &
                                     Weights_T(13,icol,ipver,lchnk,DP)*((state1%u(icol,13)-mean_preds(13,icol,lchnk,DP))/std_preds(13,icol,lchnk,DP))+ &
                                     Weights_T(14,icol,ipver,lchnk,DP)*((state1%u(icol,14)-mean_preds(14,icol,lchnk,DP))/std_preds(14,icol,lchnk,DP))+ &
                                     Weights_T(15,icol,ipver,lchnk,DP)*((state1%u(icol,15)-mean_preds(15,icol,lchnk,DP))/std_preds(15,icol,lchnk,DP))+ &
                                     Weights_T(16,icol,ipver,lchnk,DP)*((state1%u(icol,16)-mean_preds(16,icol,lchnk,DP))/std_preds(16,icol,lchnk,DP))+ &
                                     Weights_T(17,icol,ipver,lchnk,DP)*((state1%u(icol,17)-mean_preds(17,icol,lchnk,DP))/std_preds(17,icol,lchnk,DP))+ &
                                     Weights_T(18,icol,ipver,lchnk,DP)*((state1%u(icol,18)-mean_preds(18,icol,lchnk,DP))/std_preds(18,icol,lchnk,DP))+ &
                                     Weights_T(19,icol,ipver,lchnk,DP)*((state1%u(icol,19)-mean_preds(19,icol,lchnk,DP))/std_preds(19,icol,lchnk,DP))+ &
                                     Weights_T(20,icol,ipver,lchnk,DP)*((state1%u(icol,20)-mean_preds(20,icol,lchnk,DP))/std_preds(20,icol,lchnk,DP))+ &
                                     Weights_T(21,icol,ipver,lchnk,DP)*((state1%u(icol,21)-mean_preds(21,icol,lchnk,DP))/std_preds(21,icol,lchnk,DP))+ &
                                     Weights_T(22,icol,ipver,lchnk,DP)*((state1%u(icol,22)-mean_preds(22,icol,lchnk,DP))/std_preds(22,icol,lchnk,DP))+ &
                                     Weights_T(23,icol,ipver,lchnk,DP)*((state1%u(icol,23)-mean_preds(23,icol,lchnk,DP))/std_preds(23,icol,lchnk,DP))+ &
                                     Weights_T(24,icol,ipver,lchnk,DP)*((state1%u(icol,24)-mean_preds(24,icol,lchnk,DP))/std_preds(24,icol,lchnk,DP))+ &
                                     Weights_T(25,icol,ipver,lchnk,DP)*((state1%u(icol,25)-mean_preds(25,icol,lchnk,DP))/std_preds(25,icol,lchnk,DP))+ &
                                     Weights_T(26,icol,ipver,lchnk,DP)*((state1%u(icol,26)-mean_preds(26,icol,lchnk,DP))/std_preds(26,icol,lchnk,DP))+ &
                                     Weights_T(27,icol,ipver,lchnk,DP)*((state1%u(icol,27)-mean_preds(27,icol,lchnk,DP))/std_preds(27,icol,lchnk,DP))+ &
                                     Weights_T(28,icol,ipver,lchnk,DP)*((state1%u(icol,28)-mean_preds(28,icol,lchnk,DP))/std_preds(28,icol,lchnk,DP))+ &
                                     Weights_T(29,icol,ipver,lchnk,DP)*((state1%u(icol,29)-mean_preds(29,icol,lchnk,DP))/std_preds(29,icol,lchnk,DP))+ &
                                     Weights_T(30,icol,ipver,lchnk,DP)*((state1%u(icol,30)-mean_preds(30,icol,lchnk,DP))/std_preds(30,icol,lchnk,DP))+ &
                                     Weights_T(31,icol,ipver,lchnk,DP)*((state1%u(icol,31)-mean_preds(31,icol,lchnk,DP))/std_preds(31,icol,lchnk,DP))+ &
                                     Weights_T(32,icol,ipver,lchnk,DP)*((state1%u(icol,32)-mean_preds(32,icol,lchnk,DP))/std_preds(32,icol,lchnk,DP)) 
                                     
                                     
     tend_vr_t = Weights_T(33,icol,ipver,lchnk,DP)*((state1%v(icol,1)-mean_preds(33,icol,lchnk,DP))/std_preds(33,icol,lchnk,DP))+ &
                                     Weights_T(34,icol,ipver,lchnk,DP)*((state1%v(icol,2)-mean_preds(34,icol,lchnk,DP))/std_preds(34,icol,lchnk,DP))+ &
                                     Weights_T(35,icol,ipver,lchnk,DP)*((state1%v(icol,3)-mean_preds(35,icol,lchnk,DP))/std_preds(35,icol,lchnk,DP))+ &
                                     Weights_T(36,icol,ipver,lchnk,DP)*((state1%v(icol,4)-mean_preds(36,icol,lchnk,DP))/std_preds(36,icol,lchnk,DP))+ &
                                     Weights_T(37,icol,ipver,lchnk,DP)*((state1%v(icol,5)-mean_preds(37,icol,lchnk,DP))/std_preds(37,icol,lchnk,DP))+ &
                                     Weights_T(38,icol,ipver,lchnk,DP)*((state1%v(icol,6)-mean_preds(38,icol,lchnk,DP))/std_preds(38,icol,lchnk,DP))+ &
                                     Weights_T(39,icol,ipver,lchnk,DP)*((state1%v(icol,7)-mean_preds(39,icol,lchnk,DP))/std_preds(39,icol,lchnk,DP))+ &
                                     Weights_T(40,icol,ipver,lchnk,DP)*((state1%v(icol,8)-mean_preds(40,icol,lchnk,DP))/std_preds(40,icol,lchnk,DP))+ &
                                     Weights_T(41,icol,ipver,lchnk,DP)*((state1%v(icol,9)-mean_preds(41,icol,lchnk,DP))/std_preds(41,icol,lchnk,DP))+ &
                                     Weights_T(42,icol,ipver,lchnk,DP)*((state1%v(icol,10)-mean_preds(42,icol,lchnk,DP))/std_preds(42,icol,lchnk,DP))+ &
                                     Weights_T(43,icol,ipver,lchnk,DP)*((state1%v(icol,11)-mean_preds(43,icol,lchnk,DP))/std_preds(43,icol,lchnk,DP))+ &
                                     Weights_T(44,icol,ipver,lchnk,DP)*((state1%v(icol,12)-mean_preds(44,icol,lchnk,DP))/std_preds(44,icol,lchnk,DP))+ &
                                     Weights_T(45,icol,ipver,lchnk,DP)*((state1%v(icol,13)-mean_preds(45,icol,lchnk,DP))/std_preds(45,icol,lchnk,DP))+ &
                                     Weights_T(46,icol,ipver,lchnk,DP)*((state1%v(icol,14)-mean_preds(46,icol,lchnk,DP))/std_preds(46,icol,lchnk,DP))+ &
                                     Weights_T(47,icol,ipver,lchnk,DP)*((state1%v(icol,15)-mean_preds(47,icol,lchnk,DP))/std_preds(47,icol,lchnk,DP))+ &
                                     Weights_T(48,icol,ipver,lchnk,DP)*((state1%v(icol,16)-mean_preds(48,icol,lchnk,DP))/std_preds(48,icol,lchnk,DP))+ &
                                     Weights_T(49,icol,ipver,lchnk,DP)*((state1%v(icol,17)-mean_preds(49,icol,lchnk,DP))/std_preds(49,icol,lchnk,DP))+ &
                                     Weights_T(50,icol,ipver,lchnk,DP)*((state1%v(icol,18)-mean_preds(50,icol,lchnk,DP))/std_preds(50,icol,lchnk,DP))+ &
                                     Weights_T(51,icol,ipver,lchnk,DP)*((state1%v(icol,19)-mean_preds(51,icol,lchnk,DP))/std_preds(51,icol,lchnk,DP))+ &
                                     Weights_T(52,icol,ipver,lchnk,DP)*((state1%v(icol,20)-mean_preds(52,icol,lchnk,DP))/std_preds(52,icol,lchnk,DP))+ &
                                     Weights_T(53,icol,ipver,lchnk,DP)*((state1%v(icol,21)-mean_preds(53,icol,lchnk,DP))/std_preds(53,icol,lchnk,DP))+ &
                                     Weights_T(54,icol,ipver,lchnk,DP)*((state1%v(icol,22)-mean_preds(54,icol,lchnk,DP))/std_preds(54,icol,lchnk,DP))+ &
                                     Weights_T(55,icol,ipver,lchnk,DP)*((state1%v(icol,23)-mean_preds(55,icol,lchnk,DP))/std_preds(55,icol,lchnk,DP))+ &
                                     Weights_T(56,icol,ipver,lchnk,DP)*((state1%v(icol,24)-mean_preds(56,icol,lchnk,DP))/std_preds(56,icol,lchnk,DP))+ &
                                     Weights_T(57,icol,ipver,lchnk,DP)*((state1%v(icol,25)-mean_preds(57,icol,lchnk,DP))/std_preds(57,icol,lchnk,DP))+ &
                                     Weights_T(58,icol,ipver,lchnk,DP)*((state1%v(icol,26)-mean_preds(58,icol,lchnk,DP))/std_preds(58,icol,lchnk,DP))+ &
                                     Weights_T(59,icol,ipver,lchnk,DP)*((state1%v(icol,27)-mean_preds(59,icol,lchnk,DP))/std_preds(59,icol,lchnk,DP))+ &
                                     Weights_T(60,icol,ipver,lchnk,DP)*((state1%v(icol,28)-mean_preds(60,icol,lchnk,DP))/std_preds(60,icol,lchnk,DP))+ &
                                     Weights_T(61,icol,ipver,lchnk,DP)*((state1%v(icol,29)-mean_preds(61,icol,lchnk,DP))/std_preds(61,icol,lchnk,DP))+ &
                                     Weights_T(62,icol,ipver,lchnk,DP)*((state1%v(icol,30)-mean_preds(62,icol,lchnk,DP))/std_preds(62,icol,lchnk,DP))+ &
                                     Weights_T(63,icol,ipver,lchnk,DP)*((state1%v(icol,31)-mean_preds(63,icol,lchnk,DP))/std_preds(63,icol,lchnk,DP))+ &
                                     Weights_T(64,icol,ipver,lchnk,DP)*((state1%v(icol,32)-mean_preds(64,icol,lchnk,DP))/std_preds(64,icol,lchnk,DP))
                                                                  
                                
      tend_tr_t = Weights_T(65,icol,ipver,lchnk,DP)*((state1%t(icol,1)-mean_preds(65,icol,lchnk,DP))/std_preds(65,icol,lchnk,DP))+ &
                                     Weights_T(66,icol,ipver,lchnk,DP)*((state1%t(icol,2)-mean_preds(66,icol,lchnk,DP))/std_preds(66,icol,lchnk,DP))+ &
                                     Weights_T(67,icol,ipver,lchnk,DP)*((state1%t(icol,3)-mean_preds(67,icol,lchnk,DP))/std_preds(67,icol,lchnk,DP))+ &
                                     Weights_T(68,icol,ipver,lchnk,DP)*((state1%t(icol,4)-mean_preds(68,icol,lchnk,DP))/std_preds(68,icol,lchnk,DP))+ &
                                     Weights_T(69,icol,ipver,lchnk,DP)*((state1%t(icol,5)-mean_preds(69,icol,lchnk,DP))/std_preds(69,icol,lchnk,DP))+ &
                                     Weights_T(70,icol,ipver,lchnk,DP)*((state1%t(icol,6)-mean_preds(70,icol,lchnk,DP))/std_preds(70,icol,lchnk,DP))+ &
                                     Weights_T(71,icol,ipver,lchnk,DP)*((state1%t(icol,7)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_T(71,icol,ipver,lchnk,DP)*((state1%t(icol,8)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_T(73,icol,ipver,lchnk,DP)*((state1%t(icol,9)-mean_preds(73,icol,lchnk,DP))/std_preds(73,icol,lchnk,DP))+ &
                                     Weights_T(74,icol,ipver,lchnk,DP)*((state1%t(icol,10)-mean_preds(74,icol,lchnk,DP))/std_preds(74,icol,lchnk,DP))+ &
                                     Weights_T(75,icol,ipver,lchnk,DP)*((state1%t(icol,11)-mean_preds(75,icol,lchnk,DP))/std_preds(75,icol,lchnk,DP))+ &
                                     Weights_T(76,icol,ipver,lchnk,DP)*((state1%t(icol,12)-mean_preds(76,icol,lchnk,DP))/std_preds(76,icol,lchnk,DP))+ &
                                     Weights_T(77,icol,ipver,lchnk,DP)*((state1%t(icol,13)-mean_preds(77,icol,lchnk,DP))/std_preds(77,icol,lchnk,DP))+ &
                                     Weights_T(78,icol,ipver,lchnk,DP)*((state1%t(icol,14)-mean_preds(78,icol,lchnk,DP))/std_preds(78,icol,lchnk,DP))+ &
                                     Weights_T(79,icol,ipver,lchnk,DP)*((state1%t(icol,15)-mean_preds(79,icol,lchnk,DP))/std_preds(79,icol,lchnk,DP))+ &
                                     Weights_T(80,icol,ipver,lchnk,DP)*((state1%t(icol,16)-mean_preds(80,icol,lchnk,DP))/std_preds(80,icol,lchnk,DP))+ &
                                     Weights_T(81,icol,ipver,lchnk,DP)*((state1%t(icol,17)-mean_preds(81,icol,lchnk,DP))/std_preds(81,icol,lchnk,DP))+ &
                                     Weights_T(82,icol,ipver,lchnk,DP)*((state1%t(icol,18)-mean_preds(82,icol,lchnk,DP))/std_preds(82,icol,lchnk,DP))+ &
                                     Weights_T(83,icol,ipver,lchnk,DP)*((state1%t(icol,19)-mean_preds(83,icol,lchnk,DP))/std_preds(83,icol,lchnk,DP))+ &
                                     Weights_T(84,icol,ipver,lchnk,DP)*((state1%t(icol,20)-mean_preds(84,icol,lchnk,DP))/std_preds(84,icol,lchnk,DP))+ &
                                     Weights_T(85,icol,ipver,lchnk,DP)*((state1%t(icol,21)-mean_preds(85,icol,lchnk,DP))/std_preds(85,icol,lchnk,DP))+ &
                                     Weights_T(86,icol,ipver,lchnk,DP)*((state1%t(icol,22)-mean_preds(86,icol,lchnk,DP))/std_preds(86,icol,lchnk,DP))+ &
                                     Weights_T(87,icol,ipver,lchnk,DP)*((state1%t(icol,23)-mean_preds(87,icol,lchnk,DP))/std_preds(87,icol,lchnk,DP))+ &
                                     Weights_T(88,icol,ipver,lchnk,DP)*((state1%t(icol,24)-mean_preds(88,icol,lchnk,DP))/std_preds(88,icol,lchnk,DP))+ &
                                     Weights_T(89,icol,ipver,lchnk,DP)*((state1%t(icol,25)-mean_preds(89,icol,lchnk,DP))/std_preds(89,icol,lchnk,DP))+ &
                                     Weights_T(90,icol,ipver,lchnk,DP)*((state1%t(icol,26)-mean_preds(90,icol,lchnk,DP))/std_preds(90,icol,lchnk,DP))+ &
                                     Weights_T(91,icol,ipver,lchnk,DP)*((state1%t(icol,27)-mean_preds(91,icol,lchnk,DP))/std_preds(91,icol,lchnk,DP))+ &
                                     Weights_T(92,icol,ipver,lchnk,DP)*((state1%t(icol,28)-mean_preds(92,icol,lchnk,DP))/std_preds(92,icol,lchnk,DP))+ &
                                     Weights_T(93,icol,ipver,lchnk,DP)*((state1%t(icol,29)-mean_preds(93,icol,lchnk,DP))/std_preds(93,icol,lchnk,DP))+ &
                                     Weights_T(94,icol,ipver,lchnk,DP)*((state1%t(icol,30)-mean_preds(94,icol,lchnk,DP))/std_preds(94,icol,lchnk,DP))+ &
                                     Weights_T(95,icol,ipver,lchnk,DP)*((state1%t(icol,31)-mean_preds(95,icol,lchnk,DP))/std_preds(95,icol,lchnk,DP))+ &
                                     Weights_T(96,icol,ipver,lchnk,DP)*((state1%t(icol,32)-mean_preds(96,icol,lchnk,DP))/std_preds(96,icol,lchnk,DP))
                                     
      tend_qr_t = Weights_T(97,icol,ipver,lchnk,DP)*((state1%q(icol,1,1)-mean_preds(97,icol,lchnk,DP))/std_preds(97,icol,lchnk,DP))+ &
                                     Weights_T(98,icol,ipver,lchnk,DP)*((state1%q(icol,2,1)-mean_preds(98,icol,lchnk,DP))/std_preds(98,icol,lchnk,DP))+ &
                                     Weights_T(99,icol,ipver,lchnk,DP)*((state1%q(icol,3,1)-mean_preds(99,icol,lchnk,DP))/std_preds(99,icol,lchnk,DP))+ &
                                     Weights_T(100,icol,ipver,lchnk,DP)*((state1%q(icol,4,1)-mean_preds(100,icol,lchnk,DP))/std_preds(100,icol,lchnk,DP))+ &
                                     Weights_T(101,icol,ipver,lchnk,DP)*((state1%q(icol,5,1)-mean_preds(101,icol,lchnk,DP))/std_preds(101,icol,lchnk,DP))+ &
                                     Weights_T(102,icol,ipver,lchnk,DP)*((state1%q(icol,6,1)-mean_preds(102,icol,lchnk,DP))/std_preds(102,icol,lchnk,DP))+ &
                                     Weights_T(103,icol,ipver,lchnk,DP)*((state1%q(icol,7,1)-mean_preds(103,icol,lchnk,DP))/std_preds(103,icol,lchnk,DP))+ &
                                     Weights_T(104,icol,ipver,lchnk,DP)*((state1%q(icol,8,1)-mean_preds(104,icol,lchnk,DP))/std_preds(104,icol,lchnk,DP))+ &
                                     Weights_T(105,icol,ipver,lchnk,DP)*((state1%q(icol,9,1)-mean_preds(105,icol,lchnk,DP))/std_preds(105,icol,lchnk,DP))+ &
                                     Weights_T(106,icol,ipver,lchnk,DP)*((state1%q(icol,10,1)-mean_preds(106,icol,lchnk,DP))/std_preds(106,icol,lchnk,DP))+ &
                                     Weights_T(107,icol,ipver,lchnk,DP)*((state1%q(icol,11,1)-mean_preds(107,icol,lchnk,DP))/std_preds(107,icol,lchnk,DP))+ &
                                     Weights_T(108,icol,ipver,lchnk,DP)*((state1%q(icol,12,1)-mean_preds(108,icol,lchnk,DP))/std_preds(108,icol,lchnk,DP))+ &
                                     Weights_T(109,icol,ipver,lchnk,DP)*((state1%q(icol,13,1)-mean_preds(109,icol,lchnk,DP))/std_preds(109,icol,lchnk,DP))+ &
                                     Weights_T(110,icol,ipver,lchnk,DP)*((state1%q(icol,14,1)-mean_preds(110,icol,lchnk,DP))/std_preds(110,icol,lchnk,DP))+ &
                                     Weights_T(111,icol,ipver,lchnk,DP)*((state1%q(icol,15,1)-mean_preds(111,icol,lchnk,DP))/std_preds(111,icol,lchnk,DP))+ &
                                     Weights_T(112,icol,ipver,lchnk,DP)*((state1%q(icol,16,1)-mean_preds(112,icol,lchnk,DP))/std_preds(112,icol,lchnk,DP))+ &
                                     Weights_T(113,icol,ipver,lchnk,DP)*((state1%q(icol,17,1)-mean_preds(113,icol,lchnk,DP))/std_preds(113,icol,lchnk,DP))+ &
                                     Weights_T(114,icol,ipver,lchnk,DP)*((state1%q(icol,18,1)-mean_preds(114,icol,lchnk,DP))/std_preds(114,icol,lchnk,DP))+ &
                                     Weights_T(115,icol,ipver,lchnk,DP)*((state1%q(icol,19,1)-mean_preds(115,icol,lchnk,DP))/std_preds(115,icol,lchnk,DP))+ &
                                     Weights_T(116,icol,ipver,lchnk,DP)*((state1%q(icol,20,1)-mean_preds(116,icol,lchnk,DP))/std_preds(116,icol,lchnk,DP))+ &
                                     Weights_T(117,icol,ipver,lchnk,DP)*((state1%q(icol,21,1)-mean_preds(117,icol,lchnk,DP))/std_preds(117,icol,lchnk,DP))+ &
                                     Weights_T(118,icol,ipver,lchnk,DP)*((state1%q(icol,22,1)-mean_preds(118,icol,lchnk,DP))/std_preds(118,icol,lchnk,DP))+ &
                                     Weights_T(119,icol,ipver,lchnk,DP)*((state1%q(icol,23,1)-mean_preds(119,icol,lchnk,DP))/std_preds(119,icol,lchnk,DP))+ &
                                     Weights_T(120,icol,ipver,lchnk,DP)*((state1%q(icol,24,1)-mean_preds(120,icol,lchnk,DP))/std_preds(120,icol,lchnk,DP))+ &
                                     Weights_T(121,icol,ipver,lchnk,DP)*((state1%q(icol,25,1)-mean_preds(121,icol,lchnk,DP))/std_preds(121,icol,lchnk,DP))+ &
                                     Weights_T(122,icol,ipver,lchnk,DP)*((state1%q(icol,26,1)-mean_preds(122,icol,lchnk,DP))/std_preds(122,icol,lchnk,DP))+ &
                                     Weights_T(123,icol,ipver,lchnk,DP)*((state1%q(icol,27,1)-mean_preds(123,icol,lchnk,DP))/std_preds(123,icol,lchnk,DP))+ &
                                     Weights_T(124,icol,ipver,lchnk,DP)*((state1%q(icol,28,1)-mean_preds(124,icol,lchnk,DP))/std_preds(124,icol,lchnk,DP))+ &
                                     Weights_T(125,icol,ipver,lchnk,DP)*((state1%q(icol,29,1)-mean_preds(125,icol,lchnk,DP))/std_preds(125,icol,lchnk,DP))+ &
                                     Weights_T(126,icol,ipver,lchnk,DP)*((state1%q(icol,30,1)-mean_preds(126,icol,lchnk,DP))/std_preds(126,icol,lchnk,DP))+ &
                                     Weights_T(127,icol,ipver,lchnk,DP)*((state1%q(icol,31,1)-mean_preds(127,icol,lchnk,DP))/std_preds(127,icol,lchnk,DP))+ &
                                     Weights_T(128,icol,ipver,lchnk,DP)*((state1%q(icol,32,1)-mean_preds(128,icol,lchnk,DP))/std_preds(128,icol,lchnk,DP))
                                     
      tend_or_t = Weights_T(129,icol,ipver,lchnk,DP)*((state1%omega(icol,1)-mean_preds(129,icol,lchnk,DP))/std_preds(129,icol,lchnk,DP))+ &
                                     Weights_T(130,icol,ipver,lchnk,DP)*((state1%omega(icol,2)-mean_preds(130,icol,lchnk,DP))/std_preds(130,icol,lchnk,DP))+ &
                                     Weights_T(131,icol,ipver,lchnk,DP)*((state1%omega(icol,3)-mean_preds(131,icol,lchnk,DP))/std_preds(131,icol,lchnk,DP))+ &
                                     Weights_T(132,icol,ipver,lchnk,DP)*((state1%omega(icol,4)-mean_preds(132,icol,lchnk,DP))/std_preds(132,icol,lchnk,DP))+ &
                                     Weights_T(133,icol,ipver,lchnk,DP)*((state1%omega(icol,5)-mean_preds(133,icol,lchnk,DP))/std_preds(133,icol,lchnk,DP))+ &
                                     Weights_T(134,icol,ipver,lchnk,DP)*((state1%omega(icol,6)-mean_preds(134,icol,lchnk,DP))/std_preds(134,icol,lchnk,DP))+ &
                                     Weights_T(135,icol,ipver,lchnk,DP)*((state1%omega(icol,7)-mean_preds(135,icol,lchnk,DP))/std_preds(135,icol,lchnk,DP))+ &
                                     Weights_T(136,icol,ipver,lchnk,DP)*((state1%omega(icol,8)-mean_preds(136,icol,lchnk,DP))/std_preds(136,icol,lchnk,DP))+ &
                                     Weights_T(137,icol,ipver,lchnk,DP)*((state1%omega(icol,9)-mean_preds(137,icol,lchnk,DP))/std_preds(137,icol,lchnk,DP))+ &
                                     Weights_T(138,icol,ipver,lchnk,DP)*((state1%omega(icol,10)-mean_preds(138,icol,lchnk,DP))/std_preds(138,icol,lchnk,DP))+ &
                                     Weights_T(139,icol,ipver,lchnk,DP)*((state1%omega(icol,11)-mean_preds(139,icol,lchnk,DP))/std_preds(139,icol,lchnk,DP))+ &
                                     Weights_T(140,icol,ipver,lchnk,DP)*((state1%omega(icol,12)-mean_preds(140,icol,lchnk,DP))/std_preds(140,icol,lchnk,DP))+ &
                                     Weights_T(141,icol,ipver,lchnk,DP)*((state1%omega(icol,13)-mean_preds(141,icol,lchnk,DP))/std_preds(141,icol,lchnk,DP))+ &
                                     Weights_T(142,icol,ipver,lchnk,DP)*((state1%omega(icol,14)-mean_preds(142,icol,lchnk,DP))/std_preds(142,icol,lchnk,DP))+ &
                                     Weights_T(143,icol,ipver,lchnk,DP)*((state1%omega(icol,15)-mean_preds(143,icol,lchnk,DP))/std_preds(143,icol,lchnk,DP))+ &
                                     Weights_T(144,icol,ipver,lchnk,DP)*((state1%omega(icol,16)-mean_preds(144,icol,lchnk,DP))/std_preds(144,icol,lchnk,DP))+ &
                                     Weights_T(145,icol,ipver,lchnk,DP)*((state1%omega(icol,17)-mean_preds(145,icol,lchnk,DP))/std_preds(145,icol,lchnk,DP))+ &
                                     Weights_T(146,icol,ipver,lchnk,DP)*((state1%omega(icol,18)-mean_preds(146,icol,lchnk,DP))/std_preds(146,icol,lchnk,DP))+ &
                                     Weights_T(147,icol,ipver,lchnk,DP)*((state1%omega(icol,19)-mean_preds(147,icol,lchnk,DP))/std_preds(147,icol,lchnk,DP))+ &
                                     Weights_T(148,icol,ipver,lchnk,DP)*((state1%omega(icol,20)-mean_preds(148,icol,lchnk,DP))/std_preds(148,icol,lchnk,DP))+ &
                                     Weights_T(149,icol,ipver,lchnk,DP)*((state1%omega(icol,21)-mean_preds(149,icol,lchnk,DP))/std_preds(149,icol,lchnk,DP))+ &
                                     Weights_T(150,icol,ipver,lchnk,DP)*((state1%omega(icol,22)-mean_preds(150,icol,lchnk,DP))/std_preds(150,icol,lchnk,DP))+ &
                                     Weights_T(151,icol,ipver,lchnk,DP)*((state1%omega(icol,23)-mean_preds(151,icol,lchnk,DP))/std_preds(151,icol,lchnk,DP))+ &
                                     Weights_T(152,icol,ipver,lchnk,DP)*((state1%omega(icol,24)-mean_preds(152,icol,lchnk,DP))/std_preds(152,icol,lchnk,DP))+ &
                                     Weights_T(153,icol,ipver,lchnk,DP)*((state1%omega(icol,25)-mean_preds(153,icol,lchnk,DP))/std_preds(153,icol,lchnk,DP))+ &
                                     Weights_T(154,icol,ipver,lchnk,DP)*((state1%omega(icol,26)-mean_preds(154,icol,lchnk,DP))/std_preds(154,icol,lchnk,DP))+ &
                                     Weights_T(155,icol,ipver,lchnk,DP)*((state1%omega(icol,27)-mean_preds(155,icol,lchnk,DP))/std_preds(155,icol,lchnk,DP))+ &
                                     Weights_T(156,icol,ipver,lchnk,DP)*((state1%omega(icol,28)-mean_preds(156,icol,lchnk,DP))/std_preds(156,icol,lchnk,DP))+ &
                                     Weights_T(157,icol,ipver,lchnk,DP)*((state1%omega(icol,29)-mean_preds(157,icol,lchnk,DP))/std_preds(157,icol,lchnk,DP))+ &
                                     Weights_T(158,icol,ipver,lchnk,DP)*((state1%omega(icol,30)-mean_preds(158,icol,lchnk,DP))/std_preds(158,icol,lchnk,DP))+ &
                                     Weights_T(159,icol,ipver,lchnk,DP)*((state1%omega(icol,31)-mean_preds(159,icol,lchnk,DP))/std_preds(159,icol,lchnk,DP))+ &
                                     Weights_T(160,icol,ipver,lchnk,DP)*((state1%omega(icol,32)-mean_preds(160,icol,lchnk,DP))/std_preds(160,icol,lchnk,DP))
                                     
                                     
      tend_auxr_t = Weights_T(161,icol,ipver,lchnk,DP)*((upwp(icol,24)-mean_preds(161,icol,lchnk,DP))/std_preds(161,icol,lchnk,DP))+ &
                                     Weights_T(162,icol,ipver,lchnk,DP)*((upwp(icol,25)-mean_preds(162,icol,lchnk,DP))/std_preds(162,icol,lchnk,DP))+ &
                                     Weights_T(163,icol,ipver,lchnk,DP)*((upwp(icol,26)-mean_preds(163,icol,lchnk,DP))/std_preds(163,icol,lchnk,DP))+ &
                                     Weights_T(164,icol,ipver,lchnk,DP)*((upwp(icol,27)-mean_preds(164,icol,lchnk,DP))/std_preds(164,icol,lchnk,DP))+ &
                                     Weights_T(165,icol,ipver,lchnk,DP)*((upwp(icol,28)-mean_preds(165,icol,lchnk,DP))/std_preds(165,icol,lchnk,DP))+ &
                                     Weights_T(166,icol,ipver,lchnk,DP)*((upwp(icol,29)-mean_preds(166,icol,lchnk,DP))/std_preds(166,icol,lchnk,DP))+ &
                                     Weights_T(167,icol,ipver,lchnk,DP)*((upwp(icol,30)-mean_preds(167,icol,lchnk,DP))/std_preds(167,icol,lchnk,DP))+ &
                                     Weights_T(168,icol,ipver,lchnk,DP)*((upwp(icol,31)-mean_preds(168,icol,lchnk,DP))/std_preds(168,icol,lchnk,DP))+ &
                                     Weights_T(169,icol,ipver,lchnk,DP)*((upwp(icol,32)-mean_preds(169,icol,lchnk,DP))/std_preds(169,icol,lchnk,DP))+ &
                                     Weights_T(170,icol,ipver,lchnk,DP)*((upwp(icol,33)-mean_preds(170,icol,lchnk,DP))/std_preds(170,icol,lchnk,DP))+ &
                                     Weights_T(171,icol,ipver,lchnk,DP)*((cam_in%wsx(icol)-mean_preds(171,icol,lchnk,DP))/std_preds(171,icol,lchnk,DP))+ &
                                     Weights_T(172,icol,ipver,lchnk,DP)*((cam_in%wsy(icol)-mean_preds(172,icol,lchnk,DP))/std_preds(172,icol,lchnk,DP))+ &
                                     Weights_T(173,icol,ipver,lchnk,DP)*((flntc(icol)-mean_preds(173,icol,lchnk,DP))/std_preds(173,icol,lchnk,DP))+ &
                                     Weights_T(174,icol,ipver,lchnk,DP)*((fsntoa(icol)-mean_preds(174,icol,lchnk,DP))/std_preds(174,icol,lchnk,DP))+ &
                                     Weights_T(175,icol,ipver,lchnk,DP)*((pblh(icol)-mean_preds(175,icol,lchnk,DP))/std_preds(175,icol,lchnk,DP))+ &
                                     Weights_T(176,icol,ipver,lchnk,DP)*((cam_out%tbot(icol)-mean_preds(176,icol,lchnk,DP))/std_preds(176,icol,lchnk,DP))+ &
                                     Weights_T(177,icol,ipver,lchnk,DP)*((ustarwec(icol)-mean_preds(177,icol,lchnk,DP))/std_preds(177,icol,lchnk,DP))+ &
                                     Weights_T(178,icol,ipver,lchnk,DP)*1.5_r8*DOY2+ &
                                     Weights_T(179,icol,ipver,lchnk,DP)*1.5_r8*DOY1+ &
                                     Weights_T(180,icol,ipver,lchnk,DP)*1.5_r8*HOD2+ &
                                     Weights_T(181,icol,ipver,lchnk,DP)*1.5_r8*HOD1+ &
                                     Weights_T(182,icol,ipver,lchnk,DP)
                       
      
      
      
      ttend_val = (tend_ur_t + tend_vr_t + tend_tr_t + tend_qr_t + tend_or_t + tend_auxr_t)
      
      if (abs(ttend_val) > threshy) then
          DAMLin_Sstep(icol,ipver,lchnk) = sign(threshy, ttend_val)*cpair
      else 
          DAMLin_Sstep(icol,ipver,lchnk) = ttend_val*cpair
      end if 
      
      
      !DAMLin_Sstep(icol,ipver,lchnk) = (tend_ur_t + tend_vr_t + tend_tr_t + tend_qr_t + tend_or_t + tend_auxr_t)*cpair
      
      tend_ur_q = Weights_Q(1,icol,ipver,lchnk,DP)*((state1%u(icol,1)-mean_preds(1,icol,lchnk,DP))/std_preds(1,icol,lchnk,DP))+ &
                                     Weights_Q(2,icol,ipver,lchnk,DP)*((state1%u(icol,2)-mean_preds(2,icol,lchnk,DP))/std_preds(2,icol,lchnk,DP))+ &
                                     Weights_Q(3,icol,ipver,lchnk,DP)*((state1%u(icol,3)-mean_preds(3,icol,lchnk,DP))/std_preds(3,icol,lchnk,DP))+ &
                                     Weights_Q(4,icol,ipver,lchnk,DP)*((state1%u(icol,4)-mean_preds(4,icol,lchnk,DP))/std_preds(4,icol,lchnk,DP))+ &
                                     Weights_Q(5,icol,ipver,lchnk,DP)*((state1%u(icol,5)-mean_preds(5,icol,lchnk,DP))/std_preds(5,icol,lchnk,DP))+ &
                                     Weights_Q(6,icol,ipver,lchnk,DP)*((state1%u(icol,6)-mean_preds(6,icol,lchnk,DP))/std_preds(6,icol,lchnk,DP))+ &
                                     Weights_Q(7,icol,ipver,lchnk,DP)*((state1%u(icol,7)-mean_preds(7,icol,lchnk,DP))/std_preds(7,icol,lchnk,DP))+ &
                                     Weights_Q(8,icol,ipver,lchnk,DP)*((state1%u(icol,8)-mean_preds(8,icol,lchnk,DP))/std_preds(8,icol,lchnk,DP))+ &
                                     Weights_Q(9,icol,ipver,lchnk,DP)*((state1%u(icol,9)-mean_preds(9,icol,lchnk,DP))/std_preds(9,icol,lchnk,DP))+ &
                                     Weights_Q(10,icol,ipver,lchnk,DP)*((state1%u(icol,10)-mean_preds(10,icol,lchnk,DP))/std_preds(10,icol,lchnk,DP))+ &
                                     Weights_Q(11,icol,ipver,lchnk,DP)*((state1%u(icol,11)-mean_preds(11,icol,lchnk,DP))/std_preds(11,icol,lchnk,DP))+ &
                                     Weights_Q(12,icol,ipver,lchnk,DP)*((state1%u(icol,12)-mean_preds(12,icol,lchnk,DP))/std_preds(12,icol,lchnk,DP))+ &
                                     Weights_Q(13,icol,ipver,lchnk,DP)*((state1%u(icol,13)-mean_preds(13,icol,lchnk,DP))/std_preds(13,icol,lchnk,DP))+ &
                                     Weights_Q(14,icol,ipver,lchnk,DP)*((state1%u(icol,14)-mean_preds(14,icol,lchnk,DP))/std_preds(14,icol,lchnk,DP))+ &
                                     Weights_Q(15,icol,ipver,lchnk,DP)*((state1%u(icol,15)-mean_preds(15,icol,lchnk,DP))/std_preds(15,icol,lchnk,DP))+ &
                                     Weights_Q(16,icol,ipver,lchnk,DP)*((state1%u(icol,16)-mean_preds(16,icol,lchnk,DP))/std_preds(16,icol,lchnk,DP))+ &
                                     Weights_Q(17,icol,ipver,lchnk,DP)*((state1%u(icol,17)-mean_preds(17,icol,lchnk,DP))/std_preds(17,icol,lchnk,DP))+ &
                                     Weights_Q(18,icol,ipver,lchnk,DP)*((state1%u(icol,18)-mean_preds(18,icol,lchnk,DP))/std_preds(18,icol,lchnk,DP))+ &
                                     Weights_Q(19,icol,ipver,lchnk,DP)*((state1%u(icol,19)-mean_preds(19,icol,lchnk,DP))/std_preds(19,icol,lchnk,DP))+ &
                                     Weights_Q(20,icol,ipver,lchnk,DP)*((state1%u(icol,20)-mean_preds(20,icol,lchnk,DP))/std_preds(20,icol,lchnk,DP))+ &
                                     Weights_Q(21,icol,ipver,lchnk,DP)*((state1%u(icol,21)-mean_preds(21,icol,lchnk,DP))/std_preds(21,icol,lchnk,DP))+ &
                                     Weights_Q(22,icol,ipver,lchnk,DP)*((state1%u(icol,22)-mean_preds(22,icol,lchnk,DP))/std_preds(22,icol,lchnk,DP))+ &
                                     Weights_Q(23,icol,ipver,lchnk,DP)*((state1%u(icol,23)-mean_preds(23,icol,lchnk,DP))/std_preds(23,icol,lchnk,DP))+ &
                                     Weights_Q(24,icol,ipver,lchnk,DP)*((state1%u(icol,24)-mean_preds(24,icol,lchnk,DP))/std_preds(24,icol,lchnk,DP))+ &
                                     Weights_Q(25,icol,ipver,lchnk,DP)*((state1%u(icol,25)-mean_preds(25,icol,lchnk,DP))/std_preds(25,icol,lchnk,DP))+ &
                                     Weights_Q(26,icol,ipver,lchnk,DP)*((state1%u(icol,26)-mean_preds(26,icol,lchnk,DP))/std_preds(26,icol,lchnk,DP))+ &
                                     Weights_Q(27,icol,ipver,lchnk,DP)*((state1%u(icol,27)-mean_preds(27,icol,lchnk,DP))/std_preds(27,icol,lchnk,DP))+ &
                                     Weights_Q(28,icol,ipver,lchnk,DP)*((state1%u(icol,28)-mean_preds(28,icol,lchnk,DP))/std_preds(28,icol,lchnk,DP))+ &
                                     Weights_Q(29,icol,ipver,lchnk,DP)*((state1%u(icol,29)-mean_preds(29,icol,lchnk,DP))/std_preds(29,icol,lchnk,DP))+ &
                                     Weights_Q(30,icol,ipver,lchnk,DP)*((state1%u(icol,30)-mean_preds(30,icol,lchnk,DP))/std_preds(30,icol,lchnk,DP))+ &
                                     Weights_Q(31,icol,ipver,lchnk,DP)*((state1%u(icol,31)-mean_preds(31,icol,lchnk,DP))/std_preds(31,icol,lchnk,DP))+ &
                                     Weights_Q(32,icol,ipver,lchnk,DP)*((state1%u(icol,32)-mean_preds(32,icol,lchnk,DP))/std_preds(32,icol,lchnk,DP)) 
                                     
                                     
     tend_vr_q = Weights_Q(33,icol,ipver,lchnk,DP)*((state1%v(icol,1)-mean_preds(33,icol,lchnk,DP))/std_preds(33,icol,lchnk,DP))+ &
                                     Weights_Q(34,icol,ipver,lchnk,DP)*((state1%v(icol,2)-mean_preds(34,icol,lchnk,DP))/std_preds(34,icol,lchnk,DP))+ &
                                     Weights_Q(35,icol,ipver,lchnk,DP)*((state1%v(icol,3)-mean_preds(35,icol,lchnk,DP))/std_preds(35,icol,lchnk,DP))+ &
                                     Weights_Q(36,icol,ipver,lchnk,DP)*((state1%v(icol,4)-mean_preds(36,icol,lchnk,DP))/std_preds(36,icol,lchnk,DP))+ &
                                     Weights_Q(37,icol,ipver,lchnk,DP)*((state1%v(icol,5)-mean_preds(37,icol,lchnk,DP))/std_preds(37,icol,lchnk,DP))+ &
                                     Weights_Q(38,icol,ipver,lchnk,DP)*((state1%v(icol,6)-mean_preds(38,icol,lchnk,DP))/std_preds(38,icol,lchnk,DP))+ &
                                     Weights_Q(39,icol,ipver,lchnk,DP)*((state1%v(icol,7)-mean_preds(39,icol,lchnk,DP))/std_preds(39,icol,lchnk,DP))+ &
                                     Weights_Q(40,icol,ipver,lchnk,DP)*((state1%v(icol,8)-mean_preds(40,icol,lchnk,DP))/std_preds(40,icol,lchnk,DP))+ &
                                     Weights_Q(41,icol,ipver,lchnk,DP)*((state1%v(icol,9)-mean_preds(41,icol,lchnk,DP))/std_preds(41,icol,lchnk,DP))+ &
                                     Weights_Q(42,icol,ipver,lchnk,DP)*((state1%v(icol,10)-mean_preds(42,icol,lchnk,DP))/std_preds(42,icol,lchnk,DP))+ &
                                     Weights_Q(43,icol,ipver,lchnk,DP)*((state1%v(icol,11)-mean_preds(43,icol,lchnk,DP))/std_preds(43,icol,lchnk,DP))+ &
                                     Weights_Q(44,icol,ipver,lchnk,DP)*((state1%v(icol,12)-mean_preds(44,icol,lchnk,DP))/std_preds(44,icol,lchnk,DP))+ &
                                     Weights_Q(45,icol,ipver,lchnk,DP)*((state1%v(icol,13)-mean_preds(45,icol,lchnk,DP))/std_preds(45,icol,lchnk,DP))+ &
                                     Weights_Q(46,icol,ipver,lchnk,DP)*((state1%v(icol,14)-mean_preds(46,icol,lchnk,DP))/std_preds(46,icol,lchnk,DP))+ &
                                     Weights_Q(47,icol,ipver,lchnk,DP)*((state1%v(icol,15)-mean_preds(47,icol,lchnk,DP))/std_preds(47,icol,lchnk,DP))+ &
                                     Weights_Q(48,icol,ipver,lchnk,DP)*((state1%v(icol,16)-mean_preds(48,icol,lchnk,DP))/std_preds(48,icol,lchnk,DP))+ &
                                     Weights_Q(49,icol,ipver,lchnk,DP)*((state1%v(icol,17)-mean_preds(49,icol,lchnk,DP))/std_preds(49,icol,lchnk,DP))+ &
                                     Weights_Q(50,icol,ipver,lchnk,DP)*((state1%v(icol,18)-mean_preds(50,icol,lchnk,DP))/std_preds(50,icol,lchnk,DP))+ &
                                     Weights_Q(51,icol,ipver,lchnk,DP)*((state1%v(icol,19)-mean_preds(51,icol,lchnk,DP))/std_preds(51,icol,lchnk,DP))+ &
                                     Weights_Q(52,icol,ipver,lchnk,DP)*((state1%v(icol,20)-mean_preds(52,icol,lchnk,DP))/std_preds(52,icol,lchnk,DP))+ &
                                     Weights_Q(53,icol,ipver,lchnk,DP)*((state1%v(icol,21)-mean_preds(53,icol,lchnk,DP))/std_preds(53,icol,lchnk,DP))+ &
                                     Weights_Q(54,icol,ipver,lchnk,DP)*((state1%v(icol,22)-mean_preds(54,icol,lchnk,DP))/std_preds(54,icol,lchnk,DP))+ &
                                     Weights_Q(55,icol,ipver,lchnk,DP)*((state1%v(icol,23)-mean_preds(55,icol,lchnk,DP))/std_preds(55,icol,lchnk,DP))+ &
                                     Weights_Q(56,icol,ipver,lchnk,DP)*((state1%v(icol,24)-mean_preds(56,icol,lchnk,DP))/std_preds(56,icol,lchnk,DP))+ &
                                     Weights_Q(57,icol,ipver,lchnk,DP)*((state1%v(icol,25)-mean_preds(57,icol,lchnk,DP))/std_preds(57,icol,lchnk,DP))+ &
                                     Weights_Q(58,icol,ipver,lchnk,DP)*((state1%v(icol,26)-mean_preds(58,icol,lchnk,DP))/std_preds(58,icol,lchnk,DP))+ &
                                     Weights_Q(59,icol,ipver,lchnk,DP)*((state1%v(icol,27)-mean_preds(59,icol,lchnk,DP))/std_preds(59,icol,lchnk,DP))+ &
                                     Weights_Q(60,icol,ipver,lchnk,DP)*((state1%v(icol,28)-mean_preds(60,icol,lchnk,DP))/std_preds(60,icol,lchnk,DP))+ &
                                     Weights_Q(61,icol,ipver,lchnk,DP)*((state1%v(icol,29)-mean_preds(61,icol,lchnk,DP))/std_preds(61,icol,lchnk,DP))+ &
                                     Weights_Q(62,icol,ipver,lchnk,DP)*((state1%v(icol,30)-mean_preds(62,icol,lchnk,DP))/std_preds(62,icol,lchnk,DP))+ &
                                     Weights_Q(63,icol,ipver,lchnk,DP)*((state1%v(icol,31)-mean_preds(63,icol,lchnk,DP))/std_preds(63,icol,lchnk,DP))+ &
                                     Weights_Q(64,icol,ipver,lchnk,DP)*((state1%v(icol,32)-mean_preds(64,icol,lchnk,DP))/std_preds(64,icol,lchnk,DP))
                                                                  
                                
      tend_tr_q = Weights_Q(65,icol,ipver,lchnk,DP)*((state1%t(icol,1)-mean_preds(65,icol,lchnk,DP))/std_preds(65,icol,lchnk,DP))+ &
                                     Weights_Q(66,icol,ipver,lchnk,DP)*((state1%t(icol,2)-mean_preds(66,icol,lchnk,DP))/std_preds(66,icol,lchnk,DP))+ &
                                     Weights_Q(67,icol,ipver,lchnk,DP)*((state1%t(icol,3)-mean_preds(67,icol,lchnk,DP))/std_preds(67,icol,lchnk,DP))+ &
                                     Weights_Q(68,icol,ipver,lchnk,DP)*((state1%t(icol,4)-mean_preds(68,icol,lchnk,DP))/std_preds(68,icol,lchnk,DP))+ &
                                     Weights_Q(69,icol,ipver,lchnk,DP)*((state1%t(icol,5)-mean_preds(69,icol,lchnk,DP))/std_preds(69,icol,lchnk,DP))+ &
                                     Weights_Q(70,icol,ipver,lchnk,DP)*((state1%t(icol,6)-mean_preds(70,icol,lchnk,DP))/std_preds(70,icol,lchnk,DP))+ &
                                     Weights_Q(71,icol,ipver,lchnk,DP)*((state1%t(icol,7)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_Q(71,icol,ipver,lchnk,DP)*((state1%t(icol,8)-mean_preds(71,icol,lchnk,DP))/std_preds(71,icol,lchnk,DP))+ &
                                     Weights_Q(73,icol,ipver,lchnk,DP)*((state1%t(icol,9)-mean_preds(73,icol,lchnk,DP))/std_preds(73,icol,lchnk,DP))+ &
                                     Weights_Q(74,icol,ipver,lchnk,DP)*((state1%t(icol,10)-mean_preds(74,icol,lchnk,DP))/std_preds(74,icol,lchnk,DP))+ &
                                     Weights_Q(75,icol,ipver,lchnk,DP)*((state1%t(icol,11)-mean_preds(75,icol,lchnk,DP))/std_preds(75,icol,lchnk,DP))+ &
                                     Weights_Q(76,icol,ipver,lchnk,DP)*((state1%t(icol,12)-mean_preds(76,icol,lchnk,DP))/std_preds(76,icol,lchnk,DP))+ &
                                     Weights_Q(77,icol,ipver,lchnk,DP)*((state1%t(icol,13)-mean_preds(77,icol,lchnk,DP))/std_preds(77,icol,lchnk,DP))+ &
                                     Weights_Q(78,icol,ipver,lchnk,DP)*((state1%t(icol,14)-mean_preds(78,icol,lchnk,DP))/std_preds(78,icol,lchnk,DP))+ &
                                     Weights_Q(79,icol,ipver,lchnk,DP)*((state1%t(icol,15)-mean_preds(79,icol,lchnk,DP))/std_preds(79,icol,lchnk,DP))+ &
                                     Weights_Q(80,icol,ipver,lchnk,DP)*((state1%t(icol,16)-mean_preds(80,icol,lchnk,DP))/std_preds(80,icol,lchnk,DP))+ &
                                     Weights_Q(81,icol,ipver,lchnk,DP)*((state1%t(icol,17)-mean_preds(81,icol,lchnk,DP))/std_preds(81,icol,lchnk,DP))+ &
                                     Weights_Q(82,icol,ipver,lchnk,DP)*((state1%t(icol,18)-mean_preds(82,icol,lchnk,DP))/std_preds(82,icol,lchnk,DP))+ &
                                     Weights_Q(83,icol,ipver,lchnk,DP)*((state1%t(icol,19)-mean_preds(83,icol,lchnk,DP))/std_preds(83,icol,lchnk,DP))+ &
                                     Weights_Q(84,icol,ipver,lchnk,DP)*((state1%t(icol,20)-mean_preds(84,icol,lchnk,DP))/std_preds(84,icol,lchnk,DP))+ &
                                     Weights_Q(85,icol,ipver,lchnk,DP)*((state1%t(icol,21)-mean_preds(85,icol,lchnk,DP))/std_preds(85,icol,lchnk,DP))+ &
                                     Weights_Q(86,icol,ipver,lchnk,DP)*((state1%t(icol,22)-mean_preds(86,icol,lchnk,DP))/std_preds(86,icol,lchnk,DP))+ &
                                     Weights_Q(87,icol,ipver,lchnk,DP)*((state1%t(icol,23)-mean_preds(87,icol,lchnk,DP))/std_preds(87,icol,lchnk,DP))+ &
                                     Weights_Q(88,icol,ipver,lchnk,DP)*((state1%t(icol,24)-mean_preds(88,icol,lchnk,DP))/std_preds(88,icol,lchnk,DP))+ &
                                     Weights_Q(89,icol,ipver,lchnk,DP)*((state1%t(icol,25)-mean_preds(89,icol,lchnk,DP))/std_preds(89,icol,lchnk,DP))+ &
                                     Weights_Q(90,icol,ipver,lchnk,DP)*((state1%t(icol,26)-mean_preds(90,icol,lchnk,DP))/std_preds(90,icol,lchnk,DP))+ &
                                     Weights_Q(91,icol,ipver,lchnk,DP)*((state1%t(icol,27)-mean_preds(91,icol,lchnk,DP))/std_preds(91,icol,lchnk,DP))+ &
                                     Weights_Q(92,icol,ipver,lchnk,DP)*((state1%t(icol,28)-mean_preds(92,icol,lchnk,DP))/std_preds(92,icol,lchnk,DP))+ &
                                     Weights_Q(93,icol,ipver,lchnk,DP)*((state1%t(icol,29)-mean_preds(93,icol,lchnk,DP))/std_preds(93,icol,lchnk,DP))+ &
                                     Weights_Q(94,icol,ipver,lchnk,DP)*((state1%t(icol,30)-mean_preds(94,icol,lchnk,DP))/std_preds(94,icol,lchnk,DP))+ &
                                     Weights_Q(95,icol,ipver,lchnk,DP)*((state1%t(icol,31)-mean_preds(95,icol,lchnk,DP))/std_preds(95,icol,lchnk,DP))+ &
                                     Weights_Q(96,icol,ipver,lchnk,DP)*((state1%t(icol,32)-mean_preds(96,icol,lchnk,DP))/std_preds(96,icol,lchnk,DP))
                                     
      tend_qr_q = Weights_Q(97,icol,ipver,lchnk,DP)*((state1%q(icol,1,1)-mean_preds(97,icol,lchnk,DP))/std_preds(97,icol,lchnk,DP))+ &
                                     Weights_Q(98,icol,ipver,lchnk,DP)*((state1%q(icol,2,1)-mean_preds(98,icol,lchnk,DP))/std_preds(98,icol,lchnk,DP))+ &
                                     Weights_Q(99,icol,ipver,lchnk,DP)*((state1%q(icol,3,1)-mean_preds(99,icol,lchnk,DP))/std_preds(99,icol,lchnk,DP))+ &
                                     Weights_Q(100,icol,ipver,lchnk,DP)*((state1%q(icol,4,1)-mean_preds(100,icol,lchnk,DP))/std_preds(100,icol,lchnk,DP))+ &
                                     Weights_Q(101,icol,ipver,lchnk,DP)*((state1%q(icol,5,1)-mean_preds(101,icol,lchnk,DP))/std_preds(101,icol,lchnk,DP))+ &
                                     Weights_Q(102,icol,ipver,lchnk,DP)*((state1%q(icol,6,1)-mean_preds(102,icol,lchnk,DP))/std_preds(102,icol,lchnk,DP))+ &
                                     Weights_Q(103,icol,ipver,lchnk,DP)*((state1%q(icol,7,1)-mean_preds(103,icol,lchnk,DP))/std_preds(103,icol,lchnk,DP))+ &
                                     Weights_Q(104,icol,ipver,lchnk,DP)*((state1%q(icol,8,1)-mean_preds(104,icol,lchnk,DP))/std_preds(104,icol,lchnk,DP))+ &
                                     Weights_Q(105,icol,ipver,lchnk,DP)*((state1%q(icol,9,1)-mean_preds(105,icol,lchnk,DP))/std_preds(105,icol,lchnk,DP))+ &
                                     Weights_Q(106,icol,ipver,lchnk,DP)*((state1%q(icol,10,1)-mean_preds(106,icol,lchnk,DP))/std_preds(106,icol,lchnk,DP))+ &
                                     Weights_Q(107,icol,ipver,lchnk,DP)*((state1%q(icol,11,1)-mean_preds(107,icol,lchnk,DP))/std_preds(107,icol,lchnk,DP))+ &
                                     Weights_Q(108,icol,ipver,lchnk,DP)*((state1%q(icol,12,1)-mean_preds(108,icol,lchnk,DP))/std_preds(108,icol,lchnk,DP))+ &
                                     Weights_Q(109,icol,ipver,lchnk,DP)*((state1%q(icol,13,1)-mean_preds(109,icol,lchnk,DP))/std_preds(109,icol,lchnk,DP))+ &
                                     Weights_Q(110,icol,ipver,lchnk,DP)*((state1%q(icol,14,1)-mean_preds(110,icol,lchnk,DP))/std_preds(110,icol,lchnk,DP))+ &
                                     Weights_Q(111,icol,ipver,lchnk,DP)*((state1%q(icol,15,1)-mean_preds(111,icol,lchnk,DP))/std_preds(111,icol,lchnk,DP))+ &
                                     Weights_Q(112,icol,ipver,lchnk,DP)*((state1%q(icol,16,1)-mean_preds(112,icol,lchnk,DP))/std_preds(112,icol,lchnk,DP))+ &
                                     Weights_Q(113,icol,ipver,lchnk,DP)*((state1%q(icol,17,1)-mean_preds(113,icol,lchnk,DP))/std_preds(113,icol,lchnk,DP))+ &
                                     Weights_Q(114,icol,ipver,lchnk,DP)*((state1%q(icol,18,1)-mean_preds(114,icol,lchnk,DP))/std_preds(114,icol,lchnk,DP))+ &
                                     Weights_Q(115,icol,ipver,lchnk,DP)*((state1%q(icol,19,1)-mean_preds(115,icol,lchnk,DP))/std_preds(115,icol,lchnk,DP))+ &
                                     Weights_Q(116,icol,ipver,lchnk,DP)*((state1%q(icol,20,1)-mean_preds(116,icol,lchnk,DP))/std_preds(116,icol,lchnk,DP))+ &
                                     Weights_Q(117,icol,ipver,lchnk,DP)*((state1%q(icol,21,1)-mean_preds(117,icol,lchnk,DP))/std_preds(117,icol,lchnk,DP))+ &
                                     Weights_Q(118,icol,ipver,lchnk,DP)*((state1%q(icol,22,1)-mean_preds(118,icol,lchnk,DP))/std_preds(118,icol,lchnk,DP))+ &
                                     Weights_Q(119,icol,ipver,lchnk,DP)*((state1%q(icol,23,1)-mean_preds(119,icol,lchnk,DP))/std_preds(119,icol,lchnk,DP))+ &
                                     Weights_Q(120,icol,ipver,lchnk,DP)*((state1%q(icol,24,1)-mean_preds(120,icol,lchnk,DP))/std_preds(120,icol,lchnk,DP))+ &
                                     Weights_Q(121,icol,ipver,lchnk,DP)*((state1%q(icol,25,1)-mean_preds(121,icol,lchnk,DP))/std_preds(121,icol,lchnk,DP))+ &
                                     Weights_Q(122,icol,ipver,lchnk,DP)*((state1%q(icol,26,1)-mean_preds(122,icol,lchnk,DP))/std_preds(122,icol,lchnk,DP))+ &
                                     Weights_Q(123,icol,ipver,lchnk,DP)*((state1%q(icol,27,1)-mean_preds(123,icol,lchnk,DP))/std_preds(123,icol,lchnk,DP))+ &
                                     Weights_Q(124,icol,ipver,lchnk,DP)*((state1%q(icol,28,1)-mean_preds(124,icol,lchnk,DP))/std_preds(124,icol,lchnk,DP))+ &
                                     Weights_Q(125,icol,ipver,lchnk,DP)*((state1%q(icol,29,1)-mean_preds(125,icol,lchnk,DP))/std_preds(125,icol,lchnk,DP))+ &
                                     Weights_Q(126,icol,ipver,lchnk,DP)*((state1%q(icol,30,1)-mean_preds(126,icol,lchnk,DP))/std_preds(126,icol,lchnk,DP))+ &
                                     Weights_Q(127,icol,ipver,lchnk,DP)*((state1%q(icol,31,1)-mean_preds(127,icol,lchnk,DP))/std_preds(127,icol,lchnk,DP))+ &
                                     Weights_Q(128,icol,ipver,lchnk,DP)*((state1%q(icol,32,1)-mean_preds(128,icol,lchnk,DP))/std_preds(128,icol,lchnk,DP))
                                     
      tend_or_q = Weights_Q(129,icol,ipver,lchnk,DP)*((state1%omega(icol,1)-mean_preds(129,icol,lchnk,DP))/std_preds(129,icol,lchnk,DP))+ &
                                     Weights_Q(130,icol,ipver,lchnk,DP)*((state1%omega(icol,2)-mean_preds(130,icol,lchnk,DP))/std_preds(130,icol,lchnk,DP))+ &
                                     Weights_Q(131,icol,ipver,lchnk,DP)*((state1%omega(icol,3)-mean_preds(131,icol,lchnk,DP))/std_preds(131,icol,lchnk,DP))+ &
                                     Weights_Q(132,icol,ipver,lchnk,DP)*((state1%omega(icol,4)-mean_preds(132,icol,lchnk,DP))/std_preds(132,icol,lchnk,DP))+ &
                                     Weights_Q(133,icol,ipver,lchnk,DP)*((state1%omega(icol,5)-mean_preds(133,icol,lchnk,DP))/std_preds(133,icol,lchnk,DP))+ &
                                     Weights_Q(134,icol,ipver,lchnk,DP)*((state1%omega(icol,6)-mean_preds(134,icol,lchnk,DP))/std_preds(134,icol,lchnk,DP))+ &
                                     Weights_Q(135,icol,ipver,lchnk,DP)*((state1%omega(icol,7)-mean_preds(135,icol,lchnk,DP))/std_preds(135,icol,lchnk,DP))+ &
                                     Weights_Q(136,icol,ipver,lchnk,DP)*((state1%omega(icol,8)-mean_preds(136,icol,lchnk,DP))/std_preds(136,icol,lchnk,DP))+ &
                                     Weights_Q(137,icol,ipver,lchnk,DP)*((state1%omega(icol,9)-mean_preds(137,icol,lchnk,DP))/std_preds(137,icol,lchnk,DP))+ &
                                     Weights_Q(138,icol,ipver,lchnk,DP)*((state1%omega(icol,10)-mean_preds(138,icol,lchnk,DP))/std_preds(138,icol,lchnk,DP))+ &
                                     Weights_Q(139,icol,ipver,lchnk,DP)*((state1%omega(icol,11)-mean_preds(139,icol,lchnk,DP))/std_preds(139,icol,lchnk,DP))+ &
                                     Weights_Q(140,icol,ipver,lchnk,DP)*((state1%omega(icol,12)-mean_preds(140,icol,lchnk,DP))/std_preds(140,icol,lchnk,DP))+ &
                                     Weights_Q(141,icol,ipver,lchnk,DP)*((state1%omega(icol,13)-mean_preds(141,icol,lchnk,DP))/std_preds(141,icol,lchnk,DP))+ &
                                     Weights_Q(142,icol,ipver,lchnk,DP)*((state1%omega(icol,14)-mean_preds(142,icol,lchnk,DP))/std_preds(142,icol,lchnk,DP))+ &
                                     Weights_Q(143,icol,ipver,lchnk,DP)*((state1%omega(icol,15)-mean_preds(143,icol,lchnk,DP))/std_preds(143,icol,lchnk,DP))+ &
                                     Weights_Q(144,icol,ipver,lchnk,DP)*((state1%omega(icol,16)-mean_preds(144,icol,lchnk,DP))/std_preds(144,icol,lchnk,DP))+ &
                                     Weights_Q(145,icol,ipver,lchnk,DP)*((state1%omega(icol,17)-mean_preds(145,icol,lchnk,DP))/std_preds(145,icol,lchnk,DP))+ &
                                     Weights_Q(146,icol,ipver,lchnk,DP)*((state1%omega(icol,18)-mean_preds(146,icol,lchnk,DP))/std_preds(146,icol,lchnk,DP))+ &
                                     Weights_Q(147,icol,ipver,lchnk,DP)*((state1%omega(icol,19)-mean_preds(147,icol,lchnk,DP))/std_preds(147,icol,lchnk,DP))+ &
                                     Weights_Q(148,icol,ipver,lchnk,DP)*((state1%omega(icol,20)-mean_preds(148,icol,lchnk,DP))/std_preds(148,icol,lchnk,DP))+ &
                                     Weights_Q(149,icol,ipver,lchnk,DP)*((state1%omega(icol,21)-mean_preds(149,icol,lchnk,DP))/std_preds(149,icol,lchnk,DP))+ &
                                     Weights_Q(150,icol,ipver,lchnk,DP)*((state1%omega(icol,22)-mean_preds(150,icol,lchnk,DP))/std_preds(150,icol,lchnk,DP))+ &
                                     Weights_Q(151,icol,ipver,lchnk,DP)*((state1%omega(icol,23)-mean_preds(151,icol,lchnk,DP))/std_preds(151,icol,lchnk,DP))+ &
                                     Weights_Q(152,icol,ipver,lchnk,DP)*((state1%omega(icol,24)-mean_preds(152,icol,lchnk,DP))/std_preds(152,icol,lchnk,DP))+ &
                                     Weights_Q(153,icol,ipver,lchnk,DP)*((state1%omega(icol,25)-mean_preds(153,icol,lchnk,DP))/std_preds(153,icol,lchnk,DP))+ &
                                     Weights_Q(154,icol,ipver,lchnk,DP)*((state1%omega(icol,26)-mean_preds(154,icol,lchnk,DP))/std_preds(154,icol,lchnk,DP))+ &
                                     Weights_Q(155,icol,ipver,lchnk,DP)*((state1%omega(icol,27)-mean_preds(155,icol,lchnk,DP))/std_preds(155,icol,lchnk,DP))+ &
                                     Weights_Q(156,icol,ipver,lchnk,DP)*((state1%omega(icol,28)-mean_preds(156,icol,lchnk,DP))/std_preds(156,icol,lchnk,DP))+ &
                                     Weights_Q(157,icol,ipver,lchnk,DP)*((state1%omega(icol,29)-mean_preds(157,icol,lchnk,DP))/std_preds(157,icol,lchnk,DP))+ &
                                     Weights_Q(158,icol,ipver,lchnk,DP)*((state1%omega(icol,30)-mean_preds(158,icol,lchnk,DP))/std_preds(158,icol,lchnk,DP))+ &
                                     Weights_Q(159,icol,ipver,lchnk,DP)*((state1%omega(icol,31)-mean_preds(159,icol,lchnk,DP))/std_preds(159,icol,lchnk,DP))+ &
                                     Weights_Q(160,icol,ipver,lchnk,DP)*((state1%omega(icol,32)-mean_preds(160,icol,lchnk,DP))/std_preds(160,icol,lchnk,DP))
                                     
                                     
      tend_auxr_q = Weights_Q(161,icol,ipver,lchnk,DP)*((upwp(icol,24)-mean_preds(161,icol,lchnk,DP))/std_preds(161,icol,lchnk,DP))+ &
                                     Weights_Q(162,icol,ipver,lchnk,DP)*((upwp(icol,25)-mean_preds(162,icol,lchnk,DP))/std_preds(162,icol,lchnk,DP))+ &
                                     Weights_Q(163,icol,ipver,lchnk,DP)*((upwp(icol,26)-mean_preds(163,icol,lchnk,DP))/std_preds(163,icol,lchnk,DP))+ &
                                     Weights_Q(164,icol,ipver,lchnk,DP)*((upwp(icol,27)-mean_preds(164,icol,lchnk,DP))/std_preds(164,icol,lchnk,DP))+ &
                                     Weights_Q(165,icol,ipver,lchnk,DP)*((upwp(icol,28)-mean_preds(165,icol,lchnk,DP))/std_preds(165,icol,lchnk,DP))+ &
                                     Weights_Q(166,icol,ipver,lchnk,DP)*((upwp(icol,29)-mean_preds(166,icol,lchnk,DP))/std_preds(166,icol,lchnk,DP))+ &
                                     Weights_Q(167,icol,ipver,lchnk,DP)*((upwp(icol,30)-mean_preds(167,icol,lchnk,DP))/std_preds(167,icol,lchnk,DP))+ &
                                     Weights_Q(168,icol,ipver,lchnk,DP)*((upwp(icol,31)-mean_preds(168,icol,lchnk,DP))/std_preds(168,icol,lchnk,DP))+ &
                                     Weights_Q(169,icol,ipver,lchnk,DP)*((upwp(icol,32)-mean_preds(169,icol,lchnk,DP))/std_preds(169,icol,lchnk,DP))+ &
                                     Weights_Q(170,icol,ipver,lchnk,DP)*((upwp(icol,33)-mean_preds(170,icol,lchnk,DP))/std_preds(170,icol,lchnk,DP))+ &
                                     Weights_Q(171,icol,ipver,lchnk,DP)*((cam_in%wsx(icol)-mean_preds(171,icol,lchnk,DP))/std_preds(171,icol,lchnk,DP))+ &
                                     Weights_Q(172,icol,ipver,lchnk,DP)*((cam_in%wsy(icol)-mean_preds(172,icol,lchnk,DP))/std_preds(172,icol,lchnk,DP))+ &
                                     Weights_Q(173,icol,ipver,lchnk,DP)*((flntc(icol)-mean_preds(173,icol,lchnk,DP))/std_preds(173,icol,lchnk,DP))+ &
                                     Weights_Q(174,icol,ipver,lchnk,DP)*((fsntoa(icol)-mean_preds(174,icol,lchnk,DP))/std_preds(174,icol,lchnk,DP))+ &
                                     Weights_Q(175,icol,ipver,lchnk,DP)*((pblh(icol)-mean_preds(175,icol,lchnk,DP))/std_preds(175,icol,lchnk,DP))+ &
                                     Weights_Q(176,icol,ipver,lchnk,DP)*((cam_out%tbot(icol)-mean_preds(176,icol,lchnk,DP))/std_preds(176,icol,lchnk,DP))+ &
                                     Weights_Q(177,icol,ipver,lchnk,DP)*((ustarwec(icol)-mean_preds(177,icol,lchnk,DP))/std_preds(177,icol,lchnk,DP))+ &
                                     Weights_Q(178,icol,ipver,lchnk,DP)*1.5_r8*DOY2+ &
                                     Weights_Q(179,icol,ipver,lchnk,DP)*1.5_r8*DOY1+ &
                                     Weights_Q(180,icol,ipver,lchnk,DP)*1.5_r8*HOD2+ &
                                     Weights_Q(181,icol,ipver,lchnk,DP)*1.5_r8*HOD1+ &
                                     Weights_Q(182,icol,ipver,lchnk,DP)
                             
      if (abs((tend_ur_q + tend_vr_q + tend_tr_q + tend_qr_q + tend_or_q + tend_auxr_q)) > threshy) then
          DAMLin_Qstep(icol,ipver,lchnk) = sign(threshy, (tend_ur_q + tend_vr_q + tend_tr_q + tend_qr_q + tend_or_q + tend_auxr_q))
      else 
          DAMLin_Qstep(icol,ipver,lchnk) = (tend_ur_q + tend_vr_q + tend_tr_q + tend_qr_q + tend_or_q + tend_auxr_q)
      end if 

  end do
  end do
  
  do lchnk=begchunk,endchunk !Remove:
       DAMLin_Ustep(:ncol,:pver,lchnk)=(DAMLin_Ustep(:ncol,:pver,lchnk))*DAMLin_Utau(:ncol,:pver,lchnk) !These are all tendencies already [WEC]
       DAMLin_Vstep(:ncol,:pver,lchnk)=(DAMLin_Vstep(:ncol,:pver,lchnk))*DAMLin_Vtau(:ncol,:pver,lchnk)
       DAMLin_Sstep(:ncol,:pver,lchnk)=(DAMLin_Sstep(:ncol,:pver,lchnk))*DAMLin_Stau(:ncol,:pver,lchnk) !cant do Temp [WEC]
       DAMLin_Qstep(:ncol,:pver,lchnk)=(DAMLin_Qstep(:ncol,:pver,lchnk))*DAMLin_Qtau(:ncol,:pver,lchnk)
       !DAMLin_PSstep(:ncol,     lchnk)=(Target_stoch_PS(:ncol,lchnk))*Stochai_PStau(:ncol,lchnk) !cant do PS [WEC]
  end do
  
  
  end subroutine regress_daml_timestep_tend
  
  
  
  subroutine damlining_timestep_tend(phys_state,phys_tend) 
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the DAMLin
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'damling',lu=.true.,lv=.true.,ls=.true.,lq=lq)
   DAMLin_ON=.true. !TODO WEC
   if(DAMLin_ON) then !this is where the nudging tendency gets teleported into the model. [WEC]
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =DAMLin_Ustep(:ncol,:pver,lchnk)/700.0_r8
     phys_tend%v(:ncol,:pver)     =DAMLin_Vstep(:ncol,:pver,lchnk)/700.0_r8
     phys_tend%s(:ncol,:pver)     =DAMLin_Sstep(:ncol,:pver,lchnk)/700.0_r8
     phys_tend%q(:ncol,:pver,indw)=DAMLin_Qstep(:ncol,:pver,lchnk)/700.0_r8

     call outfld( 'DAMLin_U',phys_tend%u                ,pcols,lchnk)
     call outfld( 'DAMLin_V',phys_tend%v                ,pcols,lchnk)
     call outfld( 'DAMLin_T',phys_tend%s/cpair          ,pcols,lchnk) !MUST EDIT TO MESS WITH T->DSE [WEC]
     call outfld( 'DAMLin_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)
     
     
     !call outfld('Target_U',Target_U(:,:,lchnk),pcols,lchnk)
     !call outfld('Target_V',Target_V(:,:,lchnk),pcols,lchnk)
     !call outfld('Target_T',Target_T(:,:,lchnk),pcols,lchnk)
     !call outfld('Target_Q',Target_Q(:,:,lchnk),pcols,lchnk)
   endif
   DAMLin_ON=.False. !TODO WEC

   ! End Routine
   !------------
   return
  end subroutine !damlining_timestep_tend
  
  subroutine damlining_diurnal_timestep_tend(phys_state,phys_tend) 
   !
   ! NUDGING_TIMESTEP_TEND:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the DAMLin
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'damling',lu=.true.,lv=.true.,ls=.true.,lq=lq)
   DAMLin_ON=.true. !TODO WEC
   if(DAMLin_ON) then !this is where the nudging tendency gets teleported into the model. [WEC]
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =DAMLin_Ustep_avg(:ncol,:pver,lchnk)/100.0_r8
     phys_tend%v(:ncol,:pver)     =DAMLin_Vstep_avg(:ncol,:pver,lchnk)/100.0_r8
     phys_tend%s(:ncol,:pver)     =DAMLin_Sstep_avg(:ncol,:pver,lchnk)/100.0_r8
     phys_tend%q(:ncol,:pver,indw)=DAMLin_Qstep_avg(:ncol,:pver,lchnk)/100.0_r8

     call outfld( 'DAMLin_U_avg',phys_tend%u                ,pcols,lchnk)
     call outfld( 'DAMLin_V_avg',phys_tend%v                ,pcols,lchnk)
     call outfld( 'DAMLin_T_avg',phys_tend%s/cpair          ,pcols,lchnk) !MUST EDIT TO MESS WITH T->DSE [WEC]
     call outfld( 'DAMLin_Q_avg',phys_tend%q(1,1,indw)      ,pcols,lchnk)
     
   endif
   DAMLin_ON=.False. !TODO WEC

   ! End Routine
   !------------
   return
  end subroutine !damlining_diurnal_timestep_tend
    !================================================================
  
  
  !================================================================
  subroutine DAMLining_set_profile(rlat,rlon,DAMLin_prof,Wprof,nlev)
   !
   ! DAMLinING_SET_PROFILE: for the given lat,lon, and DAMLining_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on DAMLining strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,DAMLin_prof
   real(r8) rlat,rlon
   real(r8) Wprof(nlev)

   ! Local values
   !----------------
   integer  ilev
   real(r8) Hcoef,latx,lonx,Vmax,Vmin
   real(r8) lon_lo,lon_hi,lat_lo,lat_hi,lev_lo,lev_hi

   !---------------
   ! set coeffcient
   !---------------
   if(DAMLin_prof.eq.0) then
     ! No DAMLining
     !-------------
     Wprof(:)=0.0_r8
   elseif(DAMLin_prof.eq.1) then
     ! Uniform DAMLining
     !-----------------
     Wprof(:)=1.0_r8
   elseif(DAMLin_prof.eq.2) then
     ! Localized DAMLining with specified Heaviside window function
     !------------------------------------------------------------
     if(DAMLin_Hwin_max.le.DAMLin_Hwin_min) then
       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(DAMLin_Hwin_lo,DAMLin_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-DAMLin_Hwin_lat0
       lonx=rlon-DAMLin_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(DAMLin_Hwin_lonWidthH+lonx)/DAMLin_Hwin_lonDelta
       lon_hi=(DAMLin_Hwin_lonWidthH-lonx)/DAMLin_Hwin_lonDelta
       lat_lo=(DAMLin_Hwin_latWidthH+latx)/DAMLin_Hwin_latDelta
       lat_hi=(DAMLin_Hwin_latWidthH-latx)/DAMLin_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-DAMLin_Hwin_min)/(DAMLin_Hwin_max-DAMLin_Hwin_min)
       Hcoef=(1._r8-Hcoef)*DAMLin_Hwin_lo + Hcoef*DAMLin_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-DAMLin_Vwin_Lindex)/DAMLin_Vwin_Ldelta
       lev_hi=(DAMLin_Vwin_Hindex-float(ilev))/DAMLin_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((DAMLin_Vwin_Hindex.ge.(nlev+1)).and. &
                           (DAMLin_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(DAMLin_Vwin_lo,DAMLin_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=DAMLin_Vwin_lo + Wprof(:)*(DAMLin_Vwin_hi-DAMLin_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('DAMLining_set_profile:: Unknown DAMLin_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! DAMLining_set_profile
  
  
  
  !================================================================
  real(r8) function DAMLining_set_PSprofile(rlat,rlon,DAMLin_PSprof)
   !
   ! DAMLinING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on DAMLining strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  DAMLin_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(DAMLin_PSprof.eq.0) then
     ! No DAMLining
     !-------------
     DAMLining_set_PSprofile=0.0_r8
   elseif(DAMLin_PSprof.eq.1) then
     ! Uniform DAMLining
     !-----------------
     DAMLining_set_PSprofile=1.0_r8
   else
     call endrun('DAMLining_set_PSprofile:: Unknown DAMLin_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! DAMLining_set_PSprofile
  !================================================================
  
  
  
 !================================================================
  subroutine damlining_timestep_init(phys_state)
   ! The ONLY thing this should be doing is setting DAMLin_ON
   ! which should use the time functionality to say when it's on or not....
   ! damlining_TIMESTEP_INIT:
   !                 Check the current time and update Model/Stochaiing
   !                 arrays when necessary. Toggle the Stochaiing flag
   !                 when the time is withing the stochaiing window.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state
   use constituents ,only: cnst_get_ind
   use dycore       ,only: dycore_is
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use filenames    ,only: interpret_filename_spec
   use ESMF

   ! Arguments
   !-----------
   type(physics_state),intent(in):: phys_state(begchunk:endchunk)

   ! Local values
   !----------------
   integer Year,Month,Day,Sec
   integer YMD1,YMD2,YMD,YMD3,YMD4,YMD5,YMD6
   logical Update_Model,Update_DAMLin,Sync_Error
   logical After_Beg   ,Before_End
   integer lchnk,ncol,indw

   type(ESMF_Time)         Date1,Date2
   type(ESMF_TimeInterval) DateDiff
   integer                 DeltaT
   real(r8)                Tscale
   real(r8)                Tfrac
   real(r8)                r_gen !WEC
   real(r8)                r_gen_switch !WEC
   real(r8)                r_year !WEC
   integer                 rc
   integer                 r_gen_floor
   integer                 nn
   integer                 kk
   real(r8)                Sbar,Qbar,Wsum
   integer                 dtime
   integer :: values(1:8),kwec !!WEC
   integer, dimension(:), allocatable :: seedwec !!WEC
   real(8) :: rwec
   
   
   ! Check if Stochaiing is initialized
   !---------------------------------
   if(.not.DAMLin_Initialized) then
     call endrun('damlining_timestep_init:: damlining NOT Initialized')
   endif

   ! Get time step size
   !--------------------
   dtime = get_step_size()
   YMD=(Year*10000) + (Month*100) + Day
   !write(iulog,*) 'YMD Flat ...',YMD,Sec

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(DAMLin_Beg_Year*10000) + (DAMLin_Beg_Month*100) + DAMLin_Beg_Day
   !write(iulog,*) 'YMD1.1 ...',YMD1
   call timemgr_time_ge(YMD1,DAMLin_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(DAMLin_End_Year*10000) + (DAMLin_End_Month*100) + DAMLin_End_Day
   !write(iulog,*) 'YMD1.2 ...',YMD1
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,DAMLin_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)
                        
  if((Before_End).and.(Update_Model)) then
   write(iulog,*) 'No happenings in here'
  endif
  
  
  YMD1=(DAMLin_Next_Year*10000) + (DAMLin_Next_Month*100) + DAMLin_Next_Day
   !write(iulog,*) 'YMD1.5 ...',YMD1,DAMLin_Next_Sec
  call timemgr_time_ge(YMD1,DAMLin_Next_Sec,            &
                        YMD ,Sec           ,Update_DAMLin)
                        
  if((Before_End).and.(Update_DAMLin)) then
  
     DAMLin_Curr_Year =DAMLin_Next_Year
     DAMLin_Curr_Month=DAMLin_Next_Month
     DAMLin_Curr_Day  =DAMLin_Next_Day
     DAMLin_Curr_Sec  =DAMLin_Next_Sec
     YMD1=(DAMLin_Curr_Year*10000) + (DAMLin_Curr_Month*100) + DAMLin_Curr_Day
     
     !write(iulog,*) 'YMD1.7 ...',YMD1
     !write(iulog,*) 'YMD2.2 pre ...',YMD2,DAMLin_Next_Sec
     call timemgr_time_inc(YMD1,DAMLin_Curr_Sec,              &
                           YMD2,DAMLin_Next_Sec,DAMLin_Step,0,0)
     !write(iulog,*) 'YMD2.2 post ...',YMD2,DAMLin_Next_Sec
     DAMLin_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(DAMLin_Next_Year*10000)
     DAMLin_Next_Month=(YMD2/100)
     DAMLin_Next_Day  = YMD2-(DAMLin_Next_Month*100)
         
     DAMLin_ON=.true.
  
  else
    if(masterproc) then
        write(iulog,*) 'DAMLin: WARNING - It is not time to DAMLin. Switching '
        write(iulog,*) 'DAMLin:           damlining OFF to coast thru the gap. '
    endif
    DAMLin_ON=.false.
  endif
  
  
  
  end subroutine
  
end module damlining