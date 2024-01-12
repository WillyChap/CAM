module stochaiing
!=====================================================================
!
! Purpose: Implement Stochaiing of the model state of U,V,T,Q, and/or PS
!          toward specified values from random analyses increments.
!
! Author: Wiliiam Chapman -- Based off the Nudging Module by Patrick C.
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    stochai the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the stochaiing
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the stochaiing module is set up to accomodate stochaiing of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The stochaiing of the model toward the analyses data is controlled by
!    the 'stochaiing_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which stochaiing is applied, the strength of the stochaiing
!    tendencies, and its spatial distribution.
!
!    FORCING:
!    --------
!    Stochaiing tendencies are applied as a relaxation force between the current
!    model state values and target state values derived from the avalilable
!    analyses. The form of the target values is selected by the 'Stochai_Force_Opt'
!    option, the timescale of the forcing is determined from the given
!    'Stochai_TimeScale_Opt', and the stochaiing strength Alpha=[0.,1.] for each
!    variable is specified by the 'Stochai_Xcoef' values. Where X={U,V,T,Q,PS}
!
!           F_stochai = Alpha*((Target-Model(t_curr))/TimeScale
!
!
!    WINDOWING:
!    ----------
!    The region of applied stochaiing can be limited using Horizontal/Vertical
!    window functions that are constructed using a parameterization of the
!    Heaviside step function.
!
!    The Heaviside window function is the product of separate horizonal and vertical
!    windows that are controled via 12 parameters:
!
!        Stochai_Hwin_lat0:     Specify the horizontal center of the window in degrees.
!        Stochai_Hwin_lon0:     The longitude must be in the range [0,360] and the
!                             latitude should be [-90,+90].
!        Stochai_Hwin_latWidth: Specify the lat and lon widths of the window as positive
!        Stochai_Hwin_lonWidth: values in degrees.Setting a width to a large value (e.g. 999)
!                             renders the window a constant in that direction.
!        Stochai_Hwin_latDelta: Controls the sharpness of the window transition with a
!        Stochai_Hwin_lonDelta: length in degrees. Small non-zero values yeild a step
!                             function while a large value yeilds a smoother transition.
!        Stochai_Hwin_Invert  : A logical flag used to invert the horizontal window function
!                             to get its compliment.(e.g. to stochai outside a given window).
!
!        Stochai_Vwin_Lindex:   In the vertical, the window is specified in terms of model
!        Stochai_Vwin_Ldelta:   level indcies. The High and Low transition levels should
!        Stochai_Vwin_Hindex:   range from [0,(NLEV+1)]. The transition lengths are also
!        Stochai_Vwin_Hdelta:   specified in terms of model indices. For a window function
!                             constant in the vertical, the Low index should be set to 0,
!                             the High index should be set to (NLEV+1), and the transition
!                             lengths should be set to 0.001
!        Stochai_Vwin_Invert  : A logical flag used to invert the vertical window function
!                             to get its compliment.
!
!        EXAMPLE: For a channel window function centered at the equator and independent
!                 of the vertical (30 levels):
!                        Stochai_Hwin_lat0     = 0.         Stochai_Vwin_Lindex = 0.
!                        Stochai_Hwin_latWidth = 30.        Stochai_Vwin_Ldelta = 0.001
!                        Stochai_Hwin_latDelta = 5.0        Stochai_Vwin_Hindex = 31.
!                        Stochai_Hwin_lon0     = 180.       Stochai_Vwin_Hdelta = 0.001
!                        Stochai_Hwin_lonWidth = 999.       Stochai_Vwin_Invert = .false.
!                        Stochai_Hwin_lonDelta = 1.0
!                        Stochai_Hwin_Invert   = .false.
!
!                 If on the other hand one wanted to apply stochaiing at the poles and
!                 not at the equator, the settings would be similar but with:
!                        Stochai_Hwin_Invert = .true.
!
!    A user can preview the window resulting from a given set of namelist values before
!    running the model. Lookat_StochaiWindow.ncl is a script avalable in the tools directory
!    which will read in the values for a given namelist and display the resulting window.
!
!    The module is currently configured for only 1 window function. It can readily be
!    extended for multiple windows if the need arises.
!
!
! Input/Output Values:
!    Forcing contributions are available for history file output by
!    the names:    {'Stochai_U','Stochai_V','Stochai_T',and 'Stochai_Q'}
!    The target values that the model state is stochaid toward are available for history
!    file output via the variables:  {'Target_stoch_U','Target_stoch_V','Target_stoch_T',and 'Target_stoch_Q'}
!
!    &stochaiing_nl
!      Stochai_Model         - LOGICAL toggle to activate stochaiing.
!                              TRUE  -> Stochaiing is on.
!                              FALSE -> Stochaiing is off.                            [DEFAULT]
!
!      Stochai_Path          - CHAR path to the analyses files.
!                              (e.g. '/glade/scratch/USER/inputdata/stochaiing/ERAI-Data/')
!
!      Stochai_File_Template - CHAR Analyses filename with year, month, day, and second
!                                 values replaced by %y, %m, %d, and %s respectively.
!                              (e.g. '%y/ERAI_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc')
!
!      Stochai_Times_Per_Day - INT Number of analyses files available per day.
!                              1 --> daily analyses.
!                              4 --> 6 hourly analyses.
!                              8 --> 3 hourly.
!
!      Model_Times_Per_Day_Stochai - INT Number of times to update the model state (used for stochaiing)
!                                each day. The value is restricted to be longer than the
!                                current model timestep and shorter than the analyses
!                                timestep. As this number is increased, the stochaiing
!                                force has the form of newtonian cooling.
!                              48 --> 1800 Second timestep.
!                              96 -->  900 Second timestep.
!
!      Stochai_Beg_Year      - INT stochaiing begining year.  [1979- ]
!      Stochai_Beg_Month     - INT stochaiing begining month. [1-12]
!      Stochai_Beg_Day       - INT stochaiing begining day.   [1-31]
!      Stochai_End_Year      - INT stochaiing ending year.    [1979-]
!      Stochai_End_Month     - INT stochaiing ending month.   [1-12]
!      Stochai_End_Day       - INT stochaiing ending day.     [1-31]
!
!      Stochai_Force_Opt     - INT Index to select the stochaiing Target for a relaxation
!                                forcing of the form:
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -> NEXT-OBS: Target=Anal(t'_next)                 [DEFAULT]
!                              1 -> LINEAR:   Target=(F*Anal(t'_curr) +(1-F)*Anal(t'_next))
!                                                 F =(t'_next - t_curr )/Tdlt_Anal
!
!      Stochai_TimeScale_Opt - INT Index to select the timescale for stochaiing.
!                                where (t'==Analysis times ; t==Model Times)
!
!                              0 -->  TimeScale = 1/Tdlt_Anal                      [DEFAULT]
!                              1 -->  TimeScale = 1/(t'_next - t_curr )
!
!      Stochai_Uprof         - INT index of profile structure to use for U.  [0,1,2]
!      Stochai_Vprof         - INT index of profile structure to use for V.  [0,1,2]
!      Stochai_Tprof         - INT index of profile structure to use for T.  [0,1,2]
!      Stochai_Qprof         - INT index of profile structure to use for Q.  [0,1,2]
!      Stochai_PSprof        - INT index of profile structure to use for PS. [0,N/A]
!
!                                The spatial distribution is specified with a profile index.
!                                 Where:  0 == OFF      (No Stochaiing of this variable)
!                                         1 == CONSTANT (Spatially Uniform Stochaiing)
!                                         2 == HEAVISIDE WINDOW FUNCTION
!
!      Stochai_Ucoef         - REAL fractional stochaiing coeffcient for U.
!      Stochai_Vcoef         - REAL fractional stochaiing coeffcient for V.
!      Stochai_Tcoef         - REAL fractional stochaiing coeffcient for T.
!      Stochai_Qcoef         - REAL fractional stochaiing coeffcient for Q.
!      Stochai_PScoef        - REAL fractional stochaiing coeffcient for PS.
!
!                                 The strength of the stochaiing is specified as a fractional
!                                 coeffcient between [0,1].
!
!      Stochai_Hwin_lat0     - REAL latitudinal center of window in degrees.
!      Stochai_Hwin_lon0     - REAL longitudinal center of window in degrees.
!      Stochai_Hwin_latWidth - REAL latitudinal width of window in degrees.
!      Stochai_Hwin_lonWidth - REAL longitudinal width of window in degrees.
!      Stochai_Hwin_latDelta - REAL latitudinal transition length of window in degrees.
!      Stochai_Hwin_lonDelta - REAL longitudinal transition length of window in degrees.
!      Stochai_Hwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!      Stochai_Vwin_Lindex   - REAL LO model index of transition
!      Stochai_Vwin_Hindex   - REAL HI model index of transition
!      Stochai_Vwin_Ldelta   - REAL LO transition length
!      Stochai_Vwin_Hdelta   - REAL HI transition length
!      Stochai_Vwin_Invert   - LOGICAL FALSE= value=1 inside the specified window, 0 outside
!                                    TRUE = value=0 inside the specified window, 1 outside
!    /
!
!================
!
! TO DO:
! -----------
!    ** Implement Ps Stochaiing????
!
!=====================================================================
  ! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_step_size
  use phys_grid   ,   only:scatter_field_to_chunk
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  ! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: Stochai_Model,Stochai_ON
  public:: stochaiing_readnl
  public:: stochaiing_init
  public:: stochaiing_timestep_init
  public:: stochaiing_timestep_tend
  private::stochaiing_update_analyses_se
  private::stochaiing_update_analyses_eul
  private::stochaiing_update_analyses_fv
  private::stochaiing_set_PSprofile
  private::stochaiing_set_profile
  private::calc_DryStaticEnergy

  ! Stochaiing Parameters
  !--------------------
  logical          :: Stochai_Model       =.false.
  logical          :: Stochai_ON          =.false.
  logical          :: Stochai_Initialized =.false.
  character(len=cl):: Stochai_Path
  character(len=cs):: Stochai_File,Stochai_File_Template
  character(len=cs):: Ran_File
  integer          :: Stochai_Force_Opt
  integer          :: Stochai_TimeScale_Opt
  integer          :: Stochai_TSmode
  integer          :: Stochai_Times_Per_Day
  integer          :: Model_Times_Per_Day_Stochai
  real(r8)         :: Stochai_Ucoef,Stochai_Vcoef
  integer          :: Stochai_Uprof,Stochai_Vprof
  real(r8)         :: Stochai_Qcoef,Stochai_Tcoef
  integer          :: Stochai_Qprof,Stochai_Tprof
  real(r8)         :: Stochai_PScoef
  integer          :: Stochai_PSprof
  integer          :: Stochai_Beg_Year ,Stochai_Beg_Month
  integer          :: Stochai_Beg_Day  ,Stochai_Beg_Sec
  integer          :: Stochai_End_Year ,Stochai_End_Month
  integer          :: Stochai_End_Day  ,Stochai_End_Sec
  integer          :: Stochai_Curr_Year,Stochai_Curr_Month
  integer          :: Stochai_Curr_Day ,Stochai_Curr_Sec
  integer          :: Stochai_Next_Year,Stochai_Next_Month
  integer          :: Stochai_Next_Day ,Stochai_Next_Sec
  integer          :: Random_Next_Year,Random_Next_Month !++WEC
  integer          :: Random_Next_Day ,Random_Next_Sec !++WEC
  integer          :: Memory_Rand !++WEC
  integer          :: Stochai_Step
  integer          :: Stochai_Climo_Year  = 1000
  integer          :: Model_Curr_Year,Model_Curr_Month
  integer          :: Model_Curr_Day ,Model_Curr_Sec
  integer          :: Model_Next_Year,Model_Next_Month
  integer          :: Model_Next_Day ,Model_Next_Sec
  integer          :: Model_Step
  real(r8)         :: Stochai_Hwin_lat0
  real(r8)         :: Stochai_Hwin_latWidth
  real(r8)         :: Stochai_Hwin_latDelta
  real(r8)         :: Stochai_Hwin_lon0
  real(r8)         :: Stochai_Hwin_lonWidth
  real(r8)         :: Stochai_Hwin_lonDelta
  logical          :: Stochai_Hwin_Invert = .false.
  real(r8)         :: Stochai_Hwin_lo
  real(r8)         :: Stochai_Hwin_hi
  real(r8)         :: Stochai_Vwin_Hindex
  real(r8)         :: Stochai_Vwin_Hdelta
  real(r8)         :: Stochai_Vwin_Lindex
  real(r8)         :: Stochai_Vwin_Ldelta
  logical          :: Stochai_Vwin_Invert =.false.
  real(r8)         :: Stochai_Vwin_lo
  real(r8)         :: Stochai_Vwin_hi
  real(r8)         :: Stochai_Hwin_latWidthH
  real(r8)         :: Stochai_Hwin_lonWidthH
  real(r8)         :: Stochai_Hwin_max
  real(r8)         :: Stochai_Hwin_min
  integer          :: Stochai_Do !adding control [WEC]

  ! Stochaiing State Arrays
  !-----------------------
  integer Stochai_nlon,Stochai_nlat,Stochai_ncol,Stochai_nlev
  real(r8),allocatable::Target_stoch_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_stoch_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_stoch_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_stoch_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_stoch_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable::Target_stoch_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Model_U     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_V     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_T     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_S     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_Q     (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Model_PS    (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Utau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Vtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Stau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Qtau  (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_PStau (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: Stochai_PSstep(:,:)    !(pcols,begchunk:endchunk)

  ! Stochaiing Observation Arrays
  !-----------------------------
  integer               Stochai_NumObs
  integer,allocatable:: Stochai_ObsInd(:)
  logical ,allocatable::Stochai_File_Present(:)
  real(r8),allocatable::Nobs_U (:,:,:,:) !(pcols,pver,begchunk:endchunk,Stochai_NumObs)
  real(r8),allocatable::Nobs_V (:,:,:,:) !(pcols,pver,begchunk:endchunk,Stochai_NumObs)
  real(r8),allocatable::Nobs_T (:,:,:,:) !(pcols,pver,begchunk:endchunk,Stochai_NumObs)
  real(r8),allocatable::Nobs_Q (:,:,:,:) !(pcols,pver,begchunk:endchunk,Stochai_NumObs)
  real(r8),allocatable::Nobs_PS(:,:,:)   !(pcols,begchunk:endchunk,Stochai_NumObs)

contains
  !================================================================
  subroutine stochaiing_readnl(nlfile)
   !
   ! STOCHAIING_READNL: Initialize default values controlling the Stochaiing
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

   namelist /stochaiing_nl/ Stochai_Model,Stochai_Path,                       &
                         Stochai_File_Template,Stochai_Force_Opt,          &
                         Stochai_TimeScale_Opt,                          &
                         Stochai_Times_Per_Day,Model_Times_Per_Day_Stochai,      &
                         Stochai_Ucoef ,Stochai_Uprof,                     &
                         Stochai_Vcoef ,Stochai_Vprof,                     &
                         Stochai_Qcoef ,Stochai_Qprof,                     &
                         Stochai_Tcoef ,Stochai_Tprof,                     &
                         Stochai_PScoef,Stochai_PSprof,                    &
                         Stochai_Beg_Year,Stochai_Beg_Month,Stochai_Beg_Day, &
                         Stochai_End_Year,Stochai_End_Month,Stochai_End_Day, &
                         Stochai_Hwin_lat0,Stochai_Hwin_lon0,              &
                         Stochai_Hwin_latWidth,Stochai_Hwin_lonWidth,      &
                         Stochai_Hwin_latDelta,Stochai_Hwin_lonDelta,      &
                         Stochai_Hwin_Invert,                            &
                         Stochai_Vwin_Lindex,Stochai_Vwin_Hindex,          &
                         Stochai_Vwin_Ldelta,Stochai_Vwin_Hdelta,          &
                         Stochai_Vwin_Invert

   ! Stochaiing is NOT initialized yet, For now
   ! Stochaiing will always begin/end at midnight.
   !--------------------------------------------
   Stochai_Initialized =.false.
   Stochai_ON          =.false.
   Stochai_Beg_Sec=0
   Stochai_End_Sec=0

   ! Set Default Namelist values
   !-----------------------------
   Stochai_Model         = .false.
   Stochai_Path          = './Data/YOTC_ne30np4_001/'
   Stochai_File_Template = 'YOTC_ne30np4_L30.cam2.i.%y-%m-%d-%s.nc'
   Stochai_Force_Opt     = 0
   Stochai_TimeScale_Opt = 0
   Stochai_TSmode        = 0
   Stochai_Times_Per_Day = 4
   Model_Times_Per_Day_Stochai = 4
   Stochai_Ucoef         = 0._r8
   Stochai_Vcoef         = 0._r8
   Stochai_Qcoef         = 0._r8
   Stochai_Tcoef         = 0._r8
   Stochai_PScoef        = 0._r8
   Stochai_Uprof         = 0
   Stochai_Vprof         = 0
   Stochai_Qprof         = 0
   Stochai_Tprof         = 0
   Stochai_PSprof        = 0
   Stochai_Beg_Year      = 2008
   Stochai_Beg_Month     = 5
   Stochai_Beg_Day       = 1
   Stochai_End_Year      = 2008
   Stochai_End_Month     = 9
   Stochai_End_Day       = 1
   Stochai_Hwin_lat0     = 0._r8
   Stochai_Hwin_latWidth = 9999._r8
   Stochai_Hwin_latDelta = 1.0_r8
   Stochai_Hwin_lon0     = 180._r8
   Stochai_Hwin_lonWidth = 9999._r8
   Stochai_Hwin_lonDelta = 1.0_r8
   Stochai_Hwin_Invert   = .false.
   Stochai_Hwin_lo       = 0.0_r8
   Stochai_Hwin_hi       = 1.0_r8
   Stochai_Vwin_Hindex   = float(pver+1)
   Stochai_Vwin_Hdelta   = 0.001_r8
   Stochai_Vwin_Lindex   = 0.0_r8
   Stochai_Vwin_Ldelta   = 0.001_r8
   Stochai_Vwin_Invert   = .false.
   Stochai_Vwin_lo       = 0.0_r8
   Stochai_Vwin_hi       = 1.0_r8
   Stochai_Do            = 1 !adding control to start ML [WEC]
   Memory_Rand         = 5 !++WEC

   ! Read in namelist values
   !------------------------
   if(masterproc) then
     unitn = getunit()
     open(unitn,file=trim(nlfile),status='old')
     call find_group_name(unitn,'stochaiing_nl',status=ierr)
     if(ierr.eq.0) then
       read(unitn,stochaiing_nl,iostat=ierr)
       if(ierr.ne.0) then
         call endrun('stochaiing_readnl:: ERROR reading namelist')
       endif
     endif
     close(unitn)
     call freeunit(unitn)
   endif

   ! Set hi/lo values according to the given '_Invert' parameters
   !--------------------------------------------------------------
   if(Stochai_Hwin_Invert) then
     Stochai_Hwin_lo = 1.0_r8
     Stochai_Hwin_hi = 0.0_r8
   else
     Stochai_Hwin_lo = 0.0_r8
     Stochai_Hwin_hi = 1.0_r8
   endif

   if(Stochai_Vwin_Invert) then
     Stochai_Vwin_lo = 1.0_r8
     Stochai_Vwin_hi = 0.0_r8
   else
     Stochai_Vwin_lo = 0.0_r8
     Stochai_Vwin_hi = 1.0_r8
   endif

   ! Check for valid namelist values
   !----------------------------------
   if((Stochai_Hwin_lat0.lt.-90._r8).or.(Stochai_Hwin_lat0.gt.+90._r8)) then
     write(iulog,*) 'STOCHAIING: Window lat0 must be in [-90,+90]'
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_lat0=',Stochai_Hwin_lat0
     call endrun('stochaiing_readnl:: ERROR in namelist')
   endif

   if((Stochai_Hwin_lon0.lt.0._r8).or.(Stochai_Hwin_lon0.ge.360._r8)) then
     write(iulog,*) 'STOCHAIING: Window lon0 must be in [0,+360)'
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_lon0=',Stochai_Hwin_lon0
     call endrun('stochaiing_readnl:: ERROR in namelist')
   endif

   if((Stochai_Vwin_Lindex.gt.Stochai_Vwin_Hindex)                         .or. &
      (Stochai_Vwin_Hindex.gt.float(pver+1)).or.(Stochai_Vwin_Hindex.lt.0._r8).or. &
      (Stochai_Vwin_Lindex.gt.float(pver+1)).or.(Stochai_Vwin_Lindex.lt.0._r8)   ) then
     write(iulog,*) 'STOCHAIING: Window Lindex must be in [0,pver+1]'
     write(iulog,*) 'STOCHAIING: Window Hindex must be in [0,pver+1]'
     write(iulog,*) 'STOCHAIING: Lindex must be LE than Hindex'
     write(iulog,*) 'STOCHAIING:  Stochai_Vwin_Lindex=',Stochai_Vwin_Lindex
     write(iulog,*) 'STOCHAIING:  Stochai_Vwin_Hindex=',Stochai_Vwin_Hindex
     call endrun('stochaiing_readnl:: ERROR in namelist')
   endif

   if((Stochai_Hwin_latDelta.le.0._r8).or.(Stochai_Hwin_lonDelta.le.0._r8).or. &
      (Stochai_Vwin_Hdelta  .le.0._r8).or.(Stochai_Vwin_Ldelta  .le.0._r8)    ) then
     write(iulog,*) 'STOCHAIING: Window Deltas must be positive'
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_latDelta=',Stochai_Hwin_latDelta
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_lonDelta=',Stochai_Hwin_lonDelta
     write(iulog,*) 'STOCHAIING:  Stochai_Vwin_Hdelta=',Stochai_Vwin_Hdelta
     write(iulog,*) 'STOCHAIING:  Stochai_Vwin_Ldelta=',Stochai_Vwin_Ldelta
     call endrun('stochaiing_readnl:: ERROR in namelist')

   endif

   if((Stochai_Hwin_latWidth.le.0._r8).or.(Stochai_Hwin_lonWidth.le.0._r8)) then
     write(iulog,*) 'STOCHAIING: Window widths must be positive'
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_latWidth=',Stochai_Hwin_latWidth
     write(iulog,*) 'STOCHAIING:  Stochai_Hwin_lonWidth=',Stochai_Hwin_lonWidth
     call endrun('stochaiing_readnl:: ERROR in namelist')
   endif

   ! Broadcast namelist variables
   !------------------------------
#ifdef SPMD
   call mpibcast(Stochai_Path         ,len(Stochai_Path)         ,mpichar,0,mpicom)
   call mpibcast(Stochai_File_Template,len(Stochai_File_Template),mpichar,0,mpicom)
   call mpibcast(Stochai_Model        , 1, mpilog, 0, mpicom)
   call mpibcast(Stochai_Initialized  , 1, mpilog, 0, mpicom)
   call mpibcast(Stochai_ON           , 1, mpilog, 0, mpicom)
   call mpibcast(Stochai_Force_Opt    , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_TimeScale_Opt, 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_TSmode       , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Times_Per_Day, 1, mpiint, 0, mpicom)
   call mpibcast(Model_Times_Per_Day_Stochai, 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Ucoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Tcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Qcoef        , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_PScoef       , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Uprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Vprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Tprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Qprof        , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_PSprof       , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Beg_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Beg_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Beg_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Beg_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_End_Year     , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_End_Month    , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_End_Day      , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_End_Sec      , 1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Hwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_lat0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_latWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_latDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_lon0    , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_lonWidth, 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_lonDelta, 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_Invert,   1, mpilog, 0, mpicom)
   call mpibcast(Stochai_Vwin_lo      , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_hi      , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_Hindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_Hdelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_Lindex  , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_Ldelta  , 1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Vwin_Invert,   1, mpilog, 0, mpicom)
#endif

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_readnl
  !================================================================


  !================================================================
  subroutine stochaiing_init
   !
   ! STOCHAIING_INIT: Allocate space and initialize Stochaiing values
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

   ! Get the time step size
   !------------------------
   dtime = get_step_size()

   ! Allocate Space for Stochaiing data arrays
   !-----------------------------------------
   allocate(Target_stoch_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_stoch_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_stoch_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_stoch_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_stoch_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Target_stoch_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Target_stoch_PS',pcols*((endchunk-begchunk)+1))

   allocate(Model_U(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_U',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_V(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_V',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_T(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_T',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_S(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_S',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_Q',pcols*pver*((endchunk-begchunk)+1))
   allocate(Model_PS(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Model_PS',pcols*((endchunk-begchunk)+1))

   ! Allocate Space for spatial dependence of
   ! Stochaiing Coefs and Stochaiing Forcing.
   !-------------------------------------------
   allocate(Stochai_Utau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Utau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Vtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Vtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Stau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Stau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Qtau(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Qtau',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_PStau(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_PStau',pcols*((endchunk-begchunk)+1))

   allocate(Stochai_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Ustep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Vstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Sstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_Qstep',pcols*pver*((endchunk-begchunk)+1))
   allocate(Stochai_PSstep(pcols,begchunk:endchunk),stat=istat)
   call alloc_err(istat,'stochaiing_init','Stochai_PSstep',pcols*((endchunk-begchunk)+1))

   ! Register output fields with the cam history module
   !-----------------------------------------------------
   call addfld( 'Stochai_U',(/ 'lev' /),'A','m/s/s'  ,'U Stochaiing Tendency')
   call addfld( 'Stochai_V',(/ 'lev' /),'A','m/s/s'  ,'V Stochaiing Tendency')
   call addfld( 'Stochai_T',(/ 'lev' /),'A','K/s'    ,'T Stochaiing Tendency')
   call addfld( 'Stochai_Q',(/ 'lev' /),'A','kg/kg/s','Q Stochaiing Tendency')

   if(Stochai_Do==1) then
      call addfld('Target_stoch_U',(/ 'lev' /),'A','m/s/s'    ,'U Stochaiing Target'  ) !need to edit units if we stochai ML 'm/s/s' [WEC]
      call addfld('Target_stoch_V',(/ 'lev' /),'A','m/s/s'    ,'V Stochaiing Target'  ) !need to edit units if we stochai ML 'm/s/s' [WEC]
      call addfld('Target_stoch_T',(/ 'lev' /),'A','K/s'      ,'T Stochaiing Target'  ) !need to edit units if we stochai ML 'K/s' [WEC
      call addfld('Target_stoch_Q',(/ 'lev' /),'A','kg/kg/s'  ,'Q Stochaiing Target  ') !need to edit units if we stochai ML 'kg/kg/s' [WEC]
   else
      call addfld('Target_stoch_U',(/ 'lev' /),'A','m/s'    ,'U Stochaiing Target'  )
      call addfld('Target_stoch_V',(/ 'lev' /),'A','m/s'    ,'V Stochaiing Target'  )
      call addfld('Target_stoch_T',(/ 'lev' /),'A','K'      ,'T Stochaiing Target'  )
      call addfld('Target_stoch_Q',(/ 'lev' /),'A','kg/kg'  ,'Q Stochaiing Target  ')
   endif
   !-----------------------------------------
   ! Values initialized only by masterproc
   !-----------------------------------------
   if(masterproc) then

     ! Set the Stepping intervals for Model and Stochaiing values
     ! Ensure that the Model_Step is not smaller then one timestep
     !  and not larger then the Stochai_Step.
     !--------------------------------------------------------
     Model_Step=86400/Model_Times_Per_Day_Stochai
     Stochai_Step=86400/Stochai_Times_Per_Day
     if(Model_Step.lt.dtime) then
       write(iulog,*) ' '
       write(iulog,*) 'STOCHAIING: Model_Step cannot be less than a model timestep'
       write(iulog,*) 'STOCHAIING:  Setting Model_Step=dtime , dtime=',dtime
       write(iulog,*) ' '
       Model_Step=dtime
     endif
     if(Model_Step.gt.Stochai_Step) then
       write(iulog,*) ' '
       write(iulog,*) 'STOCHAIING: Model_Step cannot be more than Stochai_Step'
       write(iulog,*) 'STOCHAIING:  Setting Model_Step=Stochai_Step, Stochai_Step=',Stochai_Step
       write(iulog,*) ' '
       Model_Step=Stochai_Step
     endif

     ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     Stochai_nlon=hdim1_d
     Stochai_nlat=hdim2_d
     Stochai_ncol=hdim1_d*hdim2_d
     Stochai_nlev=pver

     ! Check the time relative to the stochaiing window
     !------------------------------------------------
     call get_curr_date(Year,Month,Day,Sec)
     YMD=(Year*10000) + (Month*100) + Day
     YMD1=(Stochai_Beg_Year*10000) + (Stochai_Beg_Month*100) + Stochai_Beg_Day
     call timemgr_time_ge(YMD1,Stochai_Beg_Sec,         &
                          YMD ,Sec          ,After_Beg)
     YMD1=(Stochai_End_Year*10000) + (Stochai_End_Month*100) + Stochai_End_Day
     call timemgr_time_ge(YMD ,Sec          ,          &
                          YMD1,Stochai_End_Sec,Before_End)

     if((After_Beg).and.(Before_End)) then
       ! Set Time indicies so that the next call to
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Year
       Model_Next_Month=Month
       Model_Next_Day  =Day
       Model_Next_Sec  =(Sec/Model_Step)*Model_Step
       Stochai_Next_Year =Year
       Stochai_Next_Month=Month
       Stochai_Next_Day  =Day
       Stochai_Next_Sec  =(Sec/Stochai_Step)*Stochai_Step

       !++WEC
       Random_Next_Year =Year
       Random_Next_Month=Month
       Random_Next_Day  =Day
       Random_Next_Sec  =(Sec/Stochai_Step)*Stochai_Step
       !--WEC
     elseif(.not.After_Beg) then
       ! Set Time indicies to Stochaiing start,
       ! timestep_init will initialize the data arrays.
       !--------------------------------------------
       Model_Next_Year =Stochai_Beg_Year
       Model_Next_Month=Stochai_Beg_Month
       Model_Next_Day  =Stochai_Beg_Day
       Model_Next_Sec  =Stochai_Beg_Sec
       Stochai_Next_Year =Stochai_Beg_Year
       Stochai_Next_Month=Stochai_Beg_Month
       Stochai_Next_Day  =Stochai_Beg_Day
       Stochai_Next_Sec  =Stochai_Beg_Sec

       !++WEC
       Random_Next_Year =Stochai_Beg_Year
       Random_Next_Month=Stochai_Beg_Month
       Random_Next_Day  =Stochai_Beg_Day
       Random_Next_Sec  =Stochai_Beg_Sec
       !--WEC
     elseif(.not.Before_End) then
       ! Stochaiing will never occur, so switch it off
       !--------------------------------------------
       Stochai_Model=.false.
       Stochai_ON   =.false.
       write(iulog,*) ' '
       write(iulog,*) 'STOCHAIING: WARNING - Stochaiing has been requested by it will'
       write(iulog,*) 'STOCHAIING:           never occur for the given time values'
       write(iulog,*) ' '
     endif

     ! Initialize values for window function
     !----------------------------------------
     lonp= 180._r8
     lon0=   0._r8
     lonn=-180._r8
     latp=  90._r8-Stochai_Hwin_lat0
     lat0=   0._r8
     latn= -90._r8-Stochai_Hwin_lat0

     Stochai_Hwin_lonWidthH=Stochai_Hwin_lonWidth/2._r8
     Stochai_Hwin_latWidthH=Stochai_Hwin_latWidth/2._r8

     Val1_p=(1._r8+tanh((Stochai_Hwin_lonWidthH+lonp)/Stochai_Hwin_lonDelta))/2._r8
     Val2_p=(1._r8+tanh((Stochai_Hwin_lonWidthH-lonp)/Stochai_Hwin_lonDelta))/2._r8
     Val3_p=(1._r8+tanh((Stochai_Hwin_latWidthH+latp)/Stochai_Hwin_latDelta))/2._r8
     Val4_p=(1._r8+tanh((Stochai_Hwin_latWidthH-latp)/Stochai_Hwin_latDelta))/2_r8
     Val1_0=(1._r8+tanh((Stochai_Hwin_lonWidthH+lon0)/Stochai_Hwin_lonDelta))/2._r8
     Val2_0=(1._r8+tanh((Stochai_Hwin_lonWidthH-lon0)/Stochai_Hwin_lonDelta))/2._r8
     Val3_0=(1._r8+tanh((Stochai_Hwin_latWidthH+lat0)/Stochai_Hwin_latDelta))/2._r8
     Val4_0=(1._r8+tanh((Stochai_Hwin_latWidthH-lat0)/Stochai_Hwin_latDelta))/2._r8

     Val1_n=(1._r8+tanh((Stochai_Hwin_lonWidthH+lonn)/Stochai_Hwin_lonDelta))/2._r8
     Val2_n=(1._r8+tanh((Stochai_Hwin_lonWidthH-lonn)/Stochai_Hwin_lonDelta))/2._r8
     Val3_n=(1._r8+tanh((Stochai_Hwin_latWidthH+latn)/Stochai_Hwin_latDelta))/2._r8
     Val4_n=(1._r8+tanh((Stochai_Hwin_latWidthH-latn)/Stochai_Hwin_latDelta))/2._r8

     Stochai_Hwin_max=     Val1_0*Val2_0*Val3_0*Val4_0
     Stochai_Hwin_min=min((Val1_p*Val2_p*Val3_n*Val4_n), &
                        (Val1_p*Val2_p*Val3_p*Val4_p), &
                        (Val1_n*Val2_n*Val3_n*Val4_n), &
                        (Val1_n*Val2_n*Val3_p*Val4_p))

     ! Initialize number of stochaiing observation values to keep track of.
     ! Allocate and initialize observation indices
     !-----------------------------------------------------------------
     if((Stochai_Force_Opt.ge.0).and.(Stochai_Force_Opt.le.1)) then
       Stochai_NumObs=2
     else
       ! Additional Options may need OBS values at more times.
       !------------------------------------------------------
       Stochai_NumObs=2
       write(iulog,*) 'STOCHAIING: Setting Stochai_NumObs=2'
       write(iulog,*) 'STOCHAIING: WARNING: Unknown Stochai_Force_Opt=',Stochai_Force_Opt
       call endrun('STOCHAIING: Unknown Forcing Option')
     endif
     allocate(Stochai_ObsInd(Stochai_NumObs),stat=istat)
     call alloc_err(istat,'stochaiing_init','Stochai_ObsInd',Stochai_NumObs)
     allocate(Stochai_File_Present(Stochai_NumObs),stat=istat)
     call alloc_err(istat,'stochaiing_init','Stochai_File_Present',Stochai_NumObs)
     do nn=1,Stochai_NumObs
       Stochai_ObsInd(nn) = Stochai_NumObs+1-nn
     end do
     Stochai_File_Present(:)=.false.

     ! Initialization is done,
     !--------------------------
     Stochai_Initialized=.true.

     ! Check that this is a valid DYCORE model
     !------------------------------------------
     if((.not.dycore_is('UNSTRUCTURED')).and. &
        (.not.dycore_is('EUL')         ).and. &
        (.not.dycore_is('LR')          )      ) then
       call endrun('STOCHAIING IS CURRENTLY ONLY CONFIGURED FOR CAM-SE, FV, or EUL')
     endif

     ! Informational Output
     !---------------------------
     write(iulog,*) ' '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) '  MODEL STOCHAIING INITIALIZED WITH THE FOLLOWING SETTINGS: '
     write(iulog,*) '---------------------------------------------------------'
     write(iulog,*) 'STOCHAIING: Stochai_Model=',Stochai_Model
     write(iulog,*) 'STOCHAIING: Stochai_Path=',Stochai_Path
     write(iulog,*) 'STOCHAIING: Stochai_File_Template =',Stochai_File_Template
     write(iulog,*) 'STOCHAIING: Stochai_Force_Opt=',Stochai_Force_Opt
     write(iulog,*) 'STOCHAIING: Stochai_TimeScale_Opt=',Stochai_TimeScale_Opt
     write(iulog,*) 'STOCHAIING: Stochai_TSmode=',Stochai_TSmode
     write(iulog,*) 'STOCHAIING: Stochai_Times_Per_Day=',Stochai_Times_Per_Day
     write(iulog,*) 'STOCHAIING: Model_Times_Per_Day_Stochai=',Model_Times_Per_Day_Stochai
     write(iulog,*) 'STOCHAIING: Stochai_Step=',Stochai_Step
     write(iulog,*) 'STOCHAIING: Model_Step=',Model_Step
     write(iulog,*) 'STOCHAIING: Stochai_Ucoef  =',Stochai_Ucoef
     write(iulog,*) 'STOCHAIING: Stochai_Vcoef  =',Stochai_Vcoef
     write(iulog,*) 'STOCHAIING: Stochai_Qcoef  =',Stochai_Qcoef
     write(iulog,*) 'STOCHAIING: Stochai_Tcoef  =',Stochai_Tcoef
     write(iulog,*) 'STOCHAIING: Stochai_PScoef =',Stochai_PScoef
     write(iulog,*) 'STOCHAIING: Stochai_Uprof  =',Stochai_Uprof
     write(iulog,*) 'STOCHAIING: Stochai_Vprof  =',Stochai_Vprof
     write(iulog,*) 'STOCHAIING: Stochai_Qprof  =',Stochai_Qprof
     write(iulog,*) 'STOCHAIING: Stochai_Tprof  =',Stochai_Tprof
     write(iulog,*) 'STOCHAIING: Stochai_PSprof =',Stochai_PSprof
     write(iulog,*) 'STOCHAIING: Stochai_Beg_Year =',Stochai_Beg_Year
     write(iulog,*) 'STOCHAIING: Stochai_Beg_Month=',Stochai_Beg_Month
     write(iulog,*) 'STOCHAIING: Stochai_Beg_Day  =',Stochai_Beg_Day
     write(iulog,*) 'STOCHAIING: Stochai_End_Year =',Stochai_End_Year
     write(iulog,*) 'STOCHAIING: Stochai_End_Month=',Stochai_End_Month
     write(iulog,*) 'STOCHAIING: Stochai_End_Day  =',Stochai_End_Day
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lat0     =',Stochai_Hwin_lat0
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_latWidth =',Stochai_Hwin_latWidth
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_latDelta =',Stochai_Hwin_latDelta
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lon0     =',Stochai_Hwin_lon0
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lonWidth =',Stochai_Hwin_lonWidth
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lonDelta =',Stochai_Hwin_lonDelta
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_Invert   =',Stochai_Hwin_Invert
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lo       =',Stochai_Hwin_lo
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_hi       =',Stochai_Hwin_hi
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_Hindex   =',Stochai_Vwin_Hindex
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_Hdelta   =',Stochai_Vwin_Hdelta
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_Lindex   =',Stochai_Vwin_Lindex
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_Ldelta   =',Stochai_Vwin_Ldelta
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_Invert   =',Stochai_Vwin_Invert
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_lo       =',Stochai_Vwin_lo
     write(iulog,*) 'STOCHAIING: Stochai_Vwin_hi       =',Stochai_Vwin_hi
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_latWidthH=',Stochai_Hwin_latWidthH
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_lonWidthH=',Stochai_Hwin_lonWidthH
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_max      =',Stochai_Hwin_max
     write(iulog,*) 'STOCHAIING: Stochai_Hwin_min      =',Stochai_Hwin_min
     write(iulog,*) 'STOCHAIING: Stochai_Initialized   =',Stochai_Initialized
     write(iulog,*) ' '
     write(iulog,*) 'STOCHAIING: Stochai_NumObs=',Stochai_NumObs
     write(iulog,*) ' '

   endif ! (masterproc) then

   ! Broadcast other variables that have changed
   !---------------------------------------------
#ifdef SPMD
   call mpibcast(Model_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Step          ,            1, mpir8 , 0, mpicom)
   call mpibcast(Model_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Model_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Next_Year     ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Next_Month    ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Next_Day      ,            1, mpiint, 0, mpicom)
   call mpibcast(Random_Next_Sec     ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Year    ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Month   ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Random_Next_Day     ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Memory_Rand         ,            1, mpiint, 0, mpicom) !++WEC
   call mpibcast(Stochai_Next_Sec      ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Model         ,            1, mpilog, 0, mpicom)
   call mpibcast(Stochai_ON            ,            1, mpilog, 0, mpicom)
   call mpibcast(Stochai_Initialized   ,            1, mpilog, 0, mpicom)
   call mpibcast(Stochai_ncol          ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_nlev          ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_nlon          ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_nlat          ,            1, mpiint, 0, mpicom)
   call mpibcast(Stochai_Hwin_max      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_min      ,            1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_lonWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_Hwin_latWidthH,            1, mpir8 , 0, mpicom)
   call mpibcast(Stochai_NumObs        ,            1, mpiint, 0, mpicom)
#endif

   ! All non-masterproc processes also need to allocate space
   ! before the broadcast of Stochai_NumObs dependent data.
   !------------------------------------------------------------
   if(.not.masterproc) then
     allocate(Stochai_ObsInd(Stochai_NumObs),stat=istat)
     call alloc_err(istat,'stochaiing_init','Stochai_ObsInd',Stochai_NumObs)
     allocate(Stochai_File_Present(Stochai_NumObs),stat=istat)
     call alloc_err(istat,'stochaiing_init','Stochai_File_Present',Stochai_NumObs)
   endif
#ifdef SPMD
   call mpibcast(Stochai_ObsInd        , Stochai_NumObs, mpiint, 0, mpicom)
   call mpibcast(Stochai_File_Present  , Stochai_NumObs, mpilog, 0, mpicom)
#endif

   ! Allocate Space for Stochaiing observation arrays, initialize with 0's
   !---------------------------------------------------------------------
   allocate(Nobs_U(pcols,pver,begchunk:endchunk,Stochai_NumObs),stat=istat)
   call alloc_err(istat,'stochaiing_init','Nobs_U',pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)
   allocate(Nobs_V(pcols,pver,begchunk:endchunk,Stochai_NumObs),stat=istat)
   call alloc_err(istat,'stochaiing_init','Nobs_V',pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)
   allocate(Nobs_T(pcols,pver,begchunk:endchunk,Stochai_NumObs),stat=istat)
   call alloc_err(istat,'stochaiing_init','Nobs_T',pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)
   allocate(Nobs_Q(pcols,pver,begchunk:endchunk,Stochai_NumObs),stat=istat)
   call alloc_err(istat,'stochaiing_init','Nobs_Q',pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)
   allocate(Nobs_PS(pcols,begchunk:endchunk,Stochai_NumObs),stat=istat)
   call alloc_err(istat,'stochaiing_init','Nobs_PS',pcols*((endchunk-begchunk)+1)*Stochai_NumObs)

   Nobs_U(:pcols,:pver,begchunk:endchunk,:Stochai_NumObs)=0._r8
   Nobs_V(:pcols,:pver,begchunk:endchunk,:Stochai_NumObs)=0._r8
   Nobs_T(:pcols,:pver,begchunk:endchunk,:Stochai_NumObs)=0._r8
   Nobs_Q(:pcols,:pver,begchunk:endchunk,:Stochai_NumObs)=0._r8
   Nobs_PS(:pcols     ,begchunk:endchunk,:Stochai_NumObs)=0._r8

!!DIAG
   if(masterproc) then
     write(iulog,*) 'STOCHAIING: stochaiing_init() OBS arrays allocated and initialized'
     write(iulog,*) 'STOCHAIING: stochaiing_init() SIZE#',(9*pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)
     write(iulog,*) 'STOCHAIING: stochaiing_init() MB:',float(8*9*pcols*pver*((endchunk-begchunk)+1)*Stochai_NumObs)/(1024._r8*1024._r8)
     write(iulog,*) 'STOCHAIING: stochaiing_init() pcols=',pcols,' pver=',pver
     write(iulog,*) 'STOCHAIING: stochaiing_init() begchunk:',begchunk,' endchunk=',endchunk
     write(iulog,*) 'STOCHAIING: stochaiing_init() chunk:',(endchunk-begchunk+1),' Stochai_NumObs=',Stochai_NumObs
     write(iulog,*) 'STOCHAIING: stochaiing_init() Stochai_ObsInd=',Stochai_ObsInd
     write(iulog,*) 'STOCHAIING: stochaiing_init() Stochai_File_Present=',Stochai_File_Present
   endif
!!DIAG

   ! Initialize the analysis filename at the NEXT time for startup.
   !---------------------------------------------------------------
   Stochai_File=interpret_filename_spec(Stochai_File_Template      , &
                                       yr_spec=Stochai_Climo_Year , &
                                      mon_spec=Stochai_Next_Month, &
                                      day_spec=Stochai_Next_Day  , &
                                      sec_spec=Stochai_Next_Sec    )
   if(masterproc) then
    write(iulog,*) 'STOCHAIING: Reading analyses:',trim(Stochai_Path)//trim(Stochai_File)
   endif

   ! Rotate Stochai_ObsInd() indices for new data, then update
   ! the Stochai observation arrays with analysis data at the
   ! NEXT==Stochai_ObsInd(1) time.
   !----------------------------------------------------------
   if(dycore_is('UNSTRUCTURED')) then
     call stochaiing_update_analyses_se (trim(Stochai_Path)//trim(Stochai_File))
   elseif(dycore_is('EUL')) then
     call stochaiing_update_analyses_eul(trim(Stochai_Path)//trim(Stochai_File))
   else !if(dycore_is('LR')) then
     call stochaiing_update_analyses_fv (trim(Stochai_Path)//trim(Stochai_File))
   endif

   ! Initialize Stochaiing Coeffcient profiles in local arrays
   ! Load zeros into stochaiing arrays
   !------------------------------------------------------
   do lchnk=begchunk,endchunk
     ncol=get_ncols_p(lchnk)
     do icol=1,ncol
       rlat=get_rlat_p(lchnk,icol)*180._r8/SHR_CONST_PI
       rlon=get_rlon_p(lchnk,icol)*180._r8/SHR_CONST_PI

       call stochaiing_set_profile(rlat,rlon,Stochai_Uprof,Wprof,pver)
       Stochai_Utau(icol,:,lchnk)=Wprof(:)
       call stochaiing_set_profile(rlat,rlon,Stochai_Vprof,Wprof,pver)
       Stochai_Vtau(icol,:,lchnk)=Wprof(:)
       call stochaiing_set_profile(rlat,rlon,Stochai_Tprof,Wprof,pver)
       Stochai_Stau(icol,:,lchnk)=Wprof(:)
       call stochaiing_set_profile(rlat,rlon,Stochai_Qprof,Wprof,pver)
       Stochai_Qtau(icol,:,lchnk)=Wprof(:)

       Stochai_PStau(icol,lchnk)=stochaiing_set_PSprofile(rlat,rlon,Stochai_PSprof)
     end do

     if(Stochai_Do==1) then
       Stochai_Utau(:ncol,:pver,lchnk) =                             &
       Stochai_Utau(:ncol,:pver,lchnk) * Stochai_Ucoef
       Stochai_Vtau(:ncol,:pver,lchnk) =                             &
       Stochai_Vtau(:ncol,:pver,lchnk) * Stochai_Vcoef
       Stochai_Stau(:ncol,:pver,lchnk) =                             &
       Stochai_Stau(:ncol,:pver,lchnk) * Stochai_Tcoef
       Stochai_Qtau(:ncol,:pver,lchnk) =                             &
       Stochai_Qtau(:ncol,:pver,lchnk) * Stochai_Qcoef
       Stochai_PStau(:ncol,lchnk)=                             &
       Stochai_PStau(:ncol,lchnk)* Stochai_PScoef
    else
      Stochai_Utau(:ncol,:pver,lchnk) =                             &
      Stochai_Utau(:ncol,:pver,lchnk) * Stochai_Ucoef/float(Stochai_Step)
      Stochai_Vtau(:ncol,:pver,lchnk) =                             &
      Stochai_Vtau(:ncol,:pver,lchnk) * Stochai_Vcoef/float(Stochai_Step)
      Stochai_Stau(:ncol,:pver,lchnk) =                             &
      Stochai_Stau(:ncol,:pver,lchnk) * Stochai_Tcoef/float(Stochai_Step)
      Stochai_Qtau(:ncol,:pver,lchnk) =                             &
      Stochai_Qtau(:ncol,:pver,lchnk) * Stochai_Qcoef/float(Stochai_Step)
      Stochai_PStau(:ncol,lchnk)=                             &
      Stochai_PStau(:ncol,lchnk)* Stochai_PScoef/float(Stochai_Step)
    endif

     Stochai_Ustep(:pcols,:pver,lchnk)=0._r8
     Stochai_Vstep(:pcols,:pver,lchnk)=0._r8
     Stochai_Sstep(:pcols,:pver,lchnk)=0._r8
     Stochai_Qstep(:pcols,:pver,lchnk)=0._r8
     Stochai_PSstep(:pcols,lchnk)=0._r8
     Target_stoch_U(:pcols,:pver,lchnk)=0._r8
     Target_stoch_V(:pcols,:pver,lchnk)=0._r8
     Target_stoch_T(:pcols,:pver,lchnk)=0._r8
     Target_stoch_S(:pcols,:pver,lchnk)=0._r8
     Target_stoch_Q(:pcols,:pver,lchnk)=0._r8
     Target_stoch_PS(:pcols,lchnk)=0._r8
   end do

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_init
  !================================================================


  !================================================================
  subroutine stochaiing_timestep_init(phys_state)
   !
   ! STOCHAIING_TIMESTEP_INIT:
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
   logical Update_Model,Update_Stochai,Sync_Error
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
   if(.not.Stochai_Initialized) then
     call endrun('stochaiing_timestep_init:: Stochaiing NOT Initialized')
   endif

   ! Get time step size
   !--------------------
   dtime = get_step_size()

   ! Get Current time
   !--------------------
   call get_curr_date(Year,Month,Day,Sec)
   YMD=(Year*10000) + (Month*100) + Day
   !write(iulog,*) 'YMD Flat ...',YMD,Sec

   !-------------------------------------------------------
   ! Determine if the current time is AFTER the begining time
   ! and if it is BEFORE the ending time.
   !-------------------------------------------------------
   YMD1=(Stochai_Beg_Year*10000) + (Stochai_Beg_Month*100) + Stochai_Beg_Day
   !write(iulog,*) 'YMD1.1 ...',YMD1
   call timemgr_time_ge(YMD1,Stochai_Beg_Sec,         &
                        YMD ,Sec          ,After_Beg)

   YMD1=(Stochai_End_Year*10000) + (Stochai_End_Month*100) + Stochai_End_Day
   !write(iulog,*) 'YMD1.2 ...',YMD1
   call timemgr_time_ge(YMD ,Sec,                    &
                        YMD1,Stochai_End_Sec,Before_End)

   !--------------------------------------------------------------
   ! When past the NEXT time, Update Model Arrays and time indices
   !--------------------------------------------------------------
   YMD1=(Model_Next_Year*10000) + (Model_Next_Month*100) + Model_Next_Day
   call timemgr_time_ge(YMD1,Model_Next_Sec,            &
                        YMD ,Sec           ,Update_Model)

   !write(iulog,*) 'YMD1.3 ...',YMD1, Update_Model

   if((Before_End).and.(Update_Model)) then
     ! Increment the Model times by the current interval
     !---------------------------------------------------
     Model_Curr_Year =Model_Next_Year
     Model_Curr_Month=Model_Next_Month
     Model_Curr_Day  =Model_Next_Day
     Model_Curr_Sec  =Model_Next_Sec
     YMD1=(Model_Curr_Year*10000) + (Model_Curr_Month*100) + Model_Curr_Day
     !write(iulog,*) 'YMD1.4 pre ...',YMD1,Model_Next_Sec
     call timemgr_time_inc(YMD1,Model_Curr_Sec,              &
                           YMD2,Model_Next_Sec,Model_Step,0,0)
     !write(iulog,*) 'YMD1.4 post ...',YMD2,Model_Next_Sec
     ! Check for Sync Error where NEXT model time after the update
     ! is before the current time. If so, reset the next model
     ! time to a Model_Step after the current time.
     !--------------------------------------------------------------
     call timemgr_time_ge(YMD2,Model_Next_Sec,            &
                          YMD ,Sec           ,Sync_Error)

     !write(iulog,*) 'YMD2.1 ...',YMD2,Model_Next_Sec
     if(Sync_Error) then
       Model_Curr_Year =Year
       Model_Curr_Month=Month
       Model_Curr_Day  =Day
       Model_Curr_Sec  =Sec
       call timemgr_time_inc(YMD ,Model_Curr_Sec,              &
                             YMD2,Model_Next_Sec,Model_Step,0,0)
       write(iulog,*) 'STOCHAIING: WARNING - Model_Time Sync ERROR... CORRECTED'
     endif
     Model_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Model_Next_Year*10000)
     Model_Next_Month=(YMD2/100)
     Model_Next_Day  = YMD2-(Model_Next_Month*100)

     ! Load values at Current into the Model arrays
     !-----------------------------------------------
     call cnst_get_ind('Q',indw)
     do lchnk=begchunk,endchunk
       ncol=phys_state(lchnk)%ncol
       Model_U(:ncol,:pver,lchnk)=phys_state(lchnk)%u(:ncol,:pver)
       Model_V(:ncol,:pver,lchnk)=phys_state(lchnk)%v(:ncol,:pver)
       Model_T(:ncol,:pver,lchnk)=phys_state(lchnk)%t(:ncol,:pver)
       Model_Q(:ncol,:pver,lchnk)=phys_state(lchnk)%q(:ncol,:pver,indw)
       Model_PS(:ncol,lchnk)=phys_state(lchnk)%ps(:ncol)
     end do

     ! Load Dry Static Energy values for Model
     !-----------------------------------------
     if(Stochai_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Model_S(:ncol,:pver,lchnk)=cpair*Model_T(:ncol,:pver,lchnk)
       end do
     elseif(Stochai_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Model_T(:,:,lchnk)  , Model_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis,  Model_PS(:,lchnk), &
                                                  Model_S(:,:,lchnk), ncol)
       end do
     endif
   endif ! ((Before_End).and.(Update_Model)) then

   !----------------------------------------------------------------
   ! When past the NEXT time, Update Stochaiing Arrays and time indices
   !----------------------------------------------------------------
   YMD1=(Stochai_Next_Year*10000) + (Stochai_Next_Month*100) + Stochai_Next_Day
   !write(iulog,*) 'YMD1.5 ...',YMD1,Stochai_Next_Sec
   call timemgr_time_ge(YMD1,Stochai_Next_Sec,            &
                        YMD ,Sec           ,Update_Stochai)

   if((Before_End).and.(Update_Stochai)) then
     ! Increment the Stochai times by the current interval
     !---------------------------------------------------
     !write(iulog,*) 'YMD1.6 ...',YMD1

     Stochai_Curr_Year =Stochai_Next_Year
     Stochai_Curr_Month=Stochai_Next_Month
     Stochai_Curr_Day  =Stochai_Next_Day
     Stochai_Curr_Sec  =Stochai_Next_Sec
     YMD1=(Stochai_Curr_Year*10000) + (Stochai_Curr_Month*100) + Stochai_Curr_Day
     !write(iulog,*) 'YMD1.7 ...',YMD1
     !write(iulog,*) 'YMD2.2 pre ...',YMD2,Stochai_Next_Sec
     call timemgr_time_inc(YMD1,Stochai_Curr_Sec,              &
                           YMD2,Stochai_Next_Sec,Stochai_Step,0,0)
     !write(iulog,*) 'YMD2.2 post ...',YMD2,Stochai_Next_Sec
     Stochai_Next_Year =(YMD2/10000)
     YMD2            = YMD2-(Stochai_Next_Year*10000)
     Stochai_Next_Month=(YMD2/100)
     Stochai_Next_Day  = YMD2-(Stochai_Next_Month*100)

     !write(iulog,*) 'Current Stochaiing Stuff',Stochai_Next_Year,Stochai_Next_Month,Stochai_Next_Day,Stochai_Next_Sec,Stochai_Step,Model_Step
     ! Set the analysis filename at the NEXT time.
     !---------------------------------------------------------------
     Stochai_File=interpret_filename_spec(Stochai_File_Template      , &
                                         yr_spec=Stochai_Climo_Year , &
                                        mon_spec=Stochai_Next_Month, &
                                        day_spec=Stochai_Next_Day  , &
                                        sec_spec=Stochai_Next_Sec    )
     !!!++WEC Grab a random file!
     if(Memory_Rand.ge.5) then
     
         call date_and_time(values=values)
         call random_seed(size=kwec)
         allocate(seedwec(1:kwec))
         seedwec(:) = values(8)
         call random_seed(put=seedwec)
        
         call random_number(r_gen_switch)
         r_gen_switch = r_gen_switch-0.5

         if(r_gen_switch.ge.0) then
             call random_number(r_gen)
             r_gen_floor = FLOOR(12*r_gen)
             !write(iulog,*) 'Stochai Next Sec',Stochai_Next_Sec
             !write(iulog,*) 'Switch positive',r_gen_switch
             call timemgr_time_inc(YMD1,Stochai_Curr_Sec,              &
                           YMD3,Random_Next_Sec,r_gen_floor*86400,0,0)
             !write(iulog,*) 'Yogurt YMD3:',YMD3,r_gen,r_gen_floor
             Random_Next_Year =(YMD3/10000)
             YMD3            = YMD3-(Random_Next_Year*10000)
             Random_Next_Month=(YMD3/100)
             Random_Next_Day  = YMD3-(Random_Next_Month*100)


             call random_number(r_year)
             Random_Next_Year = 2011+int((2019-2011 +1)*r_year)
             !write(iulog,*) 'switch positive',Random_Next_Year,Random_Next_Month,Random_Next_Day,Random_Next_Sec
             Ran_File=interpret_filename_spec(Stochai_File_Template      , &
                                         yr_spec=Random_Next_Year , &
                                        mon_spec=Random_Next_Month, &
                                        day_spec=Random_Next_Day  , &
                                        sec_spec=Stochai_Next_Sec    )
             !write(iulog,*) 'random file positive',Ran_File

         else
             call random_number(r_gen)
             r_gen_floor = FLOOR(12*r_gen)
             !write(iulog,*) 'Stochai Next Sec',Stochai_Next_Sec
             !write(iulog,*) 'Switch negative',r_gen_switch
             call timemgr_time_inc_minus(YMD1,Stochai_Curr_Sec,              &
                           YMD3,Random_Next_Sec,r_gen_floor*86400,0,0)
             !write(iulog,*) 'Yogurt YMD3:',YMD3,r_gen,r_gen_floor
             Random_Next_Year =(YMD3/10000)
             YMD3            = YMD3-(Random_Next_Year*10000)
             Random_Next_Month=(YMD3/100)
             Random_Next_Day  = YMD3-(Random_Next_Month*100)
             call random_number(r_year)

             Random_Next_Year = 2011+int((2019-2011 +1)*r_year)
             !write(iulog,*) 'switch negative',Random_Next_Year,Random_Next_Month,Random_Next_Day,Random_Next_Sec

             Ran_File=interpret_filename_spec(Stochai_File_Template      , &
                                         yr_spec=Random_Next_Year , &
                                        mon_spec=Random_Next_Month, &
                                        day_spec=Random_Next_Day  , &
                                        sec_spec=Stochai_Next_Sec    )

             !write(iulog,*) 'random file negative',Ran_File
         endif
         Memory_Rand = 1
     else
         Memory_Rand=Memory_Rand+1


         YMD5=(Random_Next_Year*10000) + (Random_Next_Month*100) + Random_Next_Day

         call timemgr_time_inc(YMD5,Stochai_Curr_Sec,              &
                           YMD6,Random_Next_Sec,Stochai_Step,0,0)
         Random_Next_Year =(YMD6/10000)
         YMD6            = YMD6-(Random_Next_Year*10000)
         Random_Next_Month=(YMD6/100)
         Random_Next_Day  = YMD6-(Random_Next_Month*100)
         Ran_File=interpret_filename_spec(Stochai_File_Template      , &
                                         yr_spec=Random_Next_Year , &
                                        mon_spec=Random_Next_Month, &
                                        day_spec=Random_Next_Day  , &
                                        sec_spec=Stochai_Next_Sec    )
         !write(iulog,*) 'Increment Memory:',Memory_Rand
         !write(iulog,*) 'Sustained Random file',Ran_File
     endif
     
     
     if(masterproc) then
      write(iulog,*) 'STOCHAIING: Reading analyses:',trim(Stochai_Path)//trim(Ran_File)
     endif
   
     ! Rotate Stochai_ObsInd() indices for new data, then update
     ! the Stochai observation arrays with analysis data at the
     ! NEXT==Stochai_ObsInd(1) time.
     !----------------------------------------------------------
     if(dycore_is('UNSTRUCTURED')) then
       call stochaiing_update_analyses_se (trim(Stochai_Path)//trim(Ran_File))
     elseif(dycore_is('EUL')) then
       call stochaiing_update_analyses_eul(trim(Stochai_Path)//trim(Ran_File))
     else !if(dycore_is('LR')) then
       call stochaiing_update_analyses_fv (trim(Stochai_Path)//trim(Ran_File))
     endif
   endif ! ((Before_End).and.(Update_Stochai)) then
   
   !!!++WEC Grab a random file!
   !----------------------------------------------------------------
   ! Toggle Stochaiing flag when the time interval is between
   ! beginning and ending times, and all of the analyses files exist.
   !----------------------------------------------------------------
   if((After_Beg).and.(Before_End)) then
     if(Stochai_Force_Opt.eq.0) then
       ! Verify that the NEXT analyses are available
       !---------------------------------------------
       Stochai_ON=Stochai_File_Present(Stochai_ObsInd(1))
     elseif(Stochai_Force_Opt.eq.1) then
       ! Verify that the CURR and NEXT analyses are available
       !-----------------------------------------------------
       Stochai_ON=(Stochai_File_Present(Stochai_ObsInd(1)).and. &
                 Stochai_File_Present(Stochai_ObsInd(2))      )
     else
       ! Verify that the ALL analyses are available
       !---------------------------------------------
       Stochai_ON=.true.
       do nn=1,Stochai_NumObs
         if(.not.Stochai_File_Present(nn)) Stochai_ON=.false.
       end do
     endif
     if(.not.Stochai_ON) then
       if(masterproc) then
         write(iulog,*) 'STOCHAIING: WARNING - analyses file NOT FOUND. Switching '
         write(iulog,*) 'STOCHAIING:           stochaiing OFF to coast thru the gap. '
       endif
     endif
   else
     Stochai_ON=.false.
   endif

   !-------------------------------------------------------
   ! HERE Implement time dependence of Stochaiing Coefs HERE
   !-------------------------------------------------------
   !---------------------------------------------------
   ! If Data arrays have changed update stepping arrays
   !---------------------------------------------------
   if((Before_End).and.((Update_Stochai).or.(Update_Model))) then

     ! Now Load the Target values for stochaiing tendencies
     !---------------------------------------------------
     if(Stochai_Force_Opt.eq.0) then
       ! Target is OBS data at NEXT time
       !----------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_stoch_U(:ncol,:pver,lchnk)=Nobs_U(:ncol,:pver,lchnk,Stochai_ObsInd(1))
         Target_stoch_V(:ncol,:pver,lchnk)=Nobs_V(:ncol,:pver,lchnk,Stochai_ObsInd(1))
         Target_stoch_T(:ncol,:pver,lchnk)=Nobs_T(:ncol,:pver,lchnk,Stochai_ObsInd(1))
         Target_stoch_Q(:ncol,:pver,lchnk)=Nobs_Q(:ncol,:pver,lchnk,Stochai_ObsInd(1))
         Target_stoch_PS(:ncol     ,lchnk)=Nobs_PS(:ncol     ,lchnk,Stochai_ObsInd(1))
       end do
     elseif(Stochai_Force_Opt.eq.1) then
       ! Target is linear interpolation of OBS data CURR<-->NEXT time
       !---------------------------------------------------------------
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Stochai_Next_Year,MM=Stochai_Next_Month, &
                               DD=Stochai_Next_Day , S=Stochai_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tfrac= float(DeltaT)/float(Stochai_Step)
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_stoch_U(:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_U(:ncol,:pver,lchnk,Stochai_ObsInd(1)) &
                                           +Tfrac *Nobs_U(:ncol,:pver,lchnk,Stochai_ObsInd(2))
         Target_stoch_V(:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_V(:ncol,:pver,lchnk,Stochai_ObsInd(1)) &
                                           +Tfrac *Nobs_V(:ncol,:pver,lchnk,Stochai_ObsInd(2))
         Target_stoch_T(:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_T(:ncol,:pver,lchnk,Stochai_ObsInd(1)) &
                                           +Tfrac *Nobs_T(:ncol,:pver,lchnk,Stochai_ObsInd(2))
         Target_stoch_Q(:ncol,:pver,lchnk)=(1._r8-Tfrac)*Nobs_Q(:ncol,:pver,lchnk,Stochai_ObsInd(1)) &
                                           +Tfrac *Nobs_Q(:ncol,:pver,lchnk,Stochai_ObsInd(2))
         Target_stoch_PS(:ncol     ,lchnk)=(1._r8-Tfrac)*Nobs_PS(:ncol     ,lchnk,Stochai_ObsInd(1)) &
                                           +Tfrac *Nobs_PS(:ncol     ,lchnk,Stochai_ObsInd(2))
       end do
     else
       write(iulog,*) 'STOCHAIING: Unknown Stochai_Force_Opt=',Stochai_Force_Opt
       call endrun('stochaiing_timestep_init:: ERROR unknown Stochaiing_Force_Opt')
     endif

     ! Now load Dry Static Energy values for Target
     !---------------------------------------------
     if(Stochai_TSmode.eq.0) then
       ! DSE tendencies from Temperature only
       !---------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         Target_stoch_S(:ncol,:pver,lchnk)=cpair*Target_stoch_T(:ncol,:pver,lchnk)
       end do
     elseif(Stochai_TSmode.eq.1) then
       ! Caluculate DSE tendencies from Temperature, Water Vapor, and Surface Pressure
       !------------------------------------------------------------------------------
       do lchnk=begchunk,endchunk
         ncol=phys_state(lchnk)%ncol
         call calc_DryStaticEnergy(Target_stoch_T(:,:,lchnk), Target_stoch_Q(:,:,lchnk), &
                                 phys_state(lchnk)%phis, Target_stoch_PS(:,lchnk), &
                                                  Target_stoch_S(:,:,lchnk), ncol)
       end do
     endif

     ! Set Tscale for the specified Forcing Option
     !-----------------------------------------------
     if(Stochai_TimeScale_Opt.eq.0) then
       Tscale=1._r8
     elseif(Stochai_TimeScale_Opt.eq.1) then
       call ESMF_TimeSet(Date1,YY=Year,MM=Month,DD=Day,S=Sec)
       call ESMF_TimeSet(Date2,YY=Stochai_Next_Year,MM=Stochai_Next_Month, &
                               DD=Stochai_Next_Day , S=Stochai_Next_Sec    )
       DateDiff =Date2-Date1
       call ESMF_TimeIntervalGet(DateDiff,S=DeltaT,rc=rc)
       Tscale=float(Stochai_Step)/float(DeltaT)
     else
       write(iulog,*) 'STOCHAIING: Unknown Stochai_TimeScale_Opt=',Stochai_TimeScale_Opt
       call endrun('stochaiing_timestep_init:: ERROR unknown Stochaiing_TimeScale_Opt')
     endif

     ! Update the stochaiing tendencies
     !--------------------------------
     Stochai_Do=1 !Defined as integer above [WEC]
     if(Stochai_Do==0) then
     do lchnk=begchunk,endchunk !Remove:
       ncol=phys_state(lchnk)%ncol
       Stochai_Ustep(:ncol,:pver,lchnk)=(  Target_stoch_U(:ncol,:pver,lchnk)      &
                                         -Model_U(:ncol,:pver,lchnk))     &
                                      *Tscale*Stochai_Utau(:ncol,:pver,lchnk)
       Stochai_Vstep(:ncol,:pver,lchnk)=(  Target_stoch_V(:ncol,:pver,lchnk)      &
                                         -Model_V(:ncol,:pver,lchnk))     &
                                      *Tscale*Stochai_Vtau(:ncol,:pver,lchnk)
       Stochai_Sstep(:ncol,:pver,lchnk)=(  Target_stoch_S(:ncol,:pver,lchnk)      &
                                         -Model_S(:ncol,:pver,lchnk))     &
                                      *Tscale*Stochai_Stau(:ncol,:pver,lchnk)
       Stochai_Qstep(:ncol,:pver,lchnk)=(  Target_stoch_Q(:ncol,:pver,lchnk)      &
                                         -Model_Q(:ncol,:pver,lchnk))     &
                                      *Tscale*Stochai_Qtau(:ncol,:pver,lchnk)
       Stochai_PSstep(:ncol,     lchnk)=(  Target_stoch_PS(:ncol,lchnk)      &
                                         -Model_PS(:ncol,lchnk))     &
                                      *Tscale*Stochai_PStau(:ncol,lchnk)
     end do
     endif

     if(Stochai_Do==1) then!ADDED BY [WEC]
     do lchnk=begchunk,endchunk !Remove:
       ncol=phys_state(lchnk)%ncol
       Stochai_Ustep(:ncol,:pver,lchnk)=(Target_stoch_U(:ncol,:pver,lchnk))*Stochai_Utau(:ncol,:pver,lchnk) !These are all tendencies already [WEC]
       Stochai_Vstep(:ncol,:pver,lchnk)=(Target_stoch_V(:ncol,:pver,lchnk))*Stochai_Vtau(:ncol,:pver,lchnk)
       Stochai_Sstep(:ncol,:pver,lchnk)=(Target_stoch_S(:ncol,:pver,lchnk))*Stochai_Stau(:ncol,:pver,lchnk) !cant do Temp [WEC]
       Stochai_Qstep(:ncol,:pver,lchnk)=(Target_stoch_Q(:ncol,:pver,lchnk))*Stochai_Qtau(:ncol,:pver,lchnk)
       Stochai_PSstep(:ncol,     lchnk)=(Target_stoch_PS(:ncol,lchnk))*Stochai_PStau(:ncol,lchnk) !cant do PS [WEC]
     end do
     endif
     !******************
     ! DIAG
     !******************
!    if(masterproc) then
!      write(iulog,*) 'PFC: Target_stoch_T(1,:pver,begchunk)=',Target_stoch_T(1,:pver,begchunk)
!      write(iulog,*) 'PFC:  Model_T(1,:pver,begchunk)=',Model_T(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Target_stoch_S(1,:pver,begchunk)=',Target_stoch_S(1,:pver,begchunk)
!      write(iulog,*) 'PFC:  Model_S(1,:pver,begchunk)=',Model_S(1,:pver,begchunk)
!      write(iulog,*) 'PFC:      Target_stoch_PS(1,begchunk)=',Target_stoch_PS(1,begchunk)
!      write(iulog,*) 'PFC:       Model_PS(1,begchunk)=',Model_PS(1,begchunk)
!      write(iulog,*) 'PFC: Stochai_Sstep(1,:pver,begchunk)=',Stochai_Sstep(1,:pver,begchunk)
!      write(iulog,*) 'PFC: Stochai_Xstep arrays updated:'
!    endif
   endif ! ((Before_End).and.((Update_Stochai).or.(Update_Model))) then

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_timestep_init
  !================================================================


  !================================================================
  subroutine stochaiing_timestep_tend(phys_state,phys_tend)
   !
   ! STOCHAIING_TIMESTEP_TEND:
   !                If Stochaiing is ON, return the Stochaiing contributions
   !                to forcing using the current contents of the Stochai
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
   call physics_ptend_init(phys_tend,phys_state%psetcols,'stochaiing',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if(Stochai_ON) then !this is where the stochaiing tendency gets teleported into the model. [WEC]
     lchnk=phys_state%lchnk
     ncol =phys_state%ncol
     phys_tend%u(:ncol,:pver)     =Stochai_Ustep(:ncol,:pver,lchnk)
     phys_tend%v(:ncol,:pver)     =Stochai_Vstep(:ncol,:pver,lchnk)
     phys_tend%s(:ncol,:pver)     =Stochai_Sstep(:ncol,:pver,lchnk)
     phys_tend%q(:ncol,:pver,indw)=Stochai_Qstep(:ncol,:pver,lchnk)

     call outfld( 'Stochai_U',phys_tend%u                ,pcols,lchnk)
     call outfld( 'Stochai_V',phys_tend%v                ,pcols,lchnk)
     call outfld( 'Stochai_T',phys_tend%s/cpair          ,pcols,lchnk) !MUST EDIT TO MESS WITH T->DSE [WEC]
     call outfld( 'Stochai_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)
     call outfld('Target_stoch_U',Target_stoch_U(:,:,lchnk),pcols,lchnk)
     call outfld('Target_stoch_V',Target_stoch_V(:,:,lchnk),pcols,lchnk)
     call outfld('Target_stoch_T',Target_stoch_T(:,:,lchnk),pcols,lchnk)
     call outfld('Target_stoch_Q',Target_stoch_Q(:,:,lchnk),pcols,lchnk)
   endif

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_timestep_tend
  !================================================================


  !================================================================
  subroutine stochaiing_update_analyses_se(anal_file)
   !
   ! STOCHAIING_UPDATE_ANALYSES_SE:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer ncol,plev,istat
   integer ncid,varid
   real(r8) Xanal(Stochai_ncol,Stochai_nlev)
   real(r8) PSanal(Stochai_ncol)
   real(r8) Lat_anal(Stochai_ncol)
   real(r8) Lon_anal(Stochai_ncol)
   integer  nn,Nindex

   ! Rotate Stochai_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Stochai_ObsInd(Stochai_NumObs)
     do nn=Stochai_NumObs,2,-1
       Stochai_ObsInd(nn)=Stochai_ObsInd(nn-1)
     end do
     Stochai_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Stochai_File_Present(Stochai_ObsInd(1)))
     write(iulog,*)'STOCHAIING: Stochai_ObsInd=',Stochai_ObsInd
     write(iulog,*)'STOCHAIING: Stochai_File_Present=',Stochai_File_Present
   endif
#ifdef SPMD
   call mpibcast(Stochai_File_Present, Stochai_NumObs, mpilog, 0, mpicom)
   call mpibcast(Stochai_ObsInd      , Stochai_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Stochai_File_Present(Stochai_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'ncol',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=ncol)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif

     if((Stochai_ncol.ne.ncol).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: stochaiing_update_analyses_se: ncol=',ncol,' Stochai_ncol=',Stochai_ncol
      write(iulog,*) 'ERROR: stochaiing_update_analyses_se: plev=',plev,' pver=',pver
      call endrun('stochaiing_update_analyses_se: analyses dimension mismatch')
     endif

     ! Read in and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_ncol,Xanal,    &
                               Nobs_U(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_ncol,Xanal,    &
                               Nobs_V(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_ncol,Xanal,    &
                               Nobs_T(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_ncol,Xanal,    &
                               Nobs_Q(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_SE')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Stochai_ncol,PSanal,           &
                               Nobs_PS(1,begchunk,Stochai_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_update_analyses_se
  !================================================================


  !================================================================
  subroutine stochaiing_update_analyses_eul(anal_file)
   !
   ! STOCHAIING_UPDATE_ANALYSES_EUL:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Stochai_nlon,Stochai_nlat,Stochai_nlev)
   real(r8) PSanal(Stochai_nlon,Stochai_nlat)
   real(r8) Lat_anal(Stochai_nlat)
   real(r8) Lon_anal(Stochai_nlon)
   real(r8) Xtrans(Stochai_nlon,Stochai_nlev,Stochai_nlat)
   integer  nn,Nindex

   ! Rotate Stochai_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Stochai_ObsInd(Stochai_NumObs)
     do nn=Stochai_NumObs,2,-1
       Stochai_ObsInd(nn)=Stochai_ObsInd(nn-1)
     end do
     Stochai_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Stochai_File_Present(Stochai_ObsInd(1)))
   endif
#ifdef SPMD
   call mpibcast(Stochai_File_Present, Stochai_NumObs, mpilog, 0, mpicom)
   call mpibcast(Stochai_ObsInd      , Stochai_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Stochai_File_Present(Stochai_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif

     if((Stochai_nlon.ne.nlon).or.(Stochai_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: stochaiing_update_analyses_eul: nlon=',nlon,' Stochai_nlon=',Stochai_nlon
      write(iulog,*) 'ERROR: stochaiing_update_analyses_eul: nlat=',nlat,' Stochai_nlat=',Stochai_nlat
      write(iulog,*) 'ERROR: stochaiing_update_analyses_eul: plev=',plev,' pver=',pver
      call endrun('stochaiing_update_analyses_eul: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_U(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_V(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_T(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_Q(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Stochai_nlon,PSanal,           &
                               Nobs_PS(1,begchunk,Stochai_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_update_analyses_eul
  !================================================================


  !================================================================
  subroutine stochaiing_update_analyses_fv(anal_file)!EDIT THIS [WEC]
   !
   ! STOCHAIING_UPDATE_ANALYSES_FV:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid ,only: pver,begchunk
   use netcdf

   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file

   ! Local values
   !-------------
   integer lev
   integer nlon,nlat,plev,istat
   integer ncid,varid
   integer ilat,ilon,ilev
   real(r8) Xanal(Stochai_nlon,Stochai_nlat,Stochai_nlev)
   real(r8) PSanal(Stochai_nlon,Stochai_nlat)
   real(r8) Lat_anal(Stochai_nlat)
   real(r8) Lon_anal(Stochai_nlon)
   real(r8) Xtrans(Stochai_nlon,Stochai_nlev,Stochai_nlat)
   integer  nn,Nindex

   ! Rotate Stochai_ObsInd() indices, then check the existence of the analyses
   ! file; broadcast the updated indices and file status to all the other MPI nodes.
   ! If the file is not there, then just return.
   !------------------------------------------------------------------------
   if(masterproc) then
     Nindex=Stochai_ObsInd(Stochai_NumObs)
     do nn=Stochai_NumObs,2,-1
       Stochai_ObsInd(nn)=Stochai_ObsInd(nn-1)
     end do
     Stochai_ObsInd(1)=Nindex
     inquire(FILE=trim(anal_file),EXIST=Stochai_File_Present(Stochai_ObsInd(1)))
     write(iulog,*)'STOCHAIING: Stochai_ObsInd=',Stochai_ObsInd
     write(iulog,*)'STOCHAIING: Stochai_File_Present=',Stochai_File_Present
   endif
#ifdef SPMD
   call mpibcast(Stochai_File_Present, Stochai_NumObs, mpilog, 0, mpicom)
   call mpibcast(Stochai_ObsInd      , Stochai_NumObs, mpiint, 0, mpicom)
#endif
   if(.not.Stochai_File_Present(Stochai_ObsInd(1))) return

   ! masterporc does all of the work here
   !-----------------------------------------
   if(masterproc) then

     ! Open the given file
     !-----------------------
     istat=nf90_open(trim(anal_file),NF90_NOWRITE,ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*)'NF90_OPEN: failed for file ',trim(anal_file)
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     ! Read in Dimensions
     !--------------------
     istat=nf90_inq_dimid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlon)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=nlat)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_dimid(ncid,'lev',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_inquire_dimension(ncid,varid,len=plev)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lon',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lon_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     istat=nf90_inq_varid(ncid,'lat',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Lat_anal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif

     if((Stochai_nlon.ne.nlon).or.(Stochai_nlat.ne.nlat).or.(plev.ne.pver)) then
      write(iulog,*) 'ERROR: stochaiing_update_analyses_fv: nlon=',nlon,' Stochai_nlon=',Stochai_nlon
      write(iulog,*) 'ERROR: stochaiing_update_analyses_fv: nlat=',nlat,' Stochai_nlat=',Stochai_nlat
      write(iulog,*) 'ERROR: stochaiing_update_analyses_fv: plev=',plev,' pver=',pver
      call endrun('stochaiing_update_analyses_fv: analyses dimension mismatch')
     endif

     ! Read in, transpose lat/lev indices,
     ! and scatter data arrays
     !----------------------------------
     istat=nf90_inq_varid(ncid,'U',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_U(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'V',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_V(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'T',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_T(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
     istat=nf90_inq_varid(ncid,'Q',varid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     istat=nf90_get_var(ncid,varid,Xanal)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_FV')
     endif
     do ilat=1,nlat
     do ilev=1,plev
     do ilon=1,nlon
       Xtrans(ilon,ilev,ilat)=Xanal(ilon,ilat,ilev)
     end do
     end do
     end do
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,Stochai_nlev,1,Stochai_nlon,Xtrans,   &
                               Nobs_Q(1,1,begchunk,Stochai_ObsInd(1)))

   if(masterproc) then
    istat=nf90_inq_varid(ncid,'PS',varid)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif
    istat=nf90_get_var(ncid,varid,PSanal)
    if(istat.ne.NF90_NOERR) then
      write(iulog,*) nf90_strerror(istat)
      call endrun ('UPDATE_ANALYSES_SE')
    endif

     ! Close the analyses file
     !-----------------------
     istat=nf90_close(ncid)
     if(istat.ne.NF90_NOERR) then
       write(iulog,*) nf90_strerror(istat)
       call endrun ('UPDATE_ANALYSES_EUL')
     endif
   endif ! (masterproc) then
   call scatter_field_to_chunk(1,1,1,Stochai_nlon,PSanal,           &
                               Nobs_PS(1,begchunk,Stochai_ObsInd(1)))

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_update_analyses_fv
  !================================================================


  !================================================================
  subroutine stochaiing_set_profile(rlat,rlon,Stochai_prof,Wprof,nlev)
   !
   ! STOCHAIING_SET_PROFILE: for the given lat,lon, and Stochaiing_prof, set
   !                      the verical profile of window coeffcients.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on stochaiing strength.
   !===============================================================

   ! Arguments
   !--------------
   integer  nlev,Stochai_prof
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
   if(Stochai_prof.eq.0) then
     ! No Stochaiing
     !-------------
     Wprof(:)=0.0_r8
   elseif(Stochai_prof.eq.1) then
     ! Uniform Stochaiing
     !-----------------
     Wprof(:)=1.0_r8
   elseif(Stochai_prof.eq.2) then
     ! Localized Stochaiing with specified Heaviside window function
     !------------------------------------------------------------
     if(Stochai_Hwin_max.le.Stochai_Hwin_min) then
       ! For a constant Horizontal window function,
       ! just set Hcoef to the maximum of Hlo/Hhi.
       !--------------------------------------------
       Hcoef=max(Stochai_Hwin_lo,Stochai_Hwin_hi)
     else
       ! get lat/lon relative to window center
       !------------------------------------------
       latx=rlat-Stochai_Hwin_lat0
       lonx=rlon-Stochai_Hwin_lon0
       if(lonx.gt. 180._r8) lonx=lonx-360._r8
       if(lonx.le.-180._r8) lonx=lonx+360._r8

       ! Calcualte RAW window value
       !-------------------------------
       lon_lo=(Stochai_Hwin_lonWidthH+lonx)/Stochai_Hwin_lonDelta
       lon_hi=(Stochai_Hwin_lonWidthH-lonx)/Stochai_Hwin_lonDelta
       lat_lo=(Stochai_Hwin_latWidthH+latx)/Stochai_Hwin_latDelta
       lat_hi=(Stochai_Hwin_latWidthH-latx)/Stochai_Hwin_latDelta
       Hcoef=((1._r8+tanh(lon_lo))/2._r8)*((1._r8+tanh(lon_hi))/2._r8) &
            *((1._r8+tanh(lat_lo))/2._r8)*((1._r8+tanh(lat_hi))/2._r8)

       ! Scale the horizontal window coef for specfied range of values.
       !--------------------------------------------------------
       Hcoef=(Hcoef-Stochai_Hwin_min)/(Stochai_Hwin_max-Stochai_Hwin_min)
       Hcoef=(1._r8-Hcoef)*Stochai_Hwin_lo + Hcoef*Stochai_Hwin_hi
     endif

     ! Load the RAW vertical window
     !------------------------------
     do ilev=1,nlev
       lev_lo=(float(ilev)-Stochai_Vwin_Lindex)/Stochai_Vwin_Ldelta
       lev_hi=(Stochai_Vwin_Hindex-float(ilev))/Stochai_Vwin_Hdelta
       Wprof(ilev)=((1._r8+tanh(lev_lo))/2._r8)*((1._r8+tanh(lev_hi))/2._r8)
     end do

     ! Scale the Window function to span the values between Vlo and Vhi:
     !-----------------------------------------------------------------
     Vmax=maxval(Wprof)
     Vmin=minval(Wprof)
     if((Vmax.le.Vmin).or.((Stochai_Vwin_Hindex.ge.(nlev+1)).and. &
                           (Stochai_Vwin_Lindex.le. 0      )     )) then
       ! For a constant Vertical window function,
       ! load maximum of Vlo/Vhi into Wprof()
       !--------------------------------------------
       Vmax=max(Stochai_Vwin_lo,Stochai_Vwin_hi)
       Wprof(:)=Vmax
     else
       ! Scale the RAW vertical window for specfied range of values.
       !--------------------------------------------------------
       Wprof(:)=(Wprof(:)-Vmin)/(Vmax-Vmin)
       Wprof(:)=Stochai_Vwin_lo + Wprof(:)*(Stochai_Vwin_hi-Stochai_Vwin_lo)
     endif

     ! The desired result is the product of the vertical profile
     ! and the horizontal window coeffcient.
     !----------------------------------------------------
     Wprof(:)=Hcoef*Wprof(:)
   else
     call endrun('stochaiing_set_profile:: Unknown Stochai_prof value')
   endif

   ! End Routine
   !------------
   return
  end subroutine ! stochaiing_set_profile
  !================================================================


  !================================================================
  real(r8) function stochaiing_set_PSprofile(rlat,rlon,Stochai_PSprof)
   !
   ! STOCHAIING_SET_PSPROFILE: for the given lat and lon set the surface
   !                      pressure profile value for the specified index.
   !                      Values range from 0. to 1. to affect spatial
   !                      variations on stochaiing strength.
   !===============================================================

   ! Arguments
   !--------------
   real(r8) rlat,rlon
   integer  Stochai_PSprof

   ! Local values
   !----------------

   !---------------
   ! set coeffcient
   !---------------
   if(Stochai_PSprof.eq.0) then
     ! No Stochaiing
     !-------------
     stochaiing_set_PSprofile=0.0_r8
   elseif(Stochai_PSprof.eq.1) then
     ! Uniform Stochaiing
     !-----------------
     stochaiing_set_PSprofile=1.0_r8
   else
     call endrun('stochaiing_set_PSprofile:: Unknown Stochai_prof value')
   endif

   ! End Routine
   !------------
   return
  end function ! stochaiing_set_PSprofile
  !================================================================


  !================================================================
  subroutine calc_DryStaticEnergy(t, q, phis, ps, dse, ncol)
   !
   ! calc_DryStaticEnergy: Given the temperature, specific humidity, surface pressure,
   !                       and surface geopotential for a chunk containing 'ncol' columns,
   !                       calculate and return the corresponding dry static energy values.
   !--------------------------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pver, pverp
   use dycore,       only: dycore_is
   use hycoef,       only: hyai, hybi, ps0, hyam, hybm
   use physconst,    only: zvir, gravit, cpair, rair
   !
   ! Input/Output arguments
   !-----------------------
   integer , intent(in) :: ncol      ! Number of columns in chunk
   real(r8), intent(in) :: t(:,:)    ! (pcols,pver) - temperature
   real(r8), intent(in) :: q(:,:)    ! (pcols,pver) - specific humidity
   real(r8), intent(in) :: ps(:)     ! (pcols)      - surface pressure
   real(r8), intent(in) :: phis(:)   ! (pcols)      - surface geopotential
   real(r8), intent(out):: dse(:,:)  ! (pcols,pver)  - dry static energy
   !
   ! Local variables
   !------------------
   logical  :: fvdyn                 ! finite volume dynamics
   integer  :: ii,kk                 ! Lon, level, level indices
   real(r8) :: tvfac                 ! Virtual temperature factor
   real(r8) :: hkk(ncol)             ! diagonal element of hydrostatic matrix
   real(r8) :: hkl(ncol)             ! off-diagonal element
   real(r8) :: pint(ncol,pverp)      ! Interface pressures
   real(r8) :: pmid(ncol,pver )      ! Midpoint pressures
   real(r8) :: zi(ncol,pverp)        ! Height above surface at interfaces
   real(r8) :: zm(ncol,pver )        ! Geopotential height at mid level

   ! Set dynamics flag
   !-------------------
   fvdyn = dycore_is ('LR')

   ! Load Pressure values and midpoint pressures
   !----------------------------------------------
   do kk=1,pverp
     do ii=1,ncol
       pint(ii,kk)=(hyai(kk)*ps0)+(hybi(kk)*ps(ii))
     end do
   end do
   do kk=1,pver
     do ii=1,ncol
       pmid(ii,kk)=(hyam(kk)*ps0)+(hybm(kk)*ps(ii))
     end do
   end do

   ! The surface height is zero by definition.
   !-------------------------------------------
   do ii = 1,ncol
     zi(ii,pverp) = 0.0_r8
   end do

   ! Compute the dry static energy, zi, zm from bottom up
   ! Note, zi(i,k) is the interface above zm(i,k)
   !---------------------------------------------------------
   do kk=pver,1,-1

     ! First set hydrostatic elements consistent with dynamics
     !--------------------------------------------------------
     if(fvdyn) then
       do ii=1,ncol
         hkl(ii)=log(pint(ii,kk+1))-log(pint(ii,kk))
         hkk(ii)=1._r8-(hkl(ii)*pint(ii,kk)/(pint(ii,kk+1)-pint(ii,kk)))
       end do
     else
       do ii=1,ncol
         hkl(ii)=(pint(ii,kk+1)-pint(ii,kk))/pmid(ii,kk)
         hkk(ii)=0.5_r8*hkl(ii)
       end do
     endif

     ! Now compute zm, zi, and dse  (WACCM-X vars rairv/zairv/cpairv not used!)
     !------------------------------------------------------------------------
     do ii=1,ncol
       tvfac=t(ii,kk)*rair*(1._r8+(zvir*q(ii,kk)))/gravit
       zm (ii,kk)=zi(ii,kk+1) + (tvfac*hkk(ii))
       zi (ii,kk)=zi(ii,kk+1) + (tvfac*hkl(ii))
       dse(ii,kk)=(t(ii,kk)*cpair) + phis(ii) + (gravit*zm(ii,kk))
     end do

   end do ! kk=pver,1,-1

   ! End Routine
   !-----------
   return
  end subroutine calc_DryStaticEnergy
  !================================================================

end module stochaiing
