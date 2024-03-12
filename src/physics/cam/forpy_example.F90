!example of a forpy_f90 scripts
module forpyex

!=====================================================================
!
! Purpose: Implement forpyexining of the model state of U,V,T,Q, and/or PS
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
!    forpyexin the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the forpyexining
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the forpyexining module is set up to accomodate forpyexining of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The forpyexining of the model toward the analyses data is controlled by
!    the 'forpyexining_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which forpyexining is applied, the strength of the forpyexining
!    tendencies, and its spatial distribution.
!

! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
  use, intrinsic :: ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all, ieee_support_halting, ieee_overflow
  use phys_grid   ,   only:scatter_field_to_chunk
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc, mstrid=>masterprocid, mpicom
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  !use MOM_coms,                  only : PE_here,num_PEs
  use forpy_mod,                 only : module_py,list,ndarray,object,tuple
  use forpy_mod,                 only : err_print
  use forpy_mod,                 only : forpy_initialize,get_sys_path,import_py,print_py
  use forpy_mod,                 only : ndarray_create,tuple_create,call_py,cast
  use forpy_mod,                 only : forpy_finalize, dict, dict_create

! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: forpyex_Model
  public:: forpyex_init
  public:: forpyex_run
  public:: forpyex_finalize


  ! forpyexining Parameters
  !--------------------
  logical          :: forpyex_Model       =.true.
  logical          :: halting_mode(5)
  type(module_py) :: pymodule
contains
   !================================================================
  subroutine forpyex_init
    !
    ! forpyexinING_INIT: Allocate space and initialize forpyexining values
    !===============================================================
    use ppgrid        ,only: pver,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec

    implicit none

    integer :: ierror
    type(list) :: my_list
    type(list) :: paths

    character(len=:), allocatable :: return_string


    if (masterproc) write(iulog,*) 'forpyex_init starting'
    if (.not. ieee_support_halting(ieee_overflow)) then
       call endrun("ieee_halting is not supported")
    endif
    call ieee_get_halting_mode(ieee_all, halting_mode)
    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = forpy_initialize()
    ierror = get_sys_path(paths)
    ierror = paths%append(".")


    ierror = import_py(pymodule,"forpyexcnn")
    call ieee_set_halting_mode(ieee_all, halting_mode)
    call paths%destroy

  end subroutine forpyex_init
  subroutine forpyex_run
    integer :: ierror
    type(tuple) :: args
    type(dict) :: kwargs
    type(object) :: return_value

    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = tuple_create(args, 3)
    ierror = args%setitem(0, 12)
    ierror = args%setitem(1, "Hi")
    ierror = args%setitem(2, .true.)

    ierror = dict_create(kwargs)
    ierror = kwargs%setitem("message", "Hello world!")

    ierror = call_py(return_value,pymodule,"forpyexcnn_run", args, kwargs)
!    ierror = cast(return_string, return_value)
!    if (masterproc) write(iulog,*)'forpyex_init ending',return_string
    call ieee_set_halting_mode(ieee_all, halting_mode)

    call args%destroy
    call kwargs%destroy
    call return_value%destroy

  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine forpyex_run ! forpyex_init
  subroutine forpyex_finalize()
    call pymodule%destroy
    call forpy_finalize
  end subroutine forpyex_finalize
  end module forpyex
