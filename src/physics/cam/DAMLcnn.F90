module CB24cnn

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

  public:: CB24cnn_Model!,DAMLin_Model_Regress,DAMLin_ON
  !public:: damlining_readnl
  public:: CB24cnn_init
  !public:: regress_daml_timestep_tend
  !public:: regress_diurnal_daml_timestep_tend
  !public:: damlining_timestep_init
  !public:: damlining_timestep_tend
  !public:: damlining_diurnal_timestep_tend
  
  !private::DAMLining_readin_weights_all
  !private::DAMLining_readin_weights
  !private::DAMLining_readin_weights_diurnal
  !private::DAMLining_readin_normdict
  !private::DAMLining_set_PSprofile
  !private::DAMLining_set_profile


  ! DAMLining Parameters
  !--------------------
  logical          :: CB24cnn_Model       =.true.


   !================================================================
  subroutine CB24cnn_init
    !
    ! DAMLinING_INIT: Allocate space and initialize DAMLining values
    !===============================================================

    use forpy_mod
    
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


    write(iulog,*)'CB24cnn_init starting'

    ierror = forpy_initialize()
    ierror = list_create(my_list)

    ierror = my_list%append(19)
    ierror = my_list%append("Hello world!")
    ierror = my_list%append(3.14d0)
    ierror = print_py(my_list)
    
    write(iulog,*)my_list

    call my_list%destroy
    call forpy_finalize



  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine ! CB24cnn_init













  