module CB24cnn
  ! Stub: FORPY-based DAML CNN replaced with no-op.
  ! CAMulator runs as an external server; no Python-Fortran bridge required.
  implicit none
  private
  public :: CB24cnn_init
  public :: CB24cnn_Model
  logical :: CB24cnn_Model = .false.
contains
  subroutine CB24cnn_init
  end subroutine CB24cnn_init
end module CB24cnn
