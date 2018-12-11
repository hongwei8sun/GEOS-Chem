module lagrange_mod

  use precision_mod
  use CMN_SIZE_mod

  implicit none
  private

  public init_lagrange
  public run_lagrange
  public cleanup_lagrange

  integer, parameter :: n_boxes_max = 500
  real(fp), allocatable :: box_lon    !first use one point to test
  real(fp), allocatable :: box_lat

  subroutine init_lagrange(am_I_Root)
    logical, intent(in) :: am_I_Root

    ! Allocate statements
    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))

    box_lon = 0.0e+0_fp
    box_lat = 0.0e+0_fp

  end subroutine init_lagrange

  subroutine run_lagrange(am_I_Root, state_met, state_chm, input_opt)
    use state_chm_mod, only : chmstate
    use state_met_mod, only : metstate
    use input_opt_mod, only ：OptInput
    
    logical, intent(in) :: am_I_Root
    TYPE(MetState), intent(in) :: state_met
    TYPE(ChmState), intent(inout) :: state_chm
    TYPE(OptInput), intent(in) :: input_opt

    real :: Pi 
    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    ! Establish pointers
    u => state_met%U   ! figure out state_met%U is based on lat/lon or model grid(i,j)
    v => state_met%V
    omeg => state_met%OMEGA
    ! Run Lagrangian advection HERE
    Pi = 3.14   
    R_earth = 6.371e+6

    Dbox_lat = (time_step*v)/(Pi*R_earth)*180            !model time_step length
    Dbox_lon = (time_step*u)/(2.0*Pi*R_earth*cos(lat*Pi/180.0))*360

    box_lat = box_lat + Dbox_lat
    box_lon = box_lon + Dbox_lon

    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)

  end subroutine run_lagrange

  subroutine cleanup_lagrange()
    if (allocated(box_lon)) deallocate(box_lon)
    if (allocated(box_lat)) deallocate(box_lat)
  end subroutine cleanup_lagrange

end module lagrange_mod
