!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Plume_Mod

  USE precision_mod
  USE PhysConstants, ONLY : PI, Re 
  USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: plume_init
  PUBLIC :: plume_run
  PUBLIC :: plume_cleanup


  integer, parameter :: n_boxes_max = 3  ! 50*40*40
  integer            :: tt

  ! The center of the cylinder
  real(fp), allocatable :: box_lon(:)
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)

  ! Consider the box as a cylinder

  ! vertical radius of the the cross-section
  real(fp), allocatable :: box_radius1(:)  
  ! horizontal radius of the the cross-section
  real(fp), allocatable :: box_radius2(:)
  ! Theta is the clockwise angle between z-axis (P) and vertical radius1
  real(fp), allocatable :: box_theta(:)  ! 0 ~ 180 degree
  real(fp), allocatable :: box_length(:)

CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE plume_init( am_I_root )

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU

    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Plume Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))

    allocate(box_radius1(n_boxes_max))
    allocate(box_radius2(n_boxes_max))
    allocate(box_theta(n_boxes_max))
    allocate(box_length(n_boxes_max))


    box_lon    = (/5.0e+0_fp,  5.1e+0_fp,  5.2e+0_fp/)
    box_lat    = (/4.0e+0_fp,  4.1e+0_fp,  4.2e+0_fp/)
    box_lev    = (/20.0e+0_fp, 20.0e+0_fp, 20.0e+0_fp/)      ! hPa

    box_radius1  = (/100.0e+0_fp,  100.0e+0_fp,  100.0e+0_fp/)     ! m
    box_radius2  = (/100.0e+0_fp,  100.0e+0_fp,  100.0e+0_fp/)     ! m
    box_theta    = (/0.0e+0_fp,    0.0e+0_fp,    0.0e+0_fp/)             ! degree
    box_length   = (/1000.0e+0_fp, 1000.0e+0_fp, 1000.0e+0_fp/)    ! m


  END SUBROUTINE plume_init

!-----------------------------------------------------------------
!=================================================================

  SUBROUTINE plume_run(am_I_Root, State_Met, Input_Opt)

    USE Input_Opt_Mod, ONLY : OptInput

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE, XMID, YMID                 
    USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR, DLAT, DLON 
    ! DLAT( IIPAR, JJPAR, LLPAR ), DLON( IIPAR, JJPAR, LLPAR )
    ! XEDGE  ( IM+1, JM,   L ), YEDGE  ( IM,   JM+1, L ), IM=IIPAR, JM=JJPAR

    logical, intent(in) :: am_I_Root
    TYPE(MetState), intent(in) :: State_Met
    !TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in) :: Input_Opt

    REAL :: Dt          ! = 600.0e+0_fp          

    integer :: i_box

    integer :: i_lon
    integer :: i_lat
    integer :: i_lev

    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: next_lon, next_lat, next_pressure

    real(fp) :: X_edge2
    real(fp) :: Y_edge2

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp), pointer :: P_BXHEIGHT(:,:,:)

    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: next_u, next_v, next_omeg

    real(fp), pointer :: P_edge(:)
    real(fp), pointer :: P_mid(:)

    real(fp) :: Dx, Dy
    real(fp) :: wind_s_shear
    real(fp) :: dbox_theta, dbox_radius1, dbox_radius2
    ! angle between travel direction and lon (-90 ~ 90 deg)
    real(fp):: box_alpha        

    real(fp), pointer :: X_edge(:)    
    real(fp), pointer :: Y_edge(:)       
    real(fp), pointer :: X_mid(:)    
    real(fp), pointer :: Y_mid(:)       

    real(fp), pointer :: P_I, R_e     

    Dt = GET_TS_DYN()

    ! Establish pointers
    ! P_I = PI
    ! R_e = Re

    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V   ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]
    P_mid => State_Met%PMID(1,1,:)  ! Pressure (w/r/t moist air) at level centers (hPa)

    P_BXHEIGHT => State_Met%BXHEIGHT

    Dx = DLON(1,1,1)
    Dy = DLAT(1,2,1)  ! DLAT(1,1,1) is half of DLAT(1,2,1) !!!
    X_edge => XEDGE(:,1,1)   ! IIPAR+1
    Y_edge => YEDGE(1,:,1)  
    X_mid => XMID(:,1,1)   ! IIPAR+1
    Y_mid => YMID(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)

    ! Run Plume advection HERE
    do i_box = 1,n_boxes_max

       ! make sure the location is not out of range
       do while (box_lat(i_box) > Y_edge(JJPAR+1))
          box_lat(i_box) = Y_edge(JJPAR+1) - ( box_lat(i_box)-Y_edge(JJPAR+1) )
       end do

       do while (box_lat(i_box) < Y_edge(1))
          box_lat(i_box) = Y_edge(1) + ( box_lat(i_box)-Y_edge(1) )
       end do

       do while (box_lon(i_box) > X_edge(IIPAR+1))
          box_lon(i_box) = box_lon(i_box) - 360.0
       end do

       do while (box_lon(i_box) < X_edge(1))
          box_lon(i_box) = box_lon(i_box) + 360.0
       end do


       curr_lon      = box_lon(i_box)
       curr_lat      = box_lat(i_box)
       curr_pressure = box_lev(i_box)        ! hPa

       ! calcualte the box_alpha
       ! Next adjacent point
       if(i_box==n_boxes_max)then
         next_lon      = box_lon(i_box-1)
         next_lat      = box_lat(i_box-1)
         next_pressure = box_lev(i_box-1)        ! hPa

       else
         next_lon      = box_lon(i_box+1)
         next_lat      = box_lat(i_box+1)
         next_pressure = box_lev(i_box+1)        ! hPa

       endif


       if((curr_lon-next_lon).NE.0.0)then
         box_alpha = ATAN( (next_lat - curr_lat) / (next_lon - curr_lon) )
       else if((next_lat-curr_lat).LT.0.0)then
         box_alpha = 0.5*PI
       else if((next_lat-curr_lat).GT.0.0)then
         box_alpha = -0.5*PI
       endif


       i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2) 
       i_lev = Find_iPLev(curr_pressure,P_edge)

       curr_u    = Interplt_wind(u,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_v    = Interplt_wind(v,   X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       curr_omeg = Interplt_wind(omeg,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


       i_lon = Find_iLonLat(next_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(next_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(next_pressure,P_edge)


       ! along the direction of airplane's movement
       !wind_l = curr_u*COS(alpha) + curr_v*SIN(alpha)
       ! horizontal direction inside the plume cross_section
       !wind_s = curr_u*SIN(alpha) + curr_v*COS(alpha)

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = wind_shear(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)


       ! dbox_theta   = wind_s_shear * (COS(box_theta(i_box))**2) * Dt
       ! dbox_theta   = ATAN( TAN( box_theta(i_box) ) + Dt*wind_s_shear ) - box_theta(i_box)

       box_theta(i_box)   = ATAN( TAN( box_theta(i_box) ) + Dt*wind_s_shear )

       if(box_theta(i_box)>=0.0)then
         dbox_radius1 = box_radius1(i_box) * wind_s_shear * SIN(box_theta(i_box)) * COS(box_theta(i_box)) * Dt
         dbox_radius2 = -1.0 * box_radius2(i_box) * wind_s_shear * SIN(box_theta(i_box)) * COS(box_theta(i_box)) * Dt
       else
         dbox_radius1 = -1.0 * box_radius1(i_box) * wind_s_shear * SIN(box_theta(i_box)) * COS(box_theta(i_box)) * Dt
         dbox_radius2 = box_radius2(i_box) * wind_s_shear * SIN(box_theta(i_box)) * COS(box_theta(i_box)) * Dt
       endif

       !box_theta(i_box)   = box_theta(i_box)   + dbox_theta
       box_radius1(i_box) = box_radius1(i_box) + dbox_radius1
       box_radius2(i_box) = box_radius2(i_box) + dbox_radius2

       write(6,*)'= plume =>', wind_s_shear, box_theta(i_box), box_radius1(i_box), box_radius2(i_box)

    end do

    ! WRITE(6,*)'-Plume(i_lev)->',i_lev,P_edge(i_lev),P_edge(i_lev+1)

    ! Everything is done, clean up pointers
    nullify(u)
    nullify(v)
    nullify(omeg)
    nullify(X_edge)
    nullify(Y_edge)
    nullify(P_edge)
    nullify(X_mid)
    nullify(Y_mid)
    nullify(P_mid)

  END SUBROUTINE plume_run

!--------------------------------------------------------------------
! functions to find location (i,j,k) of boxes 

  integer function Find_iLonLat(curr_xy,  Dxy,  XY_edge2)
    implicit none
    real(fp) :: curr_xy, Dxy, XY_edge2
    Find_iLonLat = INT( (curr_xy - (XY_edge2 - Dxy)) / Dxy )+1
    ! Notice the difference between INT(), FLOOR(), AINT()
    ! for lon: Xedge_Sec - Dx = Xedge_first
    return
  end function


  integer function Find_iPLev(curr_pressure,P_edge)
    implicit none
    real(fp) :: curr_pressure
    real(fp), pointer :: P_edge(:)
    integer :: i_lon, i_lat
    integer :: locate(1)
    locate = MINLOC(abs( P_edge(:)-curr_pressure ))

    if(P_edge(locate(1))-curr_pressure >= 0 )then
       Find_iPLev = locate(1)
    else
       Find_iPLev = locate(1) - 1
    endif

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real function Interplt_wind(wind, X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          ::  distance(2,2), Weight(2,2)
    real(fp)          ::wind_lonlat(2), wind_lonlat_lev

    ! Identify wether particle is exactly located on the grid point
    if(curr_lon==X_mid(i_lon))then
      if(curr_lat==Y_mid(i_lat))then
        if(curr_pressure==P_mid(i_lev))then

          Interplt_wind = wind(i_lon, i_lat, i_lev)
          return

        endif
      endif
    endif

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
        ii = i + init_lon - 1
        jj = j + init_lat - 1

        ! For some special circumstance:
        if(ii==0)then
          distance(i,j) = Distance_Circle(curr_lon+360.0, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
        else if(ii==(IIPAR+1))then
          distance(i,j) = Distance_Circle(curr_lon-360.0, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
        else
          distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
        endif

    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
    enddo
    enddo


    do k = 1,2
      kk = k + init_lev - 1
      wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk)  + Weight(1,2) * wind(init_lon,init_lat+1,kk)   &
                     + Weight(2,1) * wind(init_lon+1,init_lat,kk) + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
    enddo


    ! second interpolate vertically (Linear)

    wind_lonlat_lev = wind_lonlat(1)    &
                    + (wind_lonlat(2)-wind_lonlat(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev),
    !P_mid(i_lev+1), curr_pressure )

    Interplt_wind = wind_lonlat_lev

    return
  end function


  ! calculation the great-circle distance between two points on the earth surface
  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re

    Distance_Circle = Re * 2.0 * ASIN( (sin( (y1-y2)*PI/180.0 ))**2.0   & 
                      + cos(x1*PI/180.0) * cos(x2*PI/180.0) * (sin( 0.5*(x1-x2)*PI/180.0 ))**2.0 )

    return
  end function

  
  ! calculate the wind_s shear along pressure direction
  real function Wind_shear(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge,  i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: box_alpha
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: u(:,:,:), v(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:), P_edge(:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: u_lonlat(2), v_lonlat(2)
    real(fp)          :: wind_s(2),  Delt_height


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southeast of the particle or under
    ! the particle
    if(curr_lon>=X_mid(i_lon))then
      init_lon = i_lon
    else
      init_lon = i_lon - 1
    endif

    if(curr_lat>=Y_mid(i_lat))then
      init_lat = i_lat
    else
      init_lat = i_lat - 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif


    do i = 1,2
      do j = 1,2

        ii = i + init_lon - 1
        jj = j + init_lat - 1
        distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))

      enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2

      kk = k + init_lev - 1

      u_lonlat(k) =  Weight(1,1) * u(i_lon,i_lat,kk)   + Weight(1,2) * u(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * u(i_lon+1,i_lat,kk) + Weight(2,2) * u(i_lon+1,i_lat+1,kk)

      v_lonlat(k) =  Weight(1,1) * v(i_lon,i_lat,kk)   + Weight(1,2) * v(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * v(i_lon+1,i_lat,kk) + Weight(2,2) * v(i_lon+1,i_lat+1,kk)


      if(box_alpha>=0.0)then
        wind_s(k) = -1.0*u_lonlat(k)*SIN(box_alpha) + v_lonlat(k)*COS(box_alpha)
      else 
        wind_s(k) = u_lonlat(k)*SIN(box_alpha) + v_lonlat(k)*COS(box_alpha)
      endif

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
    ! Because it is the poressure center in [hPa] instead of height center in [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) + P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    Delt_height =  Pa2m( P_BXHEIGHT(init_lon,init_lat,init_lev), P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2m( P_BXHEIGHT(init_lon,init_lat,init_lev+1), P_edge(init_lev), P_edge(init_lev+1), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    wind_shear = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function


  ! transform the pressure level [Pa] to the height level [m]
  real function Pa2m(Box_height, P1, P2, Judge)  
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge
      
      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2m = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2m = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function

!-------------------------------------------------------------------


!---------------------------------------------------------------------

  subroutine plume_cleanup()

    if (allocated(box_lon))      deallocate(box_lon)
    if (allocated(box_lat))      deallocate(box_lat)
    if (allocated(box_lev))      deallocate(box_lev)
    if (allocated(box_length))   deallocate(box_length)
    if (allocated(box_radius1))  deallocate(box_radius1)
    if (allocated(box_radius2))  deallocate(box_radius2)
    if (allocated(box_theta))    deallocate(box_theta)

    WRITE(6,'(a)') '--> Plume Module Cleanup <--'

  end subroutine plume_cleanup


END MODULE Plume_Mod
