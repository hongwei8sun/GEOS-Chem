!-------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!--------------------------------------------------------------------------

! New version
!
! May 13, 2021
! (1)
! change:         if(RK_lat<0.0)then
! to:             if(box_lat<0.0)then
!
! (2)
!    ! second interpolate vertically (Linear)
!    IF(P_mid(init_lev+1)==P_mid(init_lev))THEN
!      WRITE(6,*)"*** WARNING: two same pressure level happens! ***"
!      WRITE(6,*)"init_lev, P_mid(init_lev), P_mid(init_lev+1):"
!      WRITE(6,*)init_lev, P_mid(init_lev), P_mid(init_lev+1)
!      wind_lonlat_lev = wind_lonlat(1)
!    ELSE
!      wind_lonlat_lev = wind_lonlat(1) + (wind_lonlat(2)-wind_lonlat(1)) &
!                                 / (P_mid(init_lev+1)-P_mid(init_lev)) &
!                                     * (curr_pressure-P_mid(init_lev))
!    ENDIF
!
! (3)
! Once the particle get out of stratosphere, pin the particle in that location
! without moving any longer.
!      !-----------------------------------------------------------------------
!      ! check whether the plume is above the tropopause
!      !-----------------------------------------------------------------------
!      IF(box_lev>State_Met%TROPP(i_lon,i_lat))THEN
!        box_tropp = 0
!        GOTO 400
!      ENDIF


MODULE Trajectory_Mod

  USE precision_mod
  USE ERROR_MOD
  USE ErrCode_Mod
  USE PhysConstants,   ONLY : PI, Re, g0, AIRMW, AVO
!  USE CMN_SIZE_Mod,    ONLY : IIPAR, JJPAR, LLPAR

  USE TIME_MOD,        ONLY : GET_YEAR
  USE TIME_MOD,        ONLY : GET_MONTH
  USE TIME_MOD,        ONLY : GET_DAY
  USE TIME_MOD,        ONLY : GET_HOUR
  USE TIME_MOD,        ONLY : GET_MINUTE
  USE TIME_MOD,        ONLY : GET_SECOND

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: trajectory_init
  PUBLIC :: plume_inject
!  PUBLIC :: lagrange_run
!  PUBLIC :: plume_run
!  PUBLIC :: lagrange_write_std
  PUBLIC :: trajectory_cleanup

  ! PUBLIC VARIABLES:



  PUBLIC :: use_lagrange

  integer               :: use_lagrange = 1


  integer               :: IIPAR, JJPAR, LLPAR
  integer               :: id_PASV_LA3, id_PASV_LA2, id_PASV_LA 
  integer               :: id_PASV_EU2, id_PASV_EU


  real, parameter       :: Inject_lon =   -170.1e+0_fp
  real, parameter       :: Inject_lat = -29.9e+0_fp
  real, parameter       :: Inject_lev =  50.0e+0_fp
  ! 25.0e+0_fp ! [hPa] at about 25 km
  ! 50.0e+0_fp       ! [hPa] at about 20 km


  real(fp), pointer     :: X_mid(:), Y_mid(:), P_mid(:)
  real(fp), pointer     :: P_edge(:)

  real(fp) :: mass_eu, mass_la, mass_la2

  integer               :: N_parcel   ! 131        
  integer               :: Num_inject, Num_Plume2d, Num_dissolve     
  integer               :: tt     
  ! Aircraft would release 131 aerosol parcels every time step

  ! use for plume injection
  integer               :: Total_lon, Total_lat, Total_lev

  integer, parameter    :: i_tracer  = 1
  integer, parameter    :: i_product = 2

  real(fp), parameter       :: Kchem = 1.0e-20_fp ! chemical reaction rate
  ! use Kchem = 1.0e-20_fp for >1 year simulation

  integer               :: Stop_inject ! 1: stop injecting; 0: keep injecting

  TYPE :: Plume2d_list
    integer :: IsNew ! 1: the plume is new injected
    integer :: label ! injected rank
    integer :: AboveTropp ! 1: the plume above the tropopause;
                         ! 0: the plume is below the tropopause;

    real(fp) :: LON, LAT, LEV
    real(fp) :: LIFE


    TYPE(Plume2d_list), POINTER :: next
  END TYPE


  TYPE(Plume2d_list), POINTER :: Plume2d_old, Plume2d_head


CONTAINS


!-----------------------------------------------------------------

  SUBROUTINE trajectory_init(am_I_root, Input_Opt, State_Chm, State_Grid, State_Met, RC)

    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Met_Mod, ONLY : MetState
    USE State_Chm_Mod, ONLY : ChmState

!    USE GC_GRID_MOD,   ONLY : XMID, YMID
!    USE GC_GRID_MOD,   ONLY : GET_AREA_M2 ! new

    USE State_Grid_Mod,  ONLY : GrdState

    USE TIME_MOD,        ONLY : GET_TS_DYN

    USE TIME_MOD,      ONLY : GET_YEAR
    USE TIME_MOD,      ONLY : GET_MONTH
    USE TIME_MOD,      ONLY : GET_DAY
    USE TIME_MOD,      ONLY : GET_HOUR
    USE TIME_MOD,      ONLY : GET_MINUTE
    USE TIME_MOD,      ONLY : GET_SECOND



    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    INTEGER                       :: i_box, i_slab
    INTEGER                       :: ii, jj, kk
    CHARACTER(LEN=255)            :: FILENAME, FILENAME2

    integer :: i_lon, i_lat, i_lev            !1:IIPAR


    REAL(fp) :: lon1, lat1, lon2, lat2
    REAL(fp) :: box_lon_edge, box_lat_edge        


    TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev

    REAL(fp) :: Dt


    Dt = GET_TS_DYN()

    ! use for plume injection
    Total_lev = 3
    Total_lon = 18
    Total_lat = 60*2

!    N_parcel = NINT(132 *1000/Length_init /600*Dt)
    N_parcel = 1

    Stop_inject = 0


    IIPAR = State_Grid%NX
    JJPAR = State_Grid%NY
    LLPAR = State_Grid%NZ



    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Lagrnage Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,*) 'time step=', Dt
    WRITE(6,*) 'Injected plume every time step: ', N_parcel
    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') '--------------------------------------------------------'


    ! -----------------------------------------------------------
    ! create first node (head) for 2d linked list:
    ! -----------------------------------------------------------
    ALLOCATE(Plume2d_old)
    Plume2d_old%IsNew = 1
    Plume2d_old%label = 1

    Plume2d_old%LON = Inject_lon
    Plume2d_old%LAT = Inject_lat
    Plume2d_old%LEV = Inject_lev


    Plume2d_old%LIFE = 0.0e+0_fp
    Plume2d_old%AboveTropp = 1


    ! Here assume the injection rate is 30 kg/km (=30 g/m) for H2SO4
    ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
    mass_la2 = 0.0
    mass_eu  = 0.0

    NULLIFY(Plume2d_old%next)
    Plume2d_head => Plume2d_old



    Num_Plume2d = 1

    ! use this value to set the initial latutude for injected plume
    Num_inject = 1

    Num_dissolve = 0


    X_mid  => State_Grid%XMid(:,1) ! Grid box longitude [degrees] ! XMID(:,1,1)   ! IIPAR ! new
    Y_mid  => State_Grid%YMid(1,:) ! Grid box latitude center [degree] ! YMID(1,:,1)
    P_mid  => State_Met%PMID(1,1,:)  ! Pressure at level centers (hPa)

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]


!--------------------------------------------------------
! Set the initial value of plume character:
! box_theta, box_length,
!  max/min radius for each slab
!--------------------------------------------------------

!    box_theta    = 0.0e+0_fp     ! [radian]

!    DO i_box = 1, n_box-1
!      box_lon_edge = 0.5 * ( box_lon(i_box) + box_lon(i_box+1) )
!      box_lat_edge = 0.5 * ( box_lat(i_box) + box_lat(i_box+1) )
!      box_length(i_box) = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
!                                   box_lon_edge, box_lat_edge)  ! [m]
!    ENDDO
!
!    i_box = n_box
!    box_lon_edge = box_lon(i_box) + 0.5* ( box_lon(i_box) - box_lon(i_box-1) )
!    box_lat_edge = box_lat(i_box) + 0.5* ( box_lat(i_box) - box_lat(i_box-1) )
!    box_length(i_box) = 2.0 * Distance_Circle(box_lon(i_box), box_lat(i_box), &
!                                   box_lon_edge, box_lat_edge)  ! [m]

    !--------------------------------------------------------
    ! Set the initial concentration of injected aerosols
    ! Here assume the injection rate is 30 kg/km 
    ! (same as 30 g/m) for H2SO4
    !--------------------------------------------------------
    ! State_Chm%nAdvect: the last one is PASV
!    DO i_box = 1, n_box
!      box_concnt_2D(:,:,N_species,i_box) = 0.0e+0_fp
!      box_concnt_2D(n_x_mid,n_y_mid,N_species,i_box) = &
!                                box_length(i_box)*30.0 &
!                              / (Pdx(i_box)*Pdy(i_box)*box_length(i_box))
!
!      ! From [g/m3] to [molec/cm3], 98.0 g/mol for H2SO4
!      box_concnt_2D(n_x_mid,n_y_mid,N_species,i_box) = &
!        box_concnt_2D(n_x_mid,n_y_mid,N_species,i_box) &
!          / 1.0e+6_fp / 98.0 * AVO
!
!    ENDDO


    tt   = 0



    !-----------------------------------------------------------------
    ! Output State_Met%AD(i_lon,i_lat,i_lev) into State_Met_AD.txt
    !-----------------------------------------------------------------
    OPEN( 314,      FILE='State_Met_AD.txt', STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
    !!! change xyz order ???
    Do i_lon = 1, IIPAR
    Do i_lat = 1, JJPAR
    Do i_lev = 1, LLPAR
       WRITE(314,'(x,E12.5)') State_Met%AD(i_lon,i_lat,i_lev) ![kg]
    End Do
    End Do
    End Do



    IF(use_lagrange==0)THEN
    ! instantly dissolve injected plume into Eulerian grid 

      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' You are not using Trajectory_Mod now'
      WRITE(6,'(a)') ' set variable use_lagrange = 1 to turn on lagrnage_mod  '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' '

      id_PASV_EU  = State_Chm%nAdvect-1
      id_PASV_EU2 = State_Chm%nAdvect

      ! set initial background concentration of injected aerosol as 0 in GCM
      ! output the apecies' name for double check ???
!      State_Chm%Species(:,:,:,id_PASV_EU2) = 0.0e+0_fp  ! [kg/kg]
!      State_Chm%Species(:,:,:,id_PASV_EU)  = 0.0e+0_fp  ! [kg/kg]

    ELSE

      WRITE(6,'(a)') ' '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' You are using Trajectory_Mod now'
      WRITE(6,'(a)') ' set variable use_lagrange = 0 to turn off lagrnage_mod '
      WRITE(6,'(a)') '********************************************************'
      WRITE(6,'(a)') ' '

      id_PASV_LA  = State_Chm%nAdvect-2
      id_PASV_LA2 = State_Chm%nAdvect-1
      id_PASV_LA3 = State_Chm%nAdvect

!      State_Chm%Species(:,:,:,id_PASV_LA) = 0.0e+0_fp  ! [kg/kg]
!      State_Chm%Species(:,:,:,id_PASV_LA2) = 0.0e+0_fp  ! [kg/kg]
!      State_Chm%Species(:,:,:,id_PASV_LA3)  = 0.0e+0_fp  ! [kg/kg]

      State_Chm%Species(:,:,:,id_PASV_LA) = &
        State_Chm%Species(:,:,:,id_PASV_LA)+State_Chm%Species(:,:,:,id_PASV_LA3)

      State_Chm%Species(:,:,:,id_PASV_LA3)  = 0.0e+0_fp  ! [kg/kg]

    ENDIF



    nullify(Plume2d_new)
    nullify(PLume2d)
    nullify(Plume2d_prev)


  END SUBROUTINE trajectory_init

!=================================================================

  SUBROUTINE plume_inject(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Met_Mod,   ONLY : MetState
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState

    USE TIME_MOD,        ONLY : GET_TS_DYN
    USE UnitConv_Mod,    ONLY : Convert_Spc_Units


    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(OptInput), intent(in)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure


    real(fp), pointer :: X_edge(:), Y_edge(:)
    real(fp)          :: X_edge2, Y_edge2
    real(fp)          :: Dx, Dy


    integer  :: i_box

    integer  :: i_lon, i_lat, i_lev
    real(fp) :: box_lon, box_lat, box_lev

    real(fp), pointer  :: PASV_EU
    integer            :: nAdv
    REAL(fp)           :: MW_g

    integer  :: i_advect1

    REAL(fp) :: Dt

    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg

    ErrMsg = ''


    !=======================================================================
    !
    ! use plume model
    !
    !=======================================================================
    IF(use_lagrange==1)THEN

    ! call the lagrnage_run() and plume_run() to calculate injected plume

      State_Chm%Species(:,:,:,id_PASV_LA3) = 0.0e+0_fp  ! [kg/kg]

      CALL trajectory_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
!      CALL plume_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)
      CALL trajectory_write_std( am_I_Root, RC )

    ENDIF


  END SUBROUTINE

!=================================================================

  SUBROUTINE trajectory_run(am_I_Root, State_Chm, State_Grid, State_Met, Input_Opt, RC)

    USE Input_Opt_Mod, ONLY : OptInput

    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState

    USE TIME_MOD,      ONLY : GET_TS_DYN

!    USE GC_GRID_MOD,   ONLY : XEDGE, YEDGE
              
!    USE CMN_SIZE_Mod,  ONLY : DLAT, DLON !new
    USE State_Grid_Mod,  ONLY : GrdState

    USE UnitConv_Mod,  ONLY : Convert_Spc_Units

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE
    USE TIME_MOD,        ONLY : GET_SECOND


    logical, intent(in)           :: am_I_Root
    TYPE(MetState), intent(in)    :: State_Met
    TYPE(ChmState), intent(inout) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State objectgg
    TYPE(OptInput), intent(in)    :: Input_Opt
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

    REAL(fp) :: Dt, RK_Dt(5)                  ! = 600.0e+0_fp          

    real(fp), pointer  :: PASV_EU           
    integer            :: nAdv        
    REAL(fp)           :: MW_g

    integer :: i_box, n_lon, n_lat, n_lev

    ! 1:IIPAR
    integer :: i_lon, i_lat, i_lev
    integer :: next_i_lon, next_i_lat, next_i_lev

    integer :: ii, jj, kk

    integer :: Ki

    integer :: i_species
    integer :: i_x, i_y
    integer :: i_slab
    integer :: box_tropp


    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: RK_lon, RK_lat
    real(fp) :: X_edge2, Y_edge2

    real(fp) :: dbox_lon, dbox_lat, dbox_lev
    real(fp) :: dbox_x_PS, dbox_y_PS
    real(fp) :: box_x_PS, box_y_PS

    real(fp) :: RK_Dlon(4), RK_Dlat(4), RK_Dlev(4)
    real(fp) :: RK_x_PS, RK_y_PS
    real(fp) :: RK_Dx_PS, RK_Dy_PS
    real(fp) :: RK_u(4), RK_v(4), RK_omeg(4)

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)
    real(fp), pointer :: Ptemp(:,:,:)
    real(fp), pointer :: T1(:,:,:)
    real(fp), pointer :: T2(:,:,:)

    real(fp) :: curr_T1, next_T2, ratio
    real(fp) :: V_prev, V_new

    real(fp) :: box_u, box_v, box_omeg
    real(fp) :: curr_u, curr_v, curr_omeg
    real(fp) :: curr_u_PS, curr_v_PS

    real(fp), pointer :: X_edge(:), Y_edge(:)

    real(fp) :: Dx, Dy

    real(fp) :: length0
    real(fp) :: lon1, lon2, lat1, lat2


    real(fp) :: box_lon_edge, box_lat_edge

!    real(fp) :: box_alpha
    real(fp) :: D_wind, D_x, D_y
    real(fp) :: Ly ! Lyapunov exponent [s-1]


    real(fp) :: box_lon, box_lat, box_lev
    real(fp) :: box_length, box_alpha, box_theta
    real(fp) :: Pdx, Pdy
    real(fp) :: box_Ra, box_Rb
    real(fp) :: box_extra, box_life, box_label


    real(fp)  :: start, finish

    TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev


    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND

    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg
    CHARACTER(LEN=255)     :: FILENAME, FILENAME2

    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND

    IF(mod(tt,24*60*60)==0)THEN   ! output once every day (24 hours)

      FILENAME   = 'Lagrange_xyz_' // TRIM(ADJUSTL(YEAR_C)) // '-'   &
        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'


      OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      CLOSE(261)



      FILENAME2   = 'Lagrange_lifetime_' // TRIM(ADJUSTL(YEAR_C)) // '-'   &
        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'


      OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE',  &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

      CLOSE(262)

    ENDIF ! IF(mod(tt,1440)==0)THEN



    RC     =  GC_SUCCESS
    ErrMsg = ''


    Dt = GET_TS_DYN()

    RK_Dt(1) = 0.0
    RK_Dt(2) = 0.5*Dt
    RK_Dt(3) = 0.5*Dt
    RK_Dt(4) = Dt
    RK_Dt(5) = 0.0

    u     => State_Met%U   ! figure out state_met%U is based on
                           !  lat/lon or modelgrid(i,j)
    v     => State_Met%V   ! V [m s-1]
    omeg  => State_Met%OMEGA     ! Updraft velocity [Pa/s]
    Ptemp => State_Met%THETA     ! Potential temperature [K]
    T1    => State_Met%TMPU1     ! Temperature at start of
                                 !  timestep [K]
    T2    => State_Met%TMPU2     ! Temperature at end of
                                 !  timestep [K]

    Dx = State_Grid%DX
    Dy = State_Grid%DY 
    X_edge => State_Grid%XEdge(:,1) !XEDGE(:,1,1)   ! IIPAR+1 ! new
    Y_edge => State_Grid%YEdge(1,:) !YEDGE(1,:,1)  
    ! Use second YEDGE, because sometimes YMID(2)-YMID(1) is not DLAT

    X_edge2       = X_edge(2)
    Y_edge2       = Y_edge(2)
     
    
    ! -----------------------------------------------------------
    ! add new box every time step
    ! -----------------------------------------------------------

    IF(Stop_inject==0)THEN

    DO n_lev = 1, Total_lev, 1
    DO n_lon = 1, Total_lon, 1
    DO n_lat = 1, Total_lat, 1

      i_box = Num_inject

      ALLOCATE(Plume2d_new)

      Plume2d_new%IsNew = 1
      Plume2d_new%label = i_box
      Plume2d_new%AboveTropp = 1

      Plume2d_new%LON = Inject_lon + (n_lon-1) * 20
      Plume2d_new%LAT = Inject_lat + (n_lat-1) * 0.5
!      Plume2d_new%LAT = ( -30.005e+0_fp + 0.01e+0_fp * MOD(i_box,6000) ) &
!                     * (-1.0)**FLOOR(i_box/6000.0) ! -29.995S:29.995N:0.01
      Plume2d_new%LEV = Inject_lev + (n_lev-1)*20


      Plume2d_new%LIFE = 0.0e+0_fp




!      WRITE(6,*) 'mass_la: ', Plume2d_new%CONCNT2d(n_x_mid,n_y_mid)* &
!             (Plume2d_new%DX*Plume2d_new%DY*Plume2d_new%LENGTH*1.0e+6_fp)

      NULLIFY(Plume2d_new%next)
      Plume2d_old%next => Plume2d_new
      Plume2d_old => Plume2d_old%next



      Num_Plume2d = Num_Plume2d + N_parcel

      ! use this value to set the initial latutude for injected plume
      Num_inject = Num_inject + N_parcel

    ENDDO ! DO n_lat = 1, Total_lat, 1
    ENDDO ! DO n_lon = 1, Total_lon, 1
    ENDDO ! DO n_lon = 1, Total_lev, 1

    ENDIF ! IF(Stop_inject==0)THEN

    Stop_inject = 1

    !=======================================================================
    ! for 2D plume: Run Lagrangian trajectory-track HERE
    !=======================================================================

    IF(.NOT.ASSOCIATED(Plume2d_head)) GOTO 401

    Plume2d => Plume2d_head
    i_box = 0

    DO WHILE(ASSOCIATED(Plume2d))

      i_box = i_box+1

      box_lon    = Plume2d%LON
      box_lat    = Plume2d%LAT
      box_lev    = Plume2d%LEV

      box_label  = Plume2d%label
      box_life   = Plume2d%LIFE
      box_tropp   = Plume2d%AboveTropp



      box_life = box_life + Dt


      !-----------------------------------------------------------------------
      ! For the center of the plume
      !-----------------------------------------------------------------------

      ! make sure the location is not out of range
      do while (box_lat > Y_edge(JJPAR+1))
         box_lat = Y_edge(JJPAR+1) - ( box_lat-Y_edge(JJPAR+1) )
      end do

      do while (box_lat < Y_edge(1))
         box_lat = Y_edge(1) + ( box_lat-Y_edge(1) )
      end do

      do while (box_lon > X_edge(IIPAR+1))
         box_lon = box_lon - 360.0
      end do

      do while (box_lon < X_edge(1))
         box_lon = box_lon + 360.0
      end do


      curr_lon      = box_lon
      curr_lat      = box_lat
      curr_pressure = box_lev      ! hPa
!      curr_T1       = box_T1(i_box)        ! K


      ! check this grid box is the neast one or the one located in left-bottom
      ! ???
      i_lon = Find_iLonLat(curr_lon, Dx, X_edge2)
      if(i_lon>IIPAR) i_lon=i_lon-IIPAR
      if(i_lon<1) i_lon=i_lon+IIPAR

      i_lat = Find_iLonLat(curr_lat, Dy, Y_edge2)
      if(i_lat>JJPAR) i_lat=JJPAR
      if(i_lat<1) i_lat=1

      i_lev = Find_iPLev(curr_pressure,P_edge)
      if(i_lev>LLPAR) i_lev=LLPAR
      IF(i_lev==LLPAR) WRITE(6,*) 'box_lev1:', curr_pressure, i_lev

      IF(i_lev<1) i_lev=1
      IF(i_lev==LLPAR) WRITE(6,*) 'box_lev2:', curr_pressure, i_lev

      !-----------------------------------------------------------------------
      ! check whether the plume is above the tropopause
      !-----------------------------------------------------------------------
      IF(box_lev>State_Met%TROPP(i_lon,i_lat))THEN
        box_tropp = 0
        GOTO 400
      ENDIF


      DO Ki = 1,4,1

       !------------------------------------------------------------------
       ! For vertical wind speed:
       ! pay attention for the polar region * * *
       !------------------------------------------------------------------
       if(abs(curr_lat)>Y_mid(JJPAR))then
          curr_omeg = Interplt_wind_RLL_polar(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       else
          curr_omeg = Interplt_wind_RLL(omeg, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
       endif

       RK_omeg(Ki) = curr_omeg

       ! For 
!       dbox_lev = RK_Dt(Ki) * curr_omeg / 100.0     ! Pa => hPa
!       curr_pressure = box_lev(i_box) + dbox_lev

       RK_Dlev(Ki)   = Dt * curr_omeg / 100.0     ! Pa => hPa
       curr_pressure = box_lev + RK_Dt(Ki+1) * curr_omeg / 100.0

       if(curr_pressure<P_mid(LLPAR)) &
             curr_pressure = P_mid(LLPAR) !+ ( P_mid(LLPAR) - curr_pressure )
       if(curr_pressure>P_mid(1)) &
             curr_pressure = P_mid(1) !- ( curr_pressure - P_mid(1) )


       !------------------------------------------------------------------
       ! For the region where lat<72, use Regualr Longitude-Latitude Mesh:
       !------------------------------------------------------------------
       if(abs(box_lat)<=72.0)then

         curr_u = Interplt_wind_RLL(u, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)
         curr_v = Interplt_wind_RLL(v, i_lon, i_lat, i_lev, curr_lon, &
                                                  curr_lat, curr_pressure)

         RK_u(Ki) = curr_u
         RK_v(Ki) = curr_v

         dbox_lon  = (RK_Dt(Ki+1)*curr_u) &
                        / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         dbox_lat  = (RK_Dt(Ki+1)*curr_v) / (PI*Re) * 180.0

         curr_lon  = box_lon + dbox_lon
         curr_lat  = box_lat + dbox_lat
 

         RK_Dlon(Ki) = (Dt*curr_u) &
                      / (2.0*PI*Re*COS(box_lat*PI/180.0)) * 360.0
         RK_Dlat(Ki) = (Dt*curr_v) / (PI*Re) * 180.0

       endif

       !------------------------------------------------------------------
       ! For the polar region (lat>=72), use polar sterographic
       !------------------------------------------------------------------
       if(abs(box_lat)>72.0)then

         if(abs(curr_lat)>Y_mid(JJPAR))then 
         curr_u_PS = Interplt_uv_PS_polar(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         curr_v_PS = Interplt_uv_PS_polar(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         else
         curr_u_PS = Interplt_uv_PS(1, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         curr_v_PS = Interplt_uv_PS(0, u, v, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)    
         endif

         RK_u(Ki) = curr_u_PS
         RK_v(Ki) = curr_v_PS

         dbox_x_PS = RK_Dt(Ki+1)*curr_u_PS
         dbox_y_PS = RK_Dt(Ki+1)*curr_v_PS


         RK_Dx_PS = Dt*curr_u_PS
         RK_Dy_PS = Dt*curr_v_PS


         !------------------------------------------------------------------
         ! change from (lon,lat) in RLL to (x,y) in PS: 
         !------------------------------------------------------------------
         if(box_lat<0)then
           box_x_PS = -1.0* Re* COS(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
           box_y_PS = -1.0* Re* SIN(box_lon*PI/180.0) &
                                / TAN(box_lat*PI/180.0)
         else
           box_x_PS = Re* COS(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
           box_y_PS = Re* SIN(box_lon*PI/180.0) &
                        / TAN(box_lat*PI/180.0)
         endif

         RK_x_PS  = box_x_PS + RK_Dx_PS
         RK_y_PS  = box_y_PS + RK_Dy_PS

         box_x_PS  = box_x_PS + dbox_x_PS
         box_y_PS  = box_y_PS + dbox_y_PS


         !------------------------------------------------------------------
         ! change from (x,y) in PS to (lon,lat) in RLL
         !------------------------------------------------------------------
         if(box_x_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI 
         endif
         if(box_x_PS<0.0 .and. box_y_PS<=0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI -180.0
         endif
         if(box_x_PS<0.0 .and. box_y_PS>0.0)then
           curr_lon = ATAN( box_y_PS / box_x_PS )*180.0/PI +180.0
         endif
           
         if(curr_lat<0.0)then
           curr_lat= -1* ATAN( Re/SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         else
           curr_lat= ATAN( Re / SQRT(box_x_PS**2+box_y_PS**2) ) *180.0/PI
         endif

         !------------------------------------------------------------------
         ! For 4th order Runge Kutta
         !------------------------------------------------------------------
         if(RK_x_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS<=0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI -180.0
         endif
         if(RK_x_PS<0.0 .and. RK_y_PS>0.0)then
           RK_lon = ATAN( RK_y_PS / RK_x_PS )*180.0/PI +180.0
         endif

         if(box_lat<0.0)then
           RK_lat = -1 * ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         else
           RK_lat = ATAN( Re / SQRT(RK_x_PS**2+RK_y_PS**2) ) *180.0/PI
         endif

         RK_Dlon(Ki) = RK_lon - box_lon
         RK_Dlat(Ki) = RK_lat - box_lat

       endif ! if(abs(curr_lat)>72.0)then

      ENDDO ! Ki = 1,4,1

      box_lon = box_lon + &
                (RK_Dlon(1)+2.0*RK_Dlon(2)+2.0*RK_Dlon(3)+RK_Dlon(4))/6.0
      box_lat = box_lat + &
                (RK_Dlat(1)+2.0*RK_Dlat(2)+2.0*RK_Dlat(3)+RK_Dlat(4))/6.0
      box_lev = box_lev + &
                (RK_Dlev(1)+2.0*RK_Dlev(2)+2.0*RK_Dlev(3)+RK_Dlev(4))/6.0



      box_u    = ( RK_u(1) + 2.0*RK_u(2) &
                         + 2.0*RK_u(3) + RK_u(4) ) / 6.0
      box_v    = ( RK_v(1) + 2.0*RK_v(2) &
                         + 2.0*RK_v(3) + RK_v(4) ) / 6.0
      box_omeg = ( RK_omeg(1) + 2.0*RK_omeg(2) &
                         + 2.0*RK_omeg(3) + RK_omeg(4) ) / 6.0

400 CONTINUE

      !-------------------------------------------------------------------
      ! Output box location
      !-------------------------------------------------------------------
      IF(mod(tt,24*60*60)==0)THEN   ! output once every day (24 hours)

        OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        WRITE(261,*) box_lon, box_lat, box_lev

        CLOSE(261)


        OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD',  &
          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

        WRITE(262,*) box_label, box_life, box_tropp

        CLOSE(262)

      ENDIF


      ! ----------------------------------------------------------------
      ! update     
      ! ----------------------------------------------------------------
      Plume2d%LON     = box_lon
      Plume2d%LAT     = box_lat
      Plume2d%LEV     = box_lev
      Plume2d%label   = box_label
      Plume2d%LIFE    = box_life
      Plume2d%AboveTropp    = box_tropp


      Plume2d => Plume2d%next

    ENDDO  ! DO WHILE(ASSOCIATED(Plume2d))

401 CONTINUE

!    WRITE(6,*)'=== loop for 2d plume in lagrange_run: ', i_box, Num_Plume2d



    !------------------------------------------------------------------
    ! Everything is done, clean up pointers
    !------------------------------------------------------------------


    IF(ASSOCIATED(u)) nullify(u)
    IF(ASSOCIATED(v)) nullify(v)
    IF(ASSOCIATED(omeg)) nullify(omeg)

    IF(ASSOCIATED(Ptemp)) nullify(Ptemp)
    IF(ASSOCIATED(T1)) nullify(T1)
    IF(ASSOCIATED(T2)) nullify(T2)

    IF(ASSOCIATED(X_edge)) nullify(X_edge)
    IF(ASSOCIATED(Y_edge)) nullify(Y_edge)

    IF(ASSOCIATED(PASV_EU)) nullify(PASV_EU)

    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(PLume2d)) nullify(PLume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)

  END SUBROUTINE trajectory_run


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real(fp) function Interplt_wind_RLL(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev

    ! Identify wether particle is exactly located on the grid point
    if(curr_pressure==P_mid(i_lev))then
    if(curr_lon==X_mid(i_lon))then
    if(curr_lat==Y_mid(i_lat))then

          Interplt_wind_RLL = wind(i_lon, i_lat, i_lev)
          return

    endif
    endif
    endif


    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo


    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2

       if(distance(i,j)==0)then
          Weight(:,:) = 0
          Weight(i,j) = 1
          GOTO 100
       endif

       Weight(i,j) = 1.0/distance(i,j) / SUM( 1.0/distance(:,:) )

    enddo
    enddo

 100 CONTINUE


    do k = 1,2     
        kk = k + init_lev - 1 
        if(init_lon==0)then
            wind_lonlat(k) =  Weight(1,1) * wind(IIPAR,init_lat,kk) &
                          + Weight(1,2) * wind(IIPAR,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        else if(init_lon==IIPAR)then
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(1,init_lat,kk)   &
                          + Weight(2,2) * wind(1,init_lat+1,kk)
        else
            wind_lonlat(k) =  Weight(1,1) * wind(init_lon,init_lat,kk) &
                          + Weight(1,2) * wind(init_lon,init_lat+1,kk) &
                          + Weight(2,1) * wind(init_lon+1,init_lat,kk) &
                          + Weight(2,2) * wind(init_lon+1,init_lat+1,kk)
        endif
    enddo


    ! second interpolate vertically (Linear)
    IF(P_mid(init_lev+1)==P_mid(init_lev))THEN
      WRITE(6,*)"*** WARNING: two same pressure level happens! ***"
      WRITE(6,*)"init_lev, P_mid(init_lev), P_mid(init_lev+1):"
      WRITE(6,*)init_lev, P_mid(init_lev), P_mid(init_lev+1)
      wind_lonlat_lev = wind_lonlat(1)
    ELSE
      wind_lonlat_lev = wind_lonlat(1) + (wind_lonlat(2)-wind_lonlat(1)) &
                                 / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                     * (curr_pressure-P_mid(init_lev))
    ENDIF

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev), P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL = wind_lonlat_lev


    return
  end function


! functions to interpolate vertical wind speed (w)
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The w value at polar point is the average of all surrounding grid points.

  real(fp) function Interplt_wind_RLL_polar(wind, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: wind(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(3), Weight(3)
    real(fp)          :: wind_lonlat(2), wind_lonlat_lev
    real(fp)          :: wind_polar

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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
    
    if(init_lat==0)then
      init_lat = 1
    endif

    ! For pressure level, P_mid(1) is about surface pressure
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! calculate the distance between particle and grid point
    j = 1
    do i = 1,2
      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! For some special circumstance:
      if(ii==0)then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(1), Y_mid(jj))
      else
      distance(i)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo


   if(ii==0)then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), 90.0e+0_fp)
   else if(ii==(IIPAR+1))then
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(1), 90.0e+0_fp)
   else
    distance(3)= Distance_Circle(curr_lon, curr_lat, X_mid(ii), 90.0e+0_fp)
   endif


   IF(distance(3)==0.0)THEN
       do k=1,2
         kk = k + init_lev - 1
         wind_lonlat(k) = SUM(wind(:,init_lat,kk))/IIPAR
       enddo
   ELSE
       ! Calculate the inverse distance weight
       do i=1,3
          Weight(i) = 1.0/distance(i) / sum( 1.0/distance(:) )
       enddo

       do k=1,2
         kk = k + init_lev - 1

         wind_polar = SUM(wind(:,init_lat,kk))/IIPAR      

         if(init_lon==0)then    
             wind_lonlat(k) =  Weight(1)*wind(IIPAR,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else if(init_lon==IIPAR)then
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         else
             wind_lonlat(k) =  Weight(1)*wind(init_lon,init_lat,kk)   &
                             + Weight(2)*wind(init_lon+1,init_lat,kk)   &
                             + Weight(3)*wind_polar
         endif
     enddo
   ENDIF

    ! second interpolate vertically (Linear)
    wind_lonlat_lev =   wind_lonlat(1) &
                      + (wind_lonlat(2)-wind_lonlat(1)) &
                      / (P_mid(init_lev+1)-P_mid(init_lev)) &
                      * (curr_pressure-P_mid(init_lev))

    !Line_Interplt( wind_lonlat(1), wind_lonlat(2), P_mid(i_lev),
    !P_mid(i_lev+1), curr_pressure )

    Interplt_wind_RLL_polar = wind_lonlat_lev

    return
  end function

!------------------------------------------------------------------
! functions to interpolate wind speed (u,v) 
! based on the surrounding 3 points, one of the points is the north/south polar
! point. The u/v value at polar point is the average of all surrounding grid points.

  real(fp) function Interplt_uv_PS_polar(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(3), y_PS(3)  ! the third value x_PS(3) is the polar point
    real(fp)          :: uv_PS(3,2)
    real(fp)          :: uv_polars(IIPAR)
    real(fp)          :: distance_PS(3), Weight_PS(3)
    real(fp)          :: uv_xy(2), uv_xy_lev

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
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

    ! For pressure level, P_mid(1) is about surface pressure, has biggerst
    ! value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1


    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    endif

    ! i,j means the four grid point value that surround around the particle
!    do i=1,2
!    do j=1,2

!      ii = i + init_lon - 1
!      jj = j + init_lat - 1
!
      ! For lon=180 deg:
!      if(ii==IIPAR+1)then
!        ii = 1
!      endif
!      if(ii==0)then
!        ii = IIPAR
!      endif

      ! Get the right ii and jj for interpolation at polar point:

      ! For South Polar Point:
 !     if(jj==0)then
 !       jj = jj+1
 !     endif

      ! For North Polar Point:
!      if(jj==JJPAR+1)then
!        jj = jj-1
!      endif

    if(init_lat==0)then
           jj = init_lat + 1    ! south polar
       else
           jj = init_lat
    endif
    
    do i=1,2  ! Interpolate location and wind of grid points into Polar Stereographic Plane
   
       ii = i + init_lon - 1

       ! For lon=180 deg:
       if(ii==IIPAR+1)then
          ii = 1
       endif
       if(ii==0)then
          ii = IIPAR
       endif

      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then
        x_PS(i)= Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      else
        x_PS(i)= -1.0* Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
        y_PS(i)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
      endif


       do k=1,2
          kk = k + init_lev - 1
          IF(i_uv==1)THEN ! i_ux==1 for u
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_PS(i,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                        + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF

          IF(i_uv==0)THEN ! for v
          if(Y_mid(jj)>0)then
            uv_PS(i,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                       - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          else
            uv_PS(i,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                           + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
          ENDIF
       enddo
    enddo

    ! Third grid point: south/north polar point
    x_PS(3) = 0.0
    y_PS(3) = 0.0
    
    do k=1,2
      kk = k + init_lev - 1
      IF(i_uv==1)THEN ! i_ux==1 for u
      ! interpolate all the grid points surrounding the polar point:
      do ii = 1,IIPAR
        if(Y_mid(jj)>0)then
          uv_polars(ii) = -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                         + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                         / (SIN(Y_mid(jj)*PI/180.0)**2) )
          else
            uv_polars(ii) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
          endif
        enddo

            uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF

          IF(i_uv==0)THEN ! for v
          do ii = 1,IIPAR
             if(Y_mid(jj)>0)then
               uv_polars(ii) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                             - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             else
               uv_polars(ii) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                                 + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
             endif
          enddo
             uv_PS(3,k) = SUM(uv_polars)/IIPAR
          ENDIF
    enddo



    ! calculate the distance between particle and grid point
    do i = 1,3
       distance_PS(i)= SQRT( (x_PS(i)-curr_x)**2.0 + (y_PS(i)-curr_y)**2.0 )
    enddo

    ! Calculate the inverse distance weight
    do i = 1,3
        Weight_PS(i) = 1.0/distance_PS(i) / sum( 1.0/distance_PS(:) )
    enddo


    do k = 1,2
      uv_xy(k) = Weight_PS(1) * uv_PS(1,k) &
               + Weight_PS(2) * uv_PS(2,k) &
               + Weight_PS(3) * uv_PS(3,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ (uv_xy(2)-uv_xy(1)) &
     / (P_mid(init_lev+1)-P_mid(init_lev)) * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS_polar = uv_xy_lev

    return
  end function


!------------------------------------------------------------------
! functions to interpolate wind speed (u,v,omeg) 
! based on the surrounding 4 points.

  real(fp) function Interplt_uv_PS(i_uv, u_RLL, v_RLL, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)

    implicit none

    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: curr_x, curr_y
    real(fp), pointer :: u_RLL(:,:,:), v_RLL(:,:,:)

    integer           :: i_uv
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk

    real(fp)          :: x_PS(2,2), y_PS(2,2)
    real(fp)          :: uv_PS(2,2,2)
    real(fp)          :: distance_PS(2,2), Weight_PS(2,2)
    real(fp)          :: uv_xy(2), uv_xy_lev


    ! first interpolate horizontally (Inverse Distance Weighting)

 ! identify the grid point located in the southwest of the particle or under
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

!For pressure level, P_mid(1) is about surface pressure, has biggerst value.
    if(curr_pressure<=P_mid(i_lev))then
      init_lev = i_lev
    else
      init_lev = i_lev - 1
    endif

    if(init_lev==0) init_lev = 1
    if(init_lev==LLPAR) init_lev = LLPAR-1

    
    ! change from (lon,lat) in RLL to (x,y) in PS: 
    if(curr_lat<0)then
      curr_x = -1.0* Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = -1.0* Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    else
      curr_x = Re* COS(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
      curr_y = Re* SIN(curr_lon*PI/180.0) / TAN(curr_lat*PI/180.0)
    endif

    ! i,j means the four grid point value that surround around the particle
    do i=1,2
    do j=1,2

      ii = i + init_lon - 1
      jj = j + init_lat - 1

      ! Get the right ii and jj for interpolation at polar point:
      ! For South Polar Point:
      if(jj==0)then
        jj = jj+1
      endif
      ! For North Polar Point:
      if(jj==JJPAR+1)then
        jj = jj-1
      endif

    
      ! For lon=180 deg:
      if(ii==IIPAR+1)then
        ii = 1
      endif
      if(ii==0)then
        ii = IIPAR
      endif


      ! Interpolate location and wind into Polar Stereographic Plane
      if(Y_mid(jj)>0)then

        x_PS(i,j) = Re* COS(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)  
        y_PS(i,j) = Re* SIN(X_mid(ii)*PI/180.0) / TAN(Y_mid(jj)*PI/180.0)
          
        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then ! i_ux==1 for u
             uv_PS(i,j,k)= -1.0* ( u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                         / SIN(Y_mid(jj)*PI/180.0) &
                                + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2) )
           endif

           if(i_uv==0)then ! for v
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          - v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      else

      x_PS(i,j)= -1.0* Re* COS(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)
      y_PS(i,j)= -1.0* Re* SIN(X_mid(ii)*PI/180.0) /TAN(Y_mid(jj)*PI/180.0)

        do k=1,2
           kk = k + init_lev - 1

           if(i_uv==1)then
             uv_PS(i,j,k) = u_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                          + v_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif

           if(i_uv==0)then
             uv_PS(i,j,k) = -1* u_RLL(ii,jj,kk)*COS(X_mid(ii)*PI/180.0) &
                                        / SIN(Y_mid(jj)*PI/180.0) &
                              + v_RLL(ii,jj,kk)*SIN(X_mid(ii)*PI/180.0) &
                                        / (SIN(Y_mid(jj)*PI/180.0)**2)
           endif
        enddo

      endif

    enddo
    enddo

    ! calculate the distance between particle and grid point
    do i = 1,2
    do j = 1,2
       distance_PS(i,j) = SQRT( (x_PS(i,j)-curr_x)**2.0 &
                                                + (y_PS(i,j)-curr_y)**2.0 )
    enddo
    enddo

    ! Calculate the inverse distance weight
    do i = 1,2
    do j = 1,2
        Weight_PS(i,j) = 1.0/distance_PS(i,j) / sum( 1.0/distance_PS(:,:) )
    enddo
    enddo


    do k = 1,2
      uv_xy(k) =  Weight_PS(1,1) * uv_PS(1,1,k) &
                 + Weight_PS(1,2) * uv_PS(1,2,k) &
                 + Weight_PS(2,1) * uv_PS(2,1,k) &
                 + Weight_PS(2,2) * uv_PS(2,2,k)
    enddo


    ! second interpolate vertically (Linear)

    uv_xy_lev = uv_xy(1)+ &
              (uv_xy(2)-uv_xy(1)) / (P_mid(init_lev+1)-P_mid(init_lev)) &
                                        * (curr_pressure-P_mid(init_lev))

    Interplt_uv_PS = uv_xy_lev

    return
  end function


!---------------------------------------------------------------------
!*********************************************************************
!
!    FILENAME2   = 'Plume_lifetime_seconds.txt'
!
!    OPEN( 484,      FILE=TRIM( FILENAME2   ), STATUS='OLD',  &
!          POSITION='APPEND', FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!
!--------------------------------------------------------------------
! functions to which grid cell (i,j,k) contains the box
!--------------------------------------------------------------------

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


  !-------------------------------------------------------------------
  ! calculation the great-circle distance between two points on the earth
  ! surface
  !-------------------------------------------------------------------

  real(fp) function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2      ! unit is degree
    real(fp)     :: xx1, yy1, xx2, yy2  ! unit is radian
    !real(fp) :: PI, Re

    xx1 = x1/180.0*PI
    yy1 = y1/180.0*PI
    xx2 = x2/180.0*PI
    yy2 = y2/180.0*PI

!     Distance_Circle = Re * 
!       ATAN( SQRT( (COS(y2)*SIN(x2-x1))**2 + &
!              ( COS(y1)*SIN(y2)-SIN(y1)*OCS(y2)COS((x2-x1)) )**2 ) &
!            / (SIN(y1)*SIN(y2)+COS(y1)*COS(y2)*COS(x2-x1)) )

    ! output distance is in unit of [m]
    Distance_Circle = Re * 2.0 * ASIN(SQRT( (SIN((yy1-yy2)*0.5))**2.0 &
                     +COS(yy1)*COS(yy2)*(SIN((xx1-xx2)*0.5))**2.0 ))
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure
  ! direction
  !-------------------------------------------------------------------  

  real(fp) function Wind_shear_s(u, v, P_BXHEIGHT, plume_alpha, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    real(fp)          :: plume_alpha
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: u(:,:,:), v(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: u_lonlat(2), v_lonlat(2)
    real(fp)          :: wind_s(2),  Delt_height

    ! first interpolate horizontally (Inverse Distance Weighting)

    ! identify the grid point located in the southwest of the particle or under
    ! the particle ??? 
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

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
           Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2

      kk = k + init_lev - 1

      if(init_lon==0)then
          u_lonlat(k) =  Weight(1,1) * u(IIPAR,init_lat,kk) &
                       + Weight(1,2) * u(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(IIPAR,init_lat,kk) &
                       + Weight(1,2) * v(IIPAR,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      else if(init_lon==IIPAR)then
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(1,init_lat,kk)   &
                       + Weight(2,2) * u(1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(1,init_lat,kk)   &
                       + Weight(2,2) * v(1,init_lat+1,kk)
      else
          u_lonlat(k) =  Weight(1,1) * u(init_lon,init_lat,kk) &
                       + Weight(1,2) * u(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * u(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * u(init_lon+1,init_lat+1,kk)
          v_lonlat(k) =  Weight(1,1) * v(init_lon,init_lat,kk) &
                       + Weight(1,2) * v(init_lon,init_lat+1,kk) &
                       + Weight(2,1) * v(init_lon+1,init_lat,kk) &
                       + Weight(2,2) * v(init_lon+1,init_lat+1,kk)
      endif

      wind_s(k) = u_lonlat(k)*COS(plume_alpha-0.5*PI) + v_lonlat(k)*COS(plume_alpha-PI)

    enddo


    ! second vertical shear of wind_s

    ! This code should be changed !!!
   ! Because it is the pressure center in [hPa] instead of height center in
    ! [m]
    ! Delt_height    = 0.5 * ( P_BXHEIGHT(init_lon,init_lat,init_lev) +
    ! P_BXHEIGHT(init_lon,init_lat,init_lev+1) )
    if(init_lon==0) init_lon=IIPAR
    if(init_lon==IIPAR+1) init_lon=1

    Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                          P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                          P_edge(init_lev), P_edge(init_lev+1), 0 )


    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function

!-------------------------------------------------------------------
! transform the pressure level [Pa] to the height level [m]
!-------------------------------------------------------------------

  real(fp) function Pa2meter(Box_height, P1, P2, Judge)
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge

      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2meter = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2meter = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) &
                        / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function


!-------------------------------------------------------------------
! Calculate eddy diffusivity in stratisphere. (U.Schumann, 2012)


  !-------------------------------------------------------------------  
  ! calculate the vertical shear for calculating eddy difussivity

  real(fp) function Vertical_shear(var, P_BXHEIGHT, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: var(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    integer           :: i_lon, i_lat, i_lev
    integer           :: init_lon, init_lat, init_lev
    integer           :: i, ii, j, jj, k, kk
    real(fp)          :: distance(2,2), Weight(2,2)
    real(fp)          :: var_lonlat(2)
    real(fp)          :: Delt_height

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

      ! For some special circumstance:
      if(ii==0)then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii+IIPAR), Y_mid(jj))
      else if(ii==(IIPAR+1))then
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii-IIPAR), Y_mid(jj))
      else
        distance(i,j) = &
             Distance_Circle(curr_lon, curr_lat, X_mid(ii), Y_mid(jj))
      endif

    enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2
      kk            = k + init_lev - 1
      
      IF(init_lon==0)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE IF(init_lon==IIPAR)THEN
        var_lonlat(k) =  Weight(1,1) *var(IIPAR,i_lat,kk)   &
                       + Weight(1,2) *var(IIPAR,i_lat+1,kk)   &
                       + Weight(2,1) *var(1,i_lat,kk) &
                       + Weight(2,2) *var(1,i_lat+1,kk)
      ELSE
        var_lonlat(k) =  Weight(1,1) *var(init_lon,i_lat,kk)   &
                       + Weight(1,2) *var(init_lon,i_lat+1,kk)   &
                       + Weight(2,1) *var(init_lon+1,i_lat,kk) &
                       + Weight(2,2) *var(init_lon+1,i_lat+1,kk)
      ENDIF
    enddo


    ! second vertical shear of wind
    if(init_lon==0)then
      Delt_height = Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev),    &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(IIPAR,init_lat,init_lev+1),   &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else if(init_lon==IIPAR)then
      Delt_height = Pa2meter( P_BXHEIGHT(1,init_lat,init_lev),        &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &
                 + Pa2meter( P_BXHEIGHT(1,init_lat,init_lev+1),       &
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    else
      Delt_height = Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), &
                            P_edge(init_lev), P_edge(init_lev+1), 1 ) &   
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),&
                            P_edge(init_lev), P_edge(init_lev+1), 0 )
    endif
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!===================================================================
!
!===================================================================



!  REAL(fp) FUNCTION Find_theta(theta1, concnt1_2D, Height, Ibox)
!    ! n_x_max, n_y_max, Pdx, Pdy are global variables
!
!    IMPLICIT NONE
!
!    REAL(fp)    :: theta1, concnt1_2D(n_x_max, n_y_max), Height    
!    REAL(fp)    :: R
!    REAL(fp)    :: test(100), X(100), Y(100), Conc(100)
!
!    REAL(fp)    :: angle(1)
!
!    INTEGER     :: i, Ibox
!
!      R = Height/COS(theta1) ! 2.5 is a random number
!
!
!      DO i= 1, 100
!        test(i) = PI/2 - (PI/2-theta1)*2/100*i
!        Y(i)    = R *COS( test(i) )
!        X(i)    = R *SIN( test(i) )
!        Conc(i) = Interplt_2D(X(i), Y(i), concnt1_2D, Ibox, 1)
!      ENDDO
!
!
!      angle      = test( MAXLOC(Conc(:)) )
!      Find_theta = angle(1)
!
!    return 
!
!  END FUNCTION


  REAL(fp) FUNCTION Interplt_linear(xx, x1, x2, Cx1, Cx2)

    IMPLICIT NONE

    REAL(fp)    :: xx, x1, x2, Cx1, Cx2


    Interplt_linear = ( Cx1*(x2-xx) + Cx2*(xx-x1) ) /(x2-x1)

    return

  END FUNCTION



  SUBROUTINE trajectory_write_std( am_I_Root, RC )
        
    USE m_netCDF_io_define
    USE m_netcdf_io_read
    USE m_netcdf_io_open
    USE Ncdf_Mod,            ONLY : NC_Open
    USE Ncdf_Mod,            ONLY : NC_Read_Time
    USE Ncdf_Mod,            ONLY : NC_Read_Arr
    USE Ncdf_Mod,            ONLY : NC_Create
    USE Ncdf_Mod,            ONLY : NC_Close
    USE Ncdf_Mod,            ONLY : NC_Var_Def
    USE Ncdf_Mod,            ONLY : NC_Var_Write
    USE Ncdf_Mod,            ONLY : NC_Get_RefDateTime
    USE CHARPAK_Mod,         ONLY : TRANLC
    USE JulDay_Mod,          ONLY : JulDay

    USE TIME_MOD,        ONLY : GET_TS_DYN

    USE Input_Opt_Mod, ONLY : OptInput

    USE TIME_MOD,        ONLY : GET_YEAR
    USE TIME_MOD,        ONLY : GET_MONTH
    USE TIME_MOD,        ONLY : GET_DAY
    USE TIME_MOD,        ONLY : GET_HOUR
    USE TIME_MOD,        ONLY : GET_MINUTE
    USE TIME_MOD,        ONLY : GET_SECOND

    ! Parameters for netCDF routines
    include "netcdf.inc"


    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success

    INTEGER :: i_box

    INTEGER :: YEAR
    INTEGER :: MONTH
    INTEGER :: DAY
    INTEGER :: HOUR
    INTEGER :: MINUTE
    INTEGER :: SECOND


    CHARACTER(LEN=25) :: YEAR_C
    CHARACTER(LEN=25) :: MONTH_C
    CHARACTER(LEN=25) :: DAY_C
    CHARACTER(LEN=25) :: HOUR_C
    CHARACTER(LEN=25) :: MINUTE_C
    CHARACTER(LEN=25) :: SECOND_C


    TYPE(Plume2d_list), POINTER :: Plume2d_new, PLume2d, Plume2d_prev

    REAL(fp) :: Dt

    YEAR        = GET_YEAR()
    MONTH       = GET_MONTH()
    DAY         = GET_DAY()
    HOUR        = GET_HOUR()
    MINUTE      = GET_MINUTE()
    SECOND      = GET_SECOND()

    WRITE(YEAR_C,*) YEAR
    WRITE(MONTH_C,*) MONTH
    WRITE(DAY_C,*) DAY
    WRITE(HOUR_C,*) HOUR
    WRITE(MINUTE_C,*) MINUTE
    WRITE(SECOND_C,*) SECOND


    Dt = GET_TS_DYN()

    !-------------------------------------------------------------------
    ! Output box location
    !-------------------------------------------------------------------

!    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)
!
!    FILENAME   = 'Lagrange_xyz_' // TRIM(ADJUSTL(YEAR_C)) // '-'   &
!        //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) // '-' &
!        // TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) &
!        // ':' // TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!
!    OPEN( 261,      FILE=TRIM( FILENAME   ), STATUS='REPLACE',  &
!          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!
!    WRITE(261,*) Num_Plume2d, Num_dissolve, Num_inject
!
!
!    ENDIF ! IF(mod(tt,1440)==0)THEN


    !-------------------------------------------------------------------
    ! Output plume location
    !-------------------------------------------------------------------
!    IF(mod(tt,6)==0)THEN     ! output once every hour
!    IF(mod(tt,14400000)==0)THEN   ! output once every day (24 hours)
!
!      FILENAME2   = 'Plume_concentration_molec_' // TRIM(ADJUSTL(YEAR_C)) // &
!         '-' //TRIM(ADJUSTL(MONTH_C)) // '-' // TRIM(ADJUSTL(DAY_C)) //      &
!         '-' //TRIM(ADJUSTL(HOUR_C)) // ':' // TRIM(ADJUSTL(MINUTE_C)) //    &
!         ':' //TRIM(ADJUSTL(SECOND_C)) // '.txt'
!
!      OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
!            FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!
!      WRITE(262,*)'total plume number:', n_box
!
!
!      DO i_box = 1,n_box,1
!!      IF(i_box<3)THEN !!! shw
!
!        IF(Judge_plume(i_box)==0) THEN
!          WRITE(262,*) i_box, ': dissolved'
!        ENDIF
!
!        IF(Judge_plume(i_box)==2) THEN
!          WRITE(262,*) i_box, ': 2D model, total mass [molec]:'
!          WRITE(262,*) Pdx(i_box)*Pdy(i_box)*box_length(i_box)*1.0e+6_fp &
!                * SUM(box_concnt_2D(:,:,N_species,i_box))
!        ENDIF
!
!        IF(Judge_plume(i_box)==1) THEN
!          WRITE(262,*) i_box, ': 1D model'
!          WRITE(262,*) SUM(box_concnt_1D(1:n_slab_max,N_species,i_box)) &
!                   *box_Ra(i_box)*box_Rb(i_box)*box_length(i_box)*1.0e+6_fp
!        ENDIF
!
!        WRITE(262,*)' '
!
!!     ENDIF ! IF(i_box<3)THEN
!     ENDDO ! DO i_box = 1,n_box,1
!
!
!!       DO i_box = 1, n_box
!!          WRITE(262,*) i_box, SUM(box_concnt(i_box,:,N_species)*V_slab(i_box,:))
!!       ENDDO
!
!!        WRITE(262,*) box_concnt(1,:,N_species)
!
!    ENDIF ! IF(mod(tt,144)==0)THEN
!!    ENDIF
!
    tt = tt + Dt
!

    CLOSE(261)

    IF(ASSOCIATED(Plume2d_new)) nullify(Plume2d_new)
    IF(ASSOCIATED(PLume2d)) nullify(PLume2d)
    IF(ASSOCIATED(Plume2d_prev)) nullify(Plume2d_prev)


  END SUBROUTINE trajectory_write_std

!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine trajectory_cleanup()


!!! dissolve all the alive plume here ???


!
!    if (allocated(box_lon))      deallocate(box_lon)
!    if (allocated(box_lat))      deallocate(box_lat)
!    if (allocated(box_lev))      deallocate(box_lev)
!    if (allocated(box_length))   deallocate(box_length)
!    if (allocated(box_Ra))       deallocate(box_Ra)
!    if (allocated(box_Rb))       deallocate(box_Rb)
!    if (allocated(box_theta))    deallocate(box_theta)
!
    WRITE(6,'(a)') '--> Lagrange and Plume Module Cleanup <--'
!
!
  end subroutine trajectory_cleanup


END MODULE Trajectory_Mod
