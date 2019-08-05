!------------------------------------------------------------------------------
!         GEOS-Chem Global Chemical Stratospherical Transport Model                  !
!-----------------------------------------------------------------------------!

MODULE Plume_Mod

  USE precision_mod
  USE PhysConstants, ONLY : PI, Re, g0
  USE CMN_SIZE_Mod,  ONLY : IIPAR, JJPAR, LLPAR

  IMPLICIT NONE


  ! !PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: plume_init
  PUBLIC :: plume_run
  PUBLIC :: plume_write_std
  PUBLIC :: plume_cleanup


  integer, parameter :: n_boxes_max = 3  ! 50*40*40
  integer            :: tt

  ! The center of the cylinder
  real(fp), allocatable :: box_lon(:)
  real(fp), allocatable :: box_lat(:)
  real(fp), allocatable :: box_lev(:)

  ! Consider the box as a cylinder

  ! box_radius(n_boxes_max,N_rings)
  real(fp), allocatable :: box_radiusA(:,:)   ! vertical radius
  real(fp), allocatable :: box_radiusB(:,:)   ! horizontal radius of the the cross-section


  ! Theta is the clockwise angle between z-axis (P) and vertical radiusA
  real(fp), allocatable :: box_theta(:)    ! 0 ~ 180 degree
  real(fp), allocatable :: box_length(:)

  ! D_radius should only be used at the beginning!
  real(fp), parameter :: Init_radius = 10.0e+0_fp ! 100.0e+0_fp ! 2.0e+0_fp ! 10.0e+0_fp     ! [m], the width of each ring
  real(fp), parameter :: D_radius    = 10.0e+0_fp ! 100.0e+0_fp ! 2.0e+0_fp ! 10.0e+0_fp     ! [m], the width of each ring
  integer, parameter  :: n_rings_max = 100 ! 500  ! 100          ! Degine the number of rings in one box

  ! medical concentration of each ring
  real(fp), allocatable :: box_concnt(:,:)    ! [kg/m3], box_concnt(n_boxes_max,N_rings)
  real(fp), allocatable :: box_concnt_K(:)    ! [kg/m3], box_concnt_K(N_rings_max)
  real(fp), allocatable :: RK(:,:), AA_env(:)  ! for Runge-Kutta method
  real(fp), allocatable :: eddy_h(:)  ! 
  real(fp), allocatable :: eddy_diffA(:)  ! 
  real(fp), allocatable :: eddy_diffB(:)  ! 
  real(fp), allocatable :: kA(:), kB(:)  ! 

  real(fp), allocatable :: env_amount(:)

CONTAINS


!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  SUBROUTINE plume_init( am_I_root )

    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU
    integer                       :: i_ring

    CHARACTER(LEN=255)            :: FILENAME2


    WRITE(6,'(a)') '--------------------------------------------------------'
    WRITE(6,'(a)') ' Initial Plume Module (Using Dynamic time step)'
    WRITE(6,'(a)') '--------------------------------------------------------'

    allocate(box_lon(n_boxes_max))
    allocate(box_lat(n_boxes_max))
    allocate(box_lev(n_boxes_max))

    allocate(box_radiusA(n_boxes_max,n_rings_max))
    allocate(box_radiusB(n_boxes_max,n_rings_max))

    allocate(box_theta(n_boxes_max))
    allocate(box_length(n_boxes_max))

    allocate(box_concnt(n_boxes_max,n_rings_max))
    allocate(box_concnt_K(n_rings_max))
    allocate(RK(4,n_rings_max))
    allocate(AA_env(4))

    allocate(eddy_h(n_rings_max))
    allocate(eddy_diffA(n_rings_max))
    allocate(eddy_diffB(n_rings_max))
    allocate(kA(n_rings_max))
    allocate(kB(n_rings_max))

    allocate(env_amount(n_boxes_max))


    box_lon    = (/5.0e+0_fp,  5.1e+0_fp,  5.2e+0_fp/)
    box_lat    = (/4.0e+0_fp,  4.1e+0_fp,  4.2e+0_fp/)
    box_lev    = (/20.0e+0_fp, 20.0e+0_fp, 20.0e+0_fp/)      ! hPa


!    box_radiusA(:,1)  = (/2.0e+0_fp,  2.0e+0_fp,  2.0e+0_fp/)     ! the value of the innest ring for every box
!    box_radiusB(:,1)  = (/2.0e+0_fp,  2.0e+0_fp,  2.0e+0_fp/)     ! m

    box_radiusA(:,1)  = (/10.0e+0_fp,  10.0e+0_fp,  10.0e+0_fp/)     ! the value of the innest ring for every box
    box_radiusB(:,1)  = (/10.0e+0_fp,  10.0e+0_fp,  10.0e+0_fp/)     ! m

!    box_radiusA(:,1)  = (/100.0e+0_fp,  100.0e+0_fp,  100.0e+0_fp/)     ! the
!    box_radiusB(:,1)  = (/100.0e+0_fp,  100.0e+0_fp,  100.0e+0_fp/)     ! m

    ! Set the initial value of max/min radius for each ring
    do i_ring=2,n_rings_max
      box_radiusA(:,i_ring) = box_radiusA(:,i_ring-1) + D_radius
      box_radiusB(:,i_ring) = box_radiusB(:,i_ring-1) + D_radius
    enddo

    box_theta         = (/0.0e+0_fp,    0.0e+0_fp,    0.0e+0_fp/)             ! degree
    box_length        = (/1000.0e+0_fp, 1000.0e+0_fp, 1000.0e+0_fp/)    ! m


    ! Set the initial concentration
    do i_ring=1,n_rings_max
      box_concnt(:,i_ring) = (/0.0e+0_fp,  0.0e+0_fp,  0.0e+0_fp/)   ! [kg/m3]
      box_concnt_K(i_ring) = 0.0e+0_fp
    enddo  


      ! the width of each ring is 2 meters 
!     box_concnt(1,1:250) = (/ 49.99500024999167, 49.95502024392637,                      &
!         49.875156119873004, 49.75559927079149, 49.59663583027856, 49.39864553154191,            &
!         49.16210019596277, 48.88756185966682, 48.57568054851479, 48.227191713840774,            &
!         47.84291334309548, 47.42374276128023, 46.97065314067379, 46.48468973784772,             &
!         45.96696587832591, 45.4186587104634, 44.84100475118934, 44.23529524717418,              &
!         43.60287137573596, 42.945119310392386, 42.26346517639093, 41.55936992180854,            &
!         40.83432412990554, 40.089842798347064, 39.327460110672746, 38.54872422500576,           &
!         37.75519210445118, 36.94842441294721, 36.12998049950968, 35.30141349285698,             &
!         34.464265527331314, 33.62006311985014, 32.770312716342026, 31.916496424753753,          &
!         31.06006795027254, 30.20244874690127, 29.345024397966895, 28.48914123654615,            &
!         27.636103215169694, 26.78716903252924, 25.94354952327262, 25.106405315342347,           &
!         24.276844757703973, 23.45592211973321, 22.64463606199473, 21.84392837666111,            &
!         21.05468299439559, 20.27772525316603, 19.513821423176058, 18.763678480900367,           &
!         18.02794412409878, 17.307207018663934, 16.60199726723303, 15.912787088668537,           &
!         15.239991696787666, 14.583970366097326, 13.94502767176999, 13.323414890676204,          &
!         12.719331549970189, 12.132927109503541, 11.564302764216434, 11.013513352622088,         &
!         10.48056935755489, 9.965438985491133, 9.468050310968708, 8.98829347292339,              &
!         8.526022910118579, 8.081059623266963, 7.65319345192056, 7.242185354733357,              &
!         6.84776968227266, 6.469656432164565, 6.1075334769995, 5.761068756089673,                &
!         5.429912422855146, 5.113698940313495, 4.812049117853894, 4.524572083184796,             &
!         4.2508671840497545, 3.9905258150039833, 3.743133165230175, 3.50826988404216,            &
!         3.2855136613751426, 3.0744407211885036, 2.8746272263082275, 2.6856505938086763,         &
!         2.507090720575147, 2.338531119197949, 2.1795599648239308, 2.029771054031513,            &
!         1.8887646771995985, 1.7561484062087604, 1.6315377996448017, 1.514557027970317,          &
!         1.4048394213896542, 1.3020279433574369, 1.2057755928714535, 1.11574573884832,           &
!         1.0316123900062355, 0.9530604037747433, 0.8797856378183132, 0.8114950478003603,         &
!         0.7479067350288746, 0.6887499476158524, 0.6337650387521545, 0.582703385649069,          &
!         0.5353272726296837, 0.49140974176898444, 0.4507344143832893, 0.413095286559015,         &
!         0.378296501789557, 0.3461521036590272, 0.31648577137428735, 0.2891305408037323,         &
!         0.26392851353406127, 0.24073055630619844, 0.21939599303986285,                          &
!         0.19979229150423164, 0.18179474654074906, 0.16528616159440682,                          &
!         0.15015653016263605, 0.13630271862708504, 0.1236281517937096,                           &
!         0.1120425023313622, 0.101461385168946, 0.09180605778662113, 0.08300312721784699,        &
!         0.07498426446649228, 0.06768592693702603, 0.061049089376051933,                         &
!         0.055018983730226116, 0.04954484823891421, 0.044579685999760976,                        &
!         0.04008003317157138, 0.03600573691141489, 0.03231974308149995,                          &
!         0.02898789370593113, 0.025978734107741923, 0.02326332961234793,                         &
!         0.020815091664535858, 0.01860961317202183, 0.01662451285919873,                         &
!         0.014839288389660534, 0.013235177995154574, 0.011795030331478853,                       &
!         0.010503182268224097, 0.009345344308876003, 0.008308493330363867,                       &
!         0.007380772326398169, 0.006551396836620228, 0.00581056774344261,                        &
!         0.0051493901202482135, 0.004559797818113302, 0.004034483483209832,                      &
!         0.0035668337033221663, 0.0031508689892959666, 0.0027811883055445952,                    &
!         0.0024529178728103862, 0.0021616639760634274, 0.0019034695205832695,                    &
!         0.0016747740897840246, 0.001472377269098696, 0.001293405011132706,                      &
!         0.0011352788282396743, 0.000995687609584682, 0.0008725618705716433,                     &
!         0.00076405025316125, 0.0006684981060421901, 0.0005844279837966539,                      &
!         0.0005105219140849567, 0.0004456052914330938, 0.00038863226541753536,                   &
!         0.00033867249988518743, 0.0002948991883098607, 0.0002565782184612303,                   &
!         0.00022305838724322438, 0.00019376257384503556, 0.00016817978624128203,                 &
!         0.00014585800258296413, 0.00012639773514460073, 0.00010944625024414215,                 &
!         9.469238294142673e-05, 8.186189035980977e-05, 7.071329217699776e-05,                    &
!         6.103415120867657e-05, 5.263775107644034e-05, 4.536013172440305e-05,                    &
!         3.9057447041522455e-05, 3.3603612073960644e-05, 2.8888210288571284e-05,                 &
!         2.48146340894679e-05, 2.12984343089681e-05, 1.826585670603756e-05,                      &
!         1.5652545623270556e-05, 1.34023968905797e-05, 1.1466543832735418e-05,                   &
!         9.802461850783841e-06, 8.373178515687453e-06, 7.1465774472215825e-06,                   &
!         6.094785462574568e-06, 5.1936335770144655e-06, 4.42218343259704e-06,                    &
!         3.7623116288224154e-06, 3.198345239473928e-06, 2.716742530120305e-06,                   &
!         2.305813547942261e-06, 1.9554758461631864e-06, 1.6570411354494168e-06,                  &
!         1.4030291297701508e-06, 1.1870052795086636e-06, 1.0034394648235848e-06,                 &
!         8.475830617142456e-07, 7.153620959283844e-07, 6.032844694044696e-07,                    &
!         5.083594836898615e-07, 4.2802809774955484e-07, 3.601025465352785e-07,                   &
!         3.027141141242445e-07, 2.5426800343049753e-07, 2.1340437549323007e-07,                  &
!         1.7896474701143994e-07, 1.4996303679675482e-07, 1.2556064166356478e-07,                 &
!         1.0504500132345863e-07, 8.781118131610907e-08, 7.334606398217998e-08,                   &
!         6.121479094186877e-08, 5.104914735796671e-08, 4.253761921353282e-08,                    &
!         3.5416890615392144e-08, 2.9464579363157648e-08, 2.449303624804273e-08,                  &
!         2.0344057253278993e-08, 1.688437844988782e-08, 1.4001841299704752e-08,                  &
!         1.1602131629273281e-08, 9.60600902130549e-09, 7.946955047258184e-09,                    &
!         6.569178866216583e-09, 5.425927444760903e-09, 4.4780551880986964e-09,                   &
!         3.692814270790277e-09, 3.042832552759166e-09, 2.505250782685961e-09,                    &
!         2.0609949309071472e-09, 1.6941630488153963e-09, 1.391509100983489e-09,                  &
!         1.1420088289968512e-09, 9.364949424197002e-10, 7.673508446499213e-10 /)


        ! the width of each ring is 10 meters    
        box_concnt(1,1:50) = (/49.875156119873004, 48.88756185966682                               &
               ,46.97065314067379, 44.23529524717418, 40.83432412990554, 36.94842441294721     &
               ,32.770312716342026, 28.48914123654615, 24.276844757703973, 20.27772525316603   &
               ,16.60199726723303, 13.323414890676204, 10.48056935755489, 8.081059623266963    &
               ,6.1075334769995, 4.524572083184796, 3.2855136613751426, 2.338531119197949      &
               ,1.6315377996448017, 1.11574573884832, 0.7479067350288746, 0.49140974176898444  &
               ,0.31648577137428735, 0.19979229150423164, 0.1236281517937096                   &
               ,0.07498426446649228, 0.044579685999760976,0.025978734107741923                 &
               ,0.014839288389660534, 0.008308493330363867, 0.004559797818113302               &
               ,0.0024529178728103862, 0.001293405011132706, 0.0006684981060421901             &
               ,0.00033867249988518743, 0.00016817978624128203, 8.186189035980977e-05          &
               ,3.9057447041522455e-05, 1.826585670603756e-05, 8.373178515687453e-06           &
               ,3.7623116288224154e-06, 1.6570411354494168e-06, 7.153620959283844e-07          &
               ,3.027141141242445e-07, 1.2556064166356478e-07, 5.104914735796671e-08           &
               ,2.0344057253278993e-08, 7.946955047258184e-09, 3.042832552759166e-09           &
               ,1.1420088289968512e-09/)


        ! the width of each ring is 100 meters
!        box_concnt(1,1:9) = (/38.94003915357025, 5.269961228093217,             &
!                0.09652270681138546, 0.00023925586960645045, 8.026140275928058e-08,     &
!                3.643862047909846e-12, 2.2388662208591508e-17, 1.8616815608752554e-23,  &
!                2.0950465972471985e-30/)




    env_amount = (/0.0e+0_fp, 0.0e+0_fp, 0.0e+0_fp/)


    ! Create output file
    FILENAME2   = 'Plume_theta_max_min_radius.txt'
    tt = 0

    OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='REPLACE', &
          FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

    !Do i_ring = 1, n_rings_max
    !   WRITE(262,'(I0.4,3(x,E16.5E4))') i_ring, box_theta(1), box_radiusA(1,i_ring), box_radiusB(1,i_ring)
    !End Do

    write(262,*)box_concnt(1,1:n_rings_max)

  END SUBROUTINE plume_init


!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

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

    REAL :: Dt, Dt2          ! = 600.0e+0_fp          

    integer :: i_box, i_ring

    integer :: i_lon, i_lat, i_lev
    
    integer :: t1s, Ki

    real(fp) :: curr_lon, curr_lat, curr_pressure
    real(fp) :: next_lon, next_lat, next_pressure

    real(fp) :: X_edge2, Y_edge2

    real(fp), pointer :: u(:,:,:)
    real(fp), pointer :: v(:,:,:)
    real(fp), pointer :: omeg(:,:,:)

    real(fp), pointer :: Ptemp(:,:,:)

    real(fp), pointer :: P_BXHEIGHT(:,:,:)

    real(fp) :: curr_u, curr_v, curr_omeg

    real(fp) :: curr_Ptemp, Ptemp_shear  ! potential temperature
    real(fp) :: U_shear, V_shear, UV_shear  

    real(fp), pointer :: P_edge(:)
    real(fp), pointer :: P_mid(:)

    real(fp) :: Dx, Dy
    real(fp) :: wind_s_shear
    real(fp) :: theta_previous       
    real(fp) :: dbox_theta, dbox_radiusA, dbox_radiusB
    ! angle between travel direction and lon (-90 ~ 90 deg)
    real(fp) :: box_alpha        

    real(fp), pointer :: X_edge(:)    
    real(fp), pointer :: Y_edge(:)       
    real(fp), pointer :: X_mid(:)    
    real(fp), pointer :: Y_mid(:)       

    real(fp), pointer :: P_I, R_e     

    real(fp) :: AA(n_rings_max), BB(n_rings_max), DD(n_rings_max)  ! used to calculate concentration distribution
    real(fp) :: eddy_v
    real(fp) :: Cv, Ch, Omega_N

    !real(fp) :: RK(4,n_rings_max)
    !real(fp) :: box_concnt_K(n_rings_max)

    CHARACTER(LEN=255)            :: FILENAME2

    FILENAME2   = 'Plume_theta_max_min_radius.txt'


    Dt  = GET_TS_DYN()

    u => State_Met%U   ! figure out state_met%U is based on lat/lon or modelgrid(i,j)
    v => State_Met%V   ! V [m s-1]
    omeg => State_Met%OMEGA  ! Updraft velocity [Pa/s]

    Ptemp => State_Met%THETA  ! Updraft velocity [Pa/s]

    P_edge => State_Met%PEDGE(1,1,:)  ! Wet air press @ level edges [hPa]
    P_mid => State_Met%PMID(1,1,:)  ! Pressure (w/r/t moist air) at level centers (hPa)

    P_BXHEIGHT => State_Met%BXHEIGHT  ![IIPAR,JJPAR,KKPAR]

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

       ! Judge the sign of box_alpha
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

       curr_Ptemp = Interplt_wind(Ptemp,X_mid, Y_mid, P_mid, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)


       i_lon = Find_iLonLat(next_lon, Dx, X_edge2)
       i_lat = Find_iLonLat(next_lat, Dy, Y_edge2)
       i_lev = Find_iPLev(next_pressure,P_edge)

       !!!!!!
       ! For deformation of cross-section caused by wind shear (A.D.Naiman et al., 2010):
       !!!!!!

       ! calculate the wind_s shear along pressure direction
       wind_s_shear = Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
       ! attention ***
       wind_s_shear = 0.0


       theta_previous   = box_theta(i_box)
       box_theta(i_box) = ATAN( TAN(box_theta(i_box)) + wind_s_shear*Dt )


       do i_ring=1,n_rings_max
       ! make sure use box_theta or TAN(box_theta)  ??? 
       box_radiusA(i_box,i_ring) = box_radiusA(i_box,i_ring) * (TAN(box_theta(i_box))**2+1)**0.5    / (TAN(theta_previous)**2+1)**0.5
       box_radiusB(i_box,i_ring) = box_radiusB(i_box,i_ring) * (TAN(box_theta(i_box))**2+1)**(-0.5) / (TAN(theta_previous)**2+1)**(-0.5)
       enddo


       !!!!!!
       ! For the concentration change caused by eddy diffusion:
       ! box_concnt(n_boxes_max,N_rings)
       !!!!!!

         ! Calculate vertical eddy diffusivity (U.Schumann, 2012) :
         Cv = 0.2
         Omega_N = 0.1
         Ptemp_shear = Vertical_shear(Ptemp, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev,curr_lon, curr_lat, curr_pressure)
         eddy_v = Cv * Omega_N**2 / sqrt( (Ptemp_shear*g0/curr_Ptemp) )

         ! attenttion ***
         eddy_v = 1.0


         ! Calculate horizontal eddy diffusivity:
         Ch = 0.1
         U_shear = Vertical_shear(u, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         V_shear = Vertical_shear(v, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge, i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
         UV_shear = sqrt( U_shear**2 + V_shear**2 )
         do i_ring=1,n_rings_max
          ! eddy_h(i_ring) = Ch*UV_shear*(Init_radius+(i_ring-1)*D_radius)**2
          ! attention ***
          eddy_h(i_ring) = 1.0

          eddy_diffA(i_ring) = eddy_v*cos(box_theta(i_box)) + eddy_h(i_ring)*sin(abs(box_theta(i_box))) ! a
          eddy_diffB(i_ring) = eddy_v*sin(abs(box_theta(i_box))) + eddy_h(i_ring)*cos(box_theta(i_box)) ! b
         enddo


         !---------------------------------------------------------------------------------------------------
         ! Calculate the transporting speed between different rings based on
         ! eddy diffusivity

         ! For the innest ring (i_ring = 1)
         ! kB should be rewrite in a more accurate equation !!!
!         kB(1) = eddy_diffB(1) / ( (box_radiusB(i_box,2)-0.0) / 2.0 )
!         kA(1) = eddy_diffA(1) / ( (box_radiusA(i_box,2)-0.0) / 2.0 )
         kB(1) = eddy_diffB(1) / ( box_radiusB(i_box,1) + (box_radiusB(i_box,2)-box_radiusB(i_box,1)) / 2.0 )
         kA(1) = eddy_diffA(1) / ( box_radiusA(i_box,1) + (box_radiusA(i_box,2)-box_radiusA(i_box,1)) / 2.0 )


         do i_ring = 2, n_rings_max-1

         kB(i_ring) = eddy_diffB(i_ring) / ((box_radiusB(i_box,i_ring+1)-box_radiusB(i_box,i_ring-1)) / 2.0 )
         kA(i_ring) = eddy_diffA(i_ring) / ((box_radiusA(i_box,i_ring+1)-box_radiusA(i_box,i_ring-1)) / 2.0 )

         enddo ! i_ring

         ! For outest ring (i_ring = n_rings_max)
         kB(n_rings_max) = eddy_diffB(n_rings_max)/(box_radiusB(i_box,n_rings_max)-box_radiusB(i_box,n_rings_max-1) )
         kA(n_rings_max) = eddy_diffA(n_rings_max)/(box_radiusA(i_box,n_rings_max)-box_radiusA(i_box,n_rings_max-1) )


         ! Another method to calculate kA, kB
!         do i_ring = 1, n_rings_max        

!         kB(i_ring) = eddy_diffB(i_ring) / box_radiusB(i_box,i_ring)
!         kA(i_ring) = eddy_diffA(i_ring) / box_radiusA(i_box,i_ring)

!         enddo ! i_ring


       ! --------------------------------------------------------------------------------------------------
       ! Begin to calculate the transportation between different rings:
       Dt2 = 1.0    ! 0.1     ! 1.0
       do t1s = 1, int(Dt/Dt2)

       ! Use classical Runge-Kutta method (RK4) to solve the diferential
       ! equation
       do Ki = 1,4

         do i_ring = 1, n_rings_max
           if(Ki==1)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring)
           else if(Ki==4)then
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring) + RK(3,i_ring)*Dt2 !Dt
           else
             box_concnt_K(i_ring) = box_concnt(i_box,i_ring) + RK(Ki-1,i_ring)*Dt2*0.5 !Dt
           endif
         enddo ! i_ring

 
         AA(1) = ( kB(1)*box_radiusA(i_box,1) + kA(1)*box_radiusB(i_box,1) ) * ( box_concnt_K(2) - box_concnt_K(1) )

         BB(1) = 0.0

         DD(1) = box_radiusA(i_box,1) * box_radiusB(i_box,1)

         RK(Ki,1)          = AA(1)/DD(1)


         ! For rings from 2 to (n_rings_max - 1)
         do i_ring = 2, n_rings_max-1

           AA(i_ring) = ( kB(i_ring)*box_radiusA(i_box,i_ring) + kA(i_ring)*box_radiusB(i_box,i_ring) ) * ( box_concnt_K(i_ring+1) - box_concnt_K(i_ring) )

           BB(i_ring) = ( kB(i_ring-1)*box_radiusA(i_box,i_ring-1) + kA(i_ring-1)*box_radiusB(i_box,i_ring-1) ) * ( box_concnt_K(i_ring-1) - box_concnt_K(i_ring) ) 

           DD(i_ring) = box_radiusA(i_box,i_ring)*box_radiusB(i_box,i_ring) - box_radiusA(i_box,i_ring-1)*box_radiusB(i_box,i_ring-1)

           RK(Ki,i_ring) = (AA(i_ring)+BB(i_ring))/DD(i_ring)

         enddo ! i_ring
 

         AA(n_rings_max) = ( kB(n_rings_max)*box_radiusA(i_box,n_rings_max) + kA(n_rings_max)*box_radiusB(i_box,n_rings_max) ) * ( 0 - box_concnt_K(n_rings_max) ) 

         BB(n_rings_max) = ( kB(n_rings_max-1)*box_radiusA(i_box,n_rings_max-1) + kA(n_rings_max-1)*box_radiusB(i_box,n_rings_max-1) ) * ( box_concnt_K(n_rings_max-1) - box_concnt_K(n_rings_max) ) 

         DD(n_rings_max) = box_radiusA(i_box,n_rings_max) * box_radiusB(i_box,n_rings_max) - box_radiusA(i_box,n_rings_max-1) * box_radiusB(i_box,n_rings_max-1)

         RK(Ki,n_rings_max) = ( AA(n_rings_max) + BB(n_rings_max) ) / DD(n_rings_max)

         AA_env(Ki) = -1.0 * AA(n_rings_max)


       enddo ! Ki


       do i_ring = 1,n_rings_max
         box_concnt(i_box,i_ring) = box_concnt(i_box,i_ring) + Dt2*( RK(1,i_ring)+2.0*RK(2,i_ring)+2.0*RK(3,i_ring)+RK(4,i_ring) )/6.0 ! Dt
       enddo !i_ring
         env_amount(i_box) = env_amount(i_box) + Dt2*( AA_env(1)+2.0*AA_env(2)+2.0*AA_env(3)+AA_env(4) )/6.0


       !if(i_box==1)then
       !  write(6,*) 'RK 1-4', RK(:,1)
       !  write(6,*)'= concentration 1-5 =>', box_concnt(i_box,1:5)
       !  write(6,*)'= concentration 6-10 =>', box_concnt(i_box,6:10)
       !  write(6,*)'= concentration 11-15 =>', box_concnt(i_box,11:15)
       !  write(6,*)'= total amount =>', t1s, sum( box_concnt(i_box,:)*DD(:) ) + env_amount(i_box)
       !endif

        if(i_box==1)then
!        if(MOD(t1s,10)==0)then

         OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD', &
               FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!         WRITE(262,*)'==>', t1s, box_theta(1), box_radiusA(1,1), box_radiusB(1,1), box_radiusA(1,2), box_radiusB(1,2)
!         WRITE(262,*)'kA', kA(:)
!         WRITE(262,*)'kB', kB(:)
!         WRITE(262,*)'AA', AA(:)
!         WRITE(262,*)'BB', BB(:)
!         WRITE(262,*)'DD', DD(:)
         write(262,*) box_concnt(1,1:n_rings_max)

        endif ! if(i_box==1)then
!        endif ! if(MOD(t1s,10)==0)then

       enddo ! t1s

    enddo ! i_box

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



!-------------------------------------------------------------------
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
! functions to interpolate wind speed (u,v,omeg) into a specific location (curr_lon, curr_lat, curr_pressure)
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



  !-------------------------------------------------------------------
  ! calculation the great-circle distance between two points on the earth surface

  real function Distance_Circle(x1, y1, x2, y2)
    implicit none
    real(fp)     :: x1, y1, x2, y2  ! unit is degree
    !real(fp) :: PI, Re

    Distance_Circle = Re * 2.0 * ASIN( (sin( (y1-y2)*PI/180.0 ))**2.0   & 
                      + cos(x1*PI/180.0) * cos(x2*PI/180.0) * (sin( 0.5*(x1-x2)*PI/180.0 ))**2.0 )
    return
  end function


  !-------------------------------------------------------------------  
  ! calculate the wind_s (inside a plume sross-section) shear along pressure direction

  real function Wind_shear_s(u, v, P_BXHEIGHT, box_alpha, X_mid, Y_mid, P_mid, P_edge,  i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
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
    Delt_height =  Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev), P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1), P_edge(init_lev), P_edge(init_lev+1), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    Wind_shear_s = ( wind_s(2) - wind_s(1) ) / Delt_height

    return
  end function


!-------------------------------------------------------------------
  ! transform the pressure level [Pa] to the height level [m]

  real function Pa2meter(Box_height, P1, P2, Judge) 
  ! Judge: 0 for bottom, 1 for top
  ! 
    implicit none
    real(fp) :: Box_height, P1, P2
    integer  :: Judge
      
      if(Judge==1)then
        ! calculate the height of the top half of the grid box
        Pa2meter = Box_height * ( DLOG(P2) - DLOG(0.5*(P2+P1)) ) / ( DLOG(P2) - DLOG(P1) )
      else
        ! calculate the height of the bottom half of the grid box
        Pa2meter = Box_height * ( DLOG(0.5*(P2+P1)) - DLOG(P1) ) / ( DLOG(P2) - DLOG(P1) )
      endif

    return
  end function


!-------------------------------------------------------------------
! Calculate eddy diffusivity in stratisphere. (U.Schumann, 2012)



  !-------------------------------------------------------------------  
  ! calculate the vertical shear for calculating eddy difussivity

  real function Vertical_shear(var, P_BXHEIGHT, X_mid, Y_mid, P_mid, P_edge,  i_lon, i_lat, i_lev, curr_lon, curr_lat, curr_pressure)
    implicit none
    real(fp)          :: curr_lon, curr_lat, curr_pressure
    !real(fp), pointer :: PI, Re
    real(fp), pointer :: var(:,:,:)
    real(fp), pointer :: P_BXHEIGHT(:,:,:)
    real(fp), pointer :: X_mid(:), Y_mid(:), P_mid(:), P_edge(:)
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
        distance(i,j) = Distance_Circle(curr_lon, curr_lat, X_mid(ii),Y_mid(jj))

      enddo
    enddo

    do i = 1,2
      do j = 1,2
        Weight(i,j) = 1.0/distance(i,j) / sum( 1.0/distance(:,:) )
      enddo
    enddo


    do k = 1,2
      kk            = k + init_lev - 1
      var_lonlat(k) =  Weight(1,1) * var(i_lon,i_lat,kk)   + Weight(1,2) * var(i_lon,i_lat+1,kk)   &
                       + Weight(2,1) * var(i_lon+1,i_lat,kk) + Weight(2,2) * var(i_lon+1,i_lat+1,kk)
    enddo


    ! second vertical shear of wind

    Delt_height =  Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev),P_edge(init_lev), P_edge(init_lev+1), 1 )   &
                 + Pa2meter( P_BXHEIGHT(init_lon,init_lat,init_lev+1),P_edge(init_lev+1), P_edge(init_lev+2), 0 )
    ! find the z height of each pressure level in GEOS-Chem

    Vertical_shear = ( var_lonlat(2) - var_lonlat(1) ) / Delt_height

    return
  end function


!-------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------


  SUBROUTINE plume_write_std( am_I_Root, RC )


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

    USE Input_Opt_Mod, ONLY : OptInput

    ! Parameters for netCDF routines
    include "netcdf.inc"


    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success

    REAL(f8), POINTER         :: Arr1D(:,:)
    REAL(f8), POINTER         :: Arr2D(:,:)
    REAL(f8), POINTER         :: Arr2DOld(:,:)
    REAL(f4), POINTER         :: Arr3D(:,:)
    REAL(f4), POINTER         :: Arr4D(:,:)
    REAL*8,   POINTER         :: timeVec(:)
    CHARACTER(LEN=255)        :: NcFile
    CHARACTER(LEN=255)        :: Pfx, title, Reference, Contact
    CHARACTER(LEN=255)        :: timeunit
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nTime
    INTEGER                   :: i_box, i_ring
!    INTEGER                   :: L
    LOGICAL                   :: IsOldFile

    CHARACTER(LEN=255)            :: FILENAME2

    FILENAME2   = 'Plume_theta_max_min_radius.txt'
    tt = tt +1

!    IF(mod(tt,6)==0)THEN     ! output once every hour
!    IF(mod(tt,144)==0)THEN   ! output once every day (24 hours)

!       OPEN( 262,      FILE=TRIM( FILENAME2   ), STATUS='OLD', &
!             FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

!       Do i_ring = 1, n_rings_max
!        WRITE(262,'(I0.4,3(x,E16.5E4))') i_ring, box_theta(2), box_radiusA(2,i_ring), box_radiusB(2,i_ring), box_concnt(2,i_ring)
!       End Do


!    ENDIF

  END SUBROUTINE plume_write_std

!---------------------------------------------------------------------
!*********************************************************************
!---------------------------------------------------------------------

  subroutine plume_cleanup()

    if (allocated(box_lon))      deallocate(box_lon)
    if (allocated(box_lat))      deallocate(box_lat)
    if (allocated(box_lev))      deallocate(box_lev)
    if (allocated(box_length))   deallocate(box_length)
    if (allocated(box_radiusA))  deallocate(box_radiusA)
    if (allocated(box_radiusB))  deallocate(box_radiusB)
    if (allocated(box_theta))    deallocate(box_theta)

    WRITE(6,'(a)') '--> Plume Module Cleanup <--'

  end subroutine plume_cleanup


END MODULE Plume_Mod
