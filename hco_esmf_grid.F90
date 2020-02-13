#define VERIFY_(A) if(.not.HCO_ESMF_VRFY(A,subname,__LINE__))call exit(-1)
#define ASSERT_(A) if(.not.HCO_ESMF_ASRT(A,subname,__LINE__))call exit(-1)
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_esmf_grid
!
! !DESCRIPTION: Module HCO\_ESMF\_GRID defines functions for the regridding and
!  control of the intermediate HEMCO grid in the HEMCO to CAM interface.
!\\
!\\
! !INTERFACE:
!
module hco_esmf_grid
!
! !USES:
!
    ! ESMF function wrappers
    use hco_esmf_wrappers

    ! ESMF types
    use ESMF,                     only: ESMF_Mesh, ESMF_DistGrid

    ! MPI
    use mpi,                      only: MPI_PROC_NULL, MPI_SUCCESS

    ! Floating point type.
    ! FIXME: May change to HEMCO precision later down the line
    use shr_kind_mod,             only: r8 => shr_kind_r8

    implicit none
    private
    save
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public        :: HCO_Grid_Init

! !REMARKS:
!
!  Horizontal grid specifications are based on GEOS-Chem definitions.
!  Vertical grid specifications are copied from hyai, hybi from CAM.
!
!  Notes:
!  (i)  In GEOS-Chem, level 1 is bottom-of-atmos, so 'bottom2top' lev_sequence.
!  (ii) 
!
!  ---------------------- BELOW ARE SIDE REMARKS ---------------------------------
!  This module was written by hplin on a gloomy day (as usual) in Cambridge/Boston
!  where I patiently copied the following modules and translated the terminology
!  between them:
!    GEOS-Chem:     state_grid_mod, hco_types_mod, state_met_mod, pressure_mod
!    CAM Ionosphere Interface:  edyn_mpi, edyn_geogrid, ionosphere_interface, hycoef, ref_pres
!  The result is a mash-up of CAM and GEOS-Chem coding terminology, but in the end
!  it is mostly GEOS-Chem. Conversions are done at initialization.
!
!  The vertical grid specification conversions were hand-written
!  with napkin-based calculations, so buyer beware. They need to be carefully
!  validated before final shipping.
!
!  Horizontal grid specifications were mostly borrowed from GEOS-Chem definitions.
!  We use a global ??? grid for testing at this point and they need to be
!  un-hardcoded in the future.
!
!  Thanks to Thibaud Fritz for supplying CESM-GC dev code which helped in the
!  conversions from CAM to GC speak.
!
!  -- hplin, 2/11/20
!
! !REVISION HISTORY:
!  11 Feb 2020 - H.P. Lin    - First crack.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    ! Global grid parameters.
    ! this may be refactored into some other structure that isn't global later.
    ! for now this will do (hplin, 2/11/20) -- and I am confident this will
    ! be the way for at least a few more years, because who touches working code? ;)
    integer, public, protected :: IM                 ! # of lons
    integer, public, protected :: JM                 ! # of lats
    integer, public, protected :: LM                 ! # of levs

    ! Computed parameters for compatibility with GEOS-Chem
    real(r8), public, protected:: DX                 ! Delta X           [deg long]
    real(r8), public, protected:: DY                 ! Delta X           [deg lat]

    ! Horizontal Coordinates
    real(r8), public, protected, allocatable ::    &
                                  XMid (:,:),      & ! Longitude centers [deg]
                                  XEdge(:,:),      & ! Longitude edges   [deg]
                                  YMid (:,:),      & ! Latitude  centers [deg]
                                  YEdge(:,:),      & ! Latitude  edges   [deg]
                                  YEdge_R(:,:),    & ! Latitude  edges R [rad]
                                  YSin (:,:)         ! SIN( lat edges )  [1]

    ! Shadow variables of geo-"meteorological fields" required by HEMCO
    !
    !  Hybrid Grid Coordinate Definition: (dsa, bmy, 8/27/02, 2/2/12)
    !  ============================================================================
    !
    !  The pressure at the bottom edge of grid box (I,J,L) is defined as follows:
    !     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    !  where
    !     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
    !     Ap(L)         has the same units as surface pressure [hPa]
    !     Bp(L)         is  a unitless constant given at level edges
    !
    ! Note: I don't know what to do about PEDGE, surface pressure and the like.
    !       they likely need a regrid through ESMF from State%PSDry or something.
    !       This is only available in the GridComp Run, not worrying about it here.
    !       (hplin, 2/11/20)
    real(r8), public, protected, allocatable ::    &
                                  AREA_M2(:,:),    & ! Area of grid box [m^2]
                                  Ap     (:),      & ! "hyai" Hybrid-sigma Ap value [Pa]
                                  Bp     (:)         ! "hybi" Hybrid-sigma Bp value [Pa]

    ! MPI Descriptors.
    ! Ported mostly from edyn_geogrid and edyn_mpi
    ! -- What everyone knows --
    integer, public, protected :: HCO_mpicom
    integer, public, protected :: nPET               ! Number of PETs
    integer, public, protected :: nPET_lon, nPET_lat ! # of PETs over lon, lat

    integer, public, allocatable :: HCO_petTable(:,:)! 2D table of tasks (dim'l nPET_lon+2, ..lat+2)
                                                     ! extra left and right used for halos

    ! -- Private to MPI process --
    ! Note L dimension (levs) not distributed
    integer, public, protected :: my_IM, my_JM       ! # of lons, levs in this task
    integer, public, protected :: my_IS, my_IE       ! First and last lons
    integer, public, protected :: my_JS, my_JE       ! First and last lats

    integer, public, protected :: my_ID              ! my task ID in HCO_Task
    integer, public, protected :: my_ID_I, my_ID_J   ! mytidi, mytidj coord for current task

    type HCO_Task
        integer  :: ID          ! identifier

        integer  :: ID_I        ! task coord in longitude dim'l of task table
        integer  :: ID_J        ! task coord in latitude  dim'l of task table
        integer  :: IM          ! # of lons on this task
        integer  :: JM          ! # of lats on this task
        integer  :: IS, IE      ! start and end longitude dim'l index
        integer  :: JS, JE      ! start and end latitude  dim'l index
    end type HCO_Task
    type(HCO_Task), allocatable:: HCO_Tasks(:)       ! HCO_Tasks(nPET) avail to all tasks
    

contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Grid_Init
!
! !DESCRIPTION: Subroutine HCO\_Grid\_Init initializes the HEMCO-CAM interface
!  grid descriptions and MPI distribution.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Grid_Init( IM_in, JM_in, nPET_in )
!
! !USES:
!
        ! MPI Properties from CAM
        ! Even though CAM's principle is that only spmd_utils uses MPI,
        ! ionos code uses MPI very liberally. Unfortunately we have to follow
        ! this example as spmd_utils does not provide many of the relevant
        ! functionality like the communicator split. But keep this in mind.
        use cam_logfile,        only: iulog
        use spmd_utils,         only: CAM_mpicom => mpicom, CAM_npes => npes
        use spmd_utils,         only: MPI_SUCCESS
        use spmd_utils,         only: iam
        use spmd_utils,         only: masterproc

        ! Physical constants
        use shr_const_mod,      only: pi => shr_const_pi
        use shr_const_mod,      only: Re => shr_const_rearth

        ! Grid specifications and information from CAM
        use hycoef,             only: ps0, hyai, hybi          ! Vertical specs
        use ppgrid,             only: pver                     ! # of levs

!
! !INPUT PARAMETERS:
!
        integer, intent(in)         :: IM_in, JM_in            ! # lon, lat, lev global
        integer, intent(in)         :: nPET_in                 ! # of PETs to distribute to?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  11 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter :: subname = 'HCO_Grid_Init'
        integer                     :: I, J, L, N, RC
        real(r8)                    :: SIN_N, SIN_S, PI_180

        ! MPI stuff
        integer                     :: color
        integer                     :: lons_per_task, lons_overflow
        integer                     :: lats_per_task, lats_overflow
        integer                     :: lon_beg, lon_end, lat_beg, lat_end, task_cnt

        ! Some physical constants...
        PI_180 = pi / 180.0_r8

        ! Accept external dimensions.
        IM    = IM_in
        JM    = JM_in
        nPET  = nPET_in

        !-----------------------------------------------------------------------
        ! Compute vertical grid parameters
        !-----------------------------------------------------------------------
        ! Can be directly retrieved from hyai, hybi
        ! Although they need to be flipped to be passed from CAM (from tfritz)
        !
        ! Note: In GEOS-Chem, Ap, Bp are defined in hPa, 1
        !       but in HEMCO, they are defined as    Pa, 1 (see hcoi_gc_main_mod.F90 :2813)
        ! So you have to be especially wary of the units.

        ! For now, use the CAM vertical grid verbatim
        LM = pver

        ! Ap, Bp has LM+1 edges for LM levels
        allocate(Ap(LM + 1), STAT=RC)          ! LM levels, LM+1 edges
        allocate(Bp(LM + 1), STAT=RC)
        ASSERT_(RC==0)

        ! Allocate PEDGE information

        ! G-C def: Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
        ! CAM def: Pifce(    L) = hyai(k)*ps0 + [ hybi(k) * ps ]
        !   w.r.t. ps0 = base state srfc prs; ps = ref srfc prs.
        do L = 1, (LM+1)
            Ap(L) = hyai(LM+2-L) * ps0
            Bp(L) = hybi(LM+2-L)
        enddo

        !-----------------------------------------------------------------------
        ! Compute horizontal grid parameters
        !-----------------------------------------------------------------------
        ! Notes: long range (i) goes from -180.0_r8 to +180.0_r8
        !        lat  range (j) goes from - 90.0_r8 to + 90.0_r8

        allocate(XMid (IM,   JM  ), STAT=RC)
        allocate(XEdge(IM+1, JM  ), STAT=RC)
        allocate(YMid (IM,   JM  ), STAT=RC)
        allocate(YEdge(IM,   JM+1), STAT=RC)
        allocate(YEdge_R(IM, JM+1), STAT=RC)
        allocate(YSin (IM,   JM+1), STAT=RC)
        allocate(AREA_M2(IM, JM  ), STAT=RC)
        ASSERT_(RC==0)

        ! Compute DX, DY (lon, lat)
        DX = 360.0_r8 / real(IM, r8)
        DY = 180.0_r8 / real((JM - 1), r8)

        ! Loop over horizontal grid
        ! Note: Might require special handling at poles. FIXME. (hplin, 2/11/20)
        do J = 1, JM
        do I = 1, IM
            ! Longitude centers [deg]
            XMid(I, J) = (DX * (I-1)) - 180.0_r8

            ! Latitude centers [deg]
            YMid(I, J) = (DY * (J-1)) -  90.0_r8

            ! Note half-sized polar boxes for global grid, multiply DY by 1/4 at poles
            if(J == 1) then
                YMid(I, 1) = -90.0_r8 + (0.25_r8 * DY)
            endif
            if(J == JM) then
                YMid(I, JM) = -90.0_r8 - (0.25_r8 * DY)
            endif

            ! Edges [deg] (or called corners in CAM ionos speak)
            XEdge(I, J) = XMid(I, J) - DX * 0.5_r8
            YEdge(I, J) = YMid(I, J) - DY * 0.5_r8
            YEdge_R(I, J) = (PI_180 * YEdge(I, J))
            YSin (I, J) = SIN( YEdge_R(I, J) ) ! Needed for MAP_A2A regridding

            ! Compute the LAST edges
            if(I == IM) then 
                XEdge(I+1,J) = XEdge(I, J) + DX
            endif

            ! Enforce half-sized polar boxes where northern edge of grid boxes
            ! along the SOUTH POLE to be -90 deg lat.
            if(J == 1) then
                YEdge(I, 1) = -90.0_r8
            endif

            if(J == JM) then 
                ! Northern edge of grid boxes along the north pole to be +90 deg lat
                YEdge(I,J+1) = 90.0_r8

                ! Adjust for second-to-last lat edge
                YEdge(I,J  ) = YEdge(I,J+1) - (DY * 0.5_r8)
                YEdge_R(I,J) = YEdge(I,J) * PI_180
                YSin(I, J)   = SIN( YEdge_R(I, J) )

                ! Last latitude edge [radians]
                YEdge_R(I,J+1) = YEdge(I,J+1) * PI_180
                YSin(I,J+1)    = SIN( YEdge_R(I,J+1) )
            endif
        enddo
        enddo

        ! Compute grid box areas after everything is populated...
        do J = 1, JM
        do I = 1, IM
            ! Sine of latitudes at N and S edges of grid box (I,J)
           SIN_N = SIN( YEdge_R(I,J+1) )
           SIN_S = SIN( YEdge_R(I,J  ) )

           ! Grid box surface areas [m2]
           AREA_M2(I,J) = ( DX * PI_180 ) * ( Re**2 ) * (SIN_N - SIN_S)
        enddo
        enddo

        ! Output debug information on the global grid information
        ! Copied from gc_grid_mod.F90 and pressure_mod.F
        if(masterproc) then
            write( iulog, '(a)' )
            write( iulog, '(''%%%%%%%%%%%%%%% HEMCO GRID %%%%%%%%%%%%%%%'')' )
            write( iulog, '(a)' )
            write( iulog, '(''Grid box longitude centers [degrees]: '')' )
            write( iulog, '(8(f8.3,1x))' ) ( XMid(I,1), I=1,IM )
            write( iulog, '(a)' )
            write( iulog, '(''Grid box longitude edges [degrees]: '')' )
            write( iulog, '(8(f8.3,1x))' ) ( XEdge(I,1), I=1,IM+1 )
            write( iulog, '(a)' )
            write( iulog, '(''Grid box latitude centers [degrees]: '')' )
            write( iulog, '(8(f8.3,1x))' ) ( YMid(1,J), J=1,JM )
            write( iulog, '(a)' )
            write( iulog, '(''Grid box latitude edges [degrees]: '')' )
            write( iulog, '(8(f8.3,1x))' ) ( YEdge(1,J), J=1,JM+1 )
            write( iulog, '(a)' )
            write( iulog, '(''SIN( grid box latitude edges )'')' )
            write( iulog, '(8(f8.3,1x))' ) ( YSin(1,J), J=1,JM+1 )

            write( iulog, '(a)'   ) REPEAT( '=', 79 )
            write( iulog, '(a,/)' ) 'V E R T I C A L   G R I D   S E T U P'
            write( iulog, '( ''Ap '', /, 6(f11.6,1x) )' ) AP(1:LM+1)
            write( iulog, '(a)'   )
            write( iulog, '( ''Bp '', /, 6(f11.6,1x) )' ) BP(1:LM+1)
            write( iulog, '(a)'   ) REPEAT( '=', 79 )
        endif

        !-----------------------------------------------------------------------
        ! Distribute among parallelization in MPI 1: Compute distribution
        !-----------------------------------------------------------------------

        ! edyn_geogrid uses 1-D latitude decomposition, so nPET_lon = 1
        ! and nPET_lat = JM. From tfritz this may not work well with GEOS-Chem
        ! so let's attempt some other decomposition.
        do nPET_lon = 1, IM
            nPET_lat = nPET / nPET_lon
            if( nPET_lon * nPET_lat == nPET ) then ! ok right now it will force 1xJM
                exit
            endif
        enddo

        ! Verify we have a correct decomposition
        ! Can't accept invalid decomp; also cannot accept IM, 1 (for sake of consistency)
        if( nPET_lon * nPET_lat /= nPET .or. nPET_lat == 1 ) then
            ! Fall back to same 1-D latitude decomposition like edyn_geogrid
            nPET_lon = 1
            nPET_lat = JM

            if(masterproc) then
                write(iulog,*) "HEMCO: HCO_Grid_Init failed to find a reasonable decomp."
            endif
        endif

        ! Commit to the decomposition at this point
        if(masterproc) then
            write(iulog,*) "HEMCO: HCO_Grid_Init IM, JM, LM", IM, JM, LM
            write(iulog,*) "HEMCO: nPET_lon * nPET_lat = ", nPET_lon, nPET_lat, nPET
        endif

        ! Figure out beginning and ending coordinates for each task
        ! copied from edyn_geogrid
        lons_per_task = IM / nPET_lon
        lons_overflow = MOD(IM, nPET_lon)
        lats_per_task = JM / nPET_lat
        lats_overflow = MOD(JM, nPET_lat)
        lon_beg       = 1
        lon_end       = 0
        lat_beg       = 1
        lat_end       = 0
        task_cnt      = 0
        jloop: do J = 0, nPET_lat - 1
            lat_beg = lat_end + 1
            lat_end = lat_beg + lats_per_task - 1
            if (J < lats_overflow) then
                lat_end = lat_end + 1
            endif
            lon_end = 0
            do I = 0, nPET_lon - 1
                lon_beg = lon_end + 1
                lon_end = lon_beg + lons_per_task - 1
                if (I < lons_overflow) then
                   lon_end = lon_end + 1
                endif
                task_cnt = task_cnt+1
                if (task_cnt > iam) exit jloop ! This makes this loop CPU specific
            enddo
        enddo jloop

        !-----------------------------------------------------------------------
        ! Distribute among parallelization in MPI 2: Distribute among MPI
        !-----------------------------------------------------------------------

        ! Create communicator
        ! Color may be unnecessary if using all CAM processes for CAM_mpicom
        ! but we will retain this functionality for now incase needed (hplin, 2/12/20)
        color = iam / (nPET_lat * nPET_lon)
        call MPI_comm_split(CAM_mpicom, color, iam, HCO_mpicom, RC)
        ASSERT_(RC==MPI_SUCCESS)

        ! Distribute among MPI (mp_distribute_geo in edyn_geogrid)
        ! Merged all into this huge monolithic routine..
        ! (lon_beg, lon_end, lat_beg, lat_end, 1, LM, nPET_lon, nPET_lat)
        ! (lonndx0, lonndx1, latndx0, latndx1, levndx0, levndx1, ntaski_in, ntaskj_in)

        ! Get my indices!
        my_IS = lon_beg
        my_IE = lon_end
        my_JS = lat_beg
        my_JE = lat_end

        ! Allocate task info table
        allocate(HCO_Tasks(0:nPET-1), stat=RC)

        ! Allocate 2D table of TASKS (not i j coordinates)
        allocate(HCO_petTable(-1:nPET_lon, -1:nPET_lat), stat=RC)
        ASSERT_(RC==0)

        ! 2D table of tasks communicates to MPI_PROC_NULL by default so talking
        ! to that PID has no effect in MPI comm
        HCO_petTable(:,:) = MPI_PROC_NULL

        ! Figure out ranks for the petTable, which is a table of I, J PETs
        ! with halo (hplin, 2/12/20)
        my_ID = iam
        N = 0
        do J = 0, nPET_lat-1
            do I = 0, nPET_lon-1
                HCO_petTable(I, J) = N
                if(iam == N) then
                    my_ID_I = I
                    my_ID_J = J ! Found my place in the PET table
                endif
                N = N + 1 ! move on to the next rank
            enddo

            ! Tasks are periodic in longitude (from edyn_mpi) for haloing
            ! FIXME: Check if this is true in HCO distribution. Maybe not
            HCO_petTable(-1, J) = HCO_petTable(nPET_lon-1, J)
            HCO_petTable(nPET_lon, J) = HCO_petTable(0, J)
        enddo

        ! Print some debug information on the distribution
        if( .true. ) then
            write(iulog, "('HEMCO: MPIGrid mytid=',i4,' nPET_lon,lat=',2i4,' my_ID_I,J=',2i4, &
              ' lon0,1=',2i4,' lat0,1=',2i4)") &
              my_ID,nPET_lon,nPET_lat,my_ID_I,my_ID_J,my_IS,my_IE,my_JS,my_JE

            write(iulog,"(/,'nPET=',i3,' nPET_lon=',i2,' nPET_lat=',i2,' Task Table:')") &
            nPET,nPET_lon,nPET_lat
            do J=-1,nPET_lat
                write(iulog,"('J=',i3,' HCO_petTable(:,J)=',100i3)") J,HCO_petTable(:,J)
            enddo
        endif



    end subroutine HCO_Grid_Init
!EOC
end module hco_esmf_grid