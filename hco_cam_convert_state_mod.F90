#define VERIFY_(A) if(.not.HCO_ESMF_VRFY(A,subname,__LINE__))call exit(-1)
#define ASSERT_(A) if(.not.HCO_ESMF_ASRT(A,subname,__LINE__))call exit(-1)
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_cam_convert_state_mod
!
! !DESCRIPTION: Module HCO\_CAM\_Convert\_State\_Mod handles state conversion
!  between the CAM meteorological fields and the HEMCO state.
!
!\\
!\\
! !INTERFACE:
!
module hco_cam_convert_state_mod
!
! !USES:
!
    ! ESMF function wrappers
    use hco_esmf_wrappers

    ! MPI status in CESM
    use cam_abortutils,           only: endrun      ! fatal terminator
    use spmd_utils,               only: iam, masterproc
    use cam_logfile,              only: iulog

    ! Grid information
    ! Remember - if information is pulled from hco_esmf_grid, GLOBAL INDICES MUST BE USED!
    use hco_esmf_grid,            only: my_IM, my_JM, LM
    use hco_esmf_grid,            only: my_CE
    use hco_esmf_grid,            only: my_IS, my_IE, my_JS, my_JE
    use hco_esmf_grid,            only: HCO_Grid_CAM2HCO_2D, HCO_Grid_CAM2HCO_3D
    use hco_esmf_grid,            only: HCO_Grid_HCO2CAM_2D
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs

    ! HEMCO types
    use HCO_Error_Mod,            only: sp, hp
    use HCO_State_Mod,            only: HCO_State
    use HCOX_State_Mod,           only: Ext_State

    ! Float types
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS
    use shr_kind_mod,             only: r8 => shr_kind_r8

    implicit none
    private

!
! !PUBLIC MEMBER FUNCTIONS:
!
    public       :: HCOI_Allocate_All
    public       :: CAM_GetBefore_HCOI
    public       :: CAM_RegridSet_HCOI
!
! !PRIVATE MEMBER FUNCTIONS:
!
    ! private      :: CAM_GC_ComputeMet
!
! !REMARKS:
!  The code here is structurally based on the WRF-GC coupler, but the actual
!  conversion logic is from CESM-GC.
!
!  Note that there are unique constraints here:
!  1) the CAM data structures have to be retrieved outside the GridComp in
!     HCOI_Chunk_Run phase two. Thus CAM_GetBefore_HCOI will populate the fields
!     as a memory copy from CAM state into State_CAM_* here. These are PRIVATE types.
!  2) the regridder has to run within the GridComp. Thus they are regridded using
!     another method CAM_RegridSet_HCOI and populates the PUBLIC, State_HCO_* data,
!     which is usable on the HEMCO grid.
!
!  All indices are native, and all your base are belong to us. (hplin, 12/15/20)
!
!------------------------------------------------------------------------------
!  LIST OF SUPPORTED MET FIELD CONVERSIONS
!------------------------------------------------------------------------------
!  HEMCO format
!  ----------------------------------------------------------------------------
!  AIR (airmass/AD)
!  
!  ALBD
!  F_OF_PBL
!  FROCEAN, FRSEAICE (added 8/10/22)
!  HNO3
!  NO, NO2
!  O3
!  PBLH
!  PSFC
!  QV2M (added 8/10/22)
!  SUNCOS
!  T2M
!  TK
!  TSKIN
!  USTAR (added 8/10/22)
!  U10M
!  V10M
!
!  For SeaFlux deposition:
!      GEOS-Chem        DMS/ACET    /ALD2  /MENO3/ETNO3/MOH
!      MOZART-T1        DMS/CH3COCH3/CH3CHO/--   /--   /CH3OH
!  * Some species are only available in one particular mechanism.
!  Species used for deposition are computed on the native CAM grid, and no regridding
!  is thus performed.
!
!  STUBS
!  ----------------------------------------------------------------------------
!  GWETTOP = 0
!  FRLANDIC = 0
!  FRLAKE = 0
!
!  Disabled: AIRVOL
!
! !REVISION HISTORY:
!  15 Dec 2020 - H.P. Lin    - Initial version
!  04 Feb 2021 - H.P. Lin    - Add State_GC_* for some GC-specific intermediate qtys
!  07 May 2021 - H.P. Lin    - Add state fluxes on CAM grid for deposition computation
!  10 Aug 2022 - H.P. Lin    - Update quantities QV2M, USTAR, FROCEAN, FRSEAICE
!  20 Jan 2023 - H.P. Lin    - Remove WLI for HEMCO 3.6.0+; add FRLAND.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    ! Flag for supported features
    logical                          :: feat_JValues

    ! On the CAM grid (state%psetcols, pver) (LM, my_CE)
    ! Arrays are flipped in order (k, i) for the regridder
    real(r8), pointer, public        :: State_CAM_t(:,:)
    real(r8), pointer, public        :: State_CAM_ps(:)
    real(r8), pointer, public        :: State_CAM_psdry(:)
    real(r8), pointer, public        :: State_CAM_pblh(:)

    real(r8), pointer, public        :: State_CAM_TS(:)
    real(r8), pointer, public        :: State_CAM_U10M(:)
    real(r8), pointer, public        :: State_CAM_V10M(:)
    real(r8), pointer, public        :: State_CAM_ALBD(:)
    real(r8), pointer, public        :: State_CAM_USTAR(:)

    real(r8), pointer, public        :: State_CAM_CSZA(:)

    real(r8), pointer, public        :: State_CAM_AREAM2(:)
    real(r8), pointer, public        :: State_CAM_AIRs  (:)
    real(r8), pointer, public        :: State_CAM_DELP_DRYs(:)

    ! Land fractions (converted to Olson) from CAM - 1D
    real(r8), pointer, public        :: State_CAM_FRLAND   (:)
    real(r8), pointer, public        :: State_CAM_FROCEAN  (:)
    real(r8), pointer, public        :: State_CAM_FRSEAICE (:)

    ! Chem Constituents on CAM grid
    real(r8), pointer, public        :: State_CAM_chmO3 (:,:)
    real(r8), pointer, public        :: State_CAM_chmNO (:,:)
    real(r8), pointer, public        :: State_CAM_chmNO2(:,:)
    real(r8), pointer, public        :: State_CAM_chmHNO3(:,:)

    ! For surface dep calculation, only surface needs to be copied
    real(r8), pointer, public        :: State_CAM_chmDMS(:)
    real(r8), pointer, public        :: State_CAM_chmACET(:)
    real(r8), pointer, public        :: State_CAM_chmALD2(:)
    real(r8), pointer, public        :: State_CAM_chmMENO3(:)
    real(r8), pointer, public        :: State_CAM_chmETNO3(:)
    real(r8), pointer, public        :: State_CAM_chmMOH(:)

    ! J-values from chemistry (2-D only, on surface)
    real(r8), pointer, public        :: State_CAM_JNO2(:)
    real(r8), pointer, public        :: State_CAM_JOH (:)

    ! Q at 2m [kg H2O/kg air]
    real(r8), pointer, public        :: State_CAM_QV2M(:)

    !------------------------------------------------------------------
    ! On the HEMCO grid (my_IM, my_JM, LM) or possibly LM+1
    ! HEMCO grid are set as POINTERs so it satisfies HEMCO which wants to point
    real(r8), pointer, public        :: State_GC_PSC2_DRY(:,:)  ! Dry pressure from PSC2_DRY
    real(r8), pointer, public        :: State_GC_DELP_DRY(:,:,:)! Delta dry pressure

    real(r8), pointer, public        :: State_HCO_AIR (:,:,:)   ! GC_AD, mass of dry air in grid box [kg]
    real(r8), pointer, public        :: State_HCO_AIRVOL(:,:,:) ! GC_AIRVOL, volume of grid box [m^3]

    real(r8), pointer, public        :: State_HCO_TK  (:,:,:)
    real(r8), pointer, public        :: State_HCO_PSFC(:,:)   ! Wet?
    real(r8), pointer, public        :: State_HCO_PBLH(:,:)   ! PBLH [m]

    real(r8), pointer, public        :: State_HCO_TS(:,:)
    real(r8), pointer, public        :: State_HCO_U10M(:,:)
    real(r8), pointer, public        :: State_HCO_V10M(:,:)
    real(r8), pointer, public        :: State_HCO_ALBD(:,:)
    real(r8), pointer, public        :: State_HCO_USTAR(:,:)
    real(r8), pointer, public        :: State_HCO_F_OF_PBL(:,:,:)

    real(r8), pointer, public        :: State_HCO_CSZA(:,:)

    real(r8), pointer, public        :: State_HCO_FRLAND  (:,:)
    real(r8), pointer, public        :: State_HCO_FRLANDIC(:,:)
    real(r8), pointer, public        :: State_HCO_FROCEAN (:,:)
    real(r8), pointer, public        :: State_HCO_FRSEAICE(:,:)
    real(r8), pointer, public        :: State_HCO_FRLAKE  (:,:)

    ! Chem Constituents on HEMCO grid
    real(r8), pointer, public        :: State_HCO_chmO3 (:,:,:)
    real(r8), pointer, public        :: State_HCO_chmNO (:,:,:)
    real(r8), pointer, public        :: State_HCO_chmNO2(:,:,:)
    real(r8), pointer, public        :: State_HCO_chmHNO3(:,:,:)

    ! J-values from chemistry (passed in reverse direction)
    real(r8), pointer, public        :: State_HCO_JNO2(:,:)
    real(r8), pointer, public        :: State_HCO_JOH (:,:)


    real(r8), pointer, public        :: State_HCO_QV2M(:,:)

    ! constituent indices:
    integer                          :: id_O3, id_NO, id_NO2, id_HNO3
    integer                          :: id_H2O, id_Q
    integer                          :: id_DMS, id_ACET, id_ALD2, id_MENO3, id_ETNO3, id_MOH

contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Allocate_All
!
! !DESCRIPTION: HCOI\_Allocate\_All allocates temporary met fields for use in
!  HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Allocate_All()
!
! !USES:
!
        use cam_logfile,      only: iulog
        use spmd_utils,       only: masterproc

        use HCO_ESMF_Grid,    only: AREA_M2, HCO_Grid_HCO2CAM_2D
!
! !REMARKS:
!  Fields are allocated here after initialization of the hco\_esmf\_grid.
!  They are regridded inside the gridcomp.
!  A "state conversion" to convert CAM met fields to GEOSFP format will
!  need to be coordinated with fritzt later down the road (hplin, 3/27/20)
!
!  Note that arrays in CAM format stored in the HEMCO interface are (k, i) idxd
!  for the regridder
!
! !REVISION HISTORY:
!  31 Mar 2020 - H.P. Lin    - Initial version
!  10 Aug 2022 - H.P. Lin    - Update FROCEAN, FRSEAICE for HEMCO 3.4.0+
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Allocate_All'
        integer                      :: RC                   ! ESMF return code

        ! Assume success
        RC = ESMF_SUCCESS

        ! Sanity check
        !if(masterproc) then
        !    write(iulog,*) "> HCOI_Allocate_All entering"
        !endif

        ! PBL height [m]
        ! Comes from pbuf
        allocate(State_CAM_pblh(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_PBLH(my_IM, my_JM), stat=RC)
        ASSERT_(RC==0)

        ! Grid box area [m2]
        ! Used for calculation of deposition fluxes directly on CAM grid
        allocate(State_CAM_AREAM2(my_CE), stat=RC)
        ASSERT_(RC==0)

        ! Surface grid box weight [kg]
        allocate(State_CAM_AIRs(my_CE), stat=RC)
        ASSERT_(RC==0)

        ! Surface delta grid box pressure differential [hPa]
        allocate(State_CAM_DELP_DRYs(my_CE), stat=RC)
        ASSERT_(RC==0)

        ! Surface pressure (wet) [Pa]
        allocate(State_CAM_ps(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_PSFC(my_IM, my_JM), stat=RC)
        ASSERT_(RC==0)

        ! Surface pressure (dry) [hPa]
        allocate(State_CAM_psdry(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_GC_PSC2_DRY(my_IM, my_JM), stat=RC)
        ASSERT_(RC==0)

        ! T
        allocate(State_CAM_t (LM, my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_TK(my_IM, my_JM, LM), stat=RC)
        ASSERT_(RC==0)


        ! For HEMCO extensions...
        allocate(State_CAM_TS   (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_U10M (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_V10M (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_ALBD (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_USTAR(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_CSZA (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_FRLAND(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_FROCEAN(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_FRSEAICE(my_CE), stat=RC)
        ASSERT_(RC==0)

        ! QV2M
        allocate(State_CAM_QV2M(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_QV2M(my_IM, my_JM), stat=RC)
        ASSERT_(RC==0)

        ! Constituents
        allocate(State_CAM_chmO3(LM, my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_CAM_chmNO(LM, my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_CAM_chmNO2(LM, my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_CAM_chmHNO3(LM, my_CE), stat=RC)
        ASSERT_(RC==0)

        ! Constituents for deposition flux, copy sfc only
        allocate(State_CAM_chmDMS(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_chmACET(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_chmALD2(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_chmMENO3(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_chmETNO3(my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_chmMOH(my_CE), stat=RC)
        ASSERT_(RC==0)

        ! J-values
        allocate(State_CAM_JNO2 (my_CE), stat=RC)
        ASSERT_(RC==0)
        allocate(State_CAM_JOH  (my_CE), stat=RC)
        ASSERT_(RC==0)

        ! On HEMCO grid
        allocate(State_HCO_AIR(my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_TS(my_IM, my_JM), stat=RC)
        allocate(State_HCO_U10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_V10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_ALBD(my_IM, my_JM), stat=RC)
        allocate(State_HCO_CSZA(my_IM, my_JM), stat=RC)
        allocate(State_HCO_USTAR(my_IM, my_JM), stat=RC)
        allocate(State_HCO_F_OF_PBL(my_IM, my_JM, LM), stat=RC)
        allocate(State_GC_DELP_DRY(my_IM, my_JM, LM), stat=RC)

        allocate(State_HCO_chmO3 (my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmNO (my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmNO2(my_IM, my_JM, LM), stat=RC)
        allocate(State_HCO_chmHNO3(my_IM, my_JM, LM), stat=RC)

        allocate(State_HCO_JOH (my_IM, my_JM), stat=RC)
        allocate(State_HCO_JNO2(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLAND(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLANDIC(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FROCEAN(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRSEAICE(my_IM, my_JM), stat=RC)
        allocate(State_HCO_FRLAKE(my_IM, my_JM), stat=RC)

        ! Clear values
        State_HCO_AIR(:,:,:) = 0.0_r8
        State_HCO_PBLH(:,:) = 0.0_r8
        State_HCO_PSFC(:,:) = 0.0_r8
        State_GC_PSC2_DRY(:,:) = 0.0_r8
        State_HCO_TK(:,:,:) = 0.0_r8
        State_HCO_TS(:,:) = 0.0_r8
        State_HCO_U10M(:,:) = 0.0_r8
        State_HCO_V10M(:,:) = 0.0_r8
        State_HCO_ALBD(:,:) = 0.0_r8
        State_HCO_USTAR(:,:) = 0.0_r8
        State_HCO_CSZA(:,:) = 0.0_r8
        State_HCO_F_OF_PBL(:,:,:) = 0.0_r8
        State_GC_DELP_DRY(:,:,:) = 0.0_r8

        State_HCO_chmO3(:,:,:) = 0.0_r8
        State_HCO_chmNO(:,:,:) = 0.0_r8
        State_HCO_chmNO2(:,:,:) = 0.0_r8
        State_HCO_chmHNO3(:,:,:) = 0.0_r8

        State_HCO_QV2M(:,:) = 0.0_r8

        State_CAM_chmDMS(:) = 0.0_r8
        State_CAM_chmACET(:) = 0.0_r8
        State_CAM_chmALD2(:) = 0.0_r8
        State_CAM_chmMENO3(:) = 0.0_r8
        State_CAM_chmETNO3(:) = 0.0_r8
        State_CAM_chmMOH(:) = 0.0_r8

        State_CAM_JNO2(:) = 0.0_r8
        State_CAM_JOH (:) = 0.0_r8
        State_CAM_QV2M(:) = 0.0_r8
        State_CAM_USTAR(:) = 0.0_r8

        State_HCO_JNO2(:,:) = 0.0_r8
        State_HCO_JOH(:,:) = 0.0_r8

        ! Unsupported fields in CESM are directly assigned to zero.
        ! These, if ever available, will be updated along with chemistry.F90
        ! under src/chemistry/geoschem. (hplin, 1/20/23)
        State_HCO_FRLANDIC(:,:) = 0.0_r8
        State_HCO_FRLAKE(:,:) = 0.0_r8

        ! Populate persistent values
        call HCO_Grid_HCO2CAM_2D(AREA_M2, State_CAM_AREAM2)

    end subroutine HCOI_Allocate_All
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CAM_GetBefore_HCOI
!
! !DESCRIPTION: CAM\_GetBefore\_HCOI populates the internal copy of CAM state
!  within the conversion module to prepare for regridding within the gridcomp.
!\\
!\\
! !INTERFACE:
!
    subroutine CAM_GetBefore_HCOI(cam_in, phys_state, pbuf2d, phase, HcoState, ExtState)
!
! !USES:
!
        ! Pbuf wrappers by hplin
        use hco_cam_exports,only: HCO_Export_Pbuf_QueryField

        ! Type descriptors
        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physics_buffer, only: physics_buffer_desc

        ! Physics grid
        use phys_grid,      only: get_ncols_p
        use phys_grid,      only: get_rlon_all_p, get_rlat_all_p

        use ppgrid,         only: pcols                   ! max. # of columns in chunk

        ! CAM physics buffer (some fields are here and some are in phys state)
        use physics_buffer, only: pbuf_get_chunk, pbuf_get_field
        use physics_buffer, only: pbuf_get_index

        ! Time description and zenith angle data
        use orbit,          only: zenith
        use time_manager,   only: get_curr_calday

        ! Constituent information to retrieve from physics state
        use constituents,   only: cnst_get_ind

        ! Output and mpi
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid, iam
!
! !INPUT PARAMETERS:
!
        type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
        type(physics_buffer_desc), pointer :: pbuf2d(:,:)
        integer, intent(in)                :: phase               ! 1, 2

        type(HCO_State), pointer           :: HcoState
        type(Ext_State), pointer           :: ExtState
!
! !REMARKS:
!  Note that arrays in CAM format stored in the HEMCO interface are (k, i) idxd
!  for the regridder
!
!  Also needs HEMCO and HEMCO extensions state information in order to check
!  whether we actually need to populate required meteorology fields
!
!  Note for CSZA:
!  - We do not have access to the state here, so we need to get geo data and zenith
!    using a different method, by looping through all chunks.
!
! !REVISION HISTORY:
!  16 Dec 2020 - H.P. Lin    - Initial version
!  04 Feb 2021 - H.P. Lin    - Add CSZA calculation with geographical data
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'CAM_GetBefore_HCOI'
        integer                      :: RC                   ! ESMF return code

        integer                      :: I, J, K, lchnk            ! Loop idx
        integer                      :: ncol

        type(physics_buffer_desc), pointer :: pbuf_chnk(:)        ! slice of pbuf

        ! temp 2-D data slice (pcol), points to raw
        real(r8), pointer            :: pbuf_tmp_pblh(:)
        real(r8), pointer            :: pbuf_tmp_JNO2(:), pbuf_tmp_JOH(:)

        ! pbuf indices:
        integer                      :: index_pblh, index_JNO2, index_JOH

        ! Temporary geographical indices, allocated to max size (pcols)
        ! need to use actual column # ncol = get_ncols_p to fill to my_CE, which is exact
        real(r8)                     :: lchnk_rlats(1:pcols), lchnk_rlons(1:pcols)
        real(r8)                     :: lchnk_zenith(1:pcols)

        ! Current calday for SZA
        real(r8)                     :: calday

        ! is this first timestep? skip reading pbuf from certain data if so
        logical, save                :: FIRST = .true.

        !----------------------------------------------------
        ! Assume success
        RC = ESMF_SUCCESS

        if(masterproc .and. FIRST) then
            write(iulog,*) "> CAM_GetBefore_HCOI entering"
        endif

        ! Regrid necessary physics quantities from the CAM grid to the HEMCO grid
        ! Phase 0: Prepare necessary pbuf indices to retrieve data from the buffer...

        ! TODO: This prep phase might possibly be better moved into initialization.
        ! Keeping it here for now but later can be optimized (hplin, 3/31/20)
        index_pblh = pbuf_get_index('pblh')

        ! J-values - might be -1 if not existent (i.e. not in CESM-GC)
        ! Note: Assuming they exist in pairs (if index_JNO2 exists, all do)
        ! (hplin, 3/3/21)
        index_JNO2 = pbuf_get_index('HCO_IN_JNO2', RC)
        index_JOH  = pbuf_get_index('HCO_IN_JOH', RC)
        RC = ESMF_SUCCESS ! dummy

        ! Set feature flag
        feat_JValues = index_JNO2 > 0

        if(masterproc .and. FIRST) write(iulog,*) "feat_JValues:", feat_JValues

        ! Get calday for cosza (current time, not midpoint of dt)
        calday = get_curr_calday()

        ! TODO: Move constituent index calculation (which only needs to be initialized once)
        ! to a place where it only runs once ...
        ! Setup constituent indices so their concentrations can be retrieved
        ! from state%q (MMR)
        call cnst_get_ind('O3', id_O3)
        call cnst_get_ind('NO', id_NO)
        call cnst_get_ind('NO2', id_NO2)
        call cnst_get_ind('HNO3', id_HNO3)

        ! Get constitutent index for specific humidity
        call cnst_get_ind('Q', id_Q)
        ! call cnst_get_ind('H2O', id_H2O)
        ! id_H2O not used for now, and also not present in CAM-chem. hplin, 9/9/22

        ! Retrieve optional - for deposition - constituent IDs
        call cnst_get_ind('DMS', id_DMS, abort=.False.)
        call cnst_get_ind('MENO3', id_MENO3, abort=.False.)
        call cnst_get_ind('ETNO3', id_ETNO3, abort=.False.)
        call cnst_get_ind('ACET', id_ACET, abort=.False.)
        if(id_ACET <= 0) then
            call cnst_get_ind('CH3COCH3', id_ACET, abort=.False.)
        endif
        call cnst_get_ind('ALD2', id_ALD2, abort=.False.)
        if(id_ALD2 <= 0) then
            call cnst_get_ind('CH3CHO', id_ALD2, abort=.False.)
        endif
        call cnst_get_ind('MOH', id_MOH, abort=.False.)
        if(id_MOH <= 0) then
            call cnst_get_ind('CH3OH', id_MOH, abort=.False.)
        endif

        ! Phase 1: Store the fields in hemco_interface (copy)
        I = 0
        do lchnk = begchunk, endchunk    ! loop over all chunks in the physics grid
            ncol = get_ncols_p(lchnk)    ! columns per chunk
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk) ! get pbuf for this chunk...

            ! Gather geographical information
            ! if(masterproc) write(iulog,*) "* hplin debug in lchnk = ", lchnk, ", ncol = ", ncol, " (max pcols = ", pcols, "), my_CE = ", my_CE
            call get_rlat_all_p(lchnk, ncol, lchnk_rlats)
            call get_rlon_all_p(lchnk, ncol, lchnk_rlons)

            ! Compute zenith for chunk
            ! (could also do it all at once but it would require a separate buffer to store ll...)
            ! FIXME hplin: this slicing might be a little inefficient
            call zenith(calday, lchnk_rlats(1:ncol), lchnk_rlons(1:ncol), lchnk_zenith(1:ncol), ncol)

            ! Gather data from pbuf before we proceed
            ! 2-D Fields: Write directly to pointer (does not need alloc).
            call pbuf_get_field(pbuf_chnk, index_pblh, pbuf_tmp_pblh)

            if(feat_JValues) then
                call pbuf_get_field(pbuf_chnk, index_JNO2, pbuf_tmp_JNO2)
                call pbuf_get_field(pbuf_chnk, index_JOH,  pbuf_tmp_JOH )
            endif

            ! Loop through the chunks and levs now
            ! Sanity check: Left indices should be I, and r.h.s. should be J!!!
            ! If it differs, you are likely wrong and reading out of bounds.
            ! Two hours of debugging a 1.5GPa surface pressure entry has resulted
            ! in the above comments. (hplin, 3/31/20)
            do J = 1, ncol               ! loop over columns in the chunk
                I = I + 1                ! advance one column

                ! 3-D Fields
                do K = 1, LM             !        chunk    col, lev
                    State_CAM_t(K,I) = phys_state(lchnk)%t(J,K)

                    ! FIXME: hplin, check if indices actually exist!!
                    State_CAM_chmO3(K,I) = phys_state(lchnk)%q(J,K,id_O3)
                    State_CAM_chmNO(K,I) = phys_state(lchnk)%q(J,K,id_NO)
                    State_CAM_chmNO2(K,I) = phys_state(lchnk)%q(J,K,id_NO2)
                    State_CAM_chmHNO3(K,I) = phys_state(lchnk)%q(J,K,id_HNO3)
                enddo

                ! Verify and retrieve surface fluxes for deposition, if available (hplin, 5/7/21)
                if(id_DMS  > 0) State_CAM_chmDMS (I)   = phys_state(lchnk)%q(J,LM,id_DMS )
                if(id_MENO3 > 0) State_CAM_chmMENO3(I) = phys_state(lchnk)%q(J,LM,id_MENO3)
                if(id_ETNO3 > 0) State_CAM_chmETNO3(I) = phys_state(lchnk)%q(J,LM,id_ETNO3)
                if(id_ACET > 0) State_CAM_chmACET(I)   = phys_state(lchnk)%q(J,LM,id_ACET)
                if(id_ALD2 > 0) State_CAM_chmALD2(I)   = phys_state(lchnk)%q(J,LM,id_ALD2)
                if(id_MOH  > 0) State_CAM_chmMOH (I)   = phys_state(lchnk)%q(J,LM,id_MOH )

                !----------------------------------------
                ! 2-D Fields
                !----------------------------------------

                ! DEBUG: Write to CSZA as a test for latitude to make sure we are doing correctly
                ! State_CAM_CSZA(I) = lchnk_rlats(J)

                ! QV2M [kg H2O/kg air] (MMR -- this is converted to VMR by *MWdry/18 in SeaSalt)
                ! at 2M, roughly surface ~ LM (1 is TOA)
                ! (hplin, 8/10/22)
                State_CAM_QV2M(I) = phys_state(lchnk)%q(J,LM,id_Q)

                ! PBLH [m]
                State_CAM_pblh(I) = pbuf_tmp_pblh(J)

                ! COSZA Cosine of zenith angle [1]
                State_CAM_CSZA(I) = lchnk_zenith(J)

                ! USTAR Friction velocity [m/s]
                State_CAM_USTAR(I) = cam_in(lchnk)%fv(J) * cam_in(lchnk)%landFrac(J) + cam_in(lchnk)%uStar(J) * (1.0_r8 - cam_in(lchnk)%landFrac(J))

                ! Sea level pressure [Pa] (note difference in units!!)
                State_CAM_ps(I) = phys_state(lchnk)%ps(J)

                ! Dry pressure [hPa] (Pa -> hPa, x0.01)
                State_CAM_psdry(I) = phys_state(lchnk)%psdry(J) * 0.01_r8

                ! Surface temperature [K] - use both for T2M and TSKIN for now according to CESM-GC,
                ! hplin 12/21/2020
                if(ExtState%T2M%DoUse .or. ExtState%TSKIN%DoUse) then
                    State_CAM_TS(I) = cam_in(lchnk)%TS(J)
                endif

                ! 10M E/W and N/S wind speed [m/s] (fixme: use pver?)
                if(ExtState%U10M%DoUse) then
                    State_CAM_U10M(I) = phys_state(lchnk)%U(J, pver)
                    State_CAM_V10M(I) = phys_state(lchnk)%V(J, pver)
                endif

                ! Visible surface albedo [1]
                if(ExtState%ALBD%DoUse) then
                    State_CAM_ALBD(I) = cam_in(lchnk)%asdir(J)
                endif

                ! Converted-to-Olson land fractions [1] (hplin, 8/10/22)
                if(ExtState%FRLAND%DoUse) then
                    State_CAM_FRLAND(I) = cam_in(lchnk)%landFrac(J)
                endif

                ! FRLANDIC unsupported

                if(ExtState%FROCEAN%DoUse) then
                    State_CAM_FROCEAN(I) = cam_in(lchnk)%ocnFrac(J) + cam_in(lchnk)%iceFrac(J)
                endif

                if(ExtState%FRSEAICE%DoUse) then
                    State_CAM_FRSEAICE(I) = cam_in(lchnk)%iceFrac(J)
                endif

                ! FRLAKE unsupported
                ! FRSNO unsupported

                ! J-values (from CESM-GC only, at the moment. Added to CAM-chem/mozart 5/16/21)
                if(feat_JValues .and. .not. FIRST) then
                    State_CAM_JNO2(I) = pbuf_tmp_JNO2(J)
                    State_CAM_JOH (I) = pbuf_tmp_JOH (J)
                    ! write(6,*) "hplin debug***: jno2, joh j-v", I,J, pbuf_tmp_JNO2(J), pbuf_tmp_JOH(J)
                endif
            enddo
        enddo

        ! The below are stubs and are not regridded or processed for now due to lack of data

        if(FIRST) then
            FIRST = .false.

            if(masterproc) then
                write(iulog,*) "> CAM_GetBefore_HCOI finished"
            endif
        endif

    end subroutine CAM_GetBefore_HCOI
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CAM_RegridSet_HCOI
!
! !DESCRIPTION: CAM\_GetBefore\_HCOI populates the internal copy of CAM state
!  within the conversion module to prepare for regridding within the gridcomp.
!\\
!\\
! !INTERFACE:
!
    subroutine CAM_RegridSet_HCOI(HcoState, ExtState, Phase)
!
! !USES:
!
        ! Type descriptors
        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physics_buffer, only: physics_buffer_desc

        ! Physics grid
        use phys_grid,      only: get_ncols_p

        ! CAM physics buffer (some fields are here and some are in phys state)
        use physics_buffer, only: pbuf_get_chunk, pbuf_get_field
        use physics_buffer, only: pbuf_get_index

        ! Output and mpi
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid, iam

        ! HEMCO output container state
        USE HCOX_State_Mod, only: ExtDat_Set

        ! Vertical grid specification
        use HCO_ESMF_Grid,  only: Ap, Bp, AREA_M2

        ! HEMCO utilities
        use HCO_GeoTools_Mod, only: HCO_GetSUNCOS
!
! !INPUT PARAMETERS:
!

        type(HCO_State), pointer           :: HcoState
        type(Ext_State), pointer           :: ExtState
        integer                            :: Phase

!
! !REMARKS:
!  Note that arrays in CAM format stored in the HEMCO interface are (k, i) idxd
!  for the regridder
!
!  Also the vertical does not need to be handled here, the CAM2HCO_3D will handle
!  it very nicely for you.
!
! !REVISION HISTORY:
!  16 Dec 2020 - H.P. Lin    - Initial version
!  21 Dec 2020 - H.P. Lin    - Now also sets ExtState fields appropriately
!  04 Feb 2021 - H.P. Lin    - Also compute air quantities (might have to move to separate later)
!  03 Mar 2021 - H.P. Lin    - Now also computes PBL quantities
!  03 Mar 2021 - H.P. Lin    - Now split into two phases
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'CAM_GetBefore_HCOI'
        integer                      :: RC                   ! ESMF return code

        logical, save                :: FIRST = .true.
        integer, save                :: nCalls = 0

        integer                      :: I, J, L              ! Loop index

        ! Temporary quantities needed for PBL computation
        real(r8)                     :: BLTOP, BLTHIK, DELP
        integer                      :: LTOP

        ! Physical constants from physconstants.F90 (from GEOS-Chem)
        ! TODO: Verify consistency with CAM model! (hplin, 3/3/21)
        real(r8), parameter          :: G0_100 = 100.e+0_r8 / 9.80665e+0_r8
        real(r8), parameter          :: SCALE_HEIGHT = 7600.0_r8

        ! Assume success
        RC = ESMF_SUCCESS

        nCalls = nCalls + 1

        if(masterproc .and. nCalls < 10) then
            write(iulog,*) "> CAM_RegridSet_HCOI entering", phase
        endif

        !-----------------------------------------------------------------------
        ! Regrid necessary physics quantities from the CAM grid to the HEMCO grid
        ! Phase 1: Regrid
        !-----------------------------------------------------------------------
        if(Phase == 1) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_ps,     State_HCO_PSFC   )
            call HCO_Grid_CAM2HCO_2D(State_CAM_psdry,  State_GC_PSC2_DRY)
            call HCO_Grid_CAM2HCO_2D(State_CAM_pblh,   State_HCO_PBLH   )
            call HCO_Grid_CAM2HCO_3D(State_CAM_t,      State_HCO_TK     )

            if(masterproc .and. nCalls < 10) then
                write(iulog,*) "> CAM_RegridSet_HCOI exiting phase 1"
            endif

            return
        endif

        ! Below only Phase 2...

        !-----------------------------------------------------------------------
        ! Compute air quantities (hplin, 2/4/21)
        !-----------------------------------------------------------------------

        ! Unified loop 1 (LJI)
        do L = 1, LM
        do J = 1, my_JM
        do I = 1, my_IM
            ! Calculate DELP_DRY (from pressure_mod)
            State_GC_DELP_DRY(I,J,L) = (Ap(L)   + (Bp(L)   * State_GC_PSC2_DRY(I,J))) - &
                                       (Ap(L+1) + (Bp(L+1) * State_GC_PSC2_DRY(I,J)))

            ! Calculate AD (AIR mass)
            ! Note that AREA_M2 are GLOBAL indices so you need to perform offsetting!!
            !
            ! DELP_DRY is in [hPa]. G0_100 is 100/g, converts to [Pa], and divides by [m/s2] (accel to kg)
            State_HCO_AIR(I,J,L) = State_GC_DELP_DRY(I,J,L) * G0_100 * AREA_M2(my_IS+I-1,my_JS+J-1)
        enddo
        enddo
        enddo

        ! Populate CAM information equivalent
        call HCO_Grid_HCO2CAM_2D(State_GC_DELP_DRY(:,:,1), State_CAM_DELP_DRYs)
        call HCO_Grid_HCO2CAM_2D(State_HCO_AIR(:,:,1), State_CAM_AIRs)

        !-----------------------------------------------------------------------
        ! Surface temperature [K] - use both for T2M and TSKIN for now according to CESM-GC,
        ! hplin 12/21/2020
        if(ExtState%T2M%DoUse .or. ExtState%TSKIN%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_TS, State_HCO_TS    )

            ! Point to appropriate location
            call ExtDat_Set(HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_TS)

            call ExtDat_Set(HcoState, ExtState%T2M,  'T2M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_TS)
        
            ! if(masterproc) then
            !     write(iulog,*) "HCO CAM_Convert_State: after ExtDat_Set T2M"
            ! endif
        endif

        ! 10M E/W and N/S wind speed [m/s] (fixme: use pver?)
        if(ExtState%U10M%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_U10M, State_HCO_U10M)
            call HCO_Grid_CAM2HCO_2D(State_CAM_V10M, State_HCO_V10M)

            call ExtDat_Set(HcoState, ExtState%U10M,  'U10M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_U10M)

            call ExtDat_Set(HcoState, ExtState%V10M,  'V10M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_V10M)
        endif

        ! Cos of Zenith Angle [1]
        if(ExtState%SUNCOS%DoUse) then
            ! call HCO_Grid_CAM2HCO_2D(State_CAM_CSZA, State_HCO_CSZA)

            ! Use native CSZA from HEMCO for consistency?
            call HCO_GetSUNCOS(HcoState, State_HCO_CSZA, 0, RC)

            call ExtDat_Set(HcoState, ExtState%SUNCOS,'SUNCOS_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_CSZA)
        endif

        ! Visible surface albedo [1]
        if(ExtState%ALBD%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_ALBD, State_HCO_ALBD)

            call ExtDat_Set(HcoState, ExtState%ALBD,  'ALBD_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_ALBD)
        endif

        ! Friction velocity
        if(ExtState%USTAR%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_USTAR, State_HCO_USTAR)

            call ExtDat_Set(HcoState, ExtState%USTAR, 'USTAR_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_USTAR)
        endif

        ! Air mass [kg]
        if(ExtState%AIR%DoUse) then
            ! This is computed above using GC routines for air quantities, so
            ! it does not necessitate a regrid from CAM.

            call ExtDat_Set(HcoState, ExtState%AIR,   'AIRMASS_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_AIR)
        endif

        ! Constituents [MMR]
        if(ExtState%O3%DoUse) then 
            call HCO_Grid_CAM2HCO_3D(State_CAM_chmO3, State_HCO_chmO3)

            call ExtDat_Set(HcoState, ExtState%O3,   'HEMCO_O3_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmO3)
        endif

        if(ExtState%NO2%DoUse) then 
            call HCO_Grid_CAM2HCO_3D(State_CAM_chmNO2, State_HCO_chmNO2)

            call ExtDat_Set(HcoState, ExtState%NO2,   'HEMCO_NO2_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmNO2)
        endif

        if(ExtState%NO%DoUse) then 
            call HCO_Grid_CAM2HCO_3D(State_CAM_chmNO, State_HCO_chmNO)

            call ExtDat_Set(HcoState, ExtState%NO,   'HEMCO_NO_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_chmNO)
        endif

        ! MMR of H2O at 2m [kg H2O/kg air] (2-D only)
        if(ExtState%QV2M%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_QV2M, State_HCO_QV2M)

            call ExtDat_Set(HcoState, ExtState%QV2M,  'QV2M_FOR_EMIS', &
                            RC,       FIRST,          State_HCO_QV2M)
        endif

        ! hplin 3/3/21: Disabling HNO3 for now - not actually used by HEMCO,
        ! not even in ParaNOx (code is commented out)

        ! if(ExtState%HNO3%DoUse) then 
        !     call HCO_Grid_CAM2HCO_3D(State_CAM_chmHNO3, State_HCO_chmHNO3)

        !     call ExtDat_Set(HcoState, ExtState%HNO3,   'HEMCO_HNO3_FOR_EMIS', &
        !                     RC,       FIRST,          State_HCO_chmHNO3)
        ! endif

        !-------------------------------------------------------
        ! Compute the PBL fraction (on HEMCO grid) so we avoid an extra regrid.
        ! This code largely ported from GeosCore/pbl_mix_mod.F90
        ! (hplin, 3/3/21)
        !
        ! Note that L = 1 is now surface, since in HEMCO grid the atmos
        ! has already been inverted from the CAM approach (1 as TOA)
        !
        ! Some temporary quantities are needed:
        !
        !  BLTOP    = Pressure at PBL top [hPa]
        !  BLTHIK   = PBL thickness [hPa]
        !  DELP     = Thickness of grid box (I,J,L) [hPa]
        !  LTOP     = Level of PBL top [nlev]
        !
        !
        ! Unified Loop 2 (JI, L)
        ! Again, note that loop indices I, J are 1~my_*M because the indices
        ! are not global here (we are in the HEMCO decomp at this point)
        !
        ! Remember, ONLY HCO_ESMF_GRID quantities are GLOBAL indices. Otherwise,
        ! there is a separate HEMCO index for each PE.
        !
        ! TODO: OMP parallel do default(shared) private(i, j, l, bltop, blthik, ltop, delp)

        do J = 1, my_JM
            do I = 1, my_IM
                ! use barometric law for pressure at PBL top
                BLTOP = HcoState%Grid%PEDGE%Val(I,J,1) * EXP(-State_HCO_PBLH(I, J)/SCALE_HEIGHT)

                ! PBL thickness [hPa]
                BLTHIK = HcoState%Grid%PEDGE%Val(I,J,1) - BLTOP

                ! Now, find the PBL top level
                do L = 1, LM
                    if(BLTOP > HcoState%Grid%PEDGE%Val(I,J,L+1)) then
                        LTOP = L
                        exit
                    endif
                enddo

                ! Find the fraction of grid box (I,J,L) within the PBL
                do L = 1, LM
                    DELP = HcoState%Grid%PEDGE%Val(I,J,L) - HcoState%Grid%PEDGE%Val(I,J,L+1)
                    ! ...again, PEDGE goes up to LM+1

                    if(L < LTOP) then
                        ! grid cell lies completely below the PBL top
                        State_HCO_F_OF_PBL(I,J,L) = DELP / BLTHIK
                    else if(L == LTOP) then
                        ! grid cell straddles PBL top
                        State_HCO_F_OF_PBL(I,J,L) = (HcoState%Grid%PEDGE%Val(I,J,L) - BLTOP) / BLTHIK
                    else
                        ! grid cells lies completely above the PBL top
                        State_HCO_F_OF_PBL(I,J,L) = 0.0
                    endif

                    ! write(6,*) "I,J/L", I, J, L, ": F_OF_PBL", State_HCO_F_OF_PBL(I,J,L)
                enddo
            enddo
        enddo

        if(ExtState%FRAC_OF_PBL%DoUse) then
            call ExtDat_Set(HcoState, ExtState%FRAC_OF_PBL,  'FRAC_OF_PBL_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_F_OF_PBL)
        endif

        ! J-values - if supported
        if(ExtState%JOH%DoUse .or. ExtState%JNO2%DoUse) then
            if(feat_JValues) then
                call HCO_Grid_CAM2HCO_2D(State_CAM_JOH, State_HCO_JOH)

                call ExtDat_Set(HcoState, ExtState%JOH,  'JOH_FOR_EMIS', &
                                RC,       FIRST,         State_HCO_JOH)

                call HCO_Grid_CAM2HCO_2D(State_CAM_JNO2, State_HCO_JNO2)

                call ExtDat_Set(HcoState, ExtState%JNO2, 'JNO2_FOR_EMIS', &
                                RC,       FIRST,         State_HCO_JNO2)
            else
                call endrun('hco_cam_convert_state_mod: requested J-values fields but J-value field unsupported, check chem or disable ParaNOx')
            endif
        endif

        if(ExtState%FRLAND%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_FRLAND, State_HCO_FRLAND)
            call ExtDat_Set(HcoState, ExtState%FRLAND,       'FRLAND_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLAND)
        endif

        if(ExtState%FRLANDIC%DoUse) then
            ! Unsupported - State_HCO_FRLANDIC is always zero
            call ExtDat_Set(HcoState, ExtState%FRLANDIC,     'FRLANDIC_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLANDIC)
        endif

        if(ExtState%FROCEAN%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_FROCEAN, State_HCO_FROCEAN)
            call ExtDat_Set(HcoState, ExtState%FROCEAN,      'FROCEAN_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FROCEAN)
        endif

        if(ExtState%FRSEAICE%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_FRSEAICE, State_HCO_FRSEAICE)
            call ExtDat_Set(HcoState, ExtState%FRSEAICE,     'FRSEAICE_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRSEAICE)
        endif

        if(ExtState%FRLAKE%DoUse) then
            ! Unsupported - State_HCO_FRLAKE is always zero
            call ExtDat_Set(HcoState, ExtState%FRLAKE,       'FRLAKE_FOR_EMIS', &
                            RC,       FIRST,                 State_HCO_FRLAKE)
        endif

        if(FIRST) then
            FIRST = .false.
        endif
        
        if(masterproc .and. nCalls < 10) then
            write(iulog,*) "> CAM_RegridSet_HCOI finished"
        endif

        ! Debugging code
        ! write(6,*) "HCO - PBLH(:,:), min, max", minval(State_HCO_PBLH(:,:)), maxval(State_HCO_PBLH(:,:))
        ! write(6,*) State_HCO_PBLH(:,:)

        ! write(6,*) "CAM - PBLH(:,:), min, max", minval(State_CAM_pblh(:)), maxval(State_CAM_pblh(:))
        ! write(6,*) State_CAM_pblh(:)

        ! write(6,*) "HCO - PSFC(:,:), min, max", minval(State_HCO_PSFC(:,:)), maxval(State_HCO_PSFC(:,:))
        ! write(6,*) State_HCO_PSFC(:,:)

        ! write(6,*) "CAM - ps(:,:), min, max", minval(State_CAM_ps(:)), maxval(State_CAM_ps(:))
        ! write(6,*) State_CAM_ps(:)

        ! write(6,*) "HCO - TK(:,:,1), min, max (gl)", minval(State_HCO_TK(:,:,:)), maxval(State_HCO_TK(:,:,:))
        ! write(6,*) State_HCO_TK(:,:,1)

        ! write(6,*) "CAM - TK(:,:), min, max", minval(State_CAM_t(1,:)), maxval(State_CAM_t(1,:))
        ! write(6,*) State_CAM_t(1,:)
    end subroutine CAM_RegridSet_HCOI
!EOC
end module hco_cam_convert_state_mod
