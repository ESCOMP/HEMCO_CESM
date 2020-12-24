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
    use hco_esmf_grid,            only: my_IM, my_JM, LM
    use hco_esmf_grid,            only: my_CE
    use hco_esmf_grid,            only: HCO_Grid_CAM2HCO_2D, HCO_Grid_CAM2HCO_3D
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
! !REVISION HISTORY:
!  15 Dec 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    ! On the CAM grid (state%psetcols, pver) (LM, my_CE)
    ! Arrays are flipped in order (k, i) for the regridder
    real(r8), pointer                :: State_CAM_t(:,:)
    real(r8), pointer                :: State_CAM_ps(:)
    real(r8), pointer                :: State_CAM_pblh(:)

    real(r8), pointer                :: State_CAM_TS(:)
    real(r8), pointer                :: State_CAM_U10M(:)
    real(r8), pointer                :: State_CAM_V10M(:)
    real(r8), pointer                :: State_CAM_ALBD(:)
    real(r8), pointer                :: State_CAM_LWI(:)

    ! On the HEMCO grid (my_IM, my_JM, LM) or possibly LM+1
    ! HEMCO grid are set as POINTERs so it satisfies HEMCO which wants to point
    real(r8), pointer, public        :: State_HCO_TK  (:,:,:)
    real(r8), pointer, public        :: State_HCO_PSFC(:,:)   ! Wet?
    real(r8), pointer, public        :: State_HCO_PBLH(:,:)   ! PBLH [m]

    real(r8), pointer, public        :: State_HCO_TS(:,:)
    real(r8), pointer, public        :: State_HCO_U10M(:,:)
    real(r8), pointer, public        :: State_HCO_V10M(:,:)
    real(r8), pointer, public        :: State_HCO_ALBD(:,:)
    real(r8), pointer, public        :: State_HCO_WLI(:,:)
    real(r8), pointer, public        :: State_HCO_F_OF_PBL(:,:,:)


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
        if(masterproc) then
            write(iulog,*) "> HCOI_Allocate_All entering"
        endif

        ! PBL height [m]
        ! Comes from pbuf
        allocate(State_CAM_pblh(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_PBLH(my_IM, my_JM), stat=RC)
        ASSERT_(RC==0)

        ! Surface pressure (wet) [Pa]
        allocate(State_CAM_ps(my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_PSFC(my_IM, my_JM), stat=RC)
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
        allocate(State_CAM_LWI  (my_CE), stat=RC)
        ASSERT_(RC==0)

        allocate(State_HCO_TS(my_IM, my_JM), stat=RC)
        allocate(State_HCO_U10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_V10M(my_IM, my_JM), stat=RC)
        allocate(State_HCO_ALBD(my_IM, my_JM), stat=RC)
        allocate(State_HCO_WLI(my_IM, my_JM), stat=RC)
        allocate(State_HCO_F_OF_PBL(my_IM, my_JM, LM), stat=RC)

        ! Clear values
        State_HCO_U10M(:,:) = 0.0_r8
        State_HCO_V10M(:,:) = 0.0_r8
        State_HCO_ALBD(:,:) = 0.0_r8
        State_HCO_WLI(:,:) = 0.0_r8
        State_HCO_F_OF_PBL(:,:,:) = 0.0_r8

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
! !REVISION HISTORY:
!  16 Dec 2020 - H.P. Lin    - Initial version
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
        ! pbuf indices:
        integer                      :: index_pblh

        ! Assume success
        RC = ESMF_SUCCESS

        if(masterproc) then
            write(iulog,*) "> CAM_GetBefore_HCOI entering"
        endif

        ! Regrid necessary physics quantities from the CAM grid to the HEMCO grid
        ! Phase 0: Prepare necessary pbuf indices to retrieve data from the buffer...

        ! TODO: This prep phase might possibly be better moved into initialization.
        ! Keeping it here for now but later can be optimized (hplin, 3/31/20)
        index_pblh = pbuf_get_index('pblh')


        ! Phase 1: Store the fields in hemco_interface (copy)
        I = 0
        do lchnk = begchunk, endchunk    ! loop over all chunks in the physics grid
            ncol = get_ncols_p(lchnk)    ! columns per chunk
            pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk) ! get pbuf for this chunk...

            ! Gather data from pbuf before we proceed
            ! 2-D Fields: Write directly to pointer (does not need alloc).
            call pbuf_get_field(pbuf_chnk, index_pblh, State_CAM_pblh)

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
                enddo

                !----------------------------------------
                ! 2-D Fields
                !----------------------------------------
                ! Sea level pressure [hPa]
                State_CAM_ps(I) = phys_state(lchnk)%ps(J)

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

                ! Land water index (0 && some albedo param = ocean, 1 is land, 2 is ice)
                if(ExtState%WLI%DoUse) then
                    ! assume land
                    State_CAM_LWI(I) = 1
                    if(cam_in(lchnk)%iceFrac(J) .gt. (cam_in(lchnk)%landFrac(J) + cam_in(lchnk)%ocnFrac(J))) then
                        State_CAM_LWI(I) = 2
                    endif

                    if(cam_in(lchnk)%ocnFrac(J) .gt. (cam_in(lchnk)%landFrac(J) + cam_in(lchnk)%ocnFrac(J))) then
                        State_CAM_LWI(I) = 0
                    endif
                endif
            enddo
        enddo

        ! The below are stubs and are not regridded or processed for now due to lack of data
        
        if(masterproc) then
            write(iulog,*) "> CAM_GetBefore_HCOI finished"
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
    subroutine CAM_RegridSet_HCOI(HcoState, ExtState)
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
!
! !INPUT PARAMETERS:
!

        type(HCO_State), pointer           :: HcoState
        type(Ext_State), pointer           :: ExtState

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'CAM_GetBefore_HCOI'
        integer                      :: RC                   ! ESMF return code

        logical, save                :: FIRST = .true.

        integer                      :: I, J                 ! Loop index

        ! Assume success
        RC = ESMF_SUCCESS

        if(masterproc) then
            write(iulog,*) "> CAM_RegridSet_HCOI entering"
        endif

        ! Regrid necessary physics quantities from the CAM grid to the HEMCO grid
        ! Phase 1: Regrid
        call HCO_Grid_CAM2HCO_2D(State_CAM_ps,     State_HCO_PSFC  )
        call HCO_Grid_CAM2HCO_2D(State_CAM_pblh,   State_HCO_PBLH  )
        call HCO_Grid_CAM2HCO_3D(State_CAM_t,      State_HCO_TK    )

        ! Surface temperature [K] - use both for T2M and TSKIN for now according to CESM-GC,
        ! hplin 12/21/2020
        if(ExtState%T2M%DoUse .or. ExtState%TSKIN%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_TS, State_HCO_TS    )

            ! Point to appropriate location
            if(FIRST) then
                call ExtDat_Set(HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                                RC,       FIRST,          State_HCO_TS)

                call ExtDat_Set(HcoState, ExtState%T2M,  'T2M_FOR_EMIS', &
                                RC,       FIRST,          State_HCO_TS)
            endif
        endif

        ! 10M E/W and N/S wind speed [m/s] (fixme: use pver?)
        if(ExtState%U10M%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_U10M, State_HCO_U10M)
            call HCO_Grid_CAM2HCO_2D(State_CAM_V10M, State_HCO_V10M)

            if(FIRST) then
                call ExtDat_Set(HcoState, ExtState%U10M,  'U10M_FOR_EMIS', &
                                RC,       FIRST,          State_HCO_U10M)

                call ExtDat_Set(HcoState, ExtState%V10M,  'V10M_FOR_EMIS', &
                                RC,       FIRST,          State_HCO_V10M)
            endif
        endif

        ! Visible surface albedo [1]
        if(ExtState%ALBD%DoUse) then
            call HCO_Grid_CAM2HCO_2D(State_CAM_ALBD, State_HCO_ALBD)

            if(FIRST) then
                call ExtDat_Set(HcoState, ExtState%ALBD,  'ALBD_FOR_EMIS', &
                                RC,       FIRST,          State_HCO_ALBD)
            endif
        endif

        ! Land water index (0 && some albedo param = ocean, 1 is land, 2 is ice)
        if(ExtState%WLI%DoUse) then
            ! warning: regrid integer values have unexpected consequences
            call HCO_Grid_CAM2HCO_2D(State_CAM_LWI, State_HCO_WLI)

            ! remap fractions to a rounded number ... hack for regridding integer values
            do J = 1, my_JM
                do I = 1, my_IM
                    if(State_HCO_WLI(I,J) < 0.5) then
                        State_HCO_WLI(I,J) = 0
                    else
                        if(State_HCO_WLI(I,J) < 1.5) then
                            State_HCO_WLI(I,J) = 1
                        else
                            State_HCO_WLI(I,J) = 2
                        endif
                    endif
                enddo
            enddo

            if(FIRST) then
                call ExtDat_Set(HcoState, ExtState%WLI,  'WLI_FOR_EMIS', &
                                RC,       FIRST,         State_HCO_WLI)
            endif
        endif

        ! FIXME STUB: Just set zeros for F_OF_PBL hplin 12/21/20
        if(ExtState%FRAC_OF_PBL%DoUse) then
            if(FIRST) then
                call ExtDat_Set(HcoState, ExtState%FRAC_OF_PBL,  'FRAC_OF_PBL_FOR_EMIS', &
                                RC,       FIRST,                 State_HCO_F_OF_PBL)
            endif
        endif

        if(FIRST) then
            FIRST = .false.
        endif
        
        if(masterproc) then
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
