#define VERIFY_(A) if(.not.HCO_ESMF_VRFY(A,subname,__LINE__))call exit(-1)
#define ASSERT_(A) if(.not.HCO_ESMF_ASRT(A,subname,__LINE__))call exit(-1)
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_cam_exports
!
! !DESCRIPTION: Module HCO\_CAM\_EXPORTS manages interfaces to CAM for exporting
!  fields calculated by the HEMCO emissions component. It contains use of both
!  the CAM history output (for diagnostic output of fluxes calculated and debug)
!  and the physics buffer (for internal passing to chemistry packages)
!\\
!\\
! !INTERFACE:
!
module hco_cam_exports
!
! !USES:
!
    ! ESMF function wrappers
    use hco_esmf_wrappers
    use ESMF,                     only: ESMF_SUCCESS

    use shr_kind_mod,             only: r8 => shr_kind_r8
    use physics_buffer,           only: dtype_r8

    ! Controls
    use cam_abortutils,           only: endrun      ! fatal terminator
    use cam_logfile,              only: iulog       ! output log handle

    ! HEMCO grid information (for registering with CAM)
    use hco_esmf_grid,            only: my_IS, my_IE, my_JS, my_JE, my_IM, my_JM
    use hco_esmf_grid,            only: my_ID, nPET ! mytid, ntask
    use hco_esmf_grid,            only: my_CE       ! # of CAM ncols on this task (total sum of ncols_p)
    use hco_esmf_grid,            only: IM, JM, LM, XMid, YMid ! nlat/lon/lev, glat/lon

    ! History output
    use cam_history,              only: addfld      ! Add field for history output
    use cam_history,              only: outfld      ! output field to history

    ! Physics buffer
    use physics_buffer,           only: physics_buffer_desc

    implicit none
    private
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public    :: HCO_Exports_Init
    public    :: HCO_Export_History_HCO3D
    public    :: HCO_Export_History_CAM2D
    public    :: HCO_Export_History_CAM3D

    public    :: HCO_Export_Pbuf_QueryField
    public    :: HCO_Export_Pbuf_AddField
    public    :: HCO_Export_Pbuf_CAM3D
    public    :: HCO_Export_Pbuf_CAM2D
!
! !REMARKS:
!  This module should NOT be aware of particular chemical constituents. It should
!  only do export functions and the field names must be passed in from a higher-
!  -level interface, for cleaniness.
!
!  The workflow is usually to register the grid here and hemco_interface will
!  initialize the fields (into both history and pbuf as necessary)
!  The fields are written into by hemco_interface through the run gridcomp,
!  usually after regridding back to the physics chunk so it fits in pbuf format.
!
!  History output doesn't need regrid as we register the HCO grid with cam_history
!  at initialization here. (Will it work? We will see) (No it does not, hplin 4/10/20)
!
!  For Pbuf-based output (to pass to other components):
!  - pbuf_idx_map() is a fixed-size array to store species-ID-to-pbuf-ID mapping.
!    However, not all data HEMCO reads in are species. In this case, the idx_map will not
!    be used (pass hcoID=-1), and we will attempt to look up the field through pbuf_get_field
!    as usual. This may be slow, but always reliable. pbuf_idx_map is a convenience tool.
!  - All pbuf exports are on the CAM grid. Regrid data before passing in.
!
! !REVISION HISTORY:
!  25 Feb 2020 - H.P. Lin    - Initial version
!  10 Apr 2020 - H.P. Lin    - Added pbuf functionality
!  25 Feb 2021 - H.P. Lin    - Added pbuf 2-D export and query functionality
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PUBLIC TYPES:
!
    type(physics_buffer_desc), &
               public, pointer :: hco_pbuf2d(:,:) ! Pointer to the pbuf
!
! !PRIVATE TYPES:
!
    integer                    :: pbuf_idx_map(1024)      ! Mapping of species IDX to pbuf ID
contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Exports_Init
!
! !DESCRIPTION: Initializes the exports module (after grid initialization)
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Exports_Init()
!
! !USES:
!
        use spmd_utils,               only: masterproc

        ! Register history output grid
        use cam_grid_support,         only: cam_grid_register
        use cam_grid_support,         only: horiz_coord_create
        use cam_grid_support,         only: horiz_coord_t, iMap
!
! !INPUT PARAMETERS:
!
        
!
! !REMARKS:
!  Only registers HCO grid with CAM history for now
!
! !REVISION HISTORY:
!  25 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Exports_Init'

        integer,          parameter  :: hco_decomp = 233      ! Unique within CAM
        type(horiz_coord_t), pointer :: lon_coord => null()
        type(horiz_coord_t), pointer :: lat_coord => null()
        integer(iMap),       pointer :: grid_map(:,:) => null()
        integer(iMap),       pointer :: coord_map(:) => null()

        integer                      :: I, J, ind

        !-----------------------------------------------------------------------
        ! Create "hco_grid" HEMCO Grid for CAM HISTORY export
        ! HISTORY is used for diagnostic purposes
        !-----------------------------------------------------------------------
        allocate(grid_map(4, ((my_IE - my_IS + 1) * (my_JE - my_JS + 1))))
        ind = 0
        do J = my_JS, my_JE
            do I = my_IS, my_IE
                ind = ind + 1
                grid_map(1, ind) = I
                grid_map(2, ind) = J
                grid_map(3, ind) = I
                grid_map(4, ind) = J
            enddo
        enddo

        ! FIXME: This part does not support curvilinear coords (assuming rectilinear here)
        ! (hplin, 3/2/2020)
        allocate(coord_map(my_JE - my_JS + 1))
        coord_map = (/(J, J = my_JS, my_JE)/)
        lat_coord => horiz_coord_create('YMid', '', my_JM, 'latitude',                       &
                                        'degrees_north', my_JS, my_JE, YMid(1,my_JS:my_JE),  &
                                        map = coord_map)
        deallocate(coord_map)
        nullify(coord_map)

        allocate(coord_map(my_IE - my_IS + 1))
        coord_map = (/(I, I = my_IS, my_IE)/)
        lon_coord => horiz_coord_create('XMid', '', my_IM, 'longitude',                      &
                                        'degrees_east', my_IS, my_IE, XMid(my_IS:my_IE,1),   &
                                        map = coord_map)
        deallocate(coord_map)
        nullify(coord_map)

        call cam_grid_register('hco_grid', hco_decomp, lat_coord, lon_coord,                 &
                               grid_map, unstruct = .false.)
        deallocate(grid_map)
        nullify(grid_map)

        if(masterproc) then
            write(iulog,*) ">> Registered HEMCO hco_grid in CAM for history exports"
        endif

        !-----------------------------------------------------------------------
        ! Set pbuf_idx_map to magic number indicating that this buffer
        ! is unassigned
        !-----------------------------------------------------------------------
        pbuf_idx_map(:) = -233

    end subroutine HCO_Exports_Init
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_History_HCO3D
!
! !DESCRIPTION: Writes to CAM history a 3-D field in the HEMCO array. This uses
!  the HEMCO lat-lon grid.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_History_HCO3D(fldname, array)
!
! !USES:
!
        use cam_history,              only: hist_fld_active, outfld
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        real(r8),         intent(in) :: array(my_IS:my_IE,my_JS:my_JE,1:LM)
!
! !REMARKS:
!  Remember fields need to be declared via addfld in CAM before history export.
!
!  For a native (physics mesh) variant, use HCO_Export_History_CAM3D.
!
!  Based off the convoluted savefld_waccm in the ionos WACCMx interface. Notably,
!  (1) outfld accepts arguments in (fname, field, idim, c, avg_subcol_field) order,
!      where field is a 2-D array (idim,*) containing field values, and c is a
!      very mythical index.
!  (2) If you are outputting to lat-lon, c is the LATITUDE index of your output, so "j",
!      and outfld needs to be called in loops over the lat index.
!  (3) ... this is to accommodate that in the physics mesh, you have only (k,i) idxes,
!      which means that the data is passed in (fname, field, pcols, lchnk), field(i, k)
!      where pcols is the number of columns in the mesh, and lchnk is a loop index over
!      begchunk, endchunk (ppgrid). i is 1, ncol from ncol = get_ncols_p(lchnk).
!
!  At the time this code was written (3/3/2020) I absolutely understand none of the way
!  the physics mesh data is written, hence the rant above.
!
!  OK now I get it. See the note below in the CAM3D variant.
!
!  The below code has nothing to do with the rant above,
!  as HCO_Export_History_HCO3D is operating on the HEMCO lat-lon grid.
!
! !REVISION HISTORY:
!  25 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_History_HCO3D'
        integer                      :: I, J, K
        real(r8)                     :: tmpfld_ik(my_IS:my_IE, 1:LM)     ! lon-lev by lat

        if(.not. hist_fld_active(fldname)) then
            ! This routine is ALWAYS called but may fail silently if this history field
            ! is not to be outputted.
            return
        endif

        ! Not the most efficient; can probably use slicing. This will do for now (hplin, 3/3/20)
        do J = my_JS, my_JE
            do I = my_IS, my_IE
                do K = 1, LM
                    tmpfld_ik(I, K) = array(I, J, K)
                enddo
            enddo
            call outfld(fldname, tmpfld_ik, my_IE - my_IS + 1, J)         ! By lat convert to lon glob idx
        enddo
    end subroutine HCO_Export_History_HCO3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_History_CAM3D
!
! !DESCRIPTION: Writes to CAM history a 3-D field in the CAM array format. This
!  uses the CAM physics mesh (physgrid).
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_History_CAM3D(fldname, array)
!
! !USES:
!
        use cam_history,              only: hist_fld_active, outfld
        use ppgrid,                   only: pcols, pver
        use phys_grid,                only: begchunk, endchunk, get_ncols_p

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        real(r8),         intent(in) :: array(1:LM, 1:my_CE)
!
! !REMARKS:
!  Remember fields need to be declared via addfld in CAM before history export.
!  See rant above.
!
!  my_CE is the sum of get_ncols_p(lchnk) over begchunk, endchunk, called blksize
!  in the ionos code. It is the TOTAL number of columns on this PET.
!
!  The columns on this PET are divided into "chunks", lchnk = begchunk, endchunk.
!  The chunks each contain (up to) pcols each, specific number is from get_ncols_p.
!
!  This means that while the physics array is sized (1:LM, 1:my_CE) = (pver, blksize)
!  when they are written back, you have to account for putting them back into chunks
!  and writing using format outfld(..., array(1:pcols, 1:LM), pcols, lchnk)
!  where the data is sized 1:pcols, filled to 1:get_ncols_p(lchnk) and rest zeroed,
!  and the data is ordered in (i, k) called with pcols, lchnk as dim'ls.
!
!  
!
! !REVISION HISTORY:
!  03 Mar 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_History_CAM3D'
        integer                      :: lchnk, ncol
        integer                      :: I, K, J

        real(r8)                     :: tmpfld_ik(pcols, pver)   ! Temporary array for per-column data (i, k)

        if(.not. hist_fld_active(trim(fldname))) then
            ! This routine is ALWAYS called but may fail silently if this history field
            ! is not to be outputted.
            ! if(masterproc) write(iulog,*) "HCO_Export_History_CAM3D: Cannot export", trim(fldname), " (not active)"
            return
        endif

        ! For all chunks on this PET
        J = 0
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            ! For all columns in each chunk, organize the data
            do I = 1, ncol
                J = J + 1   ! Advance one column in the physics mesh array
                do K = 1, pver
                    tmpfld_ik(I, K) = array(K, J)
                enddo
            enddo

            ! Write to outfld chunk by chunk
            ! write(6,*) "hco_cam_exports before writing ncol, lchnk, ", ncol, lchnk
            call outfld(trim(fldname), tmpfld_ik(:ncol, :), ncol, lchnk)
            ! write(6,*) "hco_cam_exports writing ncol, lchnk, ", ncol, lchnk
        enddo

    end subroutine HCO_Export_History_CAM3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_History_CAM2D
!
! !DESCRIPTION: Writes to CAM history a 2-D field in the CAM array format. This
!  uses the CAM physics mesh (physgrid).
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_History_CAM2D(fldname, array)
!
! !USES:
!
        use cam_history,              only: hist_fld_active, outfld
        use ppgrid,                   only: pcols
        use phys_grid,                only: begchunk, endchunk, get_ncols_p

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        real(r8),         intent(in) :: array(1:my_CE)
!
! !REMARKS:
!  Remember fields need to be declared via addfld in CAM before history export.
!  Refer to the 3-D field export subroutine for further documentation.
!
! !REVISION HISTORY:
!  26 Feb 2021 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_History_CAM2D'
        integer                      :: lchnk, ncol
        integer                      :: I, J

        real(r8)                     :: tmpfld_i(pcols)   ! Temporary array for per-column data (i, k)

        if(.not. hist_fld_active(trim(fldname))) then
            ! This routine is ALWAYS called but may fail silently if this history field
            ! is not to be outputted.
            ! if(masterproc) write(iulog,*) "HCO_Export_History_CAM2D: Cannot export", trim(fldname), " (not active)"
            return
        endif

        ! For all chunks on this PET
        J = 0
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            ! For all columns in each chunk, organize the data
            do I = 1, ncol
                J = J + 1   ! Advance one column in the physics mesh array
                tmpfld_i(I) = array(J)
            enddo

            ! Write to outfld chunk by chunk
            ! write(6,*) "hco_cam_exports before writing ncol, lchnk, ", ncol, lchnk
            call outfld(trim(fldname), tmpfld_i(:ncol), ncol, lchnk)
            ! write(6,*) "hco_cam_exports writing ncol, lchnk, ", ncol, lchnk
        enddo

    end subroutine HCO_Export_History_CAM2D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_Pbuf_QueryField
!
! !DESCRIPTION: Query the physics buffer for whether a field exists
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_Pbuf_QueryField(fldname, dims, result, hcoID)
!
! !USES:
!
        use ppgrid,                   only: pcols, pver
        use physics_buffer,           only: pbuf_get_index

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        integer, intent(in)          :: dims                  ! 2 or 3 dimensions data?
        integer, intent(in), optional:: hcoID                 ! Species ID

        logical, intent(out)         :: result                ! Is field present?
!
! !REMARKS:
!  Persistence is 'physpkg' for now, which means the fields are alloc/dealloc at
!  the beginning/end of each physics time step.
!
!  It may be necessary to allocate some fields for 'global' to pass data from
!  chemistry back
!
! !REVISION HISTORY:
!  25 Feb 2021 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_Pbuf_QueryField'
        character(len=255)           :: fldname_ns

        integer                      :: spcID, tmpIdx
        integer                      :: RC

        ! Create field name and verify IDs
        fldname_ns = 'HCO_' // trim(fldname)
        if(present(hcoID)) then
            spcID = hcoID
        else
            spcID = -1
        endif

        ! Not found by default
        result = .false.

        ! Verify if slot occupied
        if(spcID /= -1) then
            if(pbuf_idx_map(spcID) /= -233) then
                result = .true.
            endif
        else
            tmpIdx = pbuf_get_index(fldname_ns, RC)
            if(tmpIdx >= 0) then
                result = .true.
            endif
        endif

    end subroutine HCO_Export_Pbuf_QueryField
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_Pbuf_AddField
!
! !DESCRIPTION: Adds to the physics buffer a HEMCO field for export to the chem-
!  istry. We wrap this here so we can control the persistence and add a prefix.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_Pbuf_AddField(fldname, dims, hcoID)
!
! !USES:
!
        use ppgrid,                   only: pcols, pver
        use physics_buffer,           only: pbuf_add_field
        use physics_buffer,           only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        integer, intent(in)          :: dims                  ! 2 or 3 dimensions data?
        integer, intent(in), optional:: hcoID                 ! Species ID
!
! !REMARKS:
!  Persistence is 'physpkg' for now, which means the fields are alloc/dealloc at
!  the beginning/end of each physics time step.
!
!  We add a prefix hco_ to all the fields, to prevent namespace clashing.
!
!  All fields are exported in CAM grid format, either 2D or 3D.
!
!  Now playing: Lemon
!  "Ima demo anata wa watashi no hikari" / Even now you are still my light
!
! !REVISION HISTORY:
!  10 Apr 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_Pbuf_AddField'
        character(len=255)           :: fldname_ns

        integer                      :: spcID, tmpIdx

        ! Create field name and verify IDs
        fldname_ns = 'HCO_' // trim(fldname)
        if(present(hcoID)) then
            spcID = hcoID
        else
            spcID = -1
        endif

        ! Verify if slot free
        if(spcID /= -1) then
            ASSERT_(pbuf_idx_map(spcID) == -233)
        endif

        ! Add to pbuf field
        if(dims == 2) then
            call pbuf_add_field(trim(fldname_ns), 'physpkg', dtype_r8,  &
                                (/pcols/), tmpIdx)
        elseif(dims == 3) then
            call pbuf_add_field(trim(fldname_ns), 'physpkg', dtype_r8,  &
                                (/pcols,pver/), tmpIdx)
        else
            ASSERT_(.false.)
        endif

        ! Save field to mapping
        if(spcID /= -1) then
            pbuf_idx_map(spcID) = tmpIdx
        endif

        ! Log
        if(masterproc) write(iulog,*) "Added field " // trim(fldname_ns) // " to physpkg pbuf, idx", tmpidx, " spcID", spcID

    end subroutine HCO_Export_Pbuf_AddField
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_Pbuf_CAM3D
!
! !DESCRIPTION: Adds to the physics buffer a HEMCO field for export to the chem-
!  istry. We wrap this here so we can control the persistence and add a prefix.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_Pbuf_CAM3D(fldname, hcoID, array)
!
! !USES:
!
        use ppgrid,                   only: pcols, pver
        use phys_grid,                only: begchunk, endchunk, get_ncols_p
        use physics_buffer,           only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        integer, intent(in), optional:: hcoID                 ! Species ID
        real(r8),         intent(in) :: array(1:LM, 1:my_CE)
!
! !REMARKS:
!  pcols: maximum number of columns in a chunk
!  each pbuf is independent in each chunk, so write chunk-by-chunk, column-by-column
!
!  The pointer copy of pbuf2d is right there at the top of the module.
!  It is updated at every time step by hemco_interface::HCOI_Chunk_Run, for lack of a
!  better method to propagate it to inside the gridded component.
!
! !REVISION HISTORY:
!  10 Apr 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_Pbuf_CAM3D'
        character(len=255)           :: fldname_ns

        integer                      :: spcID, tmpIdx

        integer                      :: lchnk, ncol
        integer                      :: I, K, J
        integer                      :: RC
        type(physics_buffer_desc), pointer :: pbuf_chnk(:)       ! slice of pbuf in chnk
        real(r8), pointer            :: pbuf_ik(:,:)     ! Pointer to pbuf data (/pcols,pver/)

        ! Create field name and verify IDs
        fldname_ns = 'HCO_' // trim(fldname)
        if(present(hcoID)) then
            spcID = hcoID
        else
            spcID = -1
        endif

        ! Verify if slot occupied
        if(spcID /= -1) then
            if(pbuf_idx_map(spcID) == -233) then
                if(masterproc) write(iulog,*) "HCO_Export_Pbuf_CAM3D: Field not found", spcID, fldname_ns
                return 
            endif
            tmpIdx = pbuf_idx_map(spcID)
        else
            tmpIdx = pbuf_get_index(fldname_ns, RC)
            if(tmpIdx < 0) then
                if(masterproc) write(iulog,*) "HCO_Export_Pbuf_CAM3D: Field not found", spcID, fldname_ns
                return
            endif
        endif

        ! For all chunks on this PET
        J = 0
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, lchnk)       ! slice of pbuf for this chnk

            ! Get this pointer from pbuf
            call pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_ik)

            if(.not. associated(pbuf_ik)) then
                ! Uh-oh, potentially fatal. Throw a tantrum
                write(6,*) "HCO_Export_Pbuf_CAM3D: FATAL at lchnk", lchnk, " unassoc pbuf_get_field" // trim(fldname_ns) // " idx:", tmpIdx
                ASSERT_(.false.)
            endif

            ! For all columns in each chunk, organize the data
            do I = 1, ncol
                J = J + 1   ! Advance one column in the physics mesh array
                do K = 1, pver
                    pbuf_ik(I, K) = array(K, J)
                enddo
            enddo

            ! Nullify the field to prevent dangling stuff
            nullify(pbuf_ik)
        enddo
    end subroutine HCO_Export_Pbuf_CAM3D
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_Pbuf_CAM2D
!
! !DESCRIPTION: Adds to the physics buffer a HEMCO field for export to the chem-
!  istry. We wrap this here so we can control the persistence and add a prefix.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_Pbuf_CAM2D(fldname, hcoID, array)
!
! !USES:
!
        use ppgrid,                   only: pcols, pver
        use phys_grid,                only: begchunk, endchunk, get_ncols_p
        use physics_buffer,           only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

        use spmd_utils,               only: iam, masterproc
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        integer, intent(in), optional:: hcoID                 ! Species ID
        real(r8),         intent(in) :: array(1:my_CE)
!
! !REMARKS:
!  Read remarks at the CAM3D subroutine. This is the 2-D version (1-D in CAM speak)
!
! !REVISION HISTORY:
!  25 Feb 2021 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCO_Export_Pbuf_CAM2D'
        character(len=255)           :: fldname_ns

        integer                      :: spcID, tmpIdx

        integer                      :: lchnk, ncol
        integer                      :: I, J
        integer                      :: RC
        type(physics_buffer_desc), pointer :: pbuf_chnk(:)       ! slice of pbuf in chnk
        real(r8), pointer            :: pbuf_i(:)     ! Pointer to pbuf data (/pcols/)

        ! Create field name and verify IDs
        fldname_ns = 'HCO_' // trim(fldname)
        if(present(hcoID)) then
            spcID = hcoID
        else
            spcID = -1
        endif

        ! Verify if slot occupied
        if(spcID /= -1) then
            if(pbuf_idx_map(spcID) == -233) then
                if(masterproc) write(iulog,*) "HCO_Export_Pbuf_CAM2D: Field not found", spcID, fldname_ns
                return 
            endif
            tmpIdx = pbuf_idx_map(spcID)
        else
            tmpIdx = pbuf_get_index(fldname_ns, RC)
            if(tmpIdx < 0) then
                if(masterproc) write(iulog,*) "HCO_Export_Pbuf_CAM2D: Field not found", spcID, fldname_ns
                return
            endif
        endif

        ! For all chunks on this PET
        J = 0
        do lchnk = begchunk, endchunk
            ncol = get_ncols_p(lchnk)
            pbuf_chnk => pbuf_get_chunk(hco_pbuf2d, lchnk)       ! slice of pbuf for this chnk

            ! Get this pointer from pbuf
            call pbuf_get_field(pbuf_chnk, tmpIdx, pbuf_i)

            if(.not. associated(pbuf_i)) then
                ! Uh-oh, potentially fatal. Throw a tantrum
                write(6,*) "HCO_Export_Pbuf_CAM2D: FATAL at lchnk", lchnk, " unassoc pbuf_get_field" // trim(fldname_ns) // " idx:", tmpIdx
                ASSERT_(.false.)
            endif

            ! For all columns in each chunk, organize the data
            do I = 1, ncol
                J = J + 1   ! Advance one column in the physics mesh array
                pbuf_i(I) = array(J)
            enddo

            ! Nullify the field to prevent dangling stuff
            nullify(pbuf_i)
        enddo
    end subroutine HCO_Export_Pbuf_CAM2D
!EOC
end module hco_cam_exports