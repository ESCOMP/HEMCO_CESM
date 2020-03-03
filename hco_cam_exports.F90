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
    use shr_kind_mod,             only: r8 => shr_kind_r8


    ! Controls
    use cam_abortutils,           only: endrun      ! fatal terminator
    use cam_logfile,              only: iulog       ! output log handle

    ! HEMCO grid information (for registering with CAM)
    use hco_esmf_grid,            only: my_IS, my_IE, my_JS, my_JE
    use hco_esmf_grid,            only: my_ID, nPET ! mytid, ntask
    use hco_esmf_grid,            only: IM, JM, LM, XMid, YMid ! nlat/lon/lev, glat/lon

    ! History output
    use cam_history,              only: addfld      ! Add field for history output
    use cam_history,              only: outfld      ! output field to history

    ! Physics buffer
    use physics_buffer,           only: pbuf_get_chunk, pbuf_get_field, pbuf_get_index

    implicit none
    private
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public    :: HCO_Exports_Init
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
!  at initialization here. (Will it work? We will see)
!
! !REVISION HISTORY:
!  25 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
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
        lat_coord => horiz_coord_create('YMid', '', JM, 'latitude',                          &
                                        'degrees_north', my_JS, my_JE, YMid(1,my_JS:my_JE),  &
                                        map = coord_map)
        deallocate(coord_map)
        nullify(coord_map)

        allocate(coord_map(my_IE - my_IS + 1))
        coord_map = (/(I, I = my_IS, my_IE)/)
        lon_coord => horiz_coord_create('XMid', '', IM, 'longitude',                         &
                                        'degrees_east', my_IS, my_IE, XMid(my_IS:my_IE,1),   &
                                        map = coord_map)
        deallocate(coord_map)
        nullify(coord_map)

        call cam_grid_register('hco_grid', hco_decomp, lat_coord, lon_coord,               &
                               grid_map, unstruct=.false.)
        deallocate(grid_map)
        nullify(grid_map)

        if(masterproc) then
            write(iulog,*) ">> Registered HEMCO hco_grid in CAM for outfld exports"
        endif

    end subroutine HCO_Exports_Init
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Export_History_HCO3D
!
! !DESCRIPTION: Writes to CAM history a 3-D field in the HEMCO array.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_Export_History_HCO3D(fldname, array)
!
! !USES:
!
        
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: fldname               ! Field name
        real(r8),         intent(in) :: array(my_IS:my_IE,my_JS:my_JE,1:LM)
!
! !REMARKS:
!  Remember fields need to be declared via addfld in CAM before history export.
!  Probably check via hist_fld_active?
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

    end subroutine HCO_Export_History_HCO3D
!EOC
end module hco_cam_exports