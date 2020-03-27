#define VERIFY_(A) if(.not.HCO_ESMF_VRFY(A,subname,__LINE__))call exit(-1)
#define ASSERT_(A) if(.not.HCO_ESMF_ASRT(A,subname,__LINE__))call exit(-1)
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemco_interface
!
! !DESCRIPTION: Module HEMCO\_INTERFACE is the HEMCO-CESM interface module.
!               CESM operates on chunks thus the interface is called HCOI\_Chunk.
!               Internally it uses a gridded component to interact with the CAM
!               physics grid; these functions are internal and called HCO\_GC...
!\\
!\\
! !INTERFACE:
!
module hemco_interface
!
! !USES:
!
    ! ESMF function wrappers
    use hco_esmf_wrappers

    ! HEMCO ESMF Grid helpers
    use hco_esmf_grid

    ! CAM export helpers
    use hco_cam_exports

    ! Controls
    use cam_abortutils,           only: endrun      ! fatal terminator
    use cam_logfile,              only: iulog       ! output log handle

    ! Species information 
    use chem_mods,                only: gas_pcnst   ! # of species
    use chem_mods,                only: adv_mass    ! advected mass
    use mo_tracname,              only: solsym      ! species names
    use mo_chem_utls,             only: get_spc_ndx ! IND_

    ! Grid
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs

    ! ESMF types
    use ESMF,                     only: ESMF_State, ESMF_Clock, ESMF_GridComp
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS

    use shr_kind_mod,             only: r8 => shr_kind_r8

    implicit none
    private
!
! !PRIVATE MEMBER FUNCTIONS:
!
    private :: HCO_GC_Init
    private :: HCO_GC_SetServices
    private :: HCO_GC_Run
    private :: HCO_GC_Final
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public  :: hemco_readnl
    public  :: HCOI_Chunk_Init
    public  :: HCOI_Chunk_Run
    public  :: HCOI_Chunk_Final
!
! !REVISION HISTORY:
!  29 Jan 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    type(ESMF_GridComp)  :: HCO_GridComp        ! HEMCO GridComp
    type(ESMF_State)     :: HCO_GridCompState   ! HEMCO GridComp Import/Export State
contains
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hemco_readnl
!
! !DESCRIPTION: Reads the namelist from cam/src/control/runtime_opts.
!\\
!\\
! !INTERFACE:
!
    subroutine hemco_readnl( nlfile )
!
! !USES:
!
        use namelist_utils, only: find_group_name
        use units,          only: getunit, freeunit
        use spmd_utils,     only: mpi_real8, mpi_logical, mpi_integer, mpi_character
        use spmd_utils,     only: masterproc, mpicom, masterprocid
!
! !INPUT PARAMETERS:
!
        character(len=*), intent(in) :: nlfile ! namelist file
!
! !REVISION HISTORY:
!  31 Jan 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        integer :: unitn, ierr
        character(len=*), parameter  :: subname = 'hemco_readnl'
        character(len=256)           :: hemco_config_file = 'HEMCO_Config.rc'

        namelist /hemco_nl/ hemco_config_file

        ! Read namelist on master proc
        ! ...
        if(masterproc) then
            write(iulog,*) "This is hemco_readnl - Reading HEMCO Namelist in CESM"
            unitn = getunit()
            open(unitn, file=trim(nlfile), status='old')
            call find_group_name(unitn, 'hemco_nl', status=ierr)
            if(ierr == 0) then
                read(unitn, hemco_nl, iostat=ierr)
                if(ierr /= 0) then
                    call endrun(subname // ':: ERROR reading namelist')
                endif
            endif
            close(unitn)
            call freeunit(unitn)

            write(iulog,*) "hemco_readnl: hemco config file = ", hemco_config_file
        endif

        ! Read HEMCO Configuration file then broadcast
        ! ...

        
    end subroutine hemco_readnl
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Chunk_Init
!
! !DESCRIPTION: HCOI\_Chunk\_Init is the initialization method for the CAM
!  interface to HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Chunk_Init()
!
! !USES:
!
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid

        use mpi,            only: MPI_INTEGER
        use ESMF,           only: ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
        use ESMF,           only: ESMF_GridCompCreate, ESMF_GridCompInitialize
        use ESMF,           only: ESMF_GridCompSetServices
        use ESMF,           only: ESMF_StateCreate

        use ESMF,           only: ESMF_Initialize, ESMF_LOGKIND_MULTI

        ! CAM instance information
        use cam_instance,   only: inst_index, inst_name

        ! CAM history output (to be moved somewhere later)
        use cam_history,    only: addfld, add_default
!
! !REVISION HISTORY:
!  06 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Init'
        integer                      :: RC                   ! ESMF return code
        integer                      :: HMRC                 ! HEMCO return code

        ! Gridded component properties.
        ! Note that while edyn_grid_comp initializes the GridComp directly using
        ! cam_instance's inst_name, I think this may cause namespace clashing.
        ! So I'll be prefixing this with hco_ just incase.
        character(len=32)            :: HCO_GC_InstName = ''

        ! MPI stuff
        integer                      :: localPET, PETcount
        integer                      :: iam, npes            ! My CPU, # of PETs

        ! PETs for each instance of the physics grid
        integer, allocatable         :: PETlist(:)

        type(ESMF_VM)                :: hco_esmf_vm

        if(masterproc) then
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            write(iulog,*) "HEMCO: Harmonized Emissions Component"
            write(iulog,*) "HEMCO_CAM interface version 0.1"
        endif

        !-----------------------------------------------------------------------
        ! Setup ESMF wrapper gridded component
        ! Adapted from edyn_grid_comp_init
        !-----------------------------------------------------------------------
        call ESMF_VMGetCurrent(hco_esmf_vm, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_VMGet(hco_esmf_vm, localPet=localPET, petCount=PETcount, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Get MPI communicator stuff
        call mpi_comm_size(mpicom, npes, RC)
        call mpi_comm_rank(mpicom, iam,  RC)

        ! Allocate and collect PETs for each instance of the physics grid
        allocate(PETlist(npes))
            ! sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror
        call mpi_allgather(localPET, 1, MPI_INTEGER, PETlist, 1, MPI_INTEGER, mpicom, RC)
        
        ! Create ESMF gridded component
        ! This gridded component's IRF routines are defined in hemco_interface::HCO_GC_SetServices
        ! See there for more information. (hplin, 2/6/20)
        HCO_GC_InstName = 'HCO_' // trim(inst_name)
        HCO_GridComp = ESMF_GridCompCreate(name=trim(HCO_GC_InstName), petList=PETlist, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Create a dummy import / export state.
        call ESMF_GridCompSetServices(HCO_GridComp, HCO_GC_SetServices, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        HCO_GridCompState = ESMF_StateCreate(name='HEMCO GridComp State', rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridCompInitialize(HCO_GridComp, importState=HCO_GridCompState, exportState=HCO_GridCompState, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> Initialized ESMF environment successfully! localPET, PETcount", localPET, PETcount
            write(iulog,*) "> iam, npes", iam, npes
            write(iulog,*) "> PETlist", PETlist
        endif

        !-----------------------------------------------------------------------
        ! Setup a lat-lon "HEMCO" intermediate grid
        !-----------------------------------------------------------------------

        ! TODO: For now, hardcode using 1x1 grid for testing with fv0.9x1.25.
        ! This provides a sufficient # of grid boxes acceptable for decomp
        ! on the default Cheyenne configuration. (hplin, 2/21/20)

        ! TODO: For now, # of PEs to use for HEMCO will be total # of PEs.
        ! These will all have to be specified in the HEMCO namelist later on.
        call HCO_Grid_Init (IM_in = 360, JM_in = 181, nPET_in = npes, RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> Initialized HEMCO Grid environment successfully!"
        endif

        !-----------------------------------------------------------------------
        ! Update HEMCO regrid descriptors for the first time.
        !-----------------------------------------------------------------------
        call HCO_Grid_UpdateRegrid(RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> First refresh of HEMCO Regrid descriptors"
        endif

        !-----------------------------------------------------------------------
        ! Initialize CAM export component
        !-----------------------------------------------------------------------
        call HCO_Exports_Init()

        if(masterproc) then
            write(iulog,*) "> Initialize HEMCO/CAM exports component"
        endif

        ! Test only hplin 3/3/20: add a dummy history field in CAM to test HEMCO
        ! grid is correctly reflected.
        ! Two fields are used: HCOg_TEST, which outputs in native HCO lat-lon
        ! and HCO_TEST, which outputs in the physics mesh.
        ! These two are used to test whether we've done our math correctly.
        !
        ! AvgFlag: (cam_history) A mean, B mean00z, I instant, X max, M min, S stddev
        !
        ! This call should eventually be reflected elsewhere?
        ! call addfld("HCOg_TEST", (/'lev'/), 'I', 'kg/m2/s',         &
        !             'HEMCO 3-D Emissions Species TEST on HCO Grid', &
        !             gridname="hco_grid")
        ! call add_default("HCOg_TEST", 2, 'I')     ! Make this field always ON

        call addfld("HCO_TEST", (/'lev'/), 'I', 'kg/m2/s',          &
                    'HEMCO 3-D Emissions Species TEST',             &
                    gridname="physgrid")
        call add_default("HCO_TEST", 2, 'I')      ! Make this field always ON

        !-----------------------------------------------------------------------
        ! Initialize HEMCO!
        !-----------------------------------------------------------------------
        

        if(masterproc) then

            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            ! End the splash screen
        endif

    end subroutine HCOI_Chunk_Init
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Chunk_Run
!
! !DESCRIPTION: HCOI\_Chunk\_Run is the run method for the CAM interface to HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Chunk_Run(cam_in, phys_state, pbuf2d, phase)
!
! !USES:
!
        ! Type descriptors
        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physics_buffer, only: physics_buffer_desc

        ! Output and mpi
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid, iam

        ! ESMF
        use ESMF,           only: ESMF_GridCompRun

!
! !INPUT PARAMETERS:
!
        type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
        type(physics_buffer_desc), pointer :: pbuf2d(:,:)
        integer, intent(in)                :: phase               ! 1, 2
!
! !REVISION HISTORY:
!  06 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Run'
        integer                      :: RC                        ! Return code

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Running HCOI_Chunk_Run phase", phase
        endif

        ! We only run the gridded components on Phase 2 per recommendations
        ! from Steve, but this can be easily extended.
        !
        ! In edyn_grid_comp, the ionosphere interface runs on run1 and run2
        ! and uses a global variable to control the gridded component's
        ! run stage. This is not needed for HEMCO right now but we can always
        ! implement this in the future.
        !
        ! HEMCO also has two stages, but the distinction is not necessary
        ! in the CAM interface. For more information, look at the actual
        ! run routine in the gridded component HCO_GC_Run.
        ! (hplin, 2/6/20)
        if(phase == 2) then
            ! Run the gridded component.
            write(6,*) "HEMCO_CAM inside HCOI_Chunk_Run to enter GridComp, iam", iam
            call ESMF_GridCompRun(HCO_GridComp, rc=RC)!importState=HCO_GridCompState, &
                                                !exportState=HCO_GridCompState, &
                                                !rc=RC)

            ASSERT_(RC==ESMF_SUCCESS)
        endif

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Leaving HCOI_Chunk_Run"
        endif
    end subroutine HCOI_Chunk_Run
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Chunk_Final
!
! !DESCRIPTION: HCOI\_Chunk\_Final cleans up the CAM interface to HEMCO.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Chunk_Final()
        ! Stub...
    end subroutine HCOI_Chunk_Final
    !-----------------------------------------------------------------------
    !              H E M C O   W R A P P E R   G R I D C O M P             !
    !-----------------------------------------------------------------------
    !  Below code includes internal routines used to wrap a ESMF gridded   !
    !  component around the HEMCO interface, so ESMF can handle all the    !
    !  interaction with the physics mesh.                                  !
    !                                                                      !
    !  This is largely based on edyn_grid_comp.F90 from ionosphere/waccmx  !
    !  Thanks to Steve Goldhaber for the example                           !
    !                                                                      !
    !  (hplin, 1/31/20)                                                    !
    !-----------------------------------------------------------------------
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GC_Run
!
! !DESCRIPTION: HCO\_GC\_Run is an internal method in the HEMCO gridded component
!  in CAM. It runs the main routines of HEMCO and is called by ESMF.
!  The routines inside the GridComp operate on the HEMCO grid and are responsible
!  to run regridding routines to return data into the physics mesh.
!  In short, code goes here.
!\\
!\\
! !INTERFACE:
!
    subroutine HCO_GC_Run(GC, IMPORT, EXPORT, Clock, RC)
!
! !USES:
!
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, iam

        ! Includes actual HEMCO run routines...

!
! !INPUT/OUTPUT PARAMETERS:
!
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC
!
! !REMARKS:
!  All the input/output parameters here are dummies.
!
! !REVISION HISTORY:
!  06 Feb 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*),       parameter     :: subname = 'HCO_GC_Run'

        integer                               :: I, J, K
        real(ESMF_KIND_R8)                    :: dummy(my_IS:my_IE, my_JS:my_JE, 1:LM)
        real(ESMF_KIND_R8)                    :: dummy_CAM(1:LM, 1:my_CE)

        write(6,*) "HEMCO_CAM: Inside GridComp", iam

        !-----------------------------------------------------------------------
        ! Update regridding file handles as necessary
        !-----------------------------------------------------------------------
        call HCO_Grid_UpdateRegrid(RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> Reload (if necessary) of HEMCO Regrid descriptors"
        endif

        !-----------------------------------------------------------------------
        ! Run HEMCO!
        !-----------------------------------------------------------------------
        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Inside GridComp: running HCO_GC_Run!"
        endif

        write(6,*) "HEMCO_CAM: before dummy array"

        !-----------------------------------------------------------------------
        ! Do some dummy stuff here while we don't have actual "HEMCO"
        !-----------------------------------------------------------------------
        ! This writes a striped grid to level 1 on lat-lon. Fun!
        dummy(:,:,:) = 0.0_r8
        dummy_CAM(:,:) = 0.0_r8
        do I = my_IS, my_IE
        do J = my_JS, my_JE
            if(mod(int(XMid(I,J)), 2) == 0 .and. mod(int(YMid(I,J)), 2) == 0) then
                dummy(I,J,1) = XMid(I,J)
            endif
        enddo
        enddo

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Successful creation of dummy array!"
        endif

        ! Write to history on hco_grid (hco_grid NOT working)
        ! call HCO_Export_History_HCO3D("HCOg_TEST", dummy)

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Exported to HCOg_TEST via outfld!"
        endif

        write(6,*) "HEMCO_CAM: before dummy->dummy_CAM"

        ! Regrid to CAM physics mesh!
        call HCO_Grid_HCO2CAM_3D(dummy, dummy_CAM)

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Successful dummy regrid to CAM!"
        endif

        write(6,*) "HEMCO_CAM: after dummy->dummy_CAM"

        ! Write to history on CAM mesh
        call HCO_Export_History_CAM3D("HCO_TEST", dummy_CAM)

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Export to HCO_TEST via outfld!"
        endif

        ! ... stub ...

        RC = ESMF_SUCCESS

    end subroutine HCO_GC_Run

    !---------------------------------------------------------------------
    ! Every HEMCO GridComp routine below likely will be just boilerplate
    ! that does not need extensive maintenance
    !---------------------------------------------------------------------

    ! Init routine for the Gridded Component
    ! Largely based off edyn_grid_comp::edyn_gcomp_init
    subroutine HCO_GC_Init(GC, IMPORT, EXPORT, Clock, RC)
        use spmd_utils,     only: masterproc
        use cam_logfile,    only: iulog

        ! Dummy arguments
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC

        ! Local variables
        character(len=*),       parameter     :: subname = 'HCO_GC_Init'

        ! Note hplin 2/17/20: It seems like the physics mesh is re-created in
        ! edyn_esmf through edyn_create_physmesh. It may be redundant to do
        ! the CAM_DistGrid and CAM_PhysMesh maneuvers here and do this in
        ! HCO_ESMF_Grid::HCO_Grid_ESMF_CreateCAM instead.

    end subroutine HCO_GC_Init

    ! Finalize Gridded Component
    subroutine HCO_GC_Final(GC, IMPORT, EXPORT, Clock, RC)
        ! use ESMF,         only: ESMF_MeshDestroy

        ! Dummy arguments
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC

        ! Local variables
        character(len=*),       parameter     :: subname = 'HCO_GC_Final'

        ! call ESMF_MeshDestroy(CAM_PhysMesh, rc=RC)
        ! ASSERT_(RC==ESMF_SUCCESS)
    end subroutine HCO_GC_Final

    subroutine HCO_GC_SetServices(GC, RC)
        use ESMF,         only: ESMF_GridCompSetEntryPoint
        use ESMF,         only: ESMF_METHOD_INITIALIZE
        use ESMF,         only: ESMF_METHOD_RUN
        use ESMF,         only: ESMF_METHOD_FINALIZE

        type(ESMF_GridComp)                   :: GC
        integer, intent(out)                  :: RC
        character(len=*),       parameter     :: subname = 'HCO_GC_SetServices'

        ! Set the IRF methods as follows:
        ! HEMCO Gridded Component dummy
        !   > Init:     hemco_interface::HCO_GC_Init
        !   > Run:      hemco_interface::HCO_GC_Run
        !   > Final:    hemco_interface::HCO_GC_Final
        !
        ! Note from Steve: " all the actual regridding will have to happen inside the gridded component's run method."

        call ESMF_GridCompSetEntryPoint(GC, ESMF_METHOD_INITIALIZE, userRoutine=HCO_GC_Init,  rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridCompSetEntryPoint(GC, ESMF_METHOD_RUN,        userRoutine=HCO_GC_Run,   rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_GridCompSetEntryPoint(GC, ESMF_METHOD_FINALIZE,   userRoutine=HCO_GC_Final, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)
    end subroutine HCO_GC_SetServices
!EOC
end module hemco_interface