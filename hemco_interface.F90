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
    use ESMF,                     only: ESMF_Mesh, ESMF_DistGrid
    use ESMF,                     only: ESMF_State, ESMF_Clock, ESMF_GridComp
    use ESMF,                     only: ESMF_Field, ESMF_RouteHandle
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS

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
    type(ESMF_Mesh)      :: CAM_PhysMesh        ! Copy of CAM physics mesh
    type(ESMF_DistGrid)  :: CAM_DistGrid        ! DE-local allocation descriptor DistGrid (2D)
    type(ESMF_GridComp)  :: HCO_GridComp        ! HEMCO GridComp
    type(ESMF_State)     :: HCO_GridCompState   ! HEMCO GridComp Import/Export State
    
    integer              :: col_start, col_end  ! idx of columns in this PET
    integer              :: col_total           ! # of columns in this PET
    integer              :: nlev = 0            ! # of levs in this PET
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
            write(iulog,*) "> Initialized ESMF environment successfully!"
        endif

        !-----------------------------------------------------------------------
        ! Setup a lat-lon "HEMCO" intermediate grid
        !-----------------------------------------------------------------------

        ! TODO: For now, hardcode using 2x2.5 grid for testing. (hplin, 2/12/20)
        ! TODO: For now, # of PEs to use for HEMCO will be total # of PEs.
        ! These will all have to be specified in the HEMCO namelist later on.
        call HCO_Grid_Init (IM_in = 144, JM_in = 91, nPET_in = npes, RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> Initialized HEMCO Grid environment successfully!"
        endif

        !-----------------------------------------------------------------------
        ! Initialize HEMCO!
        !-----------------------------------------------------------------------
        

        if(masterproc) then

            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            ! End the splash screen
        endif

    end subroutine HCOI_Chunk_Init

    ! Run
    subroutine HCOI_Chunk_Run(cam_in, phys_state, pbuf2d, phase)
        ! Type descriptors
        use camsrfexch,     only: cam_in_t
        use physics_types,  only: physics_state
        use physics_buffer, only: physics_buffer_desc

        ! Output and mpi
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc, mpicom, masterprocid

        ! ESMF
        use ESMF,           only: ESMF_GridCompRun

        ! Input
        type(cam_in_t),      intent(inout) :: cam_in(begchunk:endchunk)
        type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
        type(physics_buffer_desc), pointer :: pbuf2d(:,:)
        integer, intent(in)                :: phase               ! 1, 2

        ! Local variables
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
            call ESMF_GridCompRun(HCO_GridComp, importState=HCO_GridCompState, &
                                                exportState=HCO_GridCompState, &
                                                rc=RC)

            ASSERT_(RC==ESMF_SUCCESS)
        endif

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Leaving HCOI_Chunk_Run"
        endif
    end subroutine HCOI_Chunk_Run

    ! Final
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

    ! Run routine called from the ESMF GridComp. Contains the actual computation
    ! routines (secret sauce) operating on the HEMCO grid and the regridding
    ! routines to return it into the physics mesh.
    subroutine HCO_GC_Run(GC, IMPORT, EXPORT, Clock, RC)
        ! Utilities for printing out debug output
        use cam_logfile,    only: iulog
        use spmd_utils,     only: masterproc

        ! Includes actual HEMCO run routines...

        ! Dummy arguments
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC

        ! Local variables
        character(len=*),       parameter     :: subname = 'HCO_GC_Run'

        !-----------------------------------------------------------------------
        ! Update regridding file handles as necessary
        !-----------------------------------------------------------------------
        ! TODO

        !-----------------------------------------------------------------------
        ! Run HEMCO!
        !-----------------------------------------------------------------------
        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Inside GridComp: running HCO_GC_Run!"
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
        use ESMF,         only: ESMF_DistGridCreate, ESMF_MeshCreate
        use ESMF,         only: ESMF_FILEFORMAT_ESMFMESH
        use ESMF,         only: ESMF_MeshIsCreated, ESMF_MeshDestroy

        ! CAM instance and physical grid information
        use cam_instance, only: inst_name
        use phys_control, only: phys_getopts
        use phys_grid,    only: get_ncols_p, get_gcol_p

        use spmd_utils,     only: masterproc
        use cam_logfile,    only: iulog

        ! Dummy arguments
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC

        ! Local variables
        integer                               :: ncols
        integer                               :: chnk, col, dindex
        integer,                allocatable   :: decomp(:)
        character(len=256)                    :: grid_file
        character(len=*),       parameter     :: subname = 'HCO_GC_Init'

        ! Get the physics grid information
        call phys_getopts(physics_grid_out=grid_file)

        if(masterproc) then
            write(iulog,*) "physics_grid_out=", grid_file
        endif

        ! Compute local decomposition (global variable in-module)
        col_total = 0 ! Sum of columns on this PET in all chunks
        do chnk = begchunk, endchunk
            col_total = col_total + get_ncols_p(chnk)
        enddo
        allocate(decomp(col_total))
        dindex = 0
        do chnk = begchunk, endchunk
            ncols = get_ncols_p(chnk)
            do col = 1, ncols
                dindex = dindex + 1
                decomp(dindex) = get_gcol_p(chnk, col)
            enddo
        enddo

        ! Build the 2D field CAM DistGrid based on the physics decomposition
        CAM_DistGrid = ESMF_DistGridCreate(arbSeqIndexList=decomp, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Release memory if any is being taken, to avoid leakage
        if(ESMF_MeshIsCreated(CAM_PhysMesh)) then
            call ESMF_MeshDestroy(CAM_PhysMesh)
        endif

        ! Create the physics decomposition ESMF mesh
        CAM_PhysMesh = ESMF_MeshCreate(trim(grid_file), ESMF_FILEFORMAT_ESMFMESH, &
                                       elementDistGrid=CAM_DistGrid, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        ! Todo: in edyn_gcomp they send the CAM_PhysMesh to the regridder in
        !       edyn_esmf. Maybe we should do that and send to hco_esmf_grid.
        ! But I think we can avoid doing that for now? The memory reclaiming
        ! part is now written before CAM_PhysMesh is built.rt (hplin, 2/6/20)

    end subroutine HCO_GC_Init

    ! Finalize Gridded Component
    subroutine HCO_GC_Final(GC, IMPORT, EXPORT, Clock, RC)
        use ESMF,         only: ESMF_MeshDestroy

        ! Dummy arguments
        type(ESMF_GridComp)                   :: GC
        type(ESMF_State)                      :: IMPORT
        type(ESMF_State)                      :: EXPORT
        type(ESMF_Clock)                      :: Clock
        integer, intent(out)                  :: RC

        ! Local variables
        character(len=*),       parameter     :: subname = 'HCO_GC_Final'

        call ESMF_MeshDestroy(CAM_PhysMesh, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)
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