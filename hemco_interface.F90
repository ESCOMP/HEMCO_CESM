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

    ! CAM import helpers
    use hco_cam_convert_state_mod

    ! Controls
    use cam_abortutils,           only: endrun      ! fatal terminator
    use cam_logfile,              only: iulog       ! output log handle

    ! Species information 
    use constituents,             only: pcnst       ! # of species
    use constituents,             only: cnst_name   ! species names
    use constituents,             only: cnst_mw     ! advected mass
    use mo_chem_utls,             only: get_spc_ndx ! IND_ equivalent

    ! Check chemistry option
    use chemistry,                only: chem_is

    ! Grid
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs

    ! Time
    use time_manager,             only: get_curr_time, get_prev_time, get_curr_date
    use time_manager,             only: get_step_size

    ! ESMF types
    use ESMF,                     only: ESMF_State, ESMF_Clock, ESMF_GridComp
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS

    ! HEMCO types
    use HCO_Error_Mod,            only: hp          ! HEMCO precision
    use HCO_Error_Mod,            only: sp          ! HEMCO single precision used for Ptrs
    use HCO_Error_Mod,            only: HCO_SUCCESS, HCO_FAIL, HCO_VERSION
    use HCO_State_Mod,            only: HCO_State
    use HCOX_State_Mod,           only: Ext_State
    use HCO_Types_Mod,            only: ConfigObj

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

    private :: HCOI_Initialize_Pbuf
!
! !PUBLIC MEMBER FUNCTIONS:
!
    public  :: hemco_readnl
    public  :: HCOI_Chunk_Init
    public  :: HCOI_Chunk_Run
    public  :: HCOI_Chunk_Final
!
! !REMARKS:
!  This file is both the interface of HEMCO component to CAM and the manager of the
!  underlying HEMCO gridded component (HCO\_GC\_*).
!
!  The whole file is getting a little long in the tooth, though. It might be reorg-
!  anized in the future to service hemco_init/run/final which manages the gridded
!  components, which in turn call the HEMCO chunk interface (HCOI\_Chunk\_*).
!
!  For now, the HEMCO chunk interfaces (HCOI\_Chunk\_*) are called in directly from
!  CAM cam_control.F90. I am not in the mood of writing a wrapper for now, but this
!  could be abstracted in the future. (hplin, 3/29/20)
!
!  On a side note, this is the 8th day of living in the 2019-nCoV scare. I miss
!  matcha tea and would die to eat some sweets.
!
!  It is now October 29 and COVID-19 is still raging across the globe; it feels
!  like time has frozen itself since March 21, and we've all been "on one long Zoom call"
!  ever since.
!
!  All I hope for 2021, is that the world does not pull a "You can now play as Luigi" on me.
!
!  Updated 2/24/21 from T. Fritz: Now uses constituents list, since short-term species are
!  not emitted. See PR: https://github.com/jimmielin/HEMCO_CESM/pull/5
!
!  It is now May 15, 2021 and I am vaccinated. Can we visit other peoples homes now?
!  ref: xkcd.com/2454
!
! !REVISION HISTORY:
!  29 Jan 2020 - H.P. Lin    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    type(ESMF_GridComp)              :: HCO_GridComp        ! HEMCO GridComp
    type(ESMF_State)                 :: HCO_GridCompState   ! HEMCO GridComp Import/Export State

    character(len=256)               :: HcoConfigFile       ! HEMCO configuration file loc
    type(ConfigObj), pointer         :: HcoConfig => NULL()

    type(HCO_State), pointer, public :: HcoState  => NULL()
    type(Ext_State), pointer, public :: ExtState  => NULL()

    ! HEMCO internal grid parameters
    integer                          :: HcoGridIM           ! # of lons
    integer                          :: HcoGridJM           ! # of lats

    ! HEMCO configuration parameters that are set by namelist in CESM
    integer                          :: HcoFixYY            ! if > 0, force 'Emission year'

    ! Last execution times for the HEMCO component. We are assuming that time
    ! flows unidirectionally (and forwards, for now). (hplin, 3/30/20)
    integer                          :: last_HCO_day, last_HCO_second

    ! Meteorological fields used by HEMCO to be regridded to the HEMCO grid (hplin, 3/31/20)
    ! We have to store the fields because the regridding can only take place within the GridComp.
    ! Fields are allocated after the internal grid is initialized (so my_* are avail)
    !
    ! Moved to hco_cam_convert_state_mod, 12/16/20

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
        integer                      :: hemco_grid_xdim = 0
        integer                      :: hemco_grid_ydim = 0
        integer                      :: hemco_emission_year = -1

        namelist /hemco_nl/ hemco_config_file, hemco_grid_xdim, hemco_grid_ydim, hemco_emission_year

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
            write(iulog,*) "hemco_readnl: hemco internal grid dimensions will be ", hemco_grid_xdim, " x ", hemco_grid_ydim

            if(hemco_emission_year .gt. 0) then
                write(iulog,*) "hemco_readnl: hemco will force emissions year at = ", hemco_emission_year
            endif
        endif

        ! MPI Broadcast Namelist variables
        call mpi_bcast(hemco_config_file, len(hemco_config_file), mpi_character, masterprocid, mpicom, ierr)
        call mpi_bcast(hemco_grid_xdim, 1, mpi_integer, masterprocid, mpicom, ierr)
        call mpi_bcast(hemco_grid_ydim, 1, mpi_integer, masterprocid, mpicom, ierr)
        call mpi_bcast(hemco_emission_year, 1, mpi_integer, masterprocid, mpicom, ierr)

        ! Save this to the module information
        HcoConfigFile = hemco_config_file
        HcoGridIM     = hemco_grid_xdim
        HcoGridJM     = hemco_grid_ydim
        HcoFixYY      = hemco_emission_year
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
        use cam_logfile,      only: iulog
        use spmd_utils,       only: masterproc, mpicom, masterprocid
        use spmd_utils,       only: npes, iam

        use mpi,              only: MPI_INTEGER
        use ESMF,             only: ESMF_VM, ESMF_VMGetCurrent, ESMF_VMGet
        use ESMF,             only: ESMF_GridCompCreate, ESMF_GridCompInitialize
        use ESMF,             only: ESMF_GridCompSetServices
        use ESMF,             only: ESMF_StateCreate

        use ESMF,             only: ESMF_Initialize, ESMF_LOGKIND_MULTI

        ! CAM instance information
        use cam_instance,     only: inst_index, inst_name

        ! CAM history output (to be moved somewhere later)
        use cam_history,      only: addfld, add_default, horiz_only

        ! Get extfrc list check for 3-D emissions capability
        use mo_chem_utls,     only: get_extfrc_ndx

        ! HEMCO Initialization
        use HCO_Config_Mod,   only: Config_ReadFile, ConfigInit
        use HCO_Driver_Mod,   only: HCO_Init
        use HCO_Error_Mod,    only: HCO_LOGFILE_OPEN
        use HCO_LogFile_Mod,  only: HCO_Spec2Log
        use HCO_State_Mod,    only: HcoState_Init
        use HCO_Types_Mod,    only: ConfigObj
        use HCO_Types_Mod,    only: ListCont
        use HCO_VertGrid_Mod, only: HCO_VertGrid_Define

        ! HEMCO extensions initialization
        use HCOX_Driver_Mod,  only: HCOX_Init
!
! !REMARKS:
!  None currently.
!
! !REVISION HISTORY:
!  06 Feb 2020 - H.P. Lin    - Initial version
!  15 Dec 2020 - H.P. Lin    - Implement HEMCO extensions that do require met
!  23 Mar 2021 - H.P. Lin    - Now export for diagnostics too
!  12 Jan 2023 - H.P. Lin    - Optimize memory usage by only allocating 2-D/3-D as needed
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Chunk_Init'
        integer                      :: RC                   ! ESMF return code
        integer                      :: HMRC                 ! HEMCO return code

        integer                      :: N                    ! Loop idx

        ! Gridded component properties.
        ! Note that while edyn_grid_comp initializes the GridComp directly using
        ! cam_instance's inst_name, I think this may cause namespace clashing.
        ! So I'll be prefixing this with hco_ just incase.
        character(len=32)            :: HCO_GC_InstName = ''
        type(ESMF_VM)                :: hco_esmf_vm

        ! MPI stuff
        integer                      :: localPET, PETcount
        integer, allocatable         :: PETlist(:)           ! PETs for each instance of the physics grid

        ! HEMCO properties
        integer                      :: nHcoSpc

        ! Temporary string for species naming in exports
        character(len=128)           :: exportName, exportDesc
        character(len=128)           :: exportNameTmp

        ! Timing properties
        integer                      :: year, month, day, tod
        integer                      :: hour, minute, second, dt
        integer                      :: prev_day, prev_s, now_day, now_s
        integer                      :: stepsize_tmp

        ! Temporaries
        logical                      :: IsExtfrc3DEmis

        !-----------------------------------------------------------------------

        if(masterproc) then
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            write(iulog,*) "HEMCO: Harmonized Emissions Component"
            write(iulog,*) "https://doi.org/10.5194/gmd-14-5487-2021 (Lin et al., 2021)"
            write(iulog,*) "HEMCO_CESM interface version 1.1.2"
            write(iulog,*) "You are using HEMCO version ", ADJUSTL(HCO_VERSION)
            write(iulog,*) "Config File: ", HcoConfigFile
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        endif

        ! Assume success
        RC = ESMF_SUCCESS
        HMRC = HCO_SUCCESS

        !-----------------------------------------------------------------------
        ! Get time properties
        !-----------------------------------------------------------------------
        call get_prev_time(prev_day, prev_s) ! 0 0
        call get_curr_time(now_day, now_s)   ! 0 0
        call get_curr_date(year, month, day, tod) ! 2005 1 1 0

        last_HCO_day    = now_day
        last_HCO_second = now_s   ! This means first timestep is not ran

        !-----------------------------------------------------------------------
        ! Setup ESMF wrapper gridded component
        ! Adapted from edyn_grid_comp_init
        !-----------------------------------------------------------------------
        call ESMF_VMGetCurrent(hco_esmf_vm, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        call ESMF_VMGet(hco_esmf_vm, localPet=localPET, petCount=PETcount, rc=RC)
        ASSERT_(RC==ESMF_SUCCESS)

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

        !if(masterproc) then
        !    write(iulog,*) "> Initialized ESMF environment successfully! localPET, PETcount", localPET, PETcount
        !    write(iulog,*) "> iam, npes", iam, npes
        !    write(iulog,*) "> PETlist", PETlist
        !endif

        !-----------------------------------------------------------------------
        ! Setup a lat-lon "HEMCO" intermediate grid
        !-----------------------------------------------------------------------
        ! TODO: For now, # of PEs to use for HEMCO will be total # of PEs.
        ! These will all have to be specified in the HEMCO namelist later on.

        ! The number of grid dimensions can be specified in atm namelist for
        ! hemco_grid_xdim, hemco_grid_ydim. note that half-sized polar boxes are made
        ! so increase the y-dim by one.
        ! 288x201 = 0.9x1.25
        ! 144x91  = 2.0x2.5

        ! Verify that the grid is in a reasonable state
        if(HcoGridIM .le. 1 .or. HcoGridJM .le. 1) then
            call endrun("Invalid HEMCO grid parameters - too small - in &hemco namelist. Specify hemco_grid_xdim and hemco_grid_ydim as # of grid boxes in each dimension")
        endif

        if(mod(HcoGridJM, 2) .ne. 1) then
            call endrun("Invalid HEMCO grid parameters - hemco_grid_ydim needs to be odd - in &hemco namelist. This is because y-dim has half-sized polar boxes.")
        endif

        ! Initialize the HEMCO intermediate grid
        call HCO_Grid_Init (IM_in = HcoGridIM, JM_in = HcoGridJM, nPET_in = npes, RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "> Initialized HEMCO Grid environment successfully!"
            write(iulog,*) "> Global Dimensions: ", HcoGridIM, HcoGridJM, LM
            write(iulog,*) "> my_IM, my_JM, LM, my_CE", my_IM, my_JM, LM, my_CE
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
        ! Allocate HEMCO meteorological objects
        ! We are allocating globally for the whole HEMCO component here. This may
        ! clash if my_CE changes (multiple CAM instances). To be verified.
        ! Should be a easy fix regardless, simply allocate and dealloc in the run
        ! (hplin, 3/31/20)
        !-----------------------------------------------------------------------
        call HCOI_Allocate_All()

        if(masterproc) then
            write(iulog,*) "> Allocated HEMCO temporary met fields"
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
        ! HCO_TEST outputs in the physics mesh.
        ! AvgFlag: (cam_history) A mean, B mean00z, I instant, X max, M min, S stddev
        !
        ! This call should eventually be reflected elsewhere?
        call addfld("DIAG_CAM_TEST", (/'lev'/), 'I', '1',          &
                    'HEMCO Debug, PETID written on CAM',             &
                    gridname="physgrid")

        ! Enable default for debugging:
        ! call add_default("DIAG_CAM_TEST", 2, 'I')      ! Make this field always ON

        call addfld("DIAG_HCO_TEST", (/'lev'/), 'I', '1',          &
                    'HEMCO Debug Data',             &
                    gridname="physgrid")

        ! Enable default for debugging:
        ! call add_default("DIAG_HCO_TEST", 2, 'I')      ! Make this field always ON

        !-----------------------------------------------------------------------
        ! Initialize the HEMCO configuration object...
        !-----------------------------------------------------------------------
        !if(masterproc) write(iulog,*) "> Initializing HCO configuration object"

        ! We are using pcnst here, which is # of constituents.
        nHcoSpc             = pcnst          ! # of hco species?

        call ConfigInit(HcoConfig, HMRC, nModelSpecies=nHcoSpc)
        ASSERT_(HMRC==HCO_SUCCESS)

        HcoConfig%amIRoot   = masterproc
        ! HcoConfig%amIRoot   = .true. ! for debug only so verbosity is higher

        HcoConfig%MetField  = 'MERRA2'
        HcoConfig%GridRes   = ''

        !-----------------------------------------------------------------------
        ! Retrieve the species list and register exports
        !-----------------------------------------------------------------------
        ! Below we directly use nHcoSpc which corresponds to the number of constituents
        ! (nHcoSpc = pcnst). Only constituents may be advected.
        HcoConfig%nModelSpc = nHcoSpc
        HcoConfig%nModelAdv = nHcoSpc            ! # of adv spc?

        !if(masterproc) write(iulog,*) "> Initializing HCO species list!"

        do N = 1, nHcoSpc
            HcoConfig%ModelSpc(N)%ModID   = N ! model id

            HcoConfig%ModelSpc(N)%SpcName = trim(cnst_name(N)) ! only constituents can be emitted

            !----------------------------------------------
            ! Register export properties.
            !----------------------------------------------
            ! History output (this will be moved to hco_cam_exports soon hopefully)
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(N)%SpcName)

            !if(masterproc) write(iulog,*) "Exported exportName " // trim(exportName) // " to history"

            ! Physics buffer
            ! Note that _AddField will prepend HCO_, so do not add it here
            !
            ! Update hplin 1/13/23: Verify if part of extfrc_lst. If yes,
            ! then allocate as 3-D. Otherwise, this field can be allocated
            ! as 2-D. This scan is somewhat intensive as it uses get_extfrc_ndx
            ! which loops through extcnt in extfrc_lst.
            IsExtfrc3DEmis = (get_extfrc_ndx(trim(HcoConfig%ModelSpc(N)%SpcName)) .gt. 0)

            if(masterproc) write(iulog,*) "Setting up HEMCO exportName " // trim(exportName), IsExtfrc3DEmis
            if(IsExtfrc3DEmis) then
                ! 3-D emissions are supported
                HcoConfig%ModelSpc(N)%DimMax = 3

                exportDesc = "HEMCO 3-D Emissions Species " // trim(HcoConfig%ModelSpc(N)%SpcName)
                call addfld(exportName, (/'lev'/), 'I', 'kg/m2/s',          &
                            trim(exportDesc),                               &
                            gridname='physgrid')
                call HCO_Export_Pbuf_AddField(HcoConfig%ModelSpc(N)%SpcName, 3, hcoID=N)
            else
                ! 2-D emissions only
                HcoConfig%ModelSpc(N)%DimMax = 2

                exportDesc = "HEMCO 2-D Emissions Species " // trim(HcoConfig%ModelSpc(N)%SpcName)
                call addfld(exportName, horiz_only, 'I', 'kg/m2/s',          &
                            trim(exportDesc),                               &
                            gridname='physgrid')
                call HCO_Export_Pbuf_AddField(HcoConfig%ModelSpc(N)%SpcName, 2, hcoID=N)
            endif
        enddo

        !-----------------------------------------------------------------------
        ! Read HEMCO configuration file from HcoConfigFile (location in CAM namelist)
        !-----------------------------------------------------------------------
        !if(masterproc) write(iulog,*) "> Reading HEMCO configuration file..."

        ! FIXME: Not implementing "Dry-run" functionality in HEMCO_CESM. (hplin, 3/27/20)
        ! Phase: 0 = all, 1 = sett and switches only, 2 = fields only
        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 1, HMRC, IsDryRun=.false.)
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CESM: Config_ReadFile has failed (1)!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "HEMCO configuration file could not be read."
            write(iulog,*) "This may be due to misconfiguration of the"
            write(iulog,*) "hemco_config_file namelist variable, or a"
            write(iulog,*) "misformatted HEMCO configuration file."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        ! Open the log file
        if(masterproc) then
            call HCO_LOGFILE_OPEN(HcoConfig%Err, RC=HMRC)
            ASSERT_(HMRC==HCO_SUCCESS)
        endif

        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 2, HMRC, IsDryRun=.false.)
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CESM: Config_ReadFile has failed (2)!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "HEMCO configuration file could not be read."
            write(iulog,*) "This may be due to misconfiguration of the"
            write(iulog,*) "hemco_config_file namelist variable, or a"
            write(iulog,*) "misformatted HEMCO configuration file."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        !if(masterproc) write(iulog,*) "> Read HEMCO configuration file OK!"

        !-----------------------------------------------------------------------
        ! Initialize the HEMCO state object
        !-----------------------------------------------------------------------
        call HcoState_Init(HcoState, HcoConfig, nHcoSpc, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        !if(masterproc) write(iulog,*) "> Initialize HEMCO state obj OK!"

        ! Emissions, chemistry and dynamics timestep [s]
        ! Assume 0.5h until given actual time in HCO_GC_Run!

        stepsize_tmp = get_step_size()
        if(masterproc) write(iulog,*) "> HEMCO_CESM: Step size is ", stepsize_tmp

        HcoState%TS_EMIS = stepsize_tmp * 1.0
        HcoState%TS_CHEM = stepsize_tmp * 1.0
        HcoState%TS_DYN  = stepsize_tmp * 1.0

        ! Not a MAPL simulation. isESMF is deceiving.
        HcoState%Options%isESMF = .false.

        ! Deposition length scale. Used for computing dry deposition frequencies
        ! over the entire PBL or the first model layer. Hardcoded for now,
        ! should load Input_Opt%PBL_DRYDEP from GEOS-Chem-CESM (hplin, 3/29/20)
        ! !FIXME
        HcoState%Options%PBL_DRYDEP = .false.

        ! Don't support DryRun option (for now)
        HcoState%Options%IsDryRun = .false.

        !if(masterproc) write(iulog,*) "> Set basic HEMCO state obj OK!"

        !-----------------------------------------------------------------------
        ! Register HEMCO species information (HEMCO state object)
        !-----------------------------------------------------------------------
        do N = 1, nHcoSpc
            HcoState%Spc(N)%ModID         = N               ! model id
            HcoState%Spc(N)%SpcName       = trim(cnst_name(N)) ! species name
            HcoState%Spc(N)%MW_g          = cnst_mw(N)     ! mol. weight [g/mol]

            ! !!! We don't set Henry's law coefficients in HEMCO_CESM !!!
            ! they are mostly used in HCOX_SeaFlux_Mod, but HCOX are unsupported (for now)
            ! (hplin, 3/29/20)

            ! If is CESM-GC, then set Henry's law constants for SeaFlux
            ! KLUDGE by hplin: 1/3/21
            ! For defined species, hard code the Henry* values for now so we can work
            ! with SeaFlux. This fragmentation will cause issues down the road, FIXME
            if(HcoState%Spc(N)%SpcName .eq. "CH3I") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 0.20265_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 3.6e+3_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "DMS") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 4.80e-1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 3100.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "ACET") then
                ! 101.325_r8. using new henry constants
                HcoState%Spc(N)%HenryK0  = 2.74e+1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5500.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "MOH") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 2.03e+2_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5600.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            if(HcoState%Spc(N)%SpcName .eq. "ALD2") then
                ! 101.325_r8. using new henry constants
                HcoState%Spc(N)%HenryK0  = 1.30e-01_r8 * 101.325_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5900.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif
            
            if(HcoState%Spc(N)%SpcName .eq. "MENO3") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 1.1e+1_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 4700.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif
            
            if(HcoState%Spc(N)%SpcName .eq. "ETNO3") then
                ! 101.325_r8
                HcoState%Spc(N)%HenryK0  = 1.6_r8 ! [M/atm]
                HcoState%Spc(N)%HenryCR  = 5400.0_r8 ! [K]
                HcoState%Spc(N)%HenryPKA = -999e+0_r8 ! [1] (missing_r8 from species_mod)
            endif

            ! HcoState%Spc(N)%HenryK0 ! [M/atm]
            ! HcoState%Spc(N)%HenryCR ! [K]
            ! HcoState%Spc(N)%HenryPKA ! [1]

            ! Write to log too
            if(masterproc) then
                !write(iulog,*) ">> Spc", N, " = ", cnst_name(N), "MW_g", cnst_mw(N)
                call HCO_Spec2Log(HcoState, N)
            endif
        enddo

        !if(masterproc) write(iulog,*) "> Set HEMCO species info OK!"

        !-----------------------------------------------------------------------
        ! Register HEMCO Grid information
        !-----------------------------------------------------------------------
        ! Note that HEMCO running in the CAM environment is entirely MPI and
        ! we use the grid dimensions of the local PET. Thus, remember that all
        ! data and fields are sized my_* and NOT the global indices, although
        ! all PETs are aware. This is similar to ids, ide, jds, jde ... versus
        ! its, ite, jts, jte ... in WRF, but here we use my_IS, my_IE, my_JS...
        !
        ! The vertical dimension is not decomposed and follows the CAM vertical,
        ! whatever that is. This information is all abstracted and propagated within
        ! HCO_ESMF_Grid. (hplin, 3/29/20)

        HcoState%NX = my_IM
        HcoState%NY = my_JM
        HcoState%NZ = LM

        ! Pass Ap, Bp values, units [Pa], [unitless]
        ! later remove masterproc
        call HCO_VertGrid_Define(HcoState%Config,                &
                                 zGrid = HcoState%Grid%zGrid,    &
                                 nz    = HcoState%NZ,            &
                                 Ap    = Ap,                     &
                                 Bp    = Bp,                     &
                                 RC    = HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        ! Point to grid variables
        HcoState%Grid%XMID%Val         => XMid   (my_IS:my_IE  , my_JS:my_JE  )
        HcoState%Grid%YMID%Val         => YMid   (my_IS:my_IE  , my_JS:my_JE  )
        HcoState%Grid%XEdge%Val        => XEdge  (my_IS:my_IE+1, my_JS:my_JE  )
        HcoState%Grid%YEdge%Val        => YEdge  (my_IS:my_IE  , my_JS:my_JE+1)
        HcoState%Grid%YSin%Val         => YSin   (my_IS:my_IE  , my_JS:my_JE+1)
        HcoState%Grid%AREA_M2%Val      => AREA_M2(my_IS:my_IE  , my_JS:my_JE  )

        ! Debug
        ! write(6,*) "HCOI_Chunk_Init XMid, YMid(1,1)", HcoState%Grid%XMid%Val(1,1), &
        !                                               HcoState%Grid%YMid%Val(1,1), &
        !                                               HcoState%Grid%Area_m2%Val(1,1)

        !if(masterproc) write(iulog,*) "> Set HEMCO PET-local grid info OK!"

        !-----------------------------------------------------------------------
        ! Initialize HEMCO!
        !-----------------------------------------------------------------------
        call HCO_Init(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "> HEMCO initialized successfully!"

        !-----------------------------------------------------------------------
        ! Initialize HEMCO Extensions!
        !-----------------------------------------------------------------------
        call HCOX_Init(HcoState, ExtState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "> HEMCO extensions initialized successfully!"

        !-----------------------------------------------------------------------
        ! Additional exports: Verify if additional diagnostic quantities
        ! for HEMCO extensions need to be provisioned.
        ! (hplin, 3/21/21)
        !-----------------------------------------------------------------------
        ! ParaNOx: Ship NO emissions
        ! Due to the length limit this is a non-standard name, the HEMCO names are
        ! PARANOX_O3_DEPOSITION_FLUX and PARANOX_HNO3_DEPOSITION_FLUX
        if(ExtState%ParaNOx > 0) then
            write(exportnameTmp, '(a)') 'PAR_O3_DEP'
            exportName = 'HCO_' // trim(exportNameTmp)
            exportDesc = "HEMCO Deposition Flux Name " // trim(exportNameTmp)

            call addfld(exportName, horiz_only, 'I', '1',                &
                        trim(exportDesc),                                &
                        gridname='physgrid')
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            write(exportnameTmp, '(a)') 'PAR_HNO3_DEP'
            exportName = 'HCO_' // trim(exportNameTmp)
            exportDesc = "HEMCO Deposition Flux Name " // trim(exportNameTmp)

            call addfld(exportName, horiz_only, 'I', '1',                &
                        trim(exportDesc),                                &
                        gridname='physgrid')
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            if(masterproc) write(iulog,*) "> HEMCO ParaNOx extension exports (PAR_O3_DEP, PAR_HNO3_DEP) initialized successfully"
        endif

        !-----------------------------------------------------------------------
        ! Additional exports: Verify if we need to add additional exports
        ! for integration with CESM-GC. (hplin, 4/15/20)
        !-----------------------------------------------------------------------

        ! Do additional exports!
        ! Removed debug into history output because it does not seem necessary (hplin, 9/23/22)
        if(chem_is('GEOS-Chem')) then
            do N = 0, 72
               ! LANDTYPExx
               write(exportNameTmp, '(a,i2.2)') 'LANDTYPE', N
               exportName = 'HCO_' // trim(exportNameTmp)
               call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

               ! XLAIxx
               write(exportNameTmp, '(a,i2.2)') 'XLAI', N
               exportName = 'HCO_' // trim(exportNameTmp)
               call HCO_Export_Pbuf_AddField(exportNameTmp, 2)
            enddo

            ! VMR_CH3CL
            write(exportnameTmp, '(a)') 'VMR_CH3CL'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CH2CL2
            write(exportnameTmp, '(a)') 'VMR_CH2CL2'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CHCL3
            write(exportnameTmp, '(a)') 'VMR_CHCL3'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CH3BR
            write(exportnameTmp, '(a)') 'VMR_CH3BR'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CCL4
            write(exportnameTmp, '(a)') 'VMR_CCL4'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CFC113
            write(exportnameTmp, '(a)') 'VMR_CFC113'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CFC114
            write(exportnameTmp, '(a)') 'VMR_CFC114'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CFC115
            write(exportnameTmp, '(a)') 'VMR_CFC115'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CFC11
            write(exportnameTmp, '(a)') 'VMR_CFC11'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CFC12
            write(exportnameTmp, '(a)') 'VMR_CFC12'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_CH3CCL3
            write(exportnameTmp, '(a)') 'VMR_CH3CCL3'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_H1211
            write(exportnameTmp, '(a)') 'VMR_H1211'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_H1301
            write(exportnameTmp, '(a)') 'VMR_H1301'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_H2402
            write(exportnameTmp, '(a)') 'VMR_H2402'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_HCFC141b
            write(exportnameTmp, '(a)') 'VMR_HCFC141B'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_HCFC142B
            write(exportnameTmp, '(a)') 'VMR_HCFC142B'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_HCFC22
            write(exportnameTmp, '(a)') 'VMR_HCFC22'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_N2O
            write(exportnameTmp, '(a)') 'VMR_N2O'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_OCS
            write(exportnameTmp, '(a)') 'VMR_OCS'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! VMR_H2
            write(exportnameTmp, '(a)') 'VMR_H2'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! UVALBEDO
            write(exportnameTmp, '(a)') 'UV_ALBEDO'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! SURF_IODIDE
            write(exportnameTmp, '(a)') 'iodide'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! SURF_SALINITY
            write(exportnameTmp, '(a)') 'salinity'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! OMOC_DJF
            write(exportnameTmp, '(a)') 'OMOC_DJF'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! OMOC_MAM
            write(exportnameTmp, '(a)') 'OMOC_MAM'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! OMOC_JJA
            write(exportnameTmp, '(a)') 'OMOC_JJA'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            ! OMOC_SON
            write(exportnameTmp, '(a)') 'OMOC_SON'
            exportName = 'HCO_' // trim(exportNameTmp)
            call HCO_Export_Pbuf_AddField(exportNameTmp, 2)

            if ( masterproc ) then
                write(iulog,*) "> HEMCO additional exports for CESM2-GC initialized!"
            endif
        endif
    end subroutine HCOI_Chunk_Init
!EOC
!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Initialize_Pbuf
!
! !DESCRIPTION: HCOI\_Initialize\_Pbuf resets the physics buffer.
!\\
!\\
! !INTERFACE:
!
    subroutine HCOI_Initialize_Pbuf()
!
! !USES:
!
        use cam_logfile,      only: iulog
        use spmd_utils,       only: masterproc
!
! !REVISION HISTORY:
!  14 Dec 2020 - H.P. Lin    - Initial version
!  16 Jan 2023 - H.P. Lin    - Need to reset this after pbuf fields are now 3-D or 2-D.
!                              Approximately 16 hours of debugging stack corruptions
!                              were needed to realize I wrote this routine.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*), parameter  :: subname = 'HCOI_Initialize_Pbuf'
        integer                      :: RC                   ! ESMF return code

        integer                      :: spcID
        real(ESMF_KIND_R8)           :: zeroFldCAM_3D(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)           :: zeroFldCAM_2D(1:my_CE)

        ! Zero out quantities first
        zeroFldCAM_3D(:,:)   = 0.0_r8
        zeroFldCAM_2D(:  )   = 0.0_r8

        ! Reset for each species
        do spcID = 1, HcoConfig%nModelSpc
            ! Write to physics buffer (pass model name)
            if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                call HCO_Export_Pbuf_CAM3D(HcoConfig%ModelSpc(spcID)%SpcName, spcID, zeroFldCAM_3D)
            elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                call HCO_Export_Pbuf_CAM2D(HcoConfig%ModelSpc(spcID)%SpcName, spcID, zeroFldCAM_2D)
            else
                ASSERT_(.false.)
            endif
        enddo

    end subroutine HCOI_Initialize_Pbuf
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

        ! Physics grid
        use phys_grid,      only: get_ncols_p

        ! CAM physics buffer (some fields are here and some are in phys state)
        use physics_buffer, only: pbuf_get_chunk, pbuf_get_field
        use physics_buffer, only: pbuf_get_index

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

        logical, save                :: FIRST = .true.

        if(masterproc) then
            write(iulog,*) "HEMCO_CESM: Running HCOI_Chunk_Run phase", phase
        endif

        ! For phase 1, before chemistry, reset all the physics buffer contents
        ! to prevent trash data being read by other components.
        ! (hplin, 12/15/20)
        if(phase == 1) then

            ! Only need to do this once to save time
            if(FIRST) then

                ! Pass the pbuf to the hco_cam_exports component so she has it...
                hco_pbuf2d => pbuf2d

                call HCOI_Initialize_Pbuf()

                ! No longer first call
                FIRST = .false.

            endif
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
            ! Pass the pbuf to the hco_cam_exports component so she has it...
            hco_pbuf2d => pbuf2d

            ! Set fields from CAM state before run
            call CAM_GetBefore_HCOI(cam_in, phys_state, pbuf2d, phase, &
                                    HcoState, ExtState)

            ! Run the gridded component.
            call ESMF_GridCompRun(HCO_GridComp, rc=RC)!importState=HCO_GridCompState, &
                                                !exportState=HCO_GridCompState, &
                                                !rc=RC)

            ASSERT_(RC==ESMF_SUCCESS)
        endif

        if(masterproc) then
            write(iulog,*) "HEMCO_CESM: Leaving HCOI_Chunk_Run"
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
        use cam_logfile,            only: iulog
        use spmd_utils,             only: masterproc, iam

        ! HEMCO
        use HCO_Interface_Common,   only: GetHcoVal, GetHcoDiagn
        use HCO_Clock_Mod,          only: HcoClock_Set, HcoClock_Get
        use HCO_Clock_Mod,          only: HcoClock_EmissionsDone
        use HCO_Diagn_Mod,          only: HcoDiagn_AutoUpdate
        use HCO_Driver_Mod,         only: HCO_Run
        use HCO_EmisList_Mod,       only: Hco_GetPtr
        use HCO_Calc_Mod,           only: Hco_EvalFld
        use HCO_FluxArr_Mod,        only: HCO_FluxArrReset
        use HCO_GeoTools_Mod,       only: HCO_CalcVertGrid, HCO_SetPBLm

        use HCO_State_Mod,          only: HCO_GetHcoId

        ! HEMCO Extensions
        use HCOX_Driver_Mod,        only : HCOX_Run

        ! Physical constants
        use physconst,              only : mwdry


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
!  09 Mar 2023 - H.P. Lin    - Allow for FixYY in HEMCO clock for climo runs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
        character(len=*),       parameter     :: subname = 'HCO_GC_Run'

        ! Parameters
        real(r8), parameter                   :: G0_100 = 100.e+0_r8 / 9.80665e+0_r8

        integer                               :: I, J, K
        integer                               :: HI, HJ, HL

        ! Temporaries for exports
        real(ESMF_KIND_R8)                    :: TMP
        logical                               :: FND
        character(len=128)                    :: exportName, exportNameTmp
        integer                               :: spcID, N

        ! For grabbing data from HEMCO Ptrs (uses HEMCO single-precision)
        real(sp), pointer                     :: Ptr2D(:,:)
        real(sp), pointer                     :: Ptr3D(:,:,:)

        ! Temporaries used for export
        real(ESMF_KIND_R8)                    :: exportFldHco(my_IS:my_IE, my_JS:my_JE, 1:LM)
        real(ESMF_KIND_R8)                    :: exportFldCAM(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)                    :: scratchFldCAM(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)                    :: scratchFldCAM2(1:my_CE)

        ! Temporaries used for export (2-D data)
        real(ESMF_KIND_R8)                    :: exportFldHco2(my_IS:my_IE, my_JS:my_JE)
        real(ESMF_KIND_R8)                    :: exportFldCAM2(1:my_CE)
        real(ESMF_KIND_R8)                    :: exportFldHcoDep(my_IS:my_IE, my_JS:my_JE)
        real(ESMF_KIND_R8)                    :: exportFldCAMDep(1:my_CE)

        ! For debug dummies
        real(ESMF_KIND_R8)                    :: dummy_0_CAM(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)                    :: dummy_1(my_IS:my_IE, my_JS:my_JE, 1:LM)
        real(ESMF_KIND_R8)                    :: dummy_1_CAM(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)                    :: dummy_2(my_IS:my_IE, my_JS:my_JE)

        ! Timing properties
        integer                               :: year, month, day, tod
        integer                               :: hour, minute, second
        integer                               :: prev_day, prev_s, now_day, now_s
        integer                               :: tmp_currTOD

        ! HEMCO vertical grid property pointers
        ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
        ! make PEDGE be a pointer and allocate/deallocate it on each call.
        real(hp), pointer            :: BXHEIGHT(:,:,:)    ! Grid box height      [m ]
        real(hp), pointer            :: PEDGE   (:,:,:)    ! Pressure @ lvl edges [Pa]
        real(hp), pointer            :: ZSFC    (:,:  )    ! Surface geopotential [m ]


        real(hp), pointer            :: PBLM    (:,:  )    ! PBL height           [m ]
        real(hp), pointer            :: PSFC    (:,:  )    ! Surface pressure     [Pa]
        real(hp), pointer            :: TK      (:,:,:)    ! Temperature          [K ]

        ! HEMCO return code
        integer                      :: HMRC

        logical, save                :: FIRST = .True.
        logical                      :: doExport = .False.
        integer, save                :: nCalls = 0

        ! Assume success
        RC = ESMF_SUCCESS
        HMRC = HCO_SUCCESS

        nCalls = nCalls + 1

        !-----------------------------------------------------------------------
        ! Update regridding file handles as necessary
        !-----------------------------------------------------------------------
        call HCO_Grid_UpdateRegrid(RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        !if(masterproc) then
        !    write(iulog,*) "HEMCO_CESM: Reload (if necessary) of HEMCO Regrid descriptors"
        !endif

        !-----------------------------------------------------------------------
        ! Allow for year forcing in climatological runs. (hplin, 3/9/23)
        ! This is set if HcoFixYY is set in this module and is > 0
        !-----------------------------------------------------------------------
        if(HcoFixYY .gt. 0) then
            ! Override the HEMCO clock FixYY property
            HcoState%Clock%FixYY = HcoFixYY
        endif

        !-----------------------------------------------------------------------
        ! Get time properties
        !-----------------------------------------------------------------------
        call get_prev_time(prev_day, prev_s)
        call get_curr_time(now_day, now_s)
        call get_curr_date(year, month, day, tod)

        !if(masterproc) then
        !    write(iulog,*) "hco year,month,day,tod", year, month, day, tod
        !    write(iulog,*) "hco prev_day, prev_s", prev_day, prev_s
        !    write(iulog,*) "hco now_day, now_s", now_day, now_s
        !endif
        ! 2005 1 1 1800 | 0 0 | 0 1800
        ! 2005 1 1 3600 | 0 1800 | 0 3600
        ! 2005 1 1 5400 | 0 3600 | 0 5400
        ! ...

        ! Check if we have run HEMCO for this time step already. If yes can exit
        if(last_HCO_day * 86400.0 + last_HCO_second .ge. now_day * 86400.0 + now_s) then
            ! But also do not skip the first time step
            ! FIXME: Implicit assumption of time stepping size being 1800.0
            ! FIXME hplin 2/28/21

            if(masterproc) then
                write(iulog,*) "HEMCO_CESM: !! HEMCO already ran for this time, check timestep mgr", now_day, now_s, last_HCO_day, last_HCO_second
            endif

            return
        endif

        ! Compute timestep
        if(HcoState%TS_CHEM .ne. ((now_day - prev_day) * 86400.0 + now_s - prev_s)) then
            HcoState%TS_EMIS = (now_day - prev_day) * 86400.0 + now_s - prev_s
            HcoState%TS_CHEM = (now_day - prev_day) * 86400.0 + now_s - prev_s
            HcoState%TS_DYN  = (now_day - prev_day) * 86400.0 + now_s - prev_s

            if(masterproc) then
                write(iulog,*) "HEMCO_CESM: Updated HEMCO timestep to ", HcoState%TS_CHEM
            endif
        endif

        ! Compute hour, minute, second (borrowed from tfritz)
        tmp_currTOD = tod
        hour = 0
        minute = 0
        do while(tmp_currTOD >= 3600)
            tmp_currTOD = tmp_currTOD - 3600
            hour = hour + 1
        enddo

        do while(tmp_currTOD >= 60)
            tmp_currTOD = tmp_currTOD - 60
            minute = minute + 1
        enddo
        second = tmp_currTOD

        ! Update HEMCO clock
        ! using HcoClock_Set and not common SetHcoTime because we don't have DOY
        ! and we want HEMCO to do the math for us. oh well

        !if(masterproc) then
        !    write(6,*) "HEMCO_CESM: Updating HEMCO clock to set", year, month, day, hour, minute, second
        !    write(6,*) "HEMCO_CESM: Internally HEMCO is at ", HcoState%Clock%SimHour, HcoState%Clock%SimMin, HcoState%Clock%SimSec, HcoState%Clock%nSteps
        !endif

        call HCOClock_Set(HcoState, year, month, day,  &
                          hour, minute, second, IsEmisTime=.true., RC=HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        !-----------------------------------------------------------------------
        ! Continue setting up HEMCO
        !-----------------------------------------------------------------------

        ! Reset all emission and deposition values.
        call HCO_FluxArrReset(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)
        
        !-----------------------------------------------------------------------
        ! Regrid necessary meteorological quantities (Phase 1)
        ! Computes the absolute minimum (PSFC and TK) necessary for HEMCO
        ! to define its grid.
        !-----------------------------------------------------------------------
        call CAM_RegridSet_HCOI(HcoState, ExtState, Phase=1)

        if(masterproc .and. nCalls < 10) then
            write(iulog,*) "HEMCO_CESM: Finished regridding CAM met fields to HEMCO (1)"
        endif

        !-----------------------------------------------------------------------
        ! Get grid properties
        !-----------------------------------------------------------------------
        ! Calculate HEMCO vertical grid properties, e.g. PEDGE,
        ! then PHIS, BXHEIGHT, T, ..., from GridEdge_Set
        !
        ! The below conversions mostly borrowed from tfritz's CESM2-GC.
        ! HEMCO CalcVertGrid can approximate all quantities. We provide them with
        ! PSFC (surface pressure), TK (temperature)
        ! in the form of allocated pointers. The rest can be inferred from Ap, Bp
        PSFC => State_HCO_PSFC
        TK   => State_HCO_TK

        call HCO_CalcVertGrid(HcoState, PSFC, ZSFC, TK, BXHEIGHT, PEDGE, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        ! Pass boundary layer height to HEMCO (PBLm = PBL mixing height) [m]
        call HCO_SetPBLm(HcoState, PBLM=State_HCO_PBLH, &
                         DefVal=1000.0_hp, & ! default value
                         RC=HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)
        
        !-----------------------------------------------------------------------
        ! Regrid necessary meteorological quantities (Phase 2)
        ! Has to be below grid properties because vertical grid needs to be defined
        ! for quantities to be computed!
        !-----------------------------------------------------------------------
        call CAM_RegridSet_HCOI(HcoState, ExtState, Phase=2)

        if(masterproc .and. nCalls < 10) then
            write(iulog,*) "HEMCO_CESM: Finished regridding CAM met fields to HEMCO (2)"

            ! As a test... maybe we also need to flip in the vertical
            ! write(iulog,*) State_HCO_TK(1,1,:)
            ! write(iulog,*) "PSFC(1:2,:)"
            ! write(iulog,*) State_HCO_PSFC(1:2,:) 

            !write(iulog,*) "cam state%ps dump"
            !write(iulog,*) State_CAM_ps

            ! TK: 288 283 277 271 266 261 ... 250 251 252
            ! Seems like the vertical is OK for now
        endif

        !-----------------------------------------------------------------------
        ! Set HEMCO options
        !-----------------------------------------------------------------------
        ! Range of species and emission categories.
        ! Set Extension number ExtNr to 0, indicating that the core
        ! module shall be executed.
        HcoState%Options%SpcMin = 1
        HcoState%Options%SpcMax = -1
        HcoState%Options%CatMin = 1
        HcoState%Options%CatMax = -1
        HcoState%Options%ExtNr  = 0

        ! Use temporary array?
        HcoState%Options%FillBuffer = .FALSE.

        !-----------------------------------------------------------------------
        ! Run HEMCO!
        !-----------------------------------------------------------------------

        ! Run HCO core module
        ! Pass phase as argument. Phase 1 will update the emissions list,
        ! phase 2 will calculate the emissions. Emissions will be written into
        ! the corresponding flux arrays in HcoState.
        !
        ! FIXME: hplin - setting false as last timestep of simulation. maybe see
        ! if we can figure out from CAM if we are at run end and set to true
        call HCO_Run( HcoState, 1, HMRC, IsEndStep=.false. )
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CESM: HCO_Run Phase 1 has failed!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "A critical component in HEMCO failed to run."
            write(iulog,*) "This may be due to misconfiguration, or a bug."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: HCO_Run Phase 1"

        call HCO_Run( HcoState, 2, HMRC, IsEndStep=.false. )
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CESM: HCO_Run Phase 2 has failed!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "A critical component in HEMCO failed to run."
            write(iulog,*) "This may be due to misconfiguration, or a bug."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: HCO_Run Phase 2"

        !-----------------------------------------------------------------------
        ! Run HEMCO Extensions!
        !-----------------------------------------------------------------------
        call HCOX_Run(HcoState, ExtState, HMRC)
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CESM: HCOX_Run (extensions) has failed!"
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!"
            write(iulog,*) "A critical component in HEMCO failed to run."
            write(iulog,*) "This may be due to misconfiguration, or a bug."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: HCOX_Run"

        !-----------------------------------------------------------------------
        ! Update "autofill" diagnostics.
        ! Update all 'AutoFill' diagnostics. This makes sure that all
        ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
        ! AutoFill flag is specified when creating a diagnostics container
        ! (Diagn_Create).
        !-----------------------------------------------------------------------
        call HcoDiagn_AutoUpdate(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        !if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: HcoDiagn_AutoUpdate"

        !-----------------------------------------------------------------------
        ! Tell HEMCO we are done for this timestep...
        !-----------------------------------------------------------------------
        call HcoClock_EmissionsDone(HcoState%Clock, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        !if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: HcoClock_EmissionsDone"

        !-----------------------------------------------------------------------
        ! Do some testing and write emissions to the tape
        !-----------------------------------------------------------------------

        ! HCO Index boundaries
        HI = my_IE - my_IS + 1
        HJ = my_JE - my_JS + 1
        HL = LM

        ! For each species...
        do spcID = 1, HcoConfig%nModelSpc

            ! TODO: Eventually convert aerosol number emissions from mass fluxes
            ! directly rather than using scale factors for num_ax (1/12/23, hplin)

            ! Build history / pbuf field name (HCO_NO, HCO_CO, etc.)
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(spcID)%SpcName)
            doExport   = (FIRST .or. associated(HcoState%Spc(spcID)%Emis%Val))
            ! if(masterproc) write(iulog,*) "HEMCO_CESM: Begin exporting " // trim(exportName)

            ! Get HEMCO emissions flux [kg/m2/s].
            ! For performance optimization ... tap into HEMCO structure directly (ugly ugly)
            ! No need to flip vertical here. The regridder will do it for us
            if(associated(HcoState%Spc(spcID)%Emis%Val)) then
                if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                    ! Zero out quantities first
                    exportFldHco(:,:,:) = 0.0_r8
                    exportFldCAM(:,:)   = 0.0_r8

                    ! Retrieve flux from HEMCO...
                    exportFldHco(my_IS:my_IE,my_JS:my_JE,1:LM) = HcoState%Spc(spcID)%Emis%Val(1:HI,1:HJ,1:LM)

                    ! Regrid exportFldHco to CAM grid...
                    call HCO_Grid_HCO2CAM_3D(exportFldHco, exportFldCAM)
                elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                    ! Zero out quantities first
                    exportFldHco2(:,:) = 0.0_r8
                    exportFldCAM2(:)   = 0.0_r8

                    ! Retrieve flux from HEMCO...
                    ! Only surface emissions are supported here. Note, HEMCO emissions are with 1 = surface and LM = TOA. Use index 1. (hplin, 1/12/23)
                    exportFldHco2(my_IS:my_IE,my_JS:my_JE)     = HcoState%Spc(spcID)%Emis%Val(1:HI,1:HJ,1)

                    ! Regrid exportFldHco to CAM grid...
                    call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
                else
                    ASSERT_(.false.)
                endif

                !if(masterproc) write(iulog,*) "HEMCO_CESM: Retrieved from HCO " // trim(exportName)
            else
                ! No emission value. No need to run regrid, instead populate with zeros as needed
                ! Why not populate at top, you ask? Because zeroing out arrays is expensive, and
                ! we do not want to be doing twice the work.
                if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                    exportFldCAM(:,:) = 0.0_r8
                elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                    exportFldCAM2(:)  = 0.0_r8
                else
                    ASSERT_(.false.)
                endif
            endif

            ! Handle deposition flux from deposition velocity [1/s]
            ! This can be performed on the CAM grid. (hplin, 5/7/21)
            ! GEOS-Chem: (from v/v dry)
            !  Step 1)
            !              dflx(I,J,NA) = dflx(I,J,NA)                                     &
            !              + ( dep * spc(I,J,NA) / (AIRMW / ThisSpc%MW_g)  )
            ! for kg/kg dry, no conversion is needed:
            !  Step 1a)    dflx = dflx + dep * spc
            !
            !  Step 2) Convert to 1/s
            !             dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1)                        &
            !                     / State_Grid%Area_M2(I,J)
            !
            !  Note that, AD is actually DELP_DRY * G0_100 * AREA_M2, thus the final expression is
            !  just multiplied by DELP_DRY * G0_100, unit: kg/m2 (delp_dry is hPa, g0_100 is 100 Pa/hPa * s2/m --> unit = Pa*m/s2 = kg/m/s2*s2/m = kg/m2)
            !  Multiplied by 1/s, this gives kg/m2/s
            !
            ! We retrieve the concentration flux read from the convert state module
            ! (on the CAM grid) and loop through it to apply deposition fluxes.
            if(associated(HcoState%Spc(spcID)%Depv%Val)) then
                ! Check if species is available for deposition (hard-coded)
                ! FIXME: hplin 5/7/21
                if(trim(HcoConfig%ModelSpc(spcID)%SpcName) == "DMS" .or.            &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ACET" .or.           &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3COCH3" .or.       &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ALD2" .or.           &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3CHO" .or.         &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "MENO3" .or.          &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ETNO3" .or.          &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "MOH" .or.            &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3OH") then

                    ! Clear data
                    exportFldCAMDep(:) = 0.0_r8

                    ! Regrid deposition flux HCO to CAM
                    exportFldHcoDep(my_IS:my_IE,my_JS:my_JE) = HcoState%Spc(spcID)%Depv%Val(1:HI,1:HJ)
                    call HCO_Grid_HCO2CAM_2D(exportFldHcoDep, exportFldCAMDep)

                    ! dbg:
                    ! if(trim(HcoConfig%ModelSpc(spcID)%SpcName) == "DMS") then
                    !     dummy_2(:,:) = exportFldHco2(:,:)
                    ! endif
                endif

                ! Perform handling: Note species-specific State_CAM_* data
                ! Note handling is for surface (idx LM for CAM inverted-atm)
                if(trim(HcoConfig%ModelSpc(spcID)%SpcName) == "DMS") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmDMS(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmDMS(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                elseif( &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ACET" .or. &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3COCH3") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmACET(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmACET(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                elseif( &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ALD2" .or. &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3CHO") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmALD2(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmALD2(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                elseif( &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "MOH" .or. &
                   trim(HcoConfig%ModelSpc(spcID)%SpcName) == "CH3OH") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmMOH(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmMOH(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                elseif(trim(HcoConfig%ModelSpc(spcID)%SpcName) == "MENO3") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmMENO3(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmMENO3(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                elseif(trim(HcoConfig%ModelSpc(spcID)%SpcName) == "ETNO3") then
                    if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                        exportFldCAM(LM,:) = exportFldCAM(LM,:) - exportFldCAMDep(:) * State_CAM_chmETNO3(:) * State_CAM_DELP_DRYs(:) * G0_100
                    elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                        exportFldCAM2(:)   = exportFldCAM2(:)   - exportFldCAMDep(:) * State_CAM_chmETNO3(:) * State_CAM_DELP_DRYs(:) * G0_100
                    else
                        ASSERT_(.false.)
                    endif
                endif
            endif

            if(doExport) then
                if(HcoConfig%ModelSpc(spcID)%DimMax .eq. 3) then
                    ! Write to history on CAM mesh
                    call HCO_Export_History_CAM3D(exportName, exportFldCAM)

                    ! Write to physics buffer (pass model name)
                    call HCO_Export_Pbuf_CAM3D(HcoConfig%ModelSpc(spcID)%SpcName, spcID, exportFldCAM)
                elseif(HcoConfig%ModelSpc(spcID)%DimMax .eq. 2) then
                    ! Write to history on CAM mesh
                    call HCO_Export_History_CAM2D(exportName, exportFldCAM2)

                    ! Write to physics buffer (pass model name)
                    call HCO_Export_Pbuf_CAM2D(HcoConfig%ModelSpc(spcID)%SpcName, spcID, exportFldCAM2)
                else
                    ASSERT_(.false.)
                endif
            endif

        enddo


        !-----------------------------------------------------------------------
        ! Update aerosol number emissions with correct parameters
        ! for CAM-chem and CESM2-GC
        !
        ! TODO: Implement later. Now using HEMCO config file method (hplin, 5/7/21)
        !-----------------------------------------------------------------------
        ! call HCO_Calc_Aero_Emis ( ... )

        !-----------------------------------------------------------------------
        ! Handle special diagnostics for some extensions
        !-----------------------------------------------------------------------

        ! Reset pointers first. Always do this beforehand
        Ptr2D => NULL()

        ! Eventually save necessary deposition FLUXES from extensions.
        if(ExtState%ParaNOx > 0) then
            ! PAR_O3_DEP, PAR_HNO3_DEP
            exportName = 'HCO_PAR_O3_DEP'
            exportNameTmp = 'PAR_O3_DEP'
            call GetHcoDiagn(HcoState, ExtState, DiagnName='PARANOX_O3_DEPOSITION_FLUX', &
                             StopIfNotFound=.false., Ptr2D=Ptr2D, RC=HMRC)

            if(.not. associated(Ptr2D)) then
                write(6,*) "hplin debug err: cannot find paranox o3 dep flux"
            endif

            exportFldHco2(:,:) = 0.0_r8
            exportFldCAM2(:)   = 0.0_r8
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. associated(Ptr2D)))
            if(HMRC == HCO_SUCCESS .and. associated(Ptr2D)) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif

            if(doExport) then
                call HCO_Export_History_CAM2D(exportName, exportFldCAM2)
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            exportName = 'HCO_PAR_HNO3_DEP'
            exportNameTmp = 'PAR_HNO3_DEP'
            call GetHcoDiagn(HcoState, ExtState, DiagnName='PARANOX_HNO3_DEPOSITION_FLUX', &
                             StopIfNotFound=.false., Ptr2D=Ptr2D, RC=HMRC)

            exportFldHco2(:,:) = 0.0_r8
            exportFldCAM2(:)   = 0.0_r8
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. associated(Ptr2D)))
            if(HMRC == HCO_SUCCESS .and. associated(Ptr2D)) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif

            if(doExport) then
                call HCO_Export_History_CAM2D(exportName, exportFldCAM2)
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()
        endif

        !-----------------------------------------------------------------------
        ! Do we need to do additional exports for CESM-GC?
        !-----------------------------------------------------------------------
        
        if(chem_is('GEOS-Chem')) then
            if(masterproc) write(iulog,*) "HEMCO_CESM: starting exports to GEOS-Chem"
            do N = 0, 72
                ! Assume success
                HMRC = HCO_SUCCESS

                ! LANDTYPExx
                write(exportNameTmp, '(a,i2.2)') 'LANDTYPE', N
                exportName = 'HCO_' // trim(exportNameTmp)

                exportFldCAM2(:)   = 0.0_r8

                ! Grab the pointer if available
                call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
                doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
                if(HMRC == HCO_SUCCESS .and. FND) then
                    exportFldHco2(:,:) = Ptr2D(:,:) ! Have to promote precision
                    call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
                endif

                if(doExport) then
                    call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
                endif
                Ptr2D => NULL()

                ! XLAIxx
                write(exportNameTmp, '(a,i2.2)') 'XLAI', N
                exportName = 'HCO_' // trim(exportNameTmp)

                exportFldCAM2(:)   = 0.0_r8

                ! Grab the pointer if available
                call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
                doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
                if(HMRC == HCO_SUCCESS .and. FND) then
                    exportFldHco2(:,:) = Ptr2D(:,:) ! Have to promote precision
                    call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
                endif

                if(doExport) then
                    call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
                endif
                Ptr2D => NULL()
            enddo

            ! VMR_CH3CL
            write(exportNameTmp, '(a)') 'VMR_CH3CL'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CH3CL', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CH2CL2
            write(exportNameTmp, '(a)') 'VMR_CH2CL2'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CH2CL2', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CHCL3
            write(exportNameTmp, '(a)') 'VMR_CHCL3'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CHCL3', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CH3BR
            write(exportNameTmp, '(a)') 'VMR_CH3BR'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CH3BR', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CCL4
            write(exportNameTmp, '(a)') 'VMR_CCL4'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CCL4', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CFC113
            write(exportNameTmp, '(a)') 'VMR_CFC113'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CFC113', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CFC114
            write(exportNameTmp, '(a)') 'VMR_CFC114'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CFC114', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CFC115
            write(exportNameTmp, '(a)') 'VMR_CFC115'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CFC115', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CFC11
            write(exportNameTmp, '(a)') 'VMR_CFC11'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CFC11', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CFC12
            write(exportNameTmp, '(a)') 'VMR_CFC12'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CFC12', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_CH3CCL3
            write(exportNameTmp, '(a)') 'VMR_CH3CCL3'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_CH3CCL3', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_H1211
            write(exportNameTmp, '(a)') 'VMR_H1211'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_H1211', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_H1301
            write(exportNameTmp, '(a)') 'VMR_H1301'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_H1301', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_H2402
            write(exportNameTmp, '(a)') 'VMR_H2402'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_H2402', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_HCFC141B
            write(exportNameTmp, '(a)') 'VMR_HCFC141B'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_HCFC141B', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_HCFC142B
            write(exportNameTmp, '(a)') 'VMR_HCFC142B'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_HCFC142B', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_HCFC22
            write(exportNameTmp, '(a)') 'VMR_HCFC22'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_HCFC22', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_N2O
            write(exportNameTmp, '(a)') 'VMR_N2O'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_N2O', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_OCS
            write(exportNameTmp, '(a)') 'VMR_OCS'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_OCS', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! VMR_H2
            write(exportNameTmp, '(a)') 'VMR_H2'

            exportName = 'HCO_' // trim(exportNameTmp)
            exportFldCAM2(:) = 0.0_r8

            call HCO_EvalFld(HcoState, 'SfcVMR_H2', exportFldHco2, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif

            ! UVALBEDO
            ! Warning: Keep these exportNameTmp as it allows for reuse of
            ! code below. (hplin, 2/28/21)
            write(exportNameTmp, '(a)') 'UV_ALBEDO'

            ! Reusable code templating below.
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! SURF_SALINITY
            ! Note: the name is too long, reduce to HCO_salinity, HCO_iodide
            write(exportNameTmp, '(a)') 'salinity'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, 'surf_salinity', Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! SURF_IODIDE
            write(exportNameTmp, '(a)') 'iodide'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, 'surf_iodide', Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! OMOC_DJF
            write(exportNameTmp, '(a)') 'OMOC_DJF'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! OMOC_MAM
            write(exportNameTmp, '(a)') 'OMOC_MAM'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! OMOC_JJA
            write(exportNameTmp, '(a)') 'OMOC_JJA'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()

            ! OMOC_SON
            write(exportNameTmp, '(a)') 'OMOC_SON'
            exportName = 'HCO_' // trim(exportNameTmp)

            ! Reset data for safe export at first time step
            exportFldCAM2(:) = 0.0_r8

            call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC, FOUND=FND)
            doExport = (FIRST .or. (HMRC == HCO_SUCCESS .and. FND))
            if(HMRC == HCO_SUCCESS .and. FND) then
                exportFldHco2(:,:) = Ptr2D(:,:)
                call HCO_Grid_HCO2CAM_2D(exportFldHco2, exportFldCAM2)
            endif
            if(doExport) then
                call HCO_Export_Pbuf_CAM2D(exportNameTmp, -1, exportFldCAM2)
            endif
            Ptr2D => NULL()
            
            if(masterproc .and. nCalls < 10) write(iulog,*) "HEMCO_CESM: done with exports to GEOS-Chem"
        endif

        ! dummy_0_CAM(:,:) = iam * 1.0_r8
        
        ! dummy_1(:,:,:) = iam * 1.0_r8

        ! test data, but on the CAM grid... note vertical is inverted
        ! dummy_0_CAM sizes 16384x512. remember dummy_CAM is K, I idx
        ! my_CE: 512, pver: 32
        ! write(6,*) "hplin debug: sizes dummy", size(dummy_0_CAM, 1), size(dummy_0_CAM, 2), &
        !            size(State_CAM_ps, 1), my_CE, pver

        dummy_0_CAM(:,:) = 0.0_r8
        ! dummy_0_CAM(1,:) = State_CAM_TS
        ! dummy_0_CAM(2,:) = State_CAM_U10M
        ! dummy_0_CAM(3,:) = State_CAM_V10M
        ! dummy_0_CAM(4,:) = State_CAM_ALBD
        ! dummy_0_CAM(5,:) = State_CAM_LWI
        ! dummy_0_CAM(6,:) = State_CAM_ps
        ! dummy_0_CAM(7,:) = State_CAM_pblh
        ! dummy_0_CAM(8,:) = State_CAM_CSZA
        ! dummy_0_CAM(9,:) = State_CAM_psdry
        ! dummy_0_CAM(10,:) = State_CAM_chmO3(LM,:)
        ! dummy_0_CAM(11,:) = State_CAM_JNO2
        ! dummy_0_CAM(12,:) = State_CAM_JOH
        ! dummy_0_CAM(13,:) = scratchFldCAM2

        ! fill with some test data, but clean the data first!
        !dummy_1(:,:,:) = 0.0_r8
        dummy_1_CAM(:,:) = 0.0_r8
        ! dummy_1(:,:,1) = State_HCO_TS
        ! dummy_1(:,:,2) = State_HCO_U10M
        ! dummy_1(:,:,3) = State_HCO_V10M
        ! dummy_1(:,:,4) = State_HCO_ALBD
        ! dummy_1(:,:,5) = State_HCO_WLI
        ! dummy_1(:,:,6) = State_HCO_PSFC
        ! dummy_1(:,:,7) = State_HCO_PBLH
        ! dummy_1(:,:,8) = State_HCO_CSZA
        ! dummy_1(:,:,9) = State_HCO_AIR(:,:,1)
        ! dummy_1(:,:,10) = State_HCO_AIR(:,:,2)
        ! dummy_1(:,:,11) = Area_M2(my_IS:my_IE,my_JS:my_JE)
        ! dummy_1(:,:,12) = State_HCO_chmO3(:,:,1)
        ! dummy_1(:,:,13) = State_HCO_chmNO(:,:,1)

        ! dummy_1(:,:,14) = HcoState%Grid%BXHEIGHT_M%Val(:,:,1)
        ! dummy_1(:,:,15) = HcoState%Grid%BXHEIGHT_M%Val(:,:,2)
        ! dummy_1(:,:,16) = State_HCO_F_OF_PBL(:,:,1)
        ! dummy_1(:,:,17) = State_HCO_F_OF_PBL(:,:,2)
        ! dummy_1(:,:,18) = dummy_2(:,:)

        ! Regrid to CAM physics mesh!
        !call HCO_Grid_HCO2CAM_3D(dummy_1, dummy_1_CAM)

        ! Write to history on CAM mesh
        !call HCO_Export_History_CAM3D("DIAG_HCO_TEST", dummy_1_CAM)
        !call HCO_Export_History_CAM3D("DIAG_CAM_TEST", dummy_0_CAM)

        if(masterproc .and. nCalls < 10) then
            write(iulog,*) "HEMCO_CESM: Exports completed for this timestep!"
        endif

        !-----------------------------------------------------------------------
        ! Finished!
        !-----------------------------------------------------------------------
        ! Update last execution time
        last_HCO_day    = now_day
        last_HCO_second = now_s

        IF ( FIRST ) FIRST = .False.

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

        RC = ESMF_SUCCESS

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

        RC = ESMF_SUCCESS

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
