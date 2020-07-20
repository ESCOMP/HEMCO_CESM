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

    ! Species information (CESM-GC) (temporary kludge)
#if defined( EXTERNAL_FORCING )
    use chem_mods,                only: nTracers, tracerNames
    use chem_mods,                only: MWRatio
#endif

    ! Grid
    use ppgrid,                   only: pcols, pver ! Cols, verts
    use ppgrid,                   only: begchunk, endchunk ! Chunk idxs

    ! Time
    use time_manager,             only: get_curr_time, get_prev_time, get_curr_date

    ! ESMF types
    use ESMF,                     only: ESMF_State, ESMF_Clock, ESMF_GridComp
    use ESMF,                     only: ESMF_KIND_R8, ESMF_KIND_I4, ESMF_SUCCESS

    ! HEMCO types
    use HCO_Error_Mod,            only: hp          ! HEMCO precision
    use HCO_Error_Mod,            only: sp          ! HEMCO single precision used for Ptrs
    use HCO_Error_Mod,            only: HCO_SUCCESS, HCO_FAIL, HCO_VERSION
    use HCO_State_Mod,            only: HCO_State
    use HCOX_State_Mod,           only: Ext_State   ! Note: Extensions unsupported
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

    private :: HCOI_Allocate_All
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

    ! Last execution times for the HEMCO component. We are assuming that time
    ! flows unidirectionally (and forwards, for now). (hplin, 3/30/20)
    integer                          :: last_HCO_day, last_HCO_second

    ! Are we exporting CHEM_INPUTS to CESM2-GC?
    logical                          :: HCO_CESM2GCInputs

    ! Meteorological fields used by HEMCO to be regridded to the HEMCO grid (hplin, 3/31/20)
    ! We have to store the fields because the regridding can only take place within the GridComp.
    ! Fields are allocated after the internal grid is initialized (so my_* are avail)

    ! Currently in hemco_interface.F90 for ease of devel, move to a container module later

    ! On the CAM grid (state%psetcols, pver) (LM, my_CE)
    ! Arrays are flipped in order (k, i) for the regridder
    real(r8), pointer                :: State_CAM_t(:,:)
    real(r8), pointer                :: State_CAM_ps(:)
    real(r8), pointer                :: State_CAM_pblh(:)

    ! On the HEMCO grid (my_IM, my_JM, LM) or possibly LM+1
    ! HEMCO grid are set as POINTERs so it satisfies HEMCO which wants to point
    real(r8), pointer                :: State_HCO_TK  (:,:,:)
    real(r8), pointer                :: State_HCO_PSFC(:,:)   ! Wet?
    real(r8), pointer                :: State_HCO_PBLH(:,:)   ! PBLH [m]

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

        ! Save this to the module information
        HcoConfigFile = hemco_config_file
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
        use cam_history,      only: addfld, add_default

        ! HEMCO Initialization
        use HCO_Config_Mod,   only: Config_ReadFile, ConfigInit
        use HCO_Driver_Mod,   only: HCO_Init
        use HCO_Error_Mod,    only: HCO_LOGFILE_OPEN
        use HCO_LogFile_Mod,  only: HCO_Spec2Log
        use HCO_State_Mod,    only: HcoState_Init
        use HCO_Types_Mod,    only: ConfigObj
        use HCO_Types_Mod,    only: ListCont
        use HCO_VertGrid_Mod, only: HCO_VertGrid_Define

        ! HEMCO extensions are unsupported for now.
!
! !REMARKS:
!  HEMCO extensions are unsupported in the preliminary version, due to
!  complications associated with met fields and stuff.
!  A "state conversion" to convert CAM met fields to GEOSFP format will
!  need to be coordinated with fritzt later down the road (hplin, 3/27/20)
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

        ! HEMCO types
        type(ListCont), pointer      :: TmpLct

        !-----------------------------------------------------------------------

        if(masterproc) then
            write(iulog,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
            write(iulog,*) "HEMCO: Harmonized Emissions Component"
            write(iulog,*) "HEMCO_CAM interface version 0.1"
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
        call addfld("PETID_CAM_TEST", (/'lev'/), 'I', '1',          &
                    'HEMCO Debug, PETID written on CAM',             &
                    gridname="physgrid")
        call add_default("PETID_CAM_TEST", 2, 'I')      ! Make this field always ON

        call addfld("PETID_HCO_TEST", (/'lev'/), 'I', '1',          &
                    'HEMCO Debug, PETID written on HCO',             &
                    gridname="physgrid")
        call add_default("PETID_HCO_TEST", 2, 'I')      ! Make this field always ON

        !-----------------------------------------------------------------------
        ! Initialize the HEMCO configuration object...
        !-----------------------------------------------------------------------
        if(masterproc) write(iulog,*) "> Initializing HCO configuration object"

        ! Very ugly hack to recognize CESM-GC: Change TBD
#if !defined( EXTERNAL_FORCING )
        ! We are using gas_pcnst here, which is # of "gas phase" species.
        ! This might be changed down the line but we are reading chem_mods
        ! for now.
        nHcoSpc             = gas_pcnst          ! # of hco species? using gas
#else
        ! Working with CESM-GC: use nTracers. To coordinate with tmmf
        nHcoSpc             = nTracers
#endif

        call ConfigInit(HcoConfig, HMRC, nModelSpecies=nHcoSpc)
        ASSERT_(HMRC==HCO_SUCCESS)

        HcoConfig%amIRoot   = masterproc
        ! HcoConfig%amIRoot   = .true. ! for debug only so verbosity is higher
        HcoConfig%MetField  = 'MERRA2'
        HcoConfig%GridRes   = ''

        !-----------------------------------------------------------------------
        ! Retrieve the species list and register exports
        !-----------------------------------------------------------------------
        ! Below we directly use nHcoSpc for now... it may be wrong in the future though
        ! (hplin, 3/29/20)
        HcoConfig%nModelSpc = nHcoSpc
        HcoConfig%nModelAdv = nHcoSpc            ! # of adv spc?

        if(masterproc) write(iulog,*) "> Initializing HCO species list!"

        do N = 1, nHcoSpc
            HcoConfig%ModelSpc(N)%ModID   = N ! model id

            ! Ugly kludge to recognize CESM-GC (using EXTERNAL_FORCING) (hplin, 5/16/20)
#if !defined( EXTERNAL_FORCING )
            HcoConfig%ModelSpc(N)%SpcName = trim(solsym(N))
#else
            HcoConfig%ModelSpc(N)%SpcName = trim(tracerNames(N))
#endif

            !----------------------------------------------
            ! Register export properties.
            !----------------------------------------------
            ! History output (this will be moved to hco_cam_exports soon hopefully)

            ! TODO (hplin): Writes to tape 2 by default, add namelist option to force it
            ! to go in a different tape? Or abide to CAM conventions...
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(N)%SpcName)
            exportDesc = "HEMCO 3-D Emissions Species " // trim(HcoConfig%ModelSpc(N)%SpcName)
            call addfld(exportName, (/'lev'/), 'I', 'kg/m2/s',          &
                        trim(exportDesc),                               &
                        gridname='physgrid')
            call add_default(exportName, 2, 'I') ! On by default

            if(masterproc) write(iulog,*) "Exported exportName " // trim(exportName) // " to history"

            ! Physics buffer
            ! Note that _AddField will prepend HCO_, so do not add it here
            call HCO_Export_Pbuf_AddField(HcoConfig%ModelSpc(N)%SpcName, 3, hcoID=N)
        enddo

        !-----------------------------------------------------------------------
        ! Read HEMCO configuration file from HcoConfigFile (location in CAM namelist)
        !-----------------------------------------------------------------------
        if(masterproc) write(iulog,*) "> Reading HEMCO configuration file..."

        ! FIXME: Not implementing "Dry-run" functionality in HEMCO_CESM. (hplin, 3/27/20)
        ! Phase: 0 = all, 1 = sett and switches only, 2 = fields only
        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 1, HMRC, IsDryRun=.false.)
        ASSERT_(HMRC==HCO_SUCCESS)

        ! Open the log file
        if(masterproc) then
            call HCO_LOGFILE_OPEN(HcoConfig%Err, RC=HMRC)
            ASSERT_(HMRC==HCO_SUCCESS)
        endif

        call Config_ReadFile(HcoConfig%amIRoot, HcoConfig, HcoConfigFile, 2, HMRC, IsDryRun=.false.)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "> Read HEMCO configuration file OK!"

        !-----------------------------------------------------------------------
        ! Initialize the HEMCO state object
        !-----------------------------------------------------------------------
        call HcoState_Init(HcoState, HcoConfig, nHcoSpc, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "> Initialize HEMCO state obj OK!"

        ! Emissions, chemistry and dynamics timestep [s]
        ! Assume 0.5h until given actual time in HCO_GC_Run!
        HcoState%TS_EMIS = 1800.0
        HcoState%TS_CHEM = 1800.0
        HcoState%TS_DYN  = 1800.0

        ! Not a MAPL simulation. isESMF is deceiving.
        HcoState%Options%isESMF = .false.

        ! Deposition length scale. Used for computing dry deposition frequencies
        ! over the entire PBL or the first model layer. Hardcoded for now,
        ! should load Input_Opt%PBL_DRYDEP from GEOS-Chem-CESM (hplin, 3/29/20)
        ! !FIXME
        HcoState%Options%PBL_DRYDEP = .false.

        ! Don't support DryRun option (for now)
        HcoState%Options%IsDryRun = .false.

        if(masterproc) write(iulog,*) "> Set basic HEMCO state obj OK!"

        !-----------------------------------------------------------------------
        ! Register HEMCO species information (HEMCO state object)
        !-----------------------------------------------------------------------
        do N = 1, nHcoSpc
            HcoState%Spc(N)%ModID         = N               ! model id

            ! Ugly kludge to recognize CESM-GC (using EXTERNAL_FORCING) (hplin, 5/16/20)
#if !defined( EXTERNAL_FORCING )
            HcoState%Spc(N)%SpcName       = trim(solsym(N)) ! species name

            HcoState%Spc(N)%MW_g          = adv_mass(N)     ! mol. weight [g/mol]

            ! Emitted molecules per molecules of species [1]
            ! Most spc. 1.0, for the species in the list above, will be # of moles carbon
            ! per mole species.

            ! FIXME: To resolve special case
            HcoState%Spc(N)%MolecRatio    = 1.0_hp
#else
            ! CESM-GC has all necessary data, maybe,
            ! FIXME: to coordinate EmMW_g with TMMF (hplin, 5/16/20)
            HcoState%Spc(N)%SpcName       = trim(tracerNames(N)) ! species name

            HcoState%Spc(N)%MW_g          = adv_Mass(N)     ! mol. weight [g/mol]
            HcoState%Spc(N)%MolecRatio    = MWRatio(N)      ! Emitted mol per mol of spc [1]
#endif

            ! WARNING: This is the EMITTED molecular weight of species.
            ! Some hydrocarbons (e.g. ISOP) are emitted as equiv. no. of C atoms
            ! e.g. ISOP EmMW_g is 12.0.
            !
            ! Many species in GEOS-Chem behave like this, i.e.
            ! ACET, ALD2, ALK4, BCPI, BCPO, BENZ, C2H6, C3H8, EOH, ISOP, MEK
            ! MOPI, MOPO, NAP,  OCPI, OCPO, OPOA1, OPOA2, OPOG1, OPOG2,
            ! POA1, POA2, POG1, POG2, PRPE, TOLU, XYLE, APMBCBIN**, APMOCBIN**
            !
            ! WE PREFORM A MANUAL ADJUSTMENT HERE, BECAUSE WE DO NOT KNOW ABOUT
            ! THE PRESENCE OF GEOS-CHEM. WE PROBABLY HAVE TO REVISIT THIS IN THE
            ! FUTURE.
            if(trim(solsym(N)) .eq. 'ACET'  .or.  trim(solsym(N)) .eq. 'ALD2'  .or.  &
               trim(solsym(N)) .eq. 'ALK4'  .or.  trim(solsym(N)) .eq. 'BCPI'  .or.  &
               trim(solsym(N)) .eq. 'BCPO'  .or.  trim(solsym(N)) .eq. 'BENZ'  .or.  &
               trim(solsym(N)) .eq. 'C2H6'  .or.  trim(solsym(N)) .eq. 'C3H8'  .or.  &
               trim(solsym(N)) .eq. 'EOH'   .or.  trim(solsym(N)) .eq. 'ISOP'  .or.  &
               trim(solsym(N)) .eq. 'MEK'   .or.  trim(solsym(N)) .eq. 'MOPI'  .or.  &
               trim(solsym(N)) .eq. 'MOPO'  .or.  trim(solsym(N)) .eq. 'NAP'   .or.  &
               trim(solsym(N)) .eq. 'OCPI'  .or.  trim(solsym(N)) .eq. 'OCPO'  .or.  &
               trim(solsym(N)) .eq. 'OPOA1' .or.  trim(solsym(N)) .eq. 'OPOA2' .or.  &
               trim(solsym(N)) .eq. 'OPOG1' .or.  trim(solsym(N)) .eq. 'OPOG2' .or.  &
               trim(solsym(N)) .eq. 'POA1'  .or.  trim(solsym(N)) .eq. 'POA2'  .or.  &
               trim(solsym(N)) .eq. 'POG1'  .or.  trim(solsym(N)) .eq. 'POG2'  .or.  &
               trim(solsym(N)) .eq. 'PRPE'  .or.  trim(solsym(N)) .eq. 'TOLU'  .or.  &
               trim(solsym(N)) .eq. 'XYLE') then
                HcoState%Spc(N)%EmMW_g    = 12.0_hp
            else
                HcoState%Spc(N)%EmMW_g    = adv_mass(N)     ! emitted mol. weight [g/mol]
            endif

            ! !!! We don't set Henry's law coefficients in HEMCO_CESM !!!
            ! they are mostly used in HCOX_SeaFlux_Mod, but HCOX are unsupported (for now)
            ! (hplin, 3/29/20)
            ! HcoState%Spc(N)%HenryK0 ! [M/atm]
            ! HcoState%Spc(N)%HenryCR ! [K]
            ! HcoState%Spc(N)%HenryPKA ! [1]

            ! Write to log too
            if(masterproc) then
                write(iulog,*) ">> Spc", N, " = ", solsym(N), "MW_g", adv_mass(N), "EmMW_g", HcoState%Spc(N)%EmMW_g
                call HCO_Spec2Log(HcoState, N)
            endif
        enddo

        if(masterproc) write(iulog,*) "> Set HEMCO species info OK!"

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

        if(masterproc) write(iulog,*) "> Set HEMCO PET-local grid info OK!"

        !-----------------------------------------------------------------------
        ! Initialize HEMCO!
        !-----------------------------------------------------------------------
        call HCO_Init(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "> HEMCO initialized successfully!"

        !-----------------------------------------------------------------------
        ! Additional exports: Verify if we need to add additional exports
        ! for integration with CESM-GC. (hplin, 4/15/20)
        !-----------------------------------------------------------------------
        TmpLct => HcoState%ReadLists%Once
        ! Kludge: If LANDTYPE00 is available in one-time list, all exports for CAM
        ! will be initialized in the pbuf
        HCO_CESM2GCInputs = .false.
        do while(associated(TmpLct))
            if(associated(TmpLct%Dct)) then
                if(trim(TmpLct%Dct%cName) .eq. 'LANDTYPE00') then
                    HCO_CESM2GCInputs = .true.
                    exit
                endif
                TmpLct => TmpLct%NextCont
            endif
        enddo
        TmpLct => NULL()

        ! Do additional exports!
        if(HCO_CESM2GCInputs) then
            do N = 0, 72
                ! LANDTYPExx
                write(exportNameTmp, '(a,i2.2)') 'LANDTYPE', N
                exportName = 'HCO_' // trim(exportNameTmp)
                exportDesc = "HEMCO Chemistry Input Name " // trim(exportNameTmp)

                ! FIXME (hplin): Exporting as 3-D; third dimension unused, change later...
                ! Too lazy to write an Export_CAM2D
                call addfld(exportName, (/'lev'/), 'I', '1',                &
                            trim(exportDesc),                               &
                            gridname='physgrid')
                call add_default(exportName, 2, 'I') ! On by default

                ! Also pbuf
                call HCO_Export_Pbuf_AddField(exportNameTmp, 3)
                
                if(masterproc) write(iulog,*) "Exported exportName " // trim(exportName) // " to history"

                ! XLAIxx
                write(exportNameTmp, '(a,i2.2)') 'XLAI', N
                exportName = 'HCO_' // trim(exportNameTmp)
                exportDesc = "HEMCO Chemistry Input Name " // trim(exportNameTmp)

                ! FIXME (hplin): Exporting as 3-D; third dimension unused, change later...
                ! Too lazy to write an Export_CAM2D
                call addfld(exportName, (/'lev'/), 'I', '1',                &
                            trim(exportDesc),                               &
                            gridname='physgrid')
                call add_default(exportName, 2, 'I') ! On by default

                ! Also pbuf
                call HCO_Export_Pbuf_AddField(exportNameTmp, 3)
                
                if(masterproc) write(iulog,*) "Exported exportName " // trim(exportName) // " to history"
            enddo

            if(masterproc) then
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

    end subroutine HCOI_Allocate_All
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

        integer                      :: I, J, K, lchnk            ! Loop idx
        integer                      :: ncol

        type(physics_buffer_desc), pointer :: pbuf_chnk(:)        ! slice of pbuf
        ! pbuf indices:
        integer                      :: index_pblh

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
            ! Pass the pbuf to the hco_cam_exports component so she has it...
            hco_pbuf2d => pbuf2d

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

                    ! 2-D Fields
                    State_CAM_ps(I) = phys_state(lchnk)%ps(J)

                enddo
            enddo

            ! Run the gridded component.
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
        use cam_logfile,            only: iulog
        use spmd_utils,             only: masterproc, iam

        ! HEMCO
        use HCO_Interface_Common,   only: GetHcoVal, GetHcoDiagn
        use HCO_Clock_Mod,          only: HcoClock_Set, HcoClock_Get
        use HCO_Clock_Mod,          only: HcoClock_EmissionsDone
        use HCO_Diagn_Mod,          only: HcoDiagn_AutoUpdate
        use HCO_Driver_Mod,         only: HCO_Run
        use HCO_EmisList_Mod,       only: Hco_GetPtr
        use HCO_FluxArr_Mod,        only: HCO_FluxArrReset
        use HCO_GeoTools_Mod,       only: HCO_CalcVertGrid, HCO_SetPBLm

        use HCO_State_Mod,          only: HCO_GetHcoId


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

        ! For debug dummies
        real(ESMF_KIND_R8)                    :: dummy_0_CAM(1:LM, 1:my_CE)
        real(ESMF_KIND_R8)                    :: dummy_1(my_IS:my_IE, my_JS:my_JE, 1:LM)
        real(ESMF_KIND_R8)                    :: dummy_1_CAM(1:LM, 1:my_CE)

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

        ! Assume success
        RC = ESMF_SUCCESS
        HMRC = HCO_SUCCESS

        !-----------------------------------------------------------------------
        ! Update regridding file handles as necessary
        !-----------------------------------------------------------------------
        call HCO_Grid_UpdateRegrid(RC=RC)
        ASSERT_(RC==ESMF_SUCCESS)

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Reload (if necessary) of HEMCO Regrid descriptors"
        endif

        !-----------------------------------------------------------------------
        ! Get time properties
        !-----------------------------------------------------------------------
        call get_prev_time(prev_day, prev_s)
        call get_curr_time(now_day, now_s)
        call get_curr_date(year, month, day, tod)

        if(masterproc) then
            write(iulog,*) "hco year,month,day,tod", year, month, day, tod
            write(iulog,*) "hco prev_day, prev_s", prev_day, prev_s
            write(iulog,*) "hco now_day, now_s", now_day, now_s
        endif
        ! 2005 1 1 1800 | 0 0 | 0 1800
        ! 2005 1 1 3600 | 0 1800 | 0 3600
        ! 2005 1 1 5400 | 0 3600 | 0 5400
        ! ...

        ! Check if we have run HEMCO for this time step already. If yes can exit
        if(last_HCO_day * 86400.0 + last_HCO_second .ge. now_day * 86400.0 + now_s) then
            if(masterproc) then
                write(iulog,*) "HEMCO_CAM: !! HEMCO already ran for this time, check timestep mgr", now_day, now_s, last_HCO_day, last_HCO_second
            endif

            return
        endif

        ! Compute timestep
        if(HcoState%TS_CHEM .ne. ((now_day - prev_day) * 86400.0 + now_s - prev_s)) then
            HcoState%TS_EMIS = (now_day - prev_day) * 86400.0 + now_s - prev_s
            HcoState%TS_CHEM = (now_day - prev_day) * 86400.0 + now_s - prev_s
            HcoState%TS_DYN  = (now_day - prev_day) * 86400.0 + now_s - prev_s
            if(masterproc) write(iulog,*) "HEMCO_CAM: Updated HEMCO timestep to ", HcoState%TS_CHEM
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

        if(masterproc) then
            write(6,*) "HEMCO_CAM: Updating HEMCO clock to set", year, month, day, hour, minute, second
            write(6,*) "HEMCO_CAM: Internally HEMCO is at ", HcoState%Clock%SimHour, HcoState%Clock%SimMin, HcoState%Clock%SimSec, HcoState%Clock%nSteps
        endif

        call HCOClock_Set(HcoState, year, month, day,  &
                          hour, minute, second, IsEmisTime=.true., RC=HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)
        
        !-----------------------------------------------------------------------
        ! Regrid necessary meteorological quantities
        !-----------------------------------------------------------------------
        call HCO_Grid_CAM2HCO_2D(State_CAM_ps,     State_HCO_PSFC  )
        call HCO_Grid_CAM2HCO_2D(State_CAM_pblh,   State_HCO_PBLH  )
        call HCO_Grid_CAM2HCO_3D(State_CAM_t,      State_HCO_TK    )

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Finished regridding CAM met fields to HEMCO"

            ! As a test... maybe we also need to flip in the vertical
            ! write(iulog,*) State_HCO_TK(1,1,:)
            ! write(iulog,*) "PSFC(1:2,:)"
            ! write(iulog,*) State_HCO_PSFC(1:2,:) 

            !write(iulog,*) "cam state%ps dump"
            !write(iulog,*) State_CAM_ps

            ! TK: 288 283 277 271 266 261 ... 250 251 252
            ! Seems like the vertical is OK for now
        endif
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

        !-----------------------------------------------------------------------
        ! Continue setting up HEMCO
        !-----------------------------------------------------------------------

        ! Reset all emission and deposition values.
        call HCO_FluxArrReset(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

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
            write(iulog,*) "HEMCO_CAM: HCO_Run Phase 1 has failed!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "A critical component in HEMCO failed to run."
            write(iulog,*) "This may be due to misconfiguration, or a bug."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "HEMCO_CAM: HCO_Run Phase 1"

        call HCO_Run( HcoState, 2, HMRC, IsEndStep=.false. )
        if(masterproc .and. HMRC /= HCO_SUCCESS) then
            write(iulog,*) "******************************************"
            write(iulog,*) "HEMCO_CAM: HCO_Run Phase 2 has failed!    "
            write(iulog,*) "THIS ERROR ORIGINATED WITHIN HEMCO!       "
            write(iulog,*) "A critical component in HEMCO failed to run."
            write(iulog,*) "This may be due to misconfiguration, or a bug."
            write(iulog,*) "Please refer to the HEMCO.log log file in your"
            write(iulog,*) "case run directory or as configured in HEMCO_Config.rc"
            write(iulog,*) "for more information."
            write(iulog,*) "******************************************"
        endif
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "HEMCO_CAM: HCO_Run Phase 2"

        !-----------------------------------------------------------------------
        ! Update "autofill" diagnostics.
        ! Update all 'AutoFill' diagnostics. This makes sure that all
        ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
        ! AutoFill flag is specified when creating a diagnostics container
        ! (Diagn_Create).
        !-----------------------------------------------------------------------
        call HcoDiagn_AutoUpdate(HcoState, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "HEMCO_CAM: HcoDiagn_AutoUpdate"

        !-----------------------------------------------------------------------
        ! Tell HEMCO we are done for this timestep...
        !-----------------------------------------------------------------------
        call HcoClock_EmissionsDone(HcoState%Clock, HMRC)
        ASSERT_(HMRC==HCO_SUCCESS)

        if(masterproc) write(iulog,*) "HEMCO_CAM: HcoClock_EmissionsDone"

        !-----------------------------------------------------------------------
        ! Do some testing and write emissions to the tape
        !-----------------------------------------------------------------------

        ! For each species...
        do spcID = 1, HcoConfig%nModelSpc
            ! Zero out quantities first
            exportFldHco(:,:,:) = 0.0_r8
            exportFldCAM(:,:)   = 0.0_r8

            ! Build history / pbuf field name (HCO_NO, HCO_CO, etc.)
            exportName = 'HCO_' // trim(HcoConfig%ModelSpc(spcID)%SpcName)
            ! if(masterproc) write(iulog,*) "HEMCO_CAM: Begin exporting " // trim(exportName)

            ! Get HEMCO emissions flux [kg/m2/s].
            do K = 1, LM
                HL = K                ! map to HEMCO index
            do J = my_JS, my_JE
                HJ = J - my_JS + 1
            do I = my_IS, my_IE
                HI = I - my_IS + 1

                ! Maybe no need to call GetHcoVal and just directly grab the pointer
                ! instead?
                call GetHcoVal(HcoState, ExtState, spcID, HI, HJ, HL, FND, emis=TMP)
                if(FND) exportFldHco(I,J,K) = TMP
            enddo
            enddo
            enddo
            ! if(masterproc) write(iulog,*) "HEMCO_CAM: Retrieved from HCO " // trim(exportName)

            ! Regrid exportFldHco to CAM grid...
            call HCO_Grid_HCO2CAM_3D(exportFldHco, exportFldCAM)
            ! if(masterproc) write(iulog,*) "HEMCO_CAM: Regridded " // trim(exportName)

            ! Write to history on CAM mesh
            call HCO_Export_History_CAM3D(exportName, exportFldCAM)
            ! if(masterproc) write(iulog,*) "HEMCO_CAM: Exported to history " // trim(exportName)

            ! Write to physics buffer (pass model name)
            call HCO_Export_Pbuf_CAM3D(HcoConfig%ModelSpc(spcID)%SpcName, spcID, exportFldCAM)
        enddo

        ! Do we need to do additional exports for CESM-GC?
        if(HCO_CESM2GCInputs) then
            do N = 0, 72
                ! Assume success
                HMRC = HCO_SUCCESS

                ! LANDTYPExx
                write(exportNameTmp, '(a,i2.2)') 'LANDTYPE', N
                exportName = 'HCO_' // trim(exportNameTmp)

                ! FIXME (hplin): Exporting as 3-D; third dimension unused, change later...
                ! Too lazy to write an Export_CAM2D
                exportFldHco(:,:,:) = 0.0_r8
                exportFldCAM(:,:)   = 0.0_r8

                do J = my_JS, my_JE
                    HJ = J - my_JS + 1
                do I = my_IS, my_IE
                    HI = I - my_IS + 1

                    ! Grab the pointer if available
                    call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC)
                    if(HMRC == HCO_SUCCESS) exportFldHco(:,:,1) = Ptr2D ! Copy data in
                enddo
                enddo

                call HCO_Grid_HCO2CAM_3D(exportFldHco, exportFldCAM)
                call HCO_Export_History_CAM3D(exportName, exportFldCAM)
                call HCO_Export_Pbuf_CAM3D(exportNameTmp, -1, exportFldCAM)

                ! XLAIxx
                write(exportNameTmp, '(a,i2.2)') 'XLAI', N
                exportName = 'HCO_' // trim(exportNameTmp)

                ! FIXME (hplin): Exporting as 3-D; third dimension unused, change later...
                ! Too lazy to write an Export_CAM2D
                exportFldHco(:,:,:) = 0.0_r8
                exportFldCAM(:,:)   = 0.0_r8

                do J = my_JS, my_JE
                    HJ = J - my_JS + 1
                do I = my_IS, my_IE
                    HI = I - my_IS + 1

                    ! Grab the pointer if available
                    call HCO_GetPtr(HcoState, exportNameTmp, Ptr2D, HMRC)
                    if(HMRC == HCO_SUCCESS) exportFldHco(:,:,1) = Ptr2D ! Copy data in
                enddo
                enddo

                call HCO_Grid_HCO2CAM_3D(exportFldHco, exportFldCAM)
                call HCO_Export_History_CAM3D(exportName, exportFldCAM)
                call HCO_Export_Pbuf_CAM3D(exportNameTmp, -1, exportFldCAM)
            enddo
        endif

        dummy_0_CAM(:,:) = iam * 1.0_r8
        dummy_1(:,:,:) = iam * 1.0_r8

        ! Regrid to CAM physics mesh!
        call HCO_Grid_HCO2CAM_3D(dummy_1, dummy_1_CAM)

        ! Write to history on CAM mesh
        call HCO_Export_History_CAM3D("PETID_HCO_TEST", dummy_1_CAM)
        call HCO_Export_History_CAM3D("PETID_CAM_TEST", dummy_0_CAM)

        if(masterproc) then
            write(iulog,*) "HEMCO_CAM: Exports completed for this timestep!"
        endif

        !-----------------------------------------------------------------------
        ! Finished!
        !-----------------------------------------------------------------------
        ! Update last execution time
        last_HCO_day    = now_day
        last_HCO_second = now_s

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