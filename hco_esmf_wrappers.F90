!------------------------------------------------------------------------------
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_esmf_wrappers
!
! !DESCRIPTION: Module HCO\_ESMF\_WRAPPERS defines quick wrapper functions for
!  error checking and common ESMF operations, similar to "MAPL_Generic.h" in MAPL.
!\\
!\\
! !INTERFACE:
!
module hco_esmf_wrappers
!
! !USES:
!
    ! ESMF Types
    use ESMF,                     only: ESMF_SUCCESS, ESMF_FAILURE

    ! MPI status in CESM
    use cam_abortutils,           only: endrun      ! fatal terminator
    use spmd_utils,               only: iam, masterproc
    use cam_logfile,              only: iulog

    implicit none
    private

    public :: HCO_ESMF_VRFY
    public :: HCO_ESMF_ASRT

contains

    ! Note that the VERIFY and ASSERT functions here are inverse of MAPL and
    ! return 1 on success, 0 on failure.
    logical function HCO_ESMF_VRFY(A, subname, line, ECRC)
        integer, intent(in)       :: A
        character*(*), intent(in) :: subname
        integer, intent(in)       :: line
        integer, optional, intent(out) :: ECRC

        character(len=512)        :: errmsg

        HCO_ESMF_VRFY = A==ESMF_SUCCESS
        if(.not. HCO_ESMF_VRFY) then
            if(present(ECRC)) then
                print'(A40,I10)', subname, line
                ECRC = A
            endif

            write(errmsg,*) 'ABORT in HCO_ESMF_WRAPPERS. VRFY error in ', subname, ' line ', line
            if(masterproc) then
                write(iulog,*) errmsg
            endif
            call endrun(errmsg)
        endif
    end function HCO_ESMF_VRFY

    logical function HCO_ESMF_ASRT(A, subname, line, ECRC)
        logical, intent(in)       :: A
        character*(*), intent(in) :: subname
        integer, intent(in)       :: line
        integer, optional, intent(out) :: ECRC

        character(len=512)        :: errmsg

        HCO_ESMF_ASRT = A
        if(.not. HCO_ESMF_ASRT) then
            if(present(ECRC)) then
                print'(A40,I10)', subname, line
                ECRC = ESMF_FAILURE
            endif

            write(errmsg,*) 'ABORT in HCO_ESMF_WRAPPERS. ASRT error in ', subname, ' line ', line
            if(masterproc) then
                write(iulog,*) errmsg
            endif
            call endrun(errmsg)
        endif
    end function HCO_ESMF_ASRT
!EOC
end module hco_esmf_wrappers
