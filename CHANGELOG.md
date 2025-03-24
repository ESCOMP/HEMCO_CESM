# Changelog

This file documents all notable changes to the HEMCO-CESM interface since late 2022.

## [2.1.1] - TBD
### Fixed
 - Enable b4b reproducibility in all CAM <-> HCO regridding by passing srcTermProcessing=1 to ESMF_RegridFieldStore

## [2.1.0] - 2024-12-16
### Added
- Added fleximod capability to .gitmodules
- Added print target LUN for atm.log to arguments passed to HEMCO ConfigInit to direct HEMCO root prints to atm.log rathe rather than cesm.log

### Changed
- Added new default configuration file for HEMCO which reads/regrids 3D emissions hourly
- Updated debug printout and message for where to find error messages
- Updated HEMCO submodule from 3.9.0 to 3.10.1

### Fixed
- Fixed HEMCO clock to be time at start of timestep rather than end
- Fixed formatting for Ap/Bp printout
- Fixed condition for first timestep to ensure HEMCO is not run at timestep 0

## [2.0.0] - 2024-07-31
### Added
- Added git submodule capability with HEMCO as the only submodule.

### Changed
- Update README
- Run HEMCO starting in first step
- Updates to debug printout
- **Based off upstream HEMCO 3.9.0**. Updated HEMCO submodule from version 3.8.1 to 3.9.0. Refer to HEMCO release notes.
- First timestep now runs HEMCO as there is no longer issue with timestep measurement from prior versions.

### Removed
- Removed support for manage_externals.
- Removed obsolete timestep measurement as timestep information is already obtained from time_manager reliably.
- Removed support for manage_externals.
- Removed obsolete timestep measurement as timestep information is already obtained from time_manager reliably.

## [1.3.0] - 2024-04-29
### Added
- Performance timers using CIME infrastructure for initialization and run of major HEMCO_CESM and HEMCO components.

### Changed
- **Based off upstream HEMCO 3.8.1**. Refer to the HEMCO release notes.

## [1.2.2] - 2024-03-01
### Changed
- **Based off upstream HEMCO 3.8.0**. Release 3.7.x was skipped due to regression in vertical regridding capabilities. Refer to the HEMCO release notes for changes in 3.8.0 release.

## [1.2.1] - 2023-09-15
### Fixed
- **Based off upstream HEMCO 3.6.3-cesm** - release for CESM to fix `nvhpc` compiler compatibility issues.

## [1.2.0] - 2023-06-15
### Added
- Support the `hemco_data_root` property to override `$ROOT` in HEMCO configuration file in the CESM environment.
- Support the `hemco_diagn_file` property to override `DiagnFile` in HEMCO configuration file in the CESM environment.

### Changed
- Code now uses explicit `only` includes for `hemco_interface.F90`.
- Improved error handling for `HCO_Init` call to point user to correct error messages.
- Improved debug outputs throughout.

## [1.1.6] - 2023-06-08
### Fixed
- Bit-for-bit match in identical runs now ensured through `smm_pipelinedep` fixed to 16, similarly to `ionosphere/waccmx/edyn_esmf.F90`, instead of automatic tuning.

## [1.1.5] - 2023-06-01
### Changed
- **Based off upstream HEMCO 3.6.2-cesm** - out of band release for CESM NAG compiler compatibility.
- Fixes for NAG compiler.

## [1.1.4] - 2023-04-12
### Changed
- HEMCO-CESM now runs only on CAM phase 2 and clears pbuf at that phase. Phase 1 is now unused. Interface remains unchanged and is backwards compatible.

## [1.1.3] - 2023-03-20
### Fixed
- Fixed wrong definition of `ExtState%TSKIN`. Now correctly uses SST from `cam_in%sst`.

## [1.1.2] - 2023-03-08
### Added
- Support the `hemco_emission_year` namelist option to force emission cycling year (equivalent to `HEMCO_Config.rc`'s `Emission year:` option) from CESM compsets

### Changed
- Slightly adjusted and reduced debug prints

## [1.1.1] - 2023-03-02
### Added
- **Based off upstream HEMCO 3.6.2**
- CAM-chem only: Initialization of `HCO_IN_JNO2` and `HCO_IN_JOH` fields is now done in initialization and state conversion module.

## [1.1.0] - 2023-02-09
### Added
- **Based off upstream HEMCO 3.6.0**
- State conversion module now supports `FRLAND` and stubs for `FRLANDIC` and `FRLAKE`
- Now allocates 2-D physics buffer for species not in `mo_sim_dat::extfrc_lst`

### Fixed
- Fixed stack corruption in `HCOI_Allocate_All`

### Removed
- Removed support for `LWI` from upstream

## [1.0.2] - 2023-01-11
### Fixed
- Return `ESMF_SUCCESS` for dummy ESMF gridded component routines

## [1.0.1] - 2022-11-09
### Fixed
- Initialize namelist variables before reading namelist.

## [1.0.0] (HEMCO 3.5.1) - 2022-11-03
### Added
- **Based off upstream HEMCO 3.5.1**

## [1.0.0] (HEMCO 3.5.0) - 2022-10-31
### Added
- **Based off upstream HEMCO 3.5.0**
- Support for spectral-element (including regional-refinement, SE-RR) grids
- Support the `hemco_grid_xdim` and `hemco_grid_ydim` namelist options to configure resolution of internal HEMCO grid

### Fixed
- Updates for CESM 2.2 data structures (`begchunk` and `endchunk` moved to `ppgrid`)

### Removed
- SourceMods no longer present in this version and are now part of CAM merge
