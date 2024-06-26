[![License: LGPL v2](https://img.shields.io/badge/License-LGPL%20v2-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

# HEMCO_CESM: CESM/CAM interface to modular HEMCO chemistry emissions module

This is the official repository for the **HEMCO-CESM** interface, coupling the [Harmonized Emissions Component (HEMCO)](https://github.com/geoschem/HEMCO) to the [Community Earth System Model version 2 (CESM2)](https://github.com/ESCOMP/CESM) via the [Community Atmosphere Model (CAM)](https://github.com/ESCOMP/CAM).

:bulb: HEMCO documentation: https://hemco.readthedocs.io

:bulb: HEMCO on the CAM-chem Wiki: https://wiki.ucar.edu/display/camchem/HEMCO

:book: Reference: [Lin et al., 2021](https://gmd.copernicus.org/articles/14/5487/2021/gmd-14-5487-2021.html)

:envelope: Maintainer: [Haipeng Lin](https://github.com/jimmielin) (hplin@seas.harvard.edu). Please post an issue or pull request in this repository if your request is code-related.

## How to checkout
HEMCO is fully implemented as an emissions component for [GEOS-Chem](https://gmd.copernicus.org/articles/15/8669/2022/) and [CAM-chem](https://wiki.ucar.edu/display/camchem/HEMCO) chemistry within CAM for [MUSICA](https://wiki.ucar.edu/display/MUSICA/MUSICA+Home). Simply download the latest release of CESM (2.3 and above) with CAM (6.3.118 and above) to use.

This repository, `HEMCO_CESM`, is included as an external in CAM (`src/hemco`) and contains the interface for HEMCO to communicate with CAM. This repository then includes HEMCO itself as an external (`src/hemco/HEMCO`) which is model independent and shared by all models implementing HEMCO (GEOS-Chem, GEOS, CESM, WRF, etc.)

## Main components:
* `/HEMCO/src/`: Contains the source for the [HEMCO emissions component](https://github.com/geoschem/HEMCO), imported using `manage_externals` (see `Externals_HCO.cfg`)
* `/hemco_interface.F90`: Main entry point for modular HEMCO emissions module; contains dummy ESMF gridded component used for regridding between CAM and HEMCO lat-lon.
* `/hco_esmf_grid.F90`: HEMCO internal lat-lon grid computation and ESMF regrid wrapper.
* `/hco_cam_exports.F90`: Utilities for HEMCO to export its computed results via the CAM history component (diagnostic output) and the CAM physics buffer (feeding into the chemistry components).
* `/hco_cam_convert_state_mod.F90`: State conversion module between CAM meteorology and chemistry variables into HEMCO extensions format.

### Auxiliary components:
* `/hco_esmf_wrappers.F90`: Helper routines

### Configuration files
...are available in the [HEMCO_CESM_configs](https://github.com/jimmielin/HEMCO_CESM_configs) repository with [instructions on how to translate GEOS-Chem emission species to CAM-chem](https://github.com/jimmielin/HEMCO_CESM_configs/blob/master/CAM-Chem/Mapping_Process.md).

## License
```
Copyright (c) 2019-2024, Haipeng Lin and Harvard University Atmospheric Chemistry Modeling Group.
All rights reserved.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License version 2.1 as published by the Free Software Foundation;

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
```