[![License: LGPL v2](https://img.shields.io/badge/License-LGPL%20v2-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

# HEMCO_CESM: CESM/CAM interface to modular HEMCO chemistry emissions module

This is the official repository for the **HEMCO-CESM** interface, coupling the [Harmonized Emissions Component (HEMCO)](https://github.com/geoschem/HEMCO) to the [Community Earth System Model version 2 (CESM2)](https://github.com/ESCOMP/CESM) via the [Community Atmosphere Model (CAM)](https://github.com/ESCOMP/CAM).

:bulb: HEMCO documentation: https://hemco.readthedocs.io

:bulb: HEMCO on the CAM-chem Wiki: https://wiki.ucar.edu/display/camchem/HEMCO

:book: Reference: [Lin et al., 2021](https://gmd.copernicus.org/articles/14/5487/2021/gmd-14-5487-2021.html)

## How to checkout and use:
Currently HEMCO can only operate on a slightly modified version of CAM. The corresponding [pull request](https://github.com/ESCOMP/CAM/pull/560) is currently under discussion.

**Warning: Unsupported CAM development code is unsupported. You have been warned.**

The forked repository contains `[hemco]` as an external to be deployed inside CAM in `cam/src/hemco`.

## Main components:
* `/HEMCO/src/`: Contains the source for the [HEMCO emissions component](https://github.com/geoschem/HEMCO), imported using `manage_externals` (see `Externals_HCO.cfg`)
* `/hemco_interface.F90`: Main entry point for modular HEMCO emissions module; contains dummy ESMF gridded component used for regridding between CAM and HEMCO lat-lon.
* `/hco_esmf_grid.F90`: HEMCO internal lat-lon grid computation and ESMF regrid wrapper.
* `/hco_cam_exports.F90`: Utilities for HEMCO to export its computed results via the CAM history component (diagnostic output) and the CAM physics buffer (feeding into the chemistry components).
* `/hco_cam_convert_state_mod.F90`: State conversion module between CAM meteorology and chemistry variables into HEMCO extensions format.

### Auxiliary components:
* `/hco_esmf_wrappers.F90`: Helper routines

## License
```
Copyright (c) 2019-2023, Haipeng Lin and Harvard University Atmospheric Chemistry Modeling Group.
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