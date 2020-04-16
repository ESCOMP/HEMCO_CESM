[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v2-blue.svg)](https://www.gnu.org/licenses/lgpl-2.1)

# HEMCO_CESM: CESM/CAM interface to modular HEMCO chemistry emissions module

## How to checkout and use:
Currently HEMCO can only operate on a slightly modified version of CAM. Please see [jimmielin/CAM](https://github.com/jimmielin/CAM). **Warning: Unsupported CAM development code is unsupported. You have been warned.**

The forked repository contains `[hemco]` as an external to be deployed inside CAM in `cam/src/hemco`.

## Main components:
* `/HEMCO/src/`: Contains the source for the [HEMCO emissions component](https://github.com/geoschem/HEMCO), imported using `manage_externals` (see `Externals_HCO.cfg`)
* `/hemco_interface.F90`: Main entry point for modular HEMCO emissions module; contains dummy ESMF gridded component used for regridding between CAM and HEMCO lat-lon.
* `/hco_esmf_grid.F90`: HEMCO internal lat-lon grid computation and ESMF regrid wrapper.
* `/hco_cam_exports.F90`: Utilities for HEMCO to export its computed results via the CAM history component (diagnostic output) and the CAM physics buffer (feeding into the chemistry components).

### Auxiliary components:
* `/hco_esmf_wrappers.F90`: Helper routines

## License
```
Copyright (c) 2019-2020, Haipeng Lin and Harvard University Atmospheric Chemistry Modeling Group.
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