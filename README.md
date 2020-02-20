# HEMCO_CESM: CESM/CAM interface to modular HEMCO chemistry emissions module

## How to checkout and use:
Currently HEMCO can only operate on a slightly modified version of CAM. Please see [jimmielin/CAM](https://github.com/jimmielin/CAM). **Warning: Unsupported CAM development code is unsupported. You have been warned.**

The forked repository contains `[hemco]` as an external to be deployed inside CAM in `cam/src/hemco`.

## Main components:
* `/HEMCO/src/`: Contains the source for the [HEMCO emissions component](https://github.com/geoschem/HEMCO), imported using `manage_externals` (see `Externals_HCO.cfg`)
* `/hemco_interface.F90`: Main entry point for modular HEMCO emissions module; contains dummy ESMF gridded component used for regridding between CAM and HEMCO lat-lon.
* `/hco_esmf_grid.F90`: HEMCO internal lat-lon grid computation and ESMF regrid wrapper.

### Auxiliary components:
* `/hco_esmf_wrappers.F90`: Helper routines

## License
```
Copyright (c) 2019-2020, Haipeng Lin and Harvard University Atmospheric Chemistry Modeling Group. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```