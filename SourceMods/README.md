# SourceMods

These are necessary SourceMods for HEMCO_CESM to operate correctly with CAM-Chem. Put them in `SourceMods/src.cam`.

These files are renamed `.G90` so they are not compiled. You should rename them to `.F90` when using them:

```bash
find . -name "*.G90" -exec bash -c 'mv "$1" "${1%.G90}".F90' - '{}' \;
```