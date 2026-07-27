## Resubmission

This is a resubmission of `specmine`.

In this version I have:
- added missing `\value{}` tags and clarified returned objects in `.Rd` files;
- removed commented example code and updated examples to be suitable for checks;
- replaced or guarded examples that depended on optional packages during checks;
- updated examples/tests/documentation to avoid writing to the user home filespace;
- reviewed functions that write files to avoid default output paths in the package directory or `getwd()`;
- improved package documentation consistency and metadata.

## Test environments
- Local macOS
- win-builder (R-devel)

## R CMD check results
0 errors | 0 warnings | 1 note

This package was previously archived. The remaining note on win-builder is related to CRAN incoming feasibility for a resubmission of an archived package.

The package was also checked locally on macOS with `R CMD check`, which completed with 0 errors, 0 warnings, and 0 notes.

The suggested packages `cyjShiny` and `specmine.datasets` are used only for optional functionality. `specmine.datasets` is available on CRAN.

## Reverse dependencies
None.