# Supported reproduction and historical evidence
## Supported clean-clone endpoint
Install R and, separately, `ggplot2` and `patchwork` as described in [environment/README.md](environment/README.md). From the repository root:
```sh
Rscript --vanilla scripts/validation/validate_final_synthesis.R --output=../PARPi_synthesis_check
```
The output parent must exist; the output directory must not exist. The runner refuses to write inside the repository. It verifies fingerprints, copies six frozen inputs to separate script-62/script-63 stages, parses the portable scripts, changes only their working-directory expression in memory, runs final integration/plotting, and compares against frozen references. No DE, GSEA or model fitting occurs. Required dependencies are two plotting packages plus base R; installation never occurs in this runner.

Successful verification requires 46/46 columns: eleven overlapping script-62 accession fields and 35 columns across five script-63 tables. IDs/categories/missingness are exact; numeric absolute tolerance is 1e-12. The report also records CSV byte comparison and source integrity. The generated Figure 8 is validation artwork; final manuscript artwork is under `figures/publication`.

On Windows with an invalid inherited `LC_ALL=C.UTF-8`, set `$env:LC_ALL='English_India.utf8'` and `$env:LANG='en'` in PowerShell before invoking R. Set `R_LIBS_USER` to your installed user-library directory if R cannot locate the packages. This locale correction was required on the release-preparation workstation and changes no scientific expressions. A UTF-8 compatible locale and paths are required.

## Evidence and execution order
1. For the supported endpoint, use the frozen inputs as committed; no upstream steps are needed.
2. For historical inspection, acquisition precedes the per-accession script ranges in [dataset manifest](metadata/DATASET_CONTRAST_MANIFEST.csv). Scripts 01-07 include discovery and frozen hypotheses. D1-D4 scripts 23-28 are historical, not the final endpoint.
3. Historical D1-D6 reconstruction consumes the eleven committed Hallmark tables; its exact executed source and 49-column evidence are retained. D4 reconstruction consumes acquired D4 counts and committed reference/filter/metadata tables; its previously executed evidence is retained, but DE has not been rerun for this release.
4. Historical script 45 and scripts 46-49 precede D7/D8. Recovered candidates are not proven original executed sources. Final repaired 62 and 63 are invoked only through the supported isolated runner.
5. Graphical-abstract script 64 is the exact author-executed generator. Run artwork code only in an isolated project copy; it writes publication outputs. The nine final manuscript PNGs were recovered unchanged from the author-review DOCX. Their native final assembly source was not recovered. Extraction is reproducible if the privately held source manuscript is supplied; it is not bundled.

Historical scripts may install packages, select files interactively and overwrite result paths. [HISTORICAL_SIDE_EFFECTS.csv](environment/HISTORICAL_SIDE_EFFECTS.csv) locates these dependencies. They are source records, not a supported batch pipeline. Never run all numbered scripts automatically. `PARPI_PROJECT_ROOT` configures recovered candidate paths; execute only against an isolated prepared copy. Existing imported historical copies use their recorded portability adaptations.

No full raw-to-results reproduction, exact historical environment restoration, or complete historical script-62 master equivalence is claimed. See [status](REPRODUCTION_STATUS.md).
