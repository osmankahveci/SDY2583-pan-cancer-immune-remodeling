# Data access and redistribution

This repository contains analysis source only. The multiparameter
flow-cytometry and associated metadata are obtained from the publicly available
ImmPort study SDY2583 (*Blood Immunotypes*; dataset DOI:
<https://doi.org/10.21430/M3A0B9RD5T>).

## Obtain and prepare the data

1. Locate study **SDY2583** at <https://www.immport.org/>.
2. Follow the current ImmPort access and download procedure.
3. Unpack the study tabular package and the required FCS archives locally.
4. Place panel FCS files under `data/raw/<PANEL>` or configure their actual
   locations in `config/paths.R`.
5. Set `SDY2583_IMMPORT_DOWNLOAD_DIR` to the unpacked tabular package.
6. Run the metadata integration builder, or configure an existing compatible
   subject-level matrix:

```bash
Rscript R/integration/STEP0_build_subject_metadata_matrix_RECONSTRUCTED.R
```

The builder searches for the ImmPort arm-to-subject, subject, and
subject-to-flow-result tables and creates local mappings for subject ID, FCS
file, panel, disease group, age, age group, and sex. The generated matrix is
written under `data/derived/metadata/` and is ignored by Git.

A separate local cancer-subgroup/therapy annotation table may be supplied with
`SDY2583_CLINICAL_ANNOTATION_FILE`. Such participant-level material is not
redistributed.

## Files intentionally excluded

The following must remain local:

- raw, compensated, or transformed FCS files;
- downloaded ImmPort archives and tabular packages;
- participant-level clinical or demographic tables;
- participant-level feature, score, and matched-analysis datasets;
- integrated immune/clinical matrices;
- serialized `.RData` or `.rds` objects;
- generated results, figures, validation reports, session files, and logs.

These exclusions separate reusable analysis source from the local study data
and derived participant-level material. Users remain responsible for the terms,
documentation, and ethical requirements governing their ImmPort access and use.

## Reproduction command

After local configuration, the full source workflow is launched with:

```bash
Rscript scripts/01_run_all_panels.R
```

The runner does not download or redistribute data. It operates only on the
locally configured files.

## Recommended manuscript wording

> Multiparameter flow-cytometry and associated clinical data were obtained
> from the publicly available ImmPort SDY2583 Blood Immunotypes study
> (dataset DOI: 10.21430/M3A0B9RD5T). No new participants were recruited and
> no new specimens were collected for this secondary reanalysis. Analysis code
> and reproducibility documentation are available in the accompanying software
> repository; source and participant-level data are accessed directly through
> ImmPort and are not redistributed with the code.
