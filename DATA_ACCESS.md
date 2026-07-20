# Data access and redistribution

This repository contains analysis code only. The underlying multiparameter
flow-cytometry and clinical data were obtained from the publicly available
ImmPort study SDY2583 (*Blood Immunotypes*; dataset DOI:
<https://doi.org/10.21430/M3A0B9RD5T>).

## Obtain the data

1. Visit <https://www.immport.org/> and locate study **SDY2583**.
2. Follow the current ImmPort access and download procedure.
3. Unpack the study metadata and the required panel-specific FCS archives on
   the local analysis system.
4. Configure the local locations through `config/paths.R` or the environment
   variables described in `config/paths.example.R`.

## Files intentionally excluded

The following are not distributed through this repository:

- raw or compensated FCS files;
- participant-level clinical or demographic tables;
- integrated participant-level immune/clinical matrices;
- serialized `.RData` or `.rds` analysis objects;
- locally generated results, figures, or logs.

These exclusions prevent unintended redistribution of study data and keep a
clear separation between source code and controlled local analysis material.
Users are responsible for complying with the terms, documentation, and ethical
requirements that apply to their access and use of ImmPort data.

## Recommended manuscript wording

> Multiparameter flow-cytometry and associated clinical data were obtained
> from the publicly available ImmPort SDY2583 Blood Immunotypes study
> (dataset DOI: 10.21430/M3A0B9RD5T). No new participants were recruited and
> no new specimens were collected for this secondary reanalysis.
