### Trait distributions of butterfly communities

This repository holds the code and data of the study "Beyond community-weighted means: quantifying trait distributions for detecting community assembly patterns". For more information see the publication: 

xxx

The raw count and abundance butterfly data that support the findings of this study are available from the European Butterfly Monitor Scheme (https://butterfly-monitoring.net/), or directly from each participating Butterfly Monitoring Scheme (https://www.catalanbms.org/en and https://ubms.creaf.cat/en) via a signed license agreement. These datasets are, therefore, not included in this repository. Data on abundance indexes and trait distributions is included in the /data folder. All results and figures from the paper can be derived from the included data, and are already stored in the /results folder.

#### Workflow

1) Scripts with prefix 01-03 are included for completeness, but will not work out of the box as they need the raw data not included in the repository. These scripts generate phenology curves, s-index, and collated index for the butterfly species analysed.
2) Script with prefix 04 calculates the trait distributions for each trait, species, and community.
3) Script with prefix 05 plots and analyses trait distributions.

These scripts need a few R packages to run, see the relevant scripts for them. There are further auxiliary files with specific functions or supplementary information that are also included. 

