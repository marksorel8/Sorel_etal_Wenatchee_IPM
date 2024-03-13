## [*Incorporating life history diversity in an integrated population model to inform viability analysis*](https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2023-0118)

#### Mark H. Sorel, Jeffrey C. Jorgensen, Richard W. Zabel, Mark D. Scheuerell, Andrew R. Murdoch, Cory M. Kamphaus, and Sarah J. Converse

##### Please contact the first author for questions about the code or data: Mark Sorel ([marksorel8\@gmail.com](mailto:marksorel8@gmail.com){.email})

##### Secondary contact: Sarah Converse ([sconver\@usgs.gov](mailto:sconver@usgs.gov){.email})

------------------------------------------------------------------------

## Abstract

Life history diversity can significantly affect population dynamics and effects of management actions. For instance, variation in individual responses to environmental variability can reduce extirpation risk to populations, as the portfolio effect dampens temporal variability in abundance. Moreover, differences in habitat use may cause individuals to respond differently to habitat management and climate variability. To explore the role of life history diversity in population trajectories, population models need to incorporate within-population variation. Integrated population modeling (IPM) is a population modeling approach that offers several advantages for sharing information and propagating uncertainty across datasets. In this study, we developed an IPM for an endangered population of Chinook salmon (*Oncorhynchus tshawytscha*) in the Wenatchee River, Washington, USA, that accounts for diversity in juvenile life histories, spawning location, and return age. Our analysis revealed that diversity in the age of juvenile emigration from natal streams had a portfolio effect, resulting in a 20% reduction in year-to-year variability in adult abundance in population projections. Our population viability analysis suggests that management interventions may be necessary to meet recovery goals, and our model should be useful for simulating the outcomes of proposed actions.

### Table of Contents

#### [Scripts](./scripts)

Contains main script to run analyses: *PVA-results-with-hatchery.Rmd*. Contains scripts to produce figures for supplements on juvenile production (*juvenile production.Rmd*), survival(*survival.Rmd*), and covariate simulation(*Covariate simulation appendix.Rmd*).

#### [Data](./data)

Contains some of the data that the model was fit to, including sex ratios, carcass recoveries, fecundity, broodstock removals, and processed environmental covariates.

#### [Src](./src)

Contains the TMB code for the integrated population model: *Wen_spchk_IPM.cpp*.

##### [Src/functions](./src/function)

Contains wrapper functions to build the data list to feed to the model (*make_data.R*), initialize and optimize the model (*initialize_model.R*), and conduct population projections (*pop_projection.R*)

#### [Juvenile prduction module scripts](./Wenatchee-screw-traps)

Contains estimates of juvenile emigrant abundance form an observation model (*all_emigrants_estimates.rds*), which are fed to the IPM as data. Also contains environmental covariate values, and code and data used to generate the juvenile abundance estimates. See [Sorel et al. 2023, Eco. Evo.](https://github.com/Quantitative-Conservation-Lab/Sorel_etal_2023_Ecology-Evolution) for more details.

#### [Survival module scripts](./Wenatchee-survival)

Contains mark-recpature data and code for processing that data into formats for use fitting models. Contains helper functions for examining estimates of survival and return rate from the ocean. See [Sorel et al. 2023, Ecosphere](https://github.com/Quantitative-Conservation-Lab/Sorel_etal_2023_Ecosphere) for more details.

#### Required Packages and Versions Used

base_4.2.2\
dplyr_1.1.0\
forcats_1.0.0\
ggplot2_3.4.1\
graphics_4.2.2\
grDevices_4.2.2\
here_1.0.1\
lubridate_1.9.1\
methods_4.2.2\
purrr_1.0.1\
readr_2.1.2\
readxl_1.4.3\
stats_4.2.2\
stringr_1.5.0\
tibble_3.1.8\
tidyr_1.3.0\
tidyverse_2.0.0\
TMB_1.9.3\
TMBhelper_1.4.0\
utils_4.2.2\
viridisLite_0.4.1

#### Details of Article

Sorel, MH, Jorgensen, JC, Zabel, RW, Scheuerell, MD, Murdoch, AR, Kamphaus, CM, Converse, SJ. In press. Incorporating life history diversity in an integrated population model to inform viability analysis. *Canadian Journal of Fisheries and Aquatic Scciences*.

#### How to Use this Repository

The Rmarkdown files located in the 'scripts' folder are the controlling scripts for the analysis. The *PVA-results-with-hatchery.Rmd* file is the main analysis for the paper and the others are for examining aspects of the model that are not the focus of the main manuscript. The rest of the code in the repository is sources by the Rmarkdown files so that you can understand how the results reported in the manuscript were generated from the data.
