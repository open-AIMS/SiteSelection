# Site Selection #
---

Analyses to test a site selection framework for seeded coral deployments. We test the framework proposed by Humanes et al. (2025) using spat survival data from deployments occurring between 2021 and 2023.


## Description of the folder structure ##

There should be three folders: **data** (contains the raw data), **scripts** (contains two R scripts, one that contains the functions -*functions.R*- and
one that runs the analyses and generates the figures - *analyses.R*), and **outputs** (where the plots will be saved).


### Description of the data ##
The Data folder has two files: *environment.csv* and *survival.csv*, which contain the environmental factors associated to the deployment
sites and the spat survival data, respectively.

The *environment.csv* has one deployment site per row and includes the following variables:
* Lead Researcher- researcher leading the coral seeding deployment
* Reef- name of reef where the coral spat where deployed
* Site - unique identification for each deployment site
* Longitude- site longitude
* Latitude - site latitude
* Deployment.year - years in which the spat were deployed
* Deployment.date - date of deployment
* Census.date - date of spat survival census
* Deployment.duration - number of days between deployment and census
* Species - name of species deployed at the site 
* Deployment.method - type of attachment to the substratum (fixed directly, freely deployed or tethered)
* Sedimentation - sedimentation category
* Rugosity - categorical measurement of rugosity
* Slope - slope of the site obtained from bathymetry maps
* Water.velocity - estimates of mean horizontal water velocity at the seabed (in m s-1) from Callaghan et al. 2023
* Consolidated.substrate - estimate of the percentage of hard substrate at the site
* Hard.coral - estimate of the percentage of hard coral cover at the site
* Macroalgae - estimate of the percentage of macroalgal cover at the site
* Free.space - estimate of the percentage of free consolidated substrate at the site
* score.X - score assigned to each (X) environmental factor following Humanes et al. (2025)
* Reference - study from which the deployment information was obtained

The *survival.csv* has survival data for one seeding unit and species combination per row and includes the following variables:
* Reef- name of reef where the coral spat where deployed
* Site - unique identification for each deployment site
* Deployment.year - years in which the spat were deployed
* Deployment.duration - number of days between deployment and census
* Species - species for which the survival was recorded
* n_trial - number of tabs/plugs with spat from at the time of deployment
* n_surv - number of tabs/plugs with living spat at the time of census
* Lead Researcher- researcher leading the coral seeding deployment

### Reference ##
Humanes, A., Fabricius, K.E., Ferrari, R. & Ortiz, J.C. (2025) Ecological quantitative criteria for reef site prioritisation to maximise survivorship and growth of outplanted corals. Journal of Environmental Management, 392, 126585.
