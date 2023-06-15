# PaleoVegChange

# code for Stegner MA & Spanbauer TL, North American pollen records provide macro-scale ecological evidence for the Anthropocene

# to run this code: 

# 1) Within the Wang-et-al-Cores folder, unzip the Bayesian age models and move all the individual core folders into the main Wang-et-al-Cores folder.
# Age models are from (Wang et al. 2019, Bayesian ages for pollen records since the last glaciation. Sci Data 6, 176). Code to generate these age models is available at https://github.com/yuewangpaleo/BaconAgeNeotoma NOTE: some components of the BaconAgeNeotoma code is depreciated due to the update of the 'neotoma' package to 'neotoma2'. Age models are provided here with permission from Yue Wang.

# 2) Code to select sites based on sample resolution and other criteria is provided in PaleoVegChange-choose-sites.R This code relies on the depreciated 'neotoma' package and so the final output of PaleoVegChange-choose-sites.R, which is necessary for continental and regional analyses, is provided as an Rdata object, pol_dlx_2020-12-01.Rdata

# 3) for continental analysis, run PaleoVegChange-continent.R

# 4) for regional analysis, download the EPA Level I Ecoregions of North America shapefile available at https://www.epa.gov/eco-research/ecoregions-north-america then run PaleoVegChange-regional.R