# Script S1:
* Calculation of temperature variables and correlations between the first year of data and the second year of data
* Please note that raw soil temperature data are not provided in this repository. We only provide temperature measurements for two sites here to make the script run. For the full dataset either query SoilTemp database (https://www.soiltempproject.com/) or contact authors (Kryštof Chytrý: krystof.chytry@gmail.com, https://www.mountainresearch.at/microclim/).
* The output of this code is `data\\temperature_models\\temperature_variables.csv` which is provided in the repositry.

# Script S2:
* Fitting all models and analysis producing figures 1 to 5.
* In this script we perform:
  * Modelling: Fit SDMs (binomial GLM) for all species with more than 30 occurrences, excluding aggregates or those identified on species level only. And we store D-squared values.
  * Clustering of temperature variables: We apply a hierarchical clustering algorithm on all (scaled) temperature variables, consdering multiple numbers of divisions (from 2 to 10)
