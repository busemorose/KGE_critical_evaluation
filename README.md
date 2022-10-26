Technical note: When best is the enemy of good – critical evaluation of
performance criteria in hydrological models
================
Guillaume Cinkus, Naomi Mazzilli, Hervé Jourde, Andreas Wunsch, Tanja
Liesch, Nataša Ravbar, Zhao Chen, and Nico Goldscheider

Author ORCIDs:

- Guillaume Cinkus
  [0000-0002-2877-6551](https://orcid.org/0000-0002-2877-6551)
- Naomi Mazzilli
  [0000-0002-9145-5160](https://orcid.org/0000-0002-9145-5160)
- Hervé Jourde
  [0000-0001-7124-4879](https://orcid.org/0000-0001-7124-4879)
- Andreas Wunsch
  [0000-0002-0585-9549](https://orcid.org/0000-0002-0585-9549)
- Tanja Liesch
  [0000-0001-8648-5333](https://orcid.org/0000-0001-8648-5333)
- Nataša Ravbar
  [0000-0002-0160-1460](https://orcid.org/0000-0002-0160-1460)
- Zhao Chen [0000-0003-0076-7079](https://orcid.org/0000-0003-0076-7079)
- Nico Goldscheider
  [0000-0002-8428-5001](https://orcid.org/0000-0002-8428-5001)

# Description

This repository contains the following elements:

- `/script/ANN model`: ANN model code
  - Python script to run the ANN model
  - Dummy data file to illustrate the structure of input data
- `/script/Reservoir model`: Reservoir model code
  - KarstMod `.properties` file
  - R script to perform the snow routine

# Synthetic time series

# Hydrological models

## KarstMod

<img src="img/karstmod.png" id="fig-karstmod" />

Information on the KarstMod platform can be found in Mazzilli et
al. (2019b) and the user manual (Mazzilli and Bertin, 2019a). The main
workflow is:

1.  Prepare the input data according to the KarstMod format
2.  \[Optional\] Open a KarstMod `.properties` file
3.  Import the input data
4.  Define warm-up/calibration/validation periods
5.  Define a path for the output directory
6.  Run calibration

It is possible to modify the model parameters, the objective function,
the number of iterations, the maximum time, and other options. The
`Save` button allows to save the modifications into a new `.properties`
file.

- [Download KarstMod](https://sokarst.org/en/softwares-en/karstmod-en/)
- [Download KarstMod User
  manual](https://hal.archives-ouvertes.fr/hal-01832693)

## ANNs

Information on 1D-Convolutional Neural Networks can be found in **ADD
DOI of ANN-lumped study**.

Dependencies:

- [Python 3.8](https://www.python.org/)
- [Tensorflow 2.7](https://www.tensorflow.org/)
- [BayesianOptimization
  1.2](https://github.com/fmfn/BayesianOptimization)
- [Numpy 1.21](https://numpy.org/)
- [Pandas 1.4](https://pandas.pydata.org/)
- [Scipy 1.7](https://scipy.org/)
- [Scikit-learn 1.0](https://scikit-learn.org/stable/)
- [Matplotlib 3.5](https://matplotlib.org/)

# References

Chen, Z., Hartmann, A., Wagener, T., and Goldscheider, N.: Dynamics of
water fluxes and storages in an Alpine karst catchment under current and
potential future climate conditions, Hydrol. Earth Syst. Sci., 22,
3807–3823, https://doi.org/10.5194/hess-22-3807-2018, 2018.

Hock, R.: A distributed temperature-index ice- and snowmelt model
including potential direct solar radiation, J. Glaciol., 45, 101–111,
https://doi.org/10.3189/S0022143000003087, 1999.

Mazzilli, N. and Bertin, D.: KarstMod User Guide - version 2.2,
hal-01832693, 103927, 2019a.

Mazzilli, N., Guinot, V., Jourde, H., Lecoq, N., Labat, D., Arfib, B.,
Baudement, C., Danquigny, C., Soglio, L. D., and Bertin, D.: KarstMod: A
modelling platform for rainfall - discharge analysis and modelling
dedicated to karst systems, Environ. Model. Softw., 122, 103927,
https://doi.org/10.1016/j.envsoft.2017.03.015, 2019b.

Pianosi, F., Sarrazin, F., and Wagener, T.: A Matlab toolbox for Global
Sensitivity Analysis, Environmental Modelling & Software, 70, 80–85,
https://doi.org/10.1016/j.envsoft.2015.04.009, 2015.