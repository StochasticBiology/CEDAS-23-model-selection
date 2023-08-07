# CEDAS-23-model-selection
Model selection examples and exercises for CEDAS 2023 summer school

`mombf-vignette.R` is an adapted code vignette from the `mombf` library (1), demonstrating a simple example of Bayesian model selection.

`chem-data.csv` is synthetic chemical data describing the rate of product formation as a function of reactant concentration (powers of reactant concentration are also given)

`mt-data.csv` and `pt-data.csv` is bioinformatic data on mtDNA and ptDNA gene retention and associated possible predictors, taken from (2).

`slides/` contains the TeX for the lecture slides (`cedas-23.tex`) and most images therein; those under different licenses are omitted. `cedas.R` generates example plots for these slides.

# References

(1)  Rossell D, Cook JD, Telesca D, Roebuck P, Abril O, Torrens M (2022). _mombf: Bayesian
  Model Selection and Averaging for Non-Local and Local Priors_. R package version
  3.1.3, <https://CRAN.R-project.org/package=mombf>.

(2) Giannakis, K., Arrowsmith, S.J., Richards, L., Gasparini, S., Chustecki, J.M., Røyrvik, E.C. and Johnston, I.G., 2022. _Evolutionary inference across eukaryotes identifies universal features shaping organelle gene retention_. Cell Systems, 13(11), pp.874-884. <https://www.cell.com/cell-systems/fulltext/S2405-4712(22)00351-9>
