# CEDAS-23-model-selection
Model selection examples and exercises for CEDAS 2023 summer school

`cedas-exercises.pdf` is our set of tasks for the DIY practical.

`mombf-vignette.R` is an adapted code vignette from the `mombf` library (1), demonstrating a simple example of Bayesian model selection.

`chem-data.csv` is synthetic chemical data describing the rate of product formation as a function of reactant concentration (powers of reactant concentration are also given)

`mt-data.csv` and `pt-data.csv` is bioinformatic data on mtDNA and ptDNA gene retention and associated possible predictors, taken from (2).

`cedas-23.pdf` is the lecture slides; `cedas.R` generates example plots for these slides.

Possible solutions
----
`chem-rates.R` captures some of the concepts in the chemical data example.

`odna-analysis.R` captures some of the concepts in the organelle DNA example.

Demonstration examples of model selection -- from the slides -- can be found in `cedas.R`.

# References

(1)  Rossell D, Cook JD, Telesca D, Roebuck P, Abril O, Torrens M (2022). _mombf: Bayesian
  Model Selection and Averaging for Non-Local and Local Priors_. R package version
  3.1.3, <https://CRAN.R-project.org/package=mombf>.

(2) Giannakis, K., Arrowsmith, S.J., Richards, L., Gasparini, S., Chustecki, J.M., Røyrvik, E.C. and Johnston, I.G., 2022. _Evolutionary inference across eukaryotes identifies universal features shaping organelle gene retention_. Cell Systems, 13(11), pp.874-884. <https://www.cell.com/cell-systems/fulltext/S2405-4712(22)00351-9>
