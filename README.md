# TraitGMYC

This package requires the package splits, which is available via R-Forge.
```{r, warning = F, echo = F}
install.packages("paran")
install.packages("splits", repos="http://R-Forge.R-project.org")
```

Moreover, simulations require the package phybase, which is not on cran anymore. It can be installed from github.
```{r, warning = F, echo = F}
install.packages("remotes")
remotes::install_github("bomeara/phybase")
```
