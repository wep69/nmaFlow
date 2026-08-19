# nmaFlow

`nmaFlow` is an R-first integration and teaching package for Network Meta-Analysis (NMA).
Its central rule is: **make the evidence network scientifically coherent before fitting the model**.

## Design principles

- auditable arm-level and contrast-level data architecture;
- conservative network repair with quarantine and decision logs;
- explicit fine/coarse treatment nodes;
- explicit unit conversion and outcome/time/population stratification;
- preservation of multi-arm study identity;
- frequentist orchestration through `netmeta` and `NMA`;
- Bayesian NMA/NMR and IPD/ML-NMR through `multinma`;
- study-level bootstrap that keeps multi-arm studies intact;
- transitivity, inconsistency, influence and sensitivity workflows;
- dose-response/nonlinear, component, time-course, cross-design, multivariate and SEM extensions;
- publication-oriented figures, tables, league tables and reports;
- reproducible agronomic teaching simulations;
- optional interoperability with Python, Julia, SAS, Stata, Mplus, JAGS and BUGS.

Heavy statistical backends are optional dependencies. Python and Julia are never required for the core package.

## Documentation

Start with:

- `vignettes/00-foundations-to-advanced-tutorial.Rmd`: complete ecosystem tutorial;
- `vignettes/01-data-architecture-network-formation.Rmd`: detailed network-formation and repair workflow;
- `vignettes/23-network-repair-casebook.Rmd`: applied failure-mode casebook;
- `validation/LOCAL_VALIDATION.md`: local R validation procedure.

The package contains 24 pedagogical vignettes and an 80-block instructional registry available with `nma_blocks()`.

## Development documentation

`NAMESPACE` is generated from `roxygen2` annotations. After cloning/unpacking the development tree, run:

```r
roxygen2::roxygenise()
```

Then execute the test and check workflow described in `validation/LOCAL_VALIDATION.md`.
