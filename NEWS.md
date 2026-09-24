# nmaFlow 0.1.0.9000

## Development fixes

* Corrected the NMA, netdose, MBNMAtime, crossnma, nmathresh, mixmeta and multinma backend adapters identified by the executable API audit.
* Added explicit failure handling for netmeta predictions and multinma node-split fits without posterior samples.
* Added regression tests for all corrected adapters.

* Initial development release.
* Adds audit-first data architecture with quarantining, explicit node mapping,
  fine/coarse network support and network-formation diagnostics.
* Adds frequentist, Bayesian, bootstrap and advanced-backend orchestration.
* Adds agronomic simulation datasets and extended teaching vignettes.
* Adds publication-ready plotting/table interfaces and reproducibility reports.
