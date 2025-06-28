# medRCT 0.1.0

* An initial public release of medRCT package, version 0.1.0

* Supports multiple mediators, intermediate confounders, and various variable types, 
including categorical exposures, and binary or continuous mediators, 
intermediate confounders, and outcomes.

* Handles scenarios with a single mediator and/or without intermediate confounders.

* Provides a shiny app for model assessment. 

# medRCT 0.1.1

* Updated the reference for the JOSS paper.

# medRCT 0.2.0

* JOSS paper published (added badge in README)

* Fix bug when mediator and intermediate confounder are factors

* Bootstraps that return warnings or errors are now disregarded, and the number of failed bootstraps is reported in the summary output

* Add argument to specify the effect measure

* Add package logo