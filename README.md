# ntscontrol
Interactive control cards to monitor LC-HRMS/MS performance, specifically for
non-target-analysis of environmental samples

## Funding
Publicly funded project at the German Federal Institute of Hydrology with funding
from the Federal Environmental Agency (FKZ 3717222300)

## Installation
Clone/copy the repository to a local directory. 

### Requirements

* R 4.0
* Packages
  - shiny
  - tidyverse

### Desktop shortcut to shiny app
On the Windows desktop a `.bat` file can be used to create a shortcut to the app

See `example_shortcut.bat`

## Usage
Currently, only Sciex Analyst or SciexOS CDS calibration data can be read.

The file `performance_monitor_settings.R` contains the user settings. The
`instrumentName` is only for display purposes. The `calReportType` must be
either `"analyst"` or `"sciexos"`. The `CAL_REPORTS` are the path(s) to the 
directories containing the calibration report files (txt files for Analyst, 
csv files for SciexOS). Be sure to use double backslash.

Other settings are the warning levels and weighting which need to be adjusted 
for the specific system.

Open the shiny app `shiny/app.R` to view the data (example data is provided).