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
On the Windows desktop a .bat file can be used to create a shortcut to the app

See `example_shortcut.bat`


## Usage
Currently, only Sciex Analyst CDS calibration data can be read.

### Usage with Sciex Analyst software
Export CDS calibration data by opening a cds calibration file in Analyst. View
mass spectra and right-click to export peaks as text. Save the tab-delimited
file in the corresponding directory under `calibration_data`. The name of the file
is irrelevent.

Open the shiny app `performance_monitor_viewer/app.R` to view the data (example
data is provided).