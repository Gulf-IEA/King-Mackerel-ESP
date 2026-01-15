# King Mackerel Ecosystem and Socioeconomic Profile (ESP) Report

## Contributing to this report
### Adding an indicator idea / dataset

1. On the repository home page, click the projects tab and select the King-Mackerel-ESP-metric-list project.
2. Click the plus button to add an issue using the Metric Idea Form. Make sure the linked repository is Gulf-IEA/King-Mackerel-ESP. 
3. Fill out any information you have on the data needed to create the indicator, including data source, POC, etc. If you don't have all the information right away that's fine, it can be added later.

### Writing code to manipulate data into an indicator

1. Ensure you are using R version 4.3.3 on your local machine. This is the version being used by the lock file. 
2. Clone the repository onto your local machine.
3. When you first open the repository in RStudio, you will prompted to install the packages in the lockfile. Run renv::restore() to restore the project library.
4. To work on an indicator, follow the guidance in the [ESR Best Practices](https://gulf-iea.github.io/ESR-Best-Practices/) document (chapter 4: data processing). 

## Folder structure
* data 
  * unformatted - data files retrieved from collaborators that cannot be accessed from an online source
  * intermediate - any intermediate data files created from R scripts
  * formatted - final, formatted metric data
    * formatted csvs - csvs formatted in the standard format (see below)  
    * final objects - metric objects created by running final csv files through the IEAnalyzeR data_prep function 
* scripts  
  * metrics
    * automated - scripts for processing metric data that can be pulled automatically from some online source
    * confidential - scripts for processing metric data that is confidential (e.g., fishery-dependent)
    * non-automated - scripts for processing unformatted metric data from the data/unformatted folder
  * plotting - scripts used for plotting metrics and indicators
  * other - home for other assorted code scripts
* figures
  * plots - home for any metric/indicator plots created from scripts
  * images - home for any other images used in the report/other products
* synthesis
  * data - data files related to indicator synthesis 
  * scripts - scripts related to indicator synthesis
  * outputs - indicator synthesis output files (tables, plots, etc.)
* report - home for all files needed to create the ESR quarto book report (standard NOAA tech doc format)
* sandbox - home for any misc. things, works in progress, etc.
  * prelim_code - home for code during the metric idea development before code review and branch merging

## Other components of this repository
This repository also contains two issue template files (within the .github/ISSUE_TEMPLATES folder). 
* metric idea - to be filled out whenever there is an idea for a new metric that needs to be considered for inclusion in the current ESR. These issues will populate the "Indicator scoping" project
* metric production - to be filled out once a metric that is being used in the current ESR has been developed (data identified and code written to create final formatted files). These issues will populate the quarto "Methods Document"

A Github action renders and publishes the quarto book automatically.

------------------------------------------------------------------------

### Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

### License

This content was created by U.S. Government employees as part of their official duties. This content is not subject to copyright in the United States (17 U.S.C. §105) and is in the public domain within the United States of America. Additionally, copyright is waived worldwide through the CC0 1.0 Universal public domain dedication.

  
