# Updated September 7, 2022

# Changes since GRQA\_v1.2
It was discovered that due to a preprocessing error in the previous versions of GRQA some of the parameters originating from WQP were assigned the incorrect code.
The error was caused by certain source parameter codes shifting by during the creation of the corresponding code map (WQP\_code\_map.csv).
As a result, for example, parameter pH got the code for BOD5, while TEMP got the code for TSS.
Parameters affected by this error were TPP, TDP, TP, TN, TDN, POC, DOC, TOC, BOD5, pH, TSS and TEMP.
The processing error also meant that the statistics calculated along with the number of outliers in those parameters were also affected.
In GRQA\_v1.3, this error in WQ\_code\_map.csv has been fixed and all the statistics and plots have been updated to account for these changes.
In addition, the data catalog has been updated as well.

# GRQA\_v1.3 contents
GRQA\_data\_v1.3.zip contains the CSV files for processed observation data (\*\_GRQA.csv) along with corresponding metadata (observation site information, etc) for 42 different water quality parameters.

In addition to the observation data, ZIP files containing additional metadata, figures and GRQA source datasets are provided.

## GRQA\_meta
The GRQA\_meta.zip file contains the following files:
List of parameter codes in GRQA used for looping over during parallel processing (GRQA\_param\_codes.txt)
Potential duplicate observation sites per parameter (\*\_dup\_obs.csv)
Basic statistics of GRQA observation time series per parameter (GRQA\_param\_stats.csv)

## GRQA\_data\_catalog\_v1.3
An overview of all 42 parameters in the form of statistical figures is given in the data catalog document (GRQA\_data\_catalog\_v1.3.pdf).

## GRQA\_figures
The figures shown in the data catalog are also included in the file GRQA\_figures.zip for each parameter:
Map of spatial distribution of observation sites per source dataset (\*\_spatial\_dist.png)
Map of monthly availability of time series per site (\*\_availability.png)
Map of monthly continuity of time series per site (\*\_continuity.png)
Map of median observation values per site (\*\_median.png)
Temporal distribution plot of observations per source dataset (\*\_temporal\_hist.png)
Histogram of observation values per source dataset (\*\_hist.png)
Box plot of observation values per source dataset (\*\_box.png)

Additional grid plots (GRQA\_\*\_grid.png) of each of the aforementioned seven plot types showing DO, DOC, TP and TSS were created to be included in the scientific paper describing GRQA.

## GRQA\_source\_data
In addition, the file GRQA\_source\_data.zip containing the five source datasets (CESI, GEMSTAT, GLORICH, WATERBASE and WQP) along with metadata and statistics collected from both raw and processed data is given.

The raw/meta folder in each of the source data subfolders contains the following files:
Lookup tables used for harmonizing the parameters and units (\*\_code\_map.csv)
Unit files for sources with multiple units per parameter (\*\_units.csv)
Remark code lookup for GLORICH (GLORICH\_remark\_codes.csv)
US state codes used for downloading WQP temperature data (fips\_state.csv)

The processed/meta folder in each of the source\_data subfolders contains the following files:
Number of missing values in source file columns (\*\_missing\_values.csv)
Source file size and row count information (\*\_file\_info.csv)
Observation data statistics per parameter before harmonization (\*\_raw\_stats.csv)
Observation data statistics per parameter after harmonization (\*\_processed\_stats.csv)
Duplicate site IDs in source data (\*\_dup\_sites.csv)

## Code repository
Scripts used for the creation of GRQA are found on Zenodo at https://doi.org/10.5281/zenodo.5082147
