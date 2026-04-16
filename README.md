# AC-popstr-meta-analysis-2026
Script repository for Arctic cod population structure derived using a meta-analytical approach. This repository comprises the relevant script(s) and links to data files (on [*Borealis*](https://borealisdata.ca/) database) forming the basis for the manuscript titled: 

### Population structure in Arctic cod *Boreogadus saida* across its range: a meta-analysis approach 
- submitted to ICES Journal of Marine Science.

Authored by: Marie Launay, Caroline Bouchard, and Denis Roy 

#
### Script(s):
To see the main scripts used for this research, go to 'scripts' folder.

In this repository - we wrote several scripts but merged them into a single long one split by different tasks/analyses carried out in the study.
* ACmeta25.R: This is the main script that reads in the data that was collected and compiled into one large ".csv" file that can be downloaded from the Borealis database (here) once the paper is published.

* Initially the script begins by generating the substitute data for 2 studies that do not report actual *F<sub>ST<sub>* values, but do report data that can be used to derive their approximate values. It then reads in the compiled data and searches for all unique geographic positions through all studies, which it then plots using [*ggOceanmaps*](https://mikkovihtakari.github.io/ggOceanMaps/).

* Once the spatial arrangement of studies is plotted, the script then performs the main global meta-analyses using the [*metafor*](https://www.metafor-project.org/doku.php/metafor) package in R. The script also runs a funnel plot type analysis and some model diagnostics to help with model confidence.

* Finally, the script then runs through several different models in attempts to test varying hypotheses as to what study and/or environmental factor(s) help predict the variability observed in pairwise *F<sub>ST<sub>* reported in the different compiled studies.

* Annotation for the script is contained inside the script, and the data will be available once the manuscript has been approved for publication.
