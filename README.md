# MSDtracking

R and Python code to analyze single-molecule tracking data as described in ...

For tracking the SOS Plugin (http://smal.ws/wp/software/sosplugin/, Reuter, M., Zelensky, A., Smal, I., Meijering, E., van Cappellen, W.A., de Gruiter, H.M., van Belle, G.J., van Royen, M.E., Houtsmuller, A.B., Essers, J., et al. (2014). BRCA2 diffuses as oligomeric clusters with RAD51 and changes mobility after DNA damage in live cells. The Journal of Cell Biology 207, 599–613.) is used.
Track segmentation and MSD analysis is done using the DLMSS method (Arts, M., Smal, I., Paul, M.W., Wyman, C., and Meijering, E. (2019). Particle Mobility Analysis Using Deep Learning and the Moment Scaling Spectrum. Sci Rep 9, 1–10.)

Data is imported from folders with csv files using 'scripts/analysis script v2.R' and track segmentation is done using 'analyze_MLMSS.R' a adapted script that allows to do DLMSS analysis directy from RStudio, using the Reticulate package.

The package can be installed in R Studio (www.r-studio.com)

```R
install.packages("devtools")  
library(devtools)  


install_github("maartenpaul/DBD_tracking")

library(DBD_tracking)  
```

For questions please contact Maarten Paul (m.w.paul@erasmusmc.nl)
