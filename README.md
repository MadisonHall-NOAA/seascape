This is R code for a simulation project that aims to explore the effects of omitting untrawlable habitat from the NMFS survey. I try to accomplish this by: 

1) simulating several seascapes of random autocorrelated noise 
2) assigning trawlable/untrawlable definitions at the same rate they currently appear in the GOA sampling grid 
3) simulating several populations across the GOA grid. The simulated populations are based on the grid cells abundance levels of Zack Oyafuso's operating model,
then allocating the generated grid cell abundance levels across the grid to be correlated with the trawlability variable at various strengths (.10, .25, .50, .75, .90) 
4) simulating stratified random sampling under 2 scenarios (incorporating/omitting untrawlable habitats) based on area swept data from the last decade of NMFS surveys

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project
code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the
Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub
project will be governed by all applicable Federal law. Any reference to specific commercial products,
processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or
imply their endorsement, recommendation or favoring by the Department of Commerce. The Department
of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to
imply endorsement of any commercial product or activity by DOC or the United States Government.
