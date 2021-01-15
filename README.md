This is R code for a simulation project that aims to explore the effects of omitting untrawlable habitat from the NMFS survey. I try to accomplish this by: 

1) simulating several seascapes of random autocorrelated noise 
2) assigning trawlable/untrawlable definitions at the same rate they currently appear in the GOA sampling grid 
3) simulating several populations across the GOA grid. The simulated populations are based on the grid cells abundance levels of Zack Oyafuso's operating model,
then allocating the generated grid cell abundance levels across the grid to be correlated with the trawlability variable at various strengths (.10, .25, .50, .75, .90) 
4) simulating stratified random sampling under 2 scenarios (incorporating/omitting untrawlable habitats) based on area swept data from the last decade of NMFS surveys

