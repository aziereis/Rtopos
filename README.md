# How to start Rtopos
function to make ERP topos in R
Install via 

Load the devtools package and download the package
```
library(devtools)
install_github("aziereis/Rtopos")
```
(if you don't have devtools already installed, install this first: 
```
install.packages("devtools")
library(devtools)
```

## How to start with the function
First you need to shape the data in the format that can be included in the function. 
The data should have the dimensions rows = condition(-levels), columns = condition variables and electrodes. 
See also an example data set included in the package: 

```
library(Rtopos)
View(Rtopos::xdat)
```

To create a topoplot, you have several parameters to set and customize, e.g., the number of channels included, the colour palettes etc. 
Check out the help page for details:

```
?maketopoplot
```

The topoplot will be exported as a list of ggplot objects, which you can then save or export easily. 

```
tplots <- maketopoplot(Rtopos::xdat, condition = "difficulty", nrchans = 128, extralegend = F)
tplots[[1]] # to see the first condition level "high" 

ggplot2::ggsave("example_topoplot.png", tplots[[1]], width = 3, height = 4)
```

![man/example_topoplot.png](man/example_topoplot.png)
