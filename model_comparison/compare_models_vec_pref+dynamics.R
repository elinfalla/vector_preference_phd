################################################################################################################
##### COMPARISON OF VECTOR PREFERENCE AND VECTOR DYNAMICS IN MADDEN (2000) AND DONNELLY (2019) TYPE MODELS #####
################################################################################################################

### MADDEN MODEL:
### Modified version of model by Madden et al. 2000 to include vector preference - both of plant attractiveness and
### plant acceptability. Vector population dynamics (emigration birth, death) are included to emulate the Donnelly model,
### plants are represented by an SI rather than an SEIR model, vectors represented by an SI rather than SEI model. 
### Vectors cannot be born infective (as this model represents NPT viruses).

### One version of the model has a fixed value of phi (number of plants visited pre day), the other has a variable value
### of phi to account for the fact that feeding takes longer than probing, determined as defined by Cunniffe et al. 2021.

### DONNELLY MODEL:
### Modification of deterministic model by Donnelly et al. 2019 to so it tracks I, number of infected plants, rather than 
### i, proportion of infected plants. Otherwise this model is equivalent to the one in Donnelly et al. 2019.

rm(list=ls())

# packages
library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(grid)
library(gridtext)

####