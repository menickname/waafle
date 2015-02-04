#!/usr/bin/R

# Allows me to insert commands
args <- commandArgs(TRUE)
library(ggplot2)

# Load functions
source(file='plotcontig.R')

# Import dataframes
CONTIGNAME <- args[1]
ANSWERKEY <- read.table(file=args[2], sep=' ', header=T)
GENETABLE <- read.table(file=args[3], sep=' ', header=T)
NAME <- args[4]

saveplot(CONTIGNAME, GENETABLE, ANSWERKEY, NAME)

