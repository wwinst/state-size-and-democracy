rm(list=ls())
getwd()
setwd("C:/Users/Winston/Desktop/Crucial Y4 things/Small States/Data")
install.packages("DiagrammeR")
library(DiagrammeR)
DiagrammeR::grViz("
digraph {
  graph [rankdir=LR]
  node []
    A [label = 'Population Size',shape=box]
    B [label = 'Clientelism',shape=box]
    C [label = 'Military Repression',shape=box]
    Y [label = 'Democracy',shape=box]
  edge []
    A->B->Y
    A->C->Y

}
",
engine="twopi")
