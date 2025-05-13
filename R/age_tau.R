#code to estimate tau from age proportion data
tau_comp <- function(ppnMat){
  mvLogisticLL(agePropData=ppnMat)
}