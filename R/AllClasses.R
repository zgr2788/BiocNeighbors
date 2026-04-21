# Defines the BiocNeighborParam class and derivatives.

#' @export
setClass("BiocNeighborParam", contains="VIRTUAL", slots=c(distance="character"))

#' @export
setClass("ExhaustiveParam", contains="BiocNeighborParam")

#' @export
setClass("KmknnParam", contains="BiocNeighborParam")

#' @export
setClass("VptreeParam", contains="BiocNeighborParam")

#' @export
setClass("AnnoyParam", contains="BiocNeighborParam", slots=c(ntrees="integer", search.mult="numeric")) 

#' @export
setClass("HnswParam", contains="BiocNeighborParam", slots=c(nlinks="integer", ef.construction="integer", ef.search="integer")) 

#' @export
setClass("NndescentParam", contains="BiocNeighborParam", slots=c(nneighbors="integer", ntrees="integer", max.candidates="integer", niter="integer"))

#' @export
setClass("RpforestParam", contains="BiocNeighborParam", slots=c(nneighbors="integer", ntrees="integer", leaf.size="integer", max.tree.depth="integer"))

#' @export
setClass("FlannKdtreeParam", contains="BiocNeighborParam", slots=c(checks="integer"))

#' @export
setClass("BiocNeighborIndex", contains="VIRTUAL")

#' @export
setClass("BiocNeighborGenericIndex", contains=c("VIRTUAL", "BiocNeighborIndex"), slots=c(ptr="externalptr", names="ANY"))

#' @export
setClass("ExhaustiveIndex", contains="BiocNeighborGenericIndex")

#' @export
setClass("KmknnIndex", contains="BiocNeighborGenericIndex")

#' @export
setClass("VptreeIndex", contains="BiocNeighborGenericIndex")

#' @export
setClass("AnnoyIndex", contains="BiocNeighborGenericIndex")

#' @export
setClass("HnswIndex", contains="BiocNeighborGenericIndex")

#' @export
setClass("NndescentIndex", contains="BiocNeighborIndex", slots=c(index="list", data="ANY", names="ANY", param="NndescentParam"))

#' @export
setClass("RpforestIndex", contains="BiocNeighborIndex", slots=c(index="list", data="ANY", names="ANY", param="RpforestParam"))

#' @export
setClass("FlannKdtreeIndex", contains="BiocNeighborIndex", slots=c(data="ANY", names="ANY", param="FlannKdtreeParam"))
