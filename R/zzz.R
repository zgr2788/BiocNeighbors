.onLoad <- function(libname, pkgname) {
    initialize_load_index_registry()
    registerLoadGenericIndexClass("knncolle_annoy::Annoy", AnnoyIndex)
    registerLoadGenericIndexClass("knncolle::Bruteforce", ExhaustiveIndex)
    registerLoadGenericIndexClass("knncolle_hnsw::Hnsw", HnswIndex)
    registerLoadGenericIndexClass("knncolle_kmknn::Kmknn", KmknnIndex)
    registerLoadGenericIndexClass("knncolle::Vptree", VptreeIndex)
    registerLoadIndexFunction(.nndescent_algorithm_name, .load_nndescent_index)
    registerLoadIndexFunction(.rpforest_algorithm_name, .load_rpforest_index)
}
