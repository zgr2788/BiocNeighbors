.onLoad <- function(libname, pkgname) {
    initialize_load_index_registry()
    registerLoadGenericIndexClass("knncolle_annoy::Annoy", AnnoyIndex)
    registerLoadGenericIndexClass("knncolle::Bruteforce", ExhaustiveIndex)
    registerLoadGenericIndexClass("knncolle_hnsw::Hnsw", HnswIndex)
    registerLoadGenericIndexClass("knncolle_kmknn::Kmknn", KmknnIndex)
    registerLoadGenericIndexClass("knncolle::Vptree", VptreeIndex)
    registerLoadIndexFunction(.flann_kdtree_algorithm_name, .load_flann_kdtree_index)
    registerLoadIndexFunction(.flann_kmeans_algorithm_name, .load_flann_kmeans_index)
    registerLoadIndexFunction(.mlpack_lsh_algorithm_name, .load_mlpack_lsh_index)
    registerLoadIndexFunction(.nndescent_algorithm_name, .load_nndescent_index)
    registerLoadIndexFunction(.rpforest_algorithm_name, .load_rpforest_index)
}
