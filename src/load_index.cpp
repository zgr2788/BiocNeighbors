#include "BiocNeighbors.h"
#include "knncolle/knncolle.hpp"
#include "knncolle_kmknn/knncolle_kmknn.hpp"
#include "knncolle_hnsw/knncolle_hnsw.hpp"
#include "annoy.h"

//[[Rcpp::export(rng=false)]]
SEXP load_index(std::string prefix) {
    throw std::runtime_error(
        "loading generic C++ indices is not supported by the installed 'knncolle' headers"
    );
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_load_index_registry() {
    return R_NilValue;
}

//[[Rcpp::export(rng=false)]]
SEXP get_load_index_registry() {
    return R_NilValue;
}
