/**
 * \file bm_util.h
 * \brief Typical assortement of different miscelaneous functions.
 */
#ifndef UTIL_H
#define UTIL_H
#include "bm_binmat.h"

typedef std::pair<idx_t,idx_t> aux_t;


void inc_verbosity();

char get_verbosity();

typedef bool (*filter_f)(const binary_matrix& x);
/**
 * Return a matrix with the rows i of X for which keep[i] != 0
 * @param X source matrix
 * @param functor which determines whether a given row should be kept 
 */
binary_matrix filter(const binary_matrix& X, const binary_matrix& indicator);

/**
 * true if the matrix/vector is not too close to all ones or all zeros 
 */
binary_matrix filter_by_weight(const binary_matrix& X, const idx_t thres);

/**
 * Sample rows of X, with or without reposition, and store them
 * in Xsub. The number of samples is determined by the size of Xsub,
 * which mmust be already allocated.
 * @param X source matrix
 * @param Xsub destination matrix 
 * @param ppsubidx if not NULL, use previously sampled indexes; otherwise sample new ones
 */
void subsample(const binary_matrix& X, const idx_t* pidx, binary_matrix& Xsub);


/**
 * Saves an image where each column of D is transformed into a square and shown
 * as a tile.
 * @param D matrix to be rendered as a mosaic of square tiles
 * @param fname file name of image to write to
 */
void render_mosaic(const  binary_matrix& D, const char* fname, const idx_t gm = 0, const idx_t gn = 0);

/**
 * utility for sorting according to counts
 */
void counting_sort(aux_t* s, idx_t n);


/**
 * randomize a matrix using a pseudo-random Bernoulli(p) process 
 */
void randomize(binary_matrix& X, double p);

/**
 * write a binary matrix A as a PBM image of file name fname
 * @param A matrix to be saved as image
 * @param fname file name of image to write to
 */
int write_pbm(binary_matrix& A, const char* fname);

#endif

