/**
 * \file binimage.h
 * \brief Handling of binary images
 */
#ifndef BINIMAGE_H
#define BINIMAGE_H
#include <stdlib.h>
#include "bm_binmat.h"
#include "bm_intmat.h"

#define CCLIP(x,a,b) ( (x) > (a) ? ( (x) < (b) ? (x) : (b) ) : (a) )

/**
 * ways of extracting patches from an image.
 */ 
typedef enum {
    EXTRACT_EXACT, ///< only extract patches which contain true pixels, possibly leaving bordering pixels out
    EXTRACT_FULL, ///< extract patches so that whole image is covered, extrapolating border pixels as needed
} extract_t;

size_t compute_grid_size(const size_t size, const size_t width, const size_t stride, extract_t extract_type);


/**
 * Extract patches from an image into a patches matrix
 * \param I input binary image
 * \param width the width (and height) of the patches
 * \param stride the distance between adjacent patches
 * \param e extraction method 
 * \param[out] P patches matrix
 * \see extract_t
 */
void extract_patches(const binary_matrix& I,
                     const size_t width,
                     const size_t stride,
                     const extract_t e,
                     binary_matrix& P);


/**
 * Stitch patches from a patches matrix into
 * a target binary image.
 * \param M number of rows in the image
 * \param N number of columns in the image
 * \param stride distance between adjacent patches
 * \param e extraction mode
 * \param[out] I output image
 * \param[out] auxiliary accumulator (integer matrix)
 * \param[out] auxiliary normalization factors (integer matrix)
 */
void stitch_patches(const binary_matrix& P,
                    const size_t M,
                    const size_t N,
                    const size_t stride,
                    const extract_t e,
                    binary_matrix& I,
                    integer_matrix& A,
                    integer_matrix& C);
#endif
