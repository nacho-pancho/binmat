/**
 * \file extract_patches_tool.cpp
 * \brief Takes an image, a patch width, a stride, and save a matrix whose columns are the image patches
 * as columns.
 */

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include "bm_util.h"
#include "bm_binmat.h"
#include "bm_pbm.h"
#include "bm_binimage.h"

idx_t W = 16;
idx_t s = 1;

const char* iname = "data/einstein.pbm";
const char* oname = "res/einstein_patches.pbm";

void parse_args(int argc, char **argv) {
    if (argc < 4) {
        std::cerr << "USAGE: " << argv[0] << " image.pbm patch_width stride patches.pbm" << std::endl;
        std::cerr << "USAGE: " << argv[0] << " image.pbm wxxwy+sx+sy patches.pbm (imagemagick-like syntax, for non-square patches/grids" << std::endl;
        exit(1);
    }
    iname = argv[1];
    W = atoi(argv[2]);
    s = atoi(argv[3]);
    oname = argv[4];

}

int main(int argc, char **argv) {
    parse_args(argc,argv);
    //
    // input data
    //
    idx_t rows,cols;
    int res= 0;
    FILE* fimg;
    fimg = fopen(iname,"r");
    if (!fimg) { 
        std::cerr << "Error  " << res << "  while opening image " << iname  << std::endl;
        return -1;
    }
    res = read_pbm_header(fimg,rows,cols);
    if (res !=PBM_OK) {
        std::cerr << "Error " << res << " reading image header."  << std::endl;
        std::exit(1);
    }
    std::cout << "Input: " << iname << " rows=" << rows << " cols=" << cols << std::endl;
    std::cout << "Parameters: width=" << W << " stride=" << s << std::endl;
    binary_matrix I(rows,cols);
    res = read_pbm_data(fimg,I);
    if (res !=PBM_OK) {
        std::cerr << "Error " << res << " reading image data."  << std::endl;
        std::exit(1);
    }
    fclose(fimg);

    binary_matrix X;
    if ((W > rows) || (W > cols)) {
        std::cerr << "Error: patch size " << W << " does not fit in image " << std::endl;
        exit(2);
    }
    const idx_t Ny = compute_grid_size(rows,W,s,EXTRACT_EXACT);
    const idx_t Nx = compute_grid_size(cols,W,s,EXTRACT_EXACT);
    const idx_t M = W*W;
    const idx_t N = Nx*Ny;
    std::cout << "Output: " << oname << "\tM=" << M << "\tN=" << N << "\tNy=" << Ny << "\tNx=" << Nx << std::endl;
    X.allocate(N,M);
    extract_patches(I,W,s,EXTRACT_EXACT,X);
    I.destroy();
    write_pbm(X,oname);
    X.destroy();
    return 0;
}

