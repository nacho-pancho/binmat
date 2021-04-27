#include <cassert>
#include "gsl/gsl_randist.h"

#include "bm_types.h"

long random_seed = 34503498;



//----------------------------------------------------------------------------
gsl_rng* get_rng() {
    static gsl_rng* rng = 0;
    if (rng == 0) {
        rng = gsl_rng_alloc (gsl_rng_rand48);
        gsl_rng_set (rng, random_seed); // set random_seed
    }
    return rng;
}

//----------------------------------------------------------------------------

unsigned long get_uniform_unsigned_sample(unsigned maxval) {
    return (unsigned long) gsl_rng_uniform_int(get_rng(),maxval);
}

//----------------------------------------------------------------------------

bool get_bernoulli_sample(double p) {
    return gsl_ran_bernoulli (get_rng(),p);
}

//----------------------------------------------------------------------------
#include <iostream>

void get_subsample_indexes(const idx_t n, const idx_t nsub, idx_t* psubidx, bool with_repo) {
    if (with_repo) {
        for (idx_t i = 0; i < nsub; ++i) {
            psubidx[i++] = get_uniform_unsigned_sample(n);
        }
    } else {
        assert(n >= nsub);
        bool* used= new bool[n];
        for (idx_t i = 0; i < n; ++i) {
            used[i] = false;
        }
        for (idx_t i = 0; i  < nsub; ) {
            idx_t ri = get_uniform_unsigned_sample(n);
            if (!used[ri]) {
                used[ri] = true;
                psubidx[i++] = ri;
            }
        }
        free(used);
    }
}
