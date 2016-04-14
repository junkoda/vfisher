#ifndef VOLUME_H
#define VOLUME_H 1

#include "fisher_matrix.h"
#include <vector>

std::vector<Volume>* uniform_noise(const double vol, const double ngal_d, const double ngal_v, const double sigma_v);

//std::vector<Volume>* uniform_noise2(const double vol, const double ngal1, const double ngal2);

std::vector<Volume>* uniform_galaxy(const double steradian, const double ngal_d, const double ngal_v, const double verr_fraction, const double zmax, const int nbin);

std::vector<Volume>* nz_list(const char filename[], const double verr_fraction);

#endif
