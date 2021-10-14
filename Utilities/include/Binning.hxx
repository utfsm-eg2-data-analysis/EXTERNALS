#ifndef BINNING_HXX
#define BINNING_HXX

#ifndef HEADERS_HXX
#include "Headers.hxx"
#endif

const Int_t NbinsQ2 = 6;
const Int_t NbinsNu = 6;

// binning omega
const Double_t kEdgesQ2_Omega[NbinsQ2 + 1] = {1, 1.26, 1.52, 1.94, 2.63, 3.31, 4.0};
const Double_t kEdgesNu_Omega[NbinsNu + 1] = {2.2, 2.59, 2.98, 3.36, 3.67, 3.95, 4.2};

// binning eta
const Double_t kEdgesQ2_Eta[NbinsQ2 + 1] = {1, 1.33, 1.66, 2.0, 2.5, 3.3, 4.1};
const Double_t kEdgesNu_Eta[NbinsNu + 1] = {2.2, 2.53, 2.86, 3.22, 3.58, 3.87, 4.25};

#endif
