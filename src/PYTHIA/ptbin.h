/// A file with some relevant data for pythia jet spectra generation.
///
/// To accurately simulate physical processes, the pythia flag [PhaseSpace:pTHatMin]
/// should not be fixed throughout generation, but rather should be incrementally set across
/// many generations and scaled according to the cross section in generating the final spectra.
/// The file contains some data tables for performing this splitting, both at LHC and RHIC kinematics
///
/// Created by Ryan Hamilton. Adapted from data provided by:
/// Isaac Mooney (orcid id #0000-0002-3571-0550)
/// Laura Havener (orcid id #0000-0002-4743-2885)

#ifndef ptbin_h
#define ptbin_h

// Possibly todo at some time in the future, for now seems like too much...
//enum Collider {
//  LHC = 1;
//  RHIC = 2;
//};
//
//class pTHatBins {
//  Collider kinematicScale;
//  int nbins;
//
//public:
//
//
//};

struct pTBinSettings {
  float binedge_low;
  float binedge_high;
  float cross_section;
  int nevent;
};


// Bins for RHIC kinematics (sqrt_s = 200 GeV)
const int nbins_pt_RHIC = 11;
const float binedge_pt_RHIC[nbins_pt_RHIC][2] ={
  {2,3},    {3,4},    {4,5},
  {5,7},    {7,9},    {9,11},
  {11,15},  {15,20},  {20,25},
  {25,35},  {35,-1}
};
const float xsec_pt_RHIC[nbins_pt_RHIC] = {
  9.00581646,       1.461908221,    0.3544350863,
  0.1513760388,     0.02488645725,  0.005845846143,
  0.002304880181,   0.000342661835, 4.562988397e-05,
  9.738041626e-06,  5.019978175e-07
};

const int nevent_xsec_scaled_RHIC[nbins_pt_RHIC] = {
  2100295,  600300,   600300,
  300289,   300289,   300289,
  160295,   100302,   80293,
  76303,    23307
};

// Scaled to take roughly 9h on GRACE
const int nevent_runtime_scaled_RHIC[nbins_pt_RHIC] = {
  162370,  168830,  174430,
  172350,  167030,  171291,
  173990,  172860,  172860,
  167041,  131330
};

// Bins for LHC kinematics (sqrt_s >= 900 GeV)
const int nbins_pt_LHC = 20;
const float binedge_pt_LHC[nbins_pt_LHC][2] ={
  {5,7},      {7, 9},     {9, 12},    {12, 16},   {16, 21},
  {21, 28},   {28, 36},   {36, 45},   {45, 57},   {57, 70},
  {70, 85},   {85, 99},   {99, 115},  {115, 132}, {132, 150},
  {150, 169}, {169, 190}, {190, 212}, {212, 235}, {235, -1}
};

const float xsec_pt_LHC[nbins_pt_LHC] = {
  16.07,    4.608,    2.149,    0.7818,   0.2646,
  9.750e-2, 2.928e-2, 9.884e-3, 4.044e-3, 1.353e-3,
  5.298e-4, 1.882e-4, 9.220e-5, 4.291e-5, 2.091e-5,
  1.061e-5, 5.749e-6, 3.004e-6, 1.616e-6, 2.098e-6
};

const int nevent_xsec_scaled_LHC[nbins_pt_LHC] = {
  2500000, 750000, 600000, 400000, 300000,
  200000,  150000, 130000, 115000, 105000,
  97000,   92000,  90000,  88000,  84000,
  79000,   70000,  60000,  50000,  40000
};

// Scaled to take roughly 9h when run on GRACE
const int nevent_runtime_scaled_LHC[nbins_pt_LHC] = {
  107500, 107250, 109200, 101200, 103200,
  79400,  81000,  78650,  75440,  78435,
  78376,  79856,  84510,  80960,  77080,
  74840,  76100,  73400,  76500,  77800
};

int getNbins(double sqrt_s) {
  if (sqrt_s < 900) return nbins_pt_RHIC;
  return nbins_pt_LHC;
}

pTBinSettings getBinSettings(float sqrt_s, int thread, 
                             bool use_xsec_over_runtime_scaled_nevent) {
  pTBinSettings settings;
  if (sqrt_s < 900) { // RHIC
    if (thread >= nbins_pt_RHIC) {
      std::cerr << "bad thread: out of range" << std::endl;
      return settings;
    }
    
    settings.binedge_low = binedge_pt_RHIC[thread][0];
    settings.binedge_high = binedge_pt_RHIC[thread][1];
    settings.cross_section = xsec_pt_RHIC[thread];
    if (use_xsec_over_runtime_scaled_nevent)
      settings.nevent = nevent_xsec_scaled_RHIC[thread];
    else
      settings.nevent = nevent_runtime_scaled_RHIC[thread];
  } else { // LHC
    if (thread >= nbins_pt_LHC) {
      std::cerr << "bad thread: out of range" << std::endl;
      return settings;
    }
    
    settings.binedge_low = binedge_pt_LHC[thread][0];
    settings.binedge_high = binedge_pt_LHC[thread][1];
    settings.cross_section = xsec_pt_LHC[thread];
    if (use_xsec_over_runtime_scaled_nevent)
      settings.nevent = nevent_xsec_scaled_LHC[thread];
    else
      settings.nevent = nevent_runtime_scaled_LHC[thread];
  }
  return settings;
}

#endif
