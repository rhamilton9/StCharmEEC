// Header and settings for PYTHIA EEC Study
// Created by Ryan Hamilton.

#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include "TCanvas.h"
#include "TVectorT.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPave.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "../utils/root_draw_tools.h"

// Constants
const double pion_mass = 0.13957039;



// STAR Kinematics
//const Int_t nevent_perbin[11] = {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
//const Float_t binedge_pTHat[11][2]={{2,3},{3,4},{4,5},{5,7},{7,9},{9,11},{11,15},{15,20},{20,25},{25,35},{35,-1}};
//const Float_t crosssec_pTHat[11] ={
//  9.00581646, 1.461908221, 0.3544350863, 0.1513760388, 0.02488645725, 0.005845846143,
//  0.002304880181, 0.000342661835, 4.562988397e-05, 9.738041626e-06, 5.019978175e-07
//};

// ALICE Kinematics
const Int_t nevent_perbin[20] = {1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000};
const Float_t binedge_pTHat[20][2]={
  {5,7},{7,9},{9,12},{12,16},{16,21},{21,28},{28,36},{36,45},{45,57},{57,70},
  {70,85},{85,99},{99,115},{115,132},{132,150},{150,169},{169,190},{190,212},{212,235},{235,-1}};
const Float_t crosssec_pTHat[20] ={
  16.07, 4.608, 2.149, 0.7818, 0.2646, 9.750e-2, 2.928e-2, 9.884e-3, 4.044e-3, 1.353e-3,
  5.298e-4, 1.882e-4, 9.220e-5, 4.291e-5, 2.091e-5, 1.061e-5, 5.749e-6, 3.004e-6, 1.616e-6, 2.098e-6
};



const bool charged_jets = false;
const bool flavor_dependent_norm = false;
const Float_t etajet_latecut = 0.6; // 0.6 for STAR (1.0 eta range) and 0.5 for ALICE (0.9 eta range)

constexpr char algo_string[20] = "anti-k_{T}";    // String for jet clustering algorithm
constexpr char algo_string_short[5] = "a-kT";     // Shortened version

//// Key generator settings, used to access relevant file. Set by user.
//const Double_t sqrt_s = 200;                       // Beam energy in GeV
//const Double_t pTmin_jet = 20;                      // Minimum pT cut for jet clustering
//const Double_t jetRadius = 0.4;                     // Jet Radius -- parameter used by jet clustering algorithm

// Key generator settings, used to access relevant file. Set by user.
const Double_t sqrt_s = 200;                       // Beam energy in GeV
const Double_t pTmin_jet = 5;                      // Minimum pT cut for jet clustering
const Double_t jetRadius = 0.4;                     // Jet Radius -- parameter used by jet clustering algorithm

// Jet and hadron pT thresholds for jet clustering.
const Int_t nbins_pT = 200;
const Double_t pTmax_track = 20;
const Double_t pTmin_hadron = 1;
const Double_t pTmin_jetcore = 5;

// (Pseudo)Rapidity and Azimuth thresholds and bin counts
const Int_t nbins_rap = 100;
const Int_t nbins_phi = 100;
const Int_t nbins_rap_display = 200; // Larger for better jet area resolution
const Int_t nbins_phi_display = 314/2;
const Double_t max_y = 1;
const Double_t max_eta = 1;
// Note -- this is the track particle eta cut.
// Jet eta thresholds are calculated based on the jet radius to ensure
// all jet constituents are within the track thresholds.


// --------------------------------------------- Aesthetic settings

Int_t plot_colors[5] = {
  getTColorFromHex("#5790fc"),
  getTColorFromHex("#e42536"),
  getTColorFromHex("#189c20"),
  getTColorFromHex("#9c9ca1"),
  getTColorFromHex("#ec6672")
};

Int_t color_accent1[3] = {
  getTColorFromHex("#E23813"),
  getTColorFromHex("#AA3218"),
  getTColorFromHex("#E7B0A5")
};

#endif /* config_hpp */
