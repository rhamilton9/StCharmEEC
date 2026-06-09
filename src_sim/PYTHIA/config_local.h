// Kinematic and detector acceptance settings

#ifndef config_local_h
#define config_local_h

//======================================================== Settings specific to this study

// Flavor flag:
//  - 0: light quark (u, d, s)
//  - 1: charm only
//  - 2: beauty only
//  - 3: gluon only
//  - 4: all hard processes are on
const Int_t flavor_flag = 4;

//======================================================== Detector mode swap


// At ALICE Kinematics
//const Float_t crosssec_process[4] = { // Process cross sections
//  7.904e-05 + 5.522e-06,
//  2*( 2.631e-05 + 1.839e-06),
//  2.596e-05 + 1.826e-06,
//  1.347e-04};
//const double sqrt_s = 5020;           // Beam energy in GeV
//const double pTmin_hadron = .150;     // Minimum hadron pT to be clustered into a jet
//const double max_y = 0.9;
//const double max_eta = 0.9;           // Longitudinal acceptance


// RHIC (STAR)
const Float_t crosssec_process[4] = { // Process cross sections
  1.191e-06+1.213e-06,
  2*(3.888e-07 + 3.991e-07),
  3.232e-07 + 3.534e-07,
  1.532e-06};
const double sqrt_s = 200;            // Beam energy in GeV
const double pTmin_hadron = .200;     // Minimum hadron pT to be clustered into a jet
const double max_y = 1;
const double max_eta = 1;             // Longitudinal acceptance



//======================================================== Simulation Control Settings

// Toggle between two settings for number of events per pTHat bin
//   - (true) !! recommended for final analysis !!
//     cross-section scaled (event number scaled to cross section)
//     This setting is calibrated to generate roughly similar numbers
//     of jets in each pTHat bin.
//   - (false) recommended for early studies where one just needs a jet sample
//     runtime scaled (pTHat bins run for similar amounts of time)
//     Note that this option is calibrated to take ~3h and generates
//     highly varying numbers of jets in each bin.
const bool local_flag_xsec_scale = false;

// To upscale or downscale the number of events
// a factor larger than one will give more events, while a factor smaller will give fewer.
const double nevent_scale_factor = 1./50;


//======================================================== Jet and Event Generator Settings

// Kinematics



// Minimum accepted hadron energy in GeV (due to detector efficiency)
// Note that this depends on the detector and should be:
//    - 150 MeV at ALICE
//    - 200 MeV at STAR

// (Pseudo)Rapidity and Azimuth cuts, should be the detector acceptance

const bool cylinder_uses_rapidity = false;
// This flag controls the cylinder where jet reconstruction happens.
//    true = rapidity-phi plane, false = pseudorapidity-phi plane.
//    It's more physical to reconstruct in y-phi because y is the lorentz angle
//    However, y is difficult to reconstruct in experiemnt, so it's better to
//    reconstruct in eta/phi when making direct comparisons against experiment.


// FastJet3 (jet finder) settings
const double jet_radius = 0.4;                    // Jet Radius -- parameter used by jet clustering algorithm
constexpr char algo_string[20] = "anti-k_{T}";    // String for jet clustering algorithm
constexpr char algo_string_short[5] = "a-kT";     // Shortened version
// Other available algorithms, algo_string, and algo_string_short:
//        kt_algorithm            k_{T}         kT
//        cambridge_algorithm     cambridge     cam
//        antikt_algorithm        anti-k_{T}    a-kT
//        genkt_algorithm         gen-k_{T}[n]  g-kTn     (note n is specified by user)
//        Misc. cone algo       [not implemented or considered in the scope of this macro]

// Jet area cut
const double min_area = .0*3.14159*jet_radius*jet_radius;

// Jet/jet constituents and pT cuts for jet clustering.
const double pTmin_jetcore = 0;         // Minumum pT of the leading hadron in the jet
const double jetcut_minpT_full = 20;    // in GeV
const double jetcut_minpT_chrg = 20;    // in GeV
const double D0_minpT = 0;              // in GeV
// Note on general values for these parameters:
//    Generally 10 GeV is a good cut for jets to remove many softer gluon jets
//    Extrapolating to lower pT could be good to use a lower cut (5 GeV)
//    The jetcore cut ensures the jet has at least one piece of hard radiation (not a clump of soft particles)
//    This can be relaxed according to the desired analysis

// PYTHIA generator settings
const int pythia_seed = -1;         // Negative = default, 0 = time randomized, else 1-9E9 fixed seed
const int pythia_tune = 6;          // See https://pythia.org//latest-manual/Tunes.html#anchor3
const int times_allow_errors = 3;   // How many times to allow errors in an event before skipping it


//======================================================== Plot binning and cosmetic settings
const int nbins_pT = 200;
const int nbins_rap = 100;
const int nbins_phi = 100;
const int nbins_rap_display = 200; // Larger for better jet area resolution
const int nbins_phi_display = 314/2;
const double pTmax_hadron_eventdisplay = 20;    // Purely cosmetic, for upper range on summary sheet
const int nbins_pT_jet = 50;      // Jet plot settings that differ from constituent, track plots

#endif
