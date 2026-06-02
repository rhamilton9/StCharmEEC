// Configuration settings to be used by all code files
// These are parameters that control certain data input and output settings

#include <fstream>
#include <sys/stat.h>
#include <filesystem>
#include "../utils/root_draw_tools.h"

// --------------------------------------------- Dataset labels

// String for accessing data files in ../data.nosync
const char* dataset_label = "pp2012";

// Strings for plotting "STAR [strings]" on plot windows
const char* dataset_label_plotstring = "(2012)";
const char* collision_system = "pp";
float sqrt_s = 200; // GeV


// --------------------------------------------- D0 Reconstruction and Fitting

// Point in invariant mass m_{Kpi} above which the integrals
// of like-sign and unlike-sign will be scaled to match.
// Set -1 to make no scaling, 0 to make match overall normalization.
float mass_integral_background_matching_threshold = 2.5;

char  resonance_labels[3][10] = {"K*892","K*1430","D01864"};
float resonance_masses[3] = {0.892, 1.430, 1.8648}; // {K*0(892), K*0(1430), D0}
float mass_param_range[3] = {0.08, 0.15, 0.03};

// Regions to allow the fit to match the function
float fitting_region_limits[3][2] = {
  {0.65, 1.15}, // K* (892)
  {1.15, 1.60}, // K* (1430)
  {1.75, 1.95}  // D0 (1864)
};

// Fitting mode
// Options: TH1::default, likelihood,
char fitmode[30] = "TH1::default";


// --------------------------------------------- Aesthetic settings

char plot_colors[5][10] = {"#5790fc","#e42536","#189c20","#9c9ca1","#ec6672"};

Int_t color_accent1[3] = {
  getTColorFromHex("#E23813"),
  getTColorFromHex("#AA3218"),
  getTColorFromHex("#E7B0A5")
};
