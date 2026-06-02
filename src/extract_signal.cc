// Some simple code to estimate signal yields from {K*,D0} -> K/pi reconstruction
// Inputs are histograms obtained from RCF
// 
// ----------------- Changelog -----------------
//  - 4/16/2026  ::  Created (from D0 reconstruction test study code)
//  - 4/23/2026  ::  Overhauled directory interactions to automate bookkeeping
//  - 5/09/2026  ::  Got fitting roughly working
//  - 5/26/2026  ::  Improved fitting functonality and got rough yield spectra


// TODO move TObjects to Reader class and plot by referencing the objects in the class


// Inclusions
#include "config.h"

// Compiler flags
//#define include_MC

//========================================================================== Global Variables

// Global canvas pointers
TCanvas* gCanvas_solo;
//TCanvas* gCanvas_subsolo;
TCanvas* gCanvas_quad;


//========================================================================== Forward Declarations

void plotSpectra();


//========================================================================== Main

// Generate some plots from produced D0->Kpi reconstructed histograms
//
// The input to the method is the modifier on the file wich should be
//         ../data.nosync/test_pp2012_[modifier].root
void extract_signal(const char *modifier = "fullsample") {
  
  bool _debug = false;
  
  // Setup strings for directory accessing
  char directory_modifier[50];
  char filename_modifier[50];
  if (modifier[0] == '\0') {
    snprintf(directory_modifier, 50, "%s", "nullmodifier");
    snprintf(filename_modifier, 50, "%s", modifier);
  } else {
    snprintf(directory_modifier, 50, "%s", modifier);
    snprintf(filename_modifier, 50, "_%s", modifier);
  }
  
  
  
  // Check that the input file exists
  TFile* datafile = new TFile(Form("../data.nosync/test_%s%s.root",dataset_label, filename_modifier));
  if (datafile->IsZombie()) {
    std::cout << "Data file <\033[31m" << Form("../data.nosync/test_%s%s.root",dataset_label, filename_modifier) << "\033[39m>";
    std::cout << " was not found!\nCheck the given inputs and try again." << std::endl;
    return;
  }
  
  
  // Check if the output plot directory structure exists for this modifier, and create it if not
  const int n_directories = 2;
  char list_directories[n_directories][50] = {
    "reco_spectra", "diff_spectra"
  };
  const int n_subdirectories = 7;
  char list_subdirectories[n_directories][n_subdirectories][50] = {
    {"reco_D0", "reco_D0bar", "reco_D0mixed", "reco_all", "", "", ""},
    {"diff_D0", "diff_D0bar", "diff_D0mixed", "diff_all", "fit1_K892", "fit2_K1430", "fit3_D01864"}
  };
  
  struct stat metadata;
  char base_directory[100];
  snprintf(base_directory, 100, "../plots/%s", directory_modifier);
  if (stat(base_directory, &metadata) != 0) { // no base directory for this modifier!
    std::filesystem::create_directory(base_directory);
    std::cout << "Created base directory for new modifier \033[32m" << directory_modifier << "\033[39m." << std::endl;
  }
  
  for (int i_dir = 0; i_dir < n_directories; ++i_dir) {
    char c_directory[100];
    snprintf(c_directory, 100, "%s/%s", base_directory, list_directories[i_dir]);
    if (stat(c_directory, &metadata) != 0) std::filesystem::create_directory(c_directory);
    
    for (int i_subdir = 0; i_subdir < n_subdirectories; ++i_subdir) {
      if (list_subdirectories[i_dir][i_subdir][0] == '\0') continue;
      
      char c_subdirectory[150];
      snprintf(c_subdirectory, 150, "%s/%s", c_directory, list_subdirectories[i_dir][i_subdir]);
      if (stat(c_subdirectory, &metadata) != 0) std::filesystem::create_directory(c_subdirectory);
    }// End of subdirectory loop
  }// End of directory loop
  
  
  
  // Set up TCanvas objects
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  suppress_warnings_draw_tools = true;
  
  gCanvas_solo = new TCanvas("gCanvas_solo");
  gPad->SetTicks(1,1);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.12);
  
//  gCanvas_subsolo = new TCanvas("gCanvas_subsolo");
//  gPad->SetTicks(1,1);
//  gPad->SetRightMargin(0.03);
//  gPad->SetLeftMargin(0.12);
  
  gCanvas_quad = new TCanvas("gCanvas_quad");
  gCanvas_quad->SetCanvasSize(800, 600);
  gCanvas_quad->Divide(2,2);
  for (int i = 0; i < 4; ++i) {
    gCanvas_quad->cd(i+1);
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(0.005);
    gPad->SetLeftMargin(0.115);
    gPad->SetTopMargin(0.09);
    gPad->SetBottomMargin(0.08);
  }// End quad canvas setup
  
  // Colors from input hex codes
  const int n_colors = sizeof(plot_colors) / sizeof(plot_colors[0]);
  Int_t plot_TColors[n_colors];
  for (int c = 0; c < n_colors; ++c) plot_TColors[c] = getTColorFromHex(plot_colors[c]);
  
  // Misc strings/settings for drawing
  const double axis_scale_factor = 1e-3;
  char D0string[3][10] = {"D0", "D0bar","merged"};
  char daughterstring_long[3][35] = {"K^{+}#pi^{-}", "K^{-}#pi^{+}","(K^{+}#pi^{-}) + (K^{-}#pi^{+})"};
  char resonance_string_long[3][50] = {"K*(892)", "[K*_{0}+K*_{2}](1430)","D^{0}"};
  
  // Store some plots to make ratios
//  std::vector<TH1F*> background_dist_plus;
//  std::vector<TH1F*> background_dist_minus;
//  std::vector<TLatex*> pT_bin_tex;
  
  
  // Fit settings
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1500);
  
  
  // *--------------------------- K/pi Spectra
  
  // Access the D0 data histograms and relevant binnings
  TH2F* hist_2d_D0[2];
  hist_2d_D0[0] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hD0RecoMass");
  hist_2d_D0[1] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hD0barRecoMass");
  TH2F* hist_2d_SS[2];
  hist_2d_SS[0] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hPosSameSignMass");
  hist_2d_SS[1] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hNegSameSignMass");
  TH2F* hist_2d_MC[2];
  hist_2d_MC[0] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hD0TruthMass");
  hist_2d_MC[1] = datafile->Get<TDirectoryFile>("D0")->Get<TH2F>("hD0barTruthMass");
  
  std::vector<double> pt_binedge;
  const int nbins_pt = hist_2d_D0[0]->GetXaxis()->GetNbins();
  const int nbins_mass = hist_2d_D0[0]->GetYaxis()->GetNbins();
  for (int i_pt = 1; i_pt <= nbins_pt + 1; ++ i_pt)
    pt_binedge.push_back(hist_2d_D0[0]->GetXaxis()->GetBinLowEdge(i_pt));
  double massrange[2] = {
    hist_2d_D0[0]->GetYaxis()->GetBinLowEdge(1),
    hist_2d_D0[1]->GetYaxis()->GetBinLowEdge(nbins_mass+1)
  };
  
  
  
  // Specify how (if at all) pT bins should be merged
  // use 0 to specify no further additions
  const int n_final_pt_bins = 10;
  float bin_range[n_final_pt_bins][2] = {
    {0.0,-1.0}, // All data (-1 means all above)
    {0.6, 2.2},
    {0.0, 0.6},
    {0.6, 1.0},
    {1.0, 1.5},
    {1.5, 2.0},
    {2.0, 2.2},
    {2.2, 3.5},
    {3.5, 4.5},
    {4.5,-1.0}  // Note that the last pT bin includes all candidates with higher pT
    // (i.e. should be thought of as the [x, infinity) bin).
  };
  char bin_range_string[n_final_pt_bins][2][20];
  // Arrays to keep track of bin indexing
  std::vector<int>* bin_mergers[n_final_pt_bins];
  
  
  // Store the fit parameters
  // 3 resonances each for D0, D0bar and Mixed in each pT bin
  // Each fit has 6 parameters
  const int n_params = 6;
  double fit_params[n_final_pt_bins][3][3][n_params];
  double fit_params_error[n_final_pt_bins][3][3][n_params];
  
  
  // Begin loop over pT bins to store the various plots of interest
  for (int i_pt_f = 0; i_pt_f < n_final_pt_bins; ++i_pt_f) {
    
    // Find the pT bins to include from the smaller histogram
    bin_mergers[i_pt_f] = new std::vector<int>();
    float pt_binrange_actual[2] = {INT_MAX, 0};
    for (int i_ptbin_cand = 0; i_ptbin_cand < nbins_pt; ++i_ptbin_cand) {
      // Check if the bin is in the range and if it's a new bin edge
      //
      // (note that bins are assumed to be continuous, and if a value is
      // provided within the center of a bin, the larger bin range will be
      // taken to include the desired range in a binned setting.)
      //
      // The small addition of 1e-6 is to prevent errors due float rounding
      if (pt_binedge[i_ptbin_cand  ] >= bin_range[i_pt_f][1] - 1e-6 && bin_range[i_pt_f][1] > 0) continue; // Low edge too high
      if (pt_binedge[i_ptbin_cand+1] <= bin_range[i_pt_f][0] + 1e-6)                             continue; // High edge too low
      
      // Append the bin to the list
      bin_mergers[i_pt_f]->push_back(i_ptbin_cand+1);
      if (pt_binedge[i_ptbin_cand  ] < pt_binrange_actual[0]) pt_binrange_actual[0] = pt_binedge[i_ptbin_cand  ];
      if (pt_binedge[i_ptbin_cand+1] > pt_binrange_actual[1]) pt_binrange_actual[1] = pt_binedge[i_ptbin_cand+1];
    }// End of pT range finding
    
    
    // Store strings for this pT bin
    if (bin_range[i_pt_f][0] >= 0) snprintf(bin_range_string[i_pt_f][0], 20, "%.1f", pt_binrange_actual[0]);
    else                           snprintf(bin_range_string[i_pt_f][0], 20, "0.0");
    if (bin_range[i_pt_f][1] >= 0) snprintf(bin_range_string[i_pt_f][1], 20, "%.1f", pt_binrange_actual[1]);
    else                           snprintf(bin_range_string[i_pt_f][1], 20, "inf");
    
    if (_debug) {
      std::cout << "Processing bin #" << i_pt_f << " with pT range :: [";
      std::cout << bin_range_string[i_pt_f][0] << ',' << bin_range_string[i_pt_f][1] << "]\n";
    }
    
    
    // Setup file for writing objects
    TFile* spectra_file = new TFile(Form("%s/%s/spectra_histograms.root",
                                         base_directory,list_directories[0]), "recreate");
    
    
    // Initialize histograms for this pT_f bin
    TH1F* hist_1d_proj_D0[3]; // D0 candidates--opposite signed signal
    hist_1d_proj_D0[0] = new TH1F(Form("hist_proj_ptbin%s-%s_D0",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_D0[1] = new TH1F(Form("hist_proj_ptbin%s-%s_D0bar",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_D0[2] = new TH1F(Form("hist_proj_ptbin%s-%s_D0+D0bar",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    
    TH1F* hist_1d_proj_SS[3]; // Like-Sign Background
    hist_1d_proj_SS[0] = new TH1F(Form("hist_proj_ptbin%s-%s_SS(++)",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]",nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_SS[1] = new TH1F(Form("hist_proj_ptbin%s-%s_SS(--)",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]",nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_SS[2] = new TH1F(Form("hist_proj_ptbin%s-%s_SS(++)+(--)",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]",nbins_mass,
                                  massrange[0],massrange[1]);
    
    TH1F* hist_1d_proj_MC[3]; // Geant-simulated MC (i.e. reconstructable D0s from MC)
    hist_1d_proj_MC[0] = new TH1F(Form("hist_proj_ptbin%s-%s_MC_D0",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_MC[1] = new TH1F(Form("hist_proj_ptbin%s-%s_MC_D0bar",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    hist_1d_proj_MC[2] = new TH1F(Form("hist_proj_ptbin%s-%s_MC_D0+D0bar",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];dN/dm_{K#pi} [10^{3}/GeV]", nbins_mass,
                                  massrange[0],massrange[1]);
    
    
    
    // Loop over invariant mass bins, gathering integral over the rescaling range to match
    float integral_D0_hist[3] = {0,0,0};
    float integral_SS_hist[3] = {0,0,0};
    float yield_MC[3] = {0,0,0};
    for (int i_pt = 0; i_pt < bin_mergers[i_pt_f]->size(); ++i_pt) {
      if (bin_mergers[i_pt_f]->at(i_pt) == 0) continue;
      
      for (int i_mass = 1; i_mass <= nbins_mass; ++i_mass) {
        
        // D0 hist fills
        hist_1d_proj_D0[0]->AddBinContent(i_mass, hist_2d_D0[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        hist_1d_proj_D0[1]->AddBinContent(i_mass, hist_2d_D0[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        hist_1d_proj_D0[2]->AddBinContent(i_mass,(hist_2d_D0[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) +
                                                  hist_2d_D0[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) ) ); // Merged, sum of both
        
        // Same sign hist fills
        hist_1d_proj_SS[0]->AddBinContent(i_mass, hist_2d_SS[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        hist_1d_proj_SS[1]->AddBinContent(i_mass, hist_2d_SS[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        hist_1d_proj_SS[2]->AddBinContent(i_mass,(hist_2d_SS[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) +
                                                  hist_2d_SS[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) ) ); // Merged, sum of both
        
        // MC hist fills
        hist_1d_proj_MC[0]->AddBinContent(i_mass, hist_2d_MC[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        hist_1d_proj_MC[1]->AddBinContent(i_mass, hist_2d_MC[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass));
        yield_MC[0] += hist_2d_MC[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass);
        yield_MC[1] += hist_2d_MC[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass);
        yield_MC[2] +=(hist_2d_MC[0]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) +
                       hist_2d_MC[1]->GetBinContent(bin_mergers[i_pt_f]->at(i_pt), i_mass) );
        
        // Gather integral above threshold for normalizing Like Sign/Unlike Sign
        if (mass_integral_background_matching_threshold >= 0 &&
            hist_2d_D0[0]->GetYaxis()->GetBinLowEdge(i_mass) <= mass_integral_background_matching_threshold) continue;
        for (int i_hist = 0; i_hist < 3; ++i_hist) {
          integral_D0_hist[i_hist] += hist_1d_proj_D0[i_hist]->GetBinWidth(i_mass) * hist_1d_proj_D0[i_hist]->GetBinContent(i_mass);
          integral_SS_hist[i_hist] += hist_1d_proj_SS[i_hist]->GetBinWidth(i_mass) * hist_1d_proj_SS[i_hist]->GetBinContent(i_mass);
        }// End of hist loop
        
      }// End of mass bin loop
    }// End of original pT bin loop
    
    
    // Renormalize as requested by threshold
    const bool has_MC = (yield_MC[2] != 0);
    float scale_D0[3] = {1.0, 1.0, 1.0};
    float scale_SS[3] = {1.0, 1.0, 1.0};
    for (int i_hist = 0; i_hist < 3; ++i_hist) {
      // Rescale if requested
      if (mass_integral_background_matching_threshold >= 0)
        scale_SS[i_hist] = scale_SS[i_hist]*integral_D0_hist[i_hist] / integral_SS_hist[i_hist];
      
      hist_1d_proj_D0[i_hist]->Scale(axis_scale_factor * scale_D0[i_hist], "width");
      hist_1d_proj_SS[i_hist]->Scale(axis_scale_factor * scale_SS[i_hist], "width");
      hist_1d_proj_MC[i_hist]->Scale(axis_scale_factor * scale_D0[i_hist], "width");
    }// End of rescaling hists
    
    
    // Setup elements for plotting
    for (int i_hist = 0; i_hist < 3; ++i_hist) {
      hist_1d_proj_D0[i_hist]->SetLineColor(plot_TColors[0]);
      hist_1d_proj_SS[i_hist]->SetLineColor(plot_TColors[1]);
      hist_1d_proj_MC[i_hist]->SetLineColor(plot_TColors[3]);
      
      hist_1d_proj_D0[i_hist]->GetYaxis()->SetRangeUser(0, hist_1d_proj_D0[i_hist]->GetMaximum() * 1.4);
    }// End of aesthetic settings
    
    TLine* D0line = new TLine();
    D0line->SetLineColor(kGray + 2);
    D0line->SetLineStyle(3);
    
    
    // Make plots of reconstructed spectra
    for (int i_hist = 0; i_hist < 3; ++i_hist) {
      for (int bSolo = 0; bSolo < 2; ++bSolo) {
        if (bSolo == 0) gCanvas_solo->cd();
        else            gCanvas_quad->cd(i_hist+1);
        
        // Plot spectra
        hist_1d_proj_D0[i_hist]->Draw("hist");
        hist_1d_proj_SS[i_hist]->Draw("hist same");
        if (has_MC) hist_1d_proj_MC[i_hist]->Draw("hist same");
        
        // Draw resonance mass lines
        for (int i_res = 0; i_res < 3; ++i_res) {
          D0line->DrawLine(resonance_masses[i_res], 0, resonance_masses[i_res],
                           hist_1d_proj_D0[i_hist]->GetBinContent(hist_1d_proj_D0[i_hist]->FindBin(resonance_masses[i_res])));
        }// End of resonance mass line drawing
        
        // Draw a legend
        float height_base_leg = 0.7;
        if (has_MC)height_base_leg = 0.67;
        TLegend* leg = new TLegend(0.7, height_base_leg, 0.95, 0.85);
        leg->SetLineWidth(0);
        leg->AddEntry(D0line, "PDG K^{+}#pi^{-} Res.", "l");
        switch (i_hist) {
          case 0: // D0 -> K- pi+
            leg->AddEntry(hist_1d_proj_D0[i_hist], Form("D^{0}#rightarrow#pi^{+}K^{-} (#times%.2f)",scale_D0[i_hist]), "l");
            leg->AddEntry(hist_1d_proj_SS[i_hist], Form("LS #pi^{+}K^{+} (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC D^{0} (#times%.2f)",scale_D0[i_hist]), "l");
            break;
          case 1: // D0bar -> K+ pi-
            leg->AddEntry(hist_1d_proj_D0[i_hist], Form("#bar{D^{0}}#rightarrow#pi^{-}K^{+} (#times%.2f)",scale_D0[i_hist]), "l");
            leg->AddEntry(hist_1d_proj_SS[i_hist], Form("LS #pi^{-}K^{-} (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC #bar{D^{0}} (#times%.2f)",scale_D0[i_hist]), "l");
            break;
          case 2: // D0 + D0bar
            leg->AddEntry(hist_1d_proj_D0[i_hist], Form("D^{0} + #bar{D^{0}} (#times%.2f)",scale_D0[i_hist]), "l");
            leg->AddEntry(hist_1d_proj_SS[i_hist], Form("Total LS (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC (#times%.2f)",scale_D0[i_hist]), "l");
            break;
        }// end of legend setup
        leg->Draw();
        
        // Draw label text
        drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.15, 0.83);
        drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.15, 0.775);
        if (mass_integral_background_matching_threshold >= 0)
          drawText(Form("Scaled to Integral in m_{K#pi} > %.2f GeV",mass_integral_background_matching_threshold), 0.15, 0.72);
        drawText(Form("%s < #it{p}_{T}^{D^{0} Reco} < %s",
                      bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]), 0.93, height_base_leg - 0.05, true);
        
        drawText(Form("Modifier #it{%s}",modifier),
                 1.0 - gPad->GetRightMargin(),
                 1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
        
        // Export the solo canvas
        if (!bSolo) {
          gCanvas_solo->SaveAs(Form("%s/%s/%s/%s_%s_invmass_pt_%s_%s.pdf",
                                   base_directory,
                                   list_directories[0],
                                   list_subdirectories[0][i_hist],
                                   dataset_label,
                                   list_subdirectories[0][i_hist],
                                   bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]) );
        }// End of plot save
        
      }// End of solo plot loop
    }// End of plotting spectra
    
    // Export composite canvas
    gCanvas_quad->SaveAs(Form("%s/%s/%s/%s_recoD0_invmass_pt_%s_%s.pdf",
                             base_directory,
                             list_directories[0],
                             list_subdirectories[0][3],
                             dataset_label,
                             bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1] ) );
    
    // Write histograms to file
    spectra_file->cd();
    for (int i_hist = 0; i_hist < 3; ++i_hist) {
      hist_1d_proj_D0[i_hist]->Write();
      hist_1d_proj_SS[i_hist]->Write();
      hist_1d_proj_MC[i_hist]->Write();
    }
//    spectra_file->Close();
    // Done with K/pi spectra
    
    
    
    // *--------------------------- Difference Signal - background
    
    // Set up histograms
    TH1F* hist_1d_diff_SS[3];
    hist_1d_diff_SS[0] = new TH1F(Form("hist_diff_ptbin%s-%s_D0-SS",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];Difference (dN/dm_{K^{-}#pi^{+}} - dN/dm_{K^{+}#pi^{+}}) [10^{3}/GeV]",
                                  nbins_mass,massrange[0],massrange[1]);
    hist_1d_diff_SS[1] = new TH1F(Form("hist_diff_ptbin%s-%s_D0bar-SS",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];Difference (dN/dm_{K^{+}#pi^{-}} - dN/dm_{K^{-}#pi^{-}}) [10^{3}/GeV]",
                                  nbins_mass,massrange[0],massrange[1]);
    hist_1d_diff_SS[2] = new TH1F(Form("hist_diff_ptbin%s-%s_US-SS",bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                  ";m_{K#pi} [GeV];Difference (dN/dm_{K^{#pm}#pi^{#mp}} - dN/dm_{K^{#pm}#pi^{#pm}}) [10^{3}/GeV]",
                                  nbins_mass,massrange[0],massrange[1]);
    // TODO add rotation BG estimation
    
    
    // Compute differences unlike_sign - like_sign
    for (int i_mass = 1; i_mass <= nbins_mass; ++i_mass) {
      for (int i_hist = 0; i_hist < 3; ++i_hist) {
        hist_1d_diff_SS[i_hist]->SetBinContent(i_mass,(hist_1d_proj_D0[i_hist]->GetBinContent(i_mass) -
                                                       hist_1d_proj_SS[i_hist]->GetBinContent(i_mass) ) );
        
        // TODO add rotation BG estimation
        
      }// End of histogram loop
    }// End of invariant mass bin loop
    
    
    
    // Store the resonance-localized histograms for later plotting
    TH1F* diff_to_fit[3][3];
    
    // Plot Difference along with fits
    for (int i_hist = 0; i_hist < 3; ++i_hist) {
      for (int bSolo = 0; bSolo < 2; ++bSolo) {
       
        // Draw the difference spectra
        hist_1d_diff_SS[i_hist]->SetLineColor(plot_TColors[1]);
        hist_1d_diff_SS[i_hist]->GetYaxis()->SetRangeUser((hist_1d_diff_SS[i_hist]->GetMinimum() < 0)*(hist_1d_diff_SS[i_hist]->GetMinimum() * 1.05),
                                                          hist_1d_diff_SS[i_hist]->GetMaximum()*1.05);
        if (bSolo == 0) gCanvas_solo->cd();
        else            gCanvas_quad->cd(i_hist+1);
        hist_1d_diff_SS[i_hist]->Draw("hist");
        if (has_MC) hist_1d_proj_MC[i_hist]->Draw("hist same");
        
        
        
        
        // *================= Perform fitting and estimate signal
        TF1* final_fit_form;
        for (int i_fit = 0; i_fit < 3; ++i_fit) {
          // Choose which difference to fit against (SS or Rotated)
          TH1F* diff_to_fit_parent = hist_1d_diff_SS[i_hist];
          double value_at_peak = diff_to_fit_parent->GetBinContent(diff_to_fit_parent->FindBin(resonance_masses[i_fit]));
          
          
          // Construct a histogram over the limited range with the relevant data
          int binrange[2] = {0,0};
          float new_binedge[2] = {0,0};
          for (int i_bin_parent = 1; i_bin_parent <= diff_to_fit_parent->GetXaxis()->GetNbins(); ++i_bin_parent) {
            
            // Lower bin is the first to contain content in the range
            if (binrange[0] == 0 && diff_to_fit_parent->GetXaxis()->GetBinUpEdge(i_bin_parent) > fitting_region_limits[i_fit][0]) {
              binrange[0] = i_bin_parent;
              new_binedge[0] = diff_to_fit_parent->GetXaxis()->GetBinLowEdge(i_bin_parent);
              //              std::cout << "Lower range " << new_binedge[0]<< " found at parent bin #" << i_bin_parent << std::endl;
            }
            
            // Upper bin is the last to contain content in the range
            if (diff_to_fit_parent->GetXaxis()->GetBinUpEdge(i_bin_parent) >= fitting_region_limits[i_fit][1]) {
              binrange[1] = i_bin_parent;
              new_binedge[1] = diff_to_fit_parent->GetXaxis()->GetBinUpEdge(i_bin_parent);
              //              std::cout << "Upper range " << new_binedge[1] << " found at parent bin #" << i_bin_parent << std::endl;
              break;
            }
          }// End of bin range finding
          
          //          std::cout << "Setting new hist over bin range [" << new_binedge[0] << ',' << new_binedge[1] << "]." << std::endl;
          if (bSolo == 0)
            diff_to_fit[i_hist][i_fit] = new TH1F(Form("fit_%s_%s_pTbin%i_res%i",D0string[i_hist],resonance_labels[i_fit],i_pt_f, i_fit),"",
                                                  1+binrange[1]-binrange[0],new_binedge[0],new_binedge[1]);
          diff_to_fit[i_hist][i_fit]->SetLineColor(diff_to_fit_parent->GetLineColor());
          diff_to_fit[i_hist][i_fit]->SetLineStyle(diff_to_fit_parent->GetLineStyle());
          diff_to_fit[i_hist][i_fit]->GetXaxis()->SetTitle(diff_to_fit_parent->GetXaxis()->GetTitle());
          diff_to_fit[i_hist][i_fit]->GetYaxis()->SetTitle(diff_to_fit_parent->GetYaxis()->GetTitle());
          for (int i_bin_fit = binrange[0]; i_bin_fit <= binrange[1]; ++i_bin_fit) {
            diff_to_fit[i_hist][i_fit]->SetBinContent(i_bin_fit - binrange[0] + 1, diff_to_fit_parent->GetBinContent(i_bin_fit));
          }// End of new hist construction
          
          
          
          // Perform the fitting
          TFitResultPtr fit_TH1;
          if (strcmp(fitmode, "TH1::default") == 0) { // *-- Fit using TH1::Fit()
            
            // Construct the fitting form: Gaussian + 2nd order polynomial
            final_fit_form = new TF1("fit_form",
                                     "[0] * exp(-0.5*((x-[1])/[2])**2)/sqrt(2*TMath::Pi()*[2]**2) + [3] + [4]*(x-[1]) + [5]*(x-[1])**2",
                                     fitting_region_limits[i_fit][0], fitting_region_limits[i_fit][1]);
            final_fit_form->SetLineColor(kBlack);
            final_fit_form->SetLineWidth(1);
              
            final_fit_form->SetParLimits(0,0,2*value_at_peak);
            if (true) // Allow some leniency in the peak position
              final_fit_form->SetParLimits(1,
                                           resonance_masses[i_fit] - mass_param_range[i_fit],
                                           resonance_masses[i_fit] + mass_param_range[i_fit]);
            else // Fix peak poisition at the PDG mass
              final_fit_form->FixParameter(1, resonance_masses[i_fit]);
            final_fit_form->SetParLimits(2, 0, 2*mass_param_range[i_fit]);
//            final_fit_form->SetParLimits(3, -1000, 1000);
//            final_fit_form->SetParLimits(4, -1000, 1000);
//            final_fit_form->SetParLimits(5, -1000, 1000);
            final_fit_form->SetParameters(diff_to_fit[i_hist][i_fit]->GetBinContent(diff_to_fit[i_hist][i_fit]->FindBin(resonance_masses[i_fit])),
                                          resonance_masses[i_fit], mass_param_range[i_fit],
                                          0, 0, 0);
            
            fit_TH1 = diff_to_fit[i_hist][i_fit]->Fit("fit_form", "SLBQ0");
            
            // Store the fit parameters
            for (int i_param = 0; i_param < n_params; ++i_param) {
              fit_params[i_pt_f][i_hist][i_fit][i_param] = final_fit_form->GetParameter(i_param);
              fit_params_error[i_pt_f][i_hist][i_fit][i_param] = final_fit_form->GetParError(i_param);
            }
            
            // End of fitting with TH1::Default
          } else if (strcmp(fitmode, "likelihood") == 0) { // *-- Likelihood fit using ROOT::Fit::Fitter
            
            // Fit limits and data
            ROOT::Fit::DataOptions dat_options;
            dat_options.fIntegral = true;  // Set comparison to bin integral over bin center
            ROOT::Fit::DataRange dat_range(fitting_region_limits[i_fit][0], fitting_region_limits[i_fit][1]);
            ROOT::Fit::BinData bindata(dat_options, dat_range);
            ROOT::Fit::FillData(bindata, diff_to_fit[i_hist][i_fit]);
            
            // Set up the fit form
            TF1 func_gauss("fit_form",
                           "[0] * exp(-0.5*((x-[1])/[2])**2) + [3] + [4]*x + [5]*(x**2)",
                           fitting_region_limits[i_fit][0], fitting_region_limits[i_fit][1]);
            ROOT::Math::WrappedMultiTF1 fit_form(func_gauss, 1);
            ROOT::Fit::Fitter fitter;
            fitter.SetFunction(fit_form, false);
            
            // Parameter limits
            fitter.Config().ParSettings(1).SetLimits(resonance_masses[i_fit] - mass_param_range[i_fit],
                                                     resonance_masses[i_fit] - mass_param_range[i_fit]);
            fitter.Config().ParSettings(2).SetLimits(0, 1.0);
            
            // Initial values in the parameter search
            fitter.Config().ParSettings(0).SetValue(1);
            fitter.Config().ParSettings(1).SetValue(resonance_masses[i_fit]);
            fitter.Config().ParSettings(2).SetValue(mass_param_range[i_fit]);
            fitter.Config().ParSettings(3).SetValue(0);
            fitter.Config().ParSettings(4).SetValue(0);
            fitter.Config().ParSettings(5).SetValue(0);
            
            // Perform the fit
            fitter.LikelihoodFit(bindata);
            ROOT::Fit::FitResult fit_result = fitter.Result();
            
            std::vector<double> params = fit_result.Parameters();
            final_fit_form = &func_gauss;
            if (true) {
              std::cout << "Via manual fitting:" << std::endl;
              std::cout << "[0] = " << params.at(0) << ", [1] = " << params.at(1) << ", [2] = " << params.at(2) << std::endl;
            }
            
            // End of fitting with Log Likelihood fit
          } else {
            std::cout << "\033[31mError\033[39m :: Fit mode input not recognized! Check config file." << std::endl;
            break;
          }// End of fit for this resonance
          
          
          
          // Draw the fits on the diff plot panel
          if (bSolo == 0) gCanvas_solo->cd();
          else            gCanvas_quad->cd(i_hist+1);
          final_fit_form->Draw("same");
          
          
        }// Done with fitting all resonances.
        
        
        // *--- Finish plotting the composite diff plot for this pT bin and parent type {D0, D0bar, mixed}.
        
        
        // Draw resonance mass lines
        for (int i_res = 0; i_res < 3; ++i_res) {
          D0line->DrawLine(resonance_masses[i_res], 0, resonance_masses[i_res],
                           hist_1d_diff_SS[i_hist]->GetBinContent(hist_1d_diff_SS[i_hist]->FindBin(resonance_masses[i_res])));
        }// End of resonance mass line drawing
        
        // Draw zero line
        TLine* zeroline = new TLine();
        zeroline->SetLineColor(kGray+2);
        zeroline->DrawLine(massrange[0], 0, massrange[1], 0);
        
        // Draw a legend
        float height_base_leg = 0.7;
        if (has_MC)height_base_leg = 0.67;
        TLegend* leg = new TLegend(0.65, height_base_leg, 0.95, 0.85);
        leg->SetLineWidth(0);
        leg->AddEntry(D0line, "PDG K^{+}#pi^{-} Res.", "l");
        switch (i_hist) {
          case 0: // D0 -> K- pi+
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (+-) - (++) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC D^{0} (#times%.2f)",scale_D0[i_hist]), "l");
            break;
          case 1: // D0bar -> K+ pi-
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (-+) - (--) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC #bar{D^{0}} (#times%.2f)",scale_D0[i_hist]), "l");
            break;
          case 2: // D0 + D0bar
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (#pm#mp) - (#pm#pm) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC (#times%.2f)",scale_D0[i_hist]), "l");
            break;
        }// end of legend setup
        leg->Draw();
        
        // Draw label text
        drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.25, 0.83);
        drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.25, 0.775);
        if (mass_integral_background_matching_threshold >= 0)
          drawText(Form("Scaled to Integral in m_{K#pi} > %.2f GeV",mass_integral_background_matching_threshold), 0.93, 0.58, true);
        drawText(Form("%s < #it{p}_{T}^{D^{0} Reco} < %s",
                      bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]), 0.93, height_base_leg - 0.05, true);
        
        drawText(Form("Modifier #it{%s}",modifier),
                 1.0 - gPad->GetRightMargin(),
                 1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
        
        
        // Export the solo canvas
        if (!bSolo) {
          gCanvas_solo->SaveAs(Form("%s/%s/%s/%s_%s_invmass_pt_%s_%s.pdf",
                                    base_directory,
                                    list_directories[1],
                                    list_subdirectories[1][i_hist],
                                    dataset_label,
                                    list_subdirectories[1][i_hist],
                                    bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1] ) );
        }// End of plot save
      }// End of solo (full spectrum range) plotting
    }// End of plotting difference
    
    // Export composite canvas
    gCanvas_quad->SaveAs(Form("%s/%s/%s/%s_recoD0_invmass_pt_%s_%s.pdf",
                             base_directory,
                             list_directories[1],
                             list_subdirectories[1][3],
                             dataset_label,
                             bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1] ) );
    // Done with K/pi difference + fits overlay
    
    
    
    // *--------------------------- Plot fits and yields over the pT range
    
    // Also plot the fits for each resonance to compare the yields
    bool _debug_integral = true;
    for (int i_res = 0; i_res < 3; ++i_res) {
      TH1F* hist_yields_local = new TH1F(Form("yields_%s_pt_%s_%s.pdf",
                                              resonance_labels[i_res],
                                              bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]),
                                         ";;Yield from Fit", 3, 0, 3);
      
      // Find the max amplitude for this resonance
      float range[2] = {1e20, 0};
      for (int i_hist = 0; i_hist < 3; ++i_hist) {
        if (diff_to_fit[i_hist][i_res]->GetMinimum() < range[0]) range[0] = diff_to_fit[i_hist][i_res]->GetMinimum();
        if (diff_to_fit[i_hist][i_res]->GetMaximum() > range[1]) range[1] = diff_to_fit[i_hist][i_res]->GetMaximum();
      }
      
      
      for (int i_hist = 0; i_hist < 3; ++i_hist) {
        // Plot the fit for this type and resonance
        gCanvas_quad->cd(i_hist+1);
        diff_to_fit[i_hist][i_res]->GetYaxis()->SetRangeUser(1.05*range[0], 1.4*range[1]);
        diff_to_fit[i_hist][i_res]->Draw("hist");
        
        
        // Draw fit+bg, pure signal using stored parameters
        TF1* fit_form_draw = new TF1("fit_form",
                                     "[0] * exp(-0.5*((x-[1])/[2])**2)/sqrt(2*TMath::Pi()*[2]**2) + [3] + [4]*(x-[1]) + [5]*(x-[1])**2",
                                     fitting_region_limits[i_res][0], fitting_region_limits[i_res][1]);
        TF1* fit_bg_remove = new TF1("fit_form_background_removed",
                                     "[0] * exp(-0.5*((x-[1])/[2])**2)/sqrt(2*TMath::Pi()*[2]**2)",
                                     fitting_region_limits[i_res][0], fitting_region_limits[i_res][1]);
        for (int i_param = 0; i_param < n_params; ++i_param) {
          fit_form_draw->FixParameter(i_param, fit_params[i_pt_f][i_hist][i_res][i_param]);
          if (i_param > 2) continue;
          fit_bg_remove->FixParameter(i_param, fit_params[i_pt_f][i_hist][i_res][i_param]);
        }
        fit_form_draw->SetLineColor(kBlack);
        fit_form_draw->Draw("same");
        
        // Draw line for the resonance
        D0line->DrawLine(resonance_masses[i_res], 1.05*range[0], resonance_masses[i_res],
                         diff_to_fit[i_hist][i_res]->GetBinContent(diff_to_fit[i_hist][i_res]->FindBin(resonance_masses[i_res])));
        
        // Draw line for zero
        TLine* zeroline = new TLine();
        zeroline->SetLineColor(kGray+2);
        zeroline->DrawLine(fitting_region_limits[i_res][0], 0, fitting_region_limits[i_res][1], 0);
        
        // Draw the resonance subtracted away from the background
        fit_bg_remove->SetLineColor(color_accent1[1]);
        fit_bg_remove->SetLineStyle(1);
        fit_bg_remove->Draw("same");
        
        // Draw a legend
        float height_base_leg = 0.7;
        if (has_MC)height_base_leg = 0.67;
        TLegend* leg = new TLegend(0.65, height_base_leg, 0.95, 0.85);
        leg->SetLineWidth(0);
        leg->AddEntry(D0line, "PDG K^{+}#pi^{-} Res.", "l");
        switch (i_hist) {
          case 0: // D0 -> K- pi+
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (+-) - (++) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC D^{0} (#times%.2f)",scale_D0[i_hist]), "l");
            drawText(Form("%s #rightarrow %s",resonance_string_long[i_res], daughterstring_long[i_hist]), 0.15, 0.70);
            break;
          case 1: // D0bar -> K+ pi-
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (-+) - (--) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC #bar{D^{0}} (#times%.2f)",scale_D0[i_hist]), "l");
            drawText(Form("#bar{%s} #rightarrow %s",resonance_string_long[i_res], daughterstring_long[i_hist]), 0.15, 0.70);
            break;
          case 2: // D0 + D0bar
            leg->AddEntry(hist_1d_diff_SS[i_hist], Form("K#pi (#pm#mp) - (#pm#pm) (#times%.2f)",scale_SS[i_hist]), "l");
            if (has_MC) leg->AddEntry(hist_1d_proj_MC[i_hist], Form("Geant+MC (#times%.2f)",scale_D0[i_hist]), "l");
            drawText(Form("%s + #bar{%s}",resonance_string_long[i_res],resonance_string_long[i_res]), 0.15, 0.70);
            break;
        }// end of legend setup
        leg->AddEntry(fit_bg_remove, "Fit (BG removed)", "l");
        leg->Draw();
        
        // Draw label text
        drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.15, 0.83);
        drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.15, 0.775);
        drawText(Form("%s < #it{p}_{T}^{D^{0} Reco} < %s",
                      bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]), 0.93, height_base_leg - 0.05, true);
        
        drawText(Form("Modifier #it{%s}",modifier),
                 1.0 - gPad->GetRightMargin(),
                 1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
        
        // mark the yield from the fit
        hist_yields_local->GetXaxis()->SetBinLabel(i_hist+1, daughterstring_long[i_hist]);
        hist_yields_local->SetBinContent(i_hist+1, fit_params[i_pt_f][i_hist][i_res][0]);
        hist_yields_local->SetBinError(i_hist+1, fit_params[i_pt_f][i_hist][i_res][0]);
        
        
        // Check that the parameter [0] actually gives the integral
        if (_debug_integral) {
          TF1* background = new TF1("fit_form",
                                    " [3] + [4]*(x-[1]) + [5]*(x-[1])**2",
                                    fitting_region_limits[i_res][0], fitting_region_limits[i_res][1]);
          for (int i_param = 0; i_param < n_params; ++i_param)
            background->FixParameter(i_param, fit_params[i_pt_f][i_hist][i_res][i_param]);
          
          double int_full = fit_form_draw->Integral(fitting_region_limits[i_res][0], fitting_region_limits[i_res][1]);
          double int_bg = background->Integral(fitting_region_limits[i_res][0], fitting_region_limits[i_res][1]);
          
          std::cout << "ptbin [" << bin_range_string[i_pt_f][0] << ',' << bin_range_string[i_pt_f][1] << "], ";
          std::cout << "hist daughter (" << daughterstring_long[i_hist] << "), res " << resonance_string_long[i_res] << std::endl;
          std::cout << "  Genuine integral: \033[32m" << int_full - int_bg << "\033[39m,\tpar [0] = " << fit_params[i_pt_f][i_hist][i_res][0];
          std::cout << "\t(Difference \033[31m" << (int_full - int_bg) - fit_params[i_pt_f][i_hist][i_res][0] << "\033[39m).\n\n";
        }// End of integral debug
      }// End of resonance fit plotting
      
      
      gCanvas_quad->cd(4);
      hist_yields_local->GetXaxis()->SetLabelSize(0.07);
      hist_yields_local->GetXaxis()->SetLabelOffset(0.009);
      hist_yields_local->SetMarkerStyle(20);
      hist_yields_local->SetMarkerColor(kBlack);
      hist_yields_local->SetLineColor(kBlack);
      hist_yields_local->Draw("hist b p e1");
      
      
      
      
      // Draw label text
      drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.15, 0.83);
      drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.15, 0.775);
//      if (mass_integral_background_matching_threshold >= 0)
//        drawText(Form("Scaled to Integral in m_{K#pi} > %.2f GeV",mass_integral_background_matching_threshold), 0.93, 0.58, true);
      drawText(Form("%s < #it{p}_{T}^{D^{0} Reco} < %s",
                    bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]), 0.15, 0.7, false);
      
      drawText(Form("Modifier #it{%s}",modifier),
               1.0 - gPad->GetRightMargin(),
               1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
      
      gCanvas_quad->SaveAs(Form("%s/%s/%s/%s_fit%s_pt_%s_%s.pdf",
                                base_directory,
                                list_directories[1],
                                list_subdirectories[1][4+i_res],
                                dataset_label,
                                resonance_labels[i_res],
                                bin_range_string[i_pt_f][0],bin_range_string[i_pt_f][1]) );
      
      
      gCanvas_quad->cd(4);
      gPad->Clear();
    }// End of plot for fits in each local resonance
    
    
    // *--------------------------- Study background spectra: compare ++/-- with difference/ratio
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // *--------------------------- Write TObjects to File
    
    // TODO write histograms to file
    
  }// End of final pT bin loop
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // *--------------------------- Analyze candidate multiplicities from combinatorics
  
  
  
  
#if 0


  
  // Compare the ++/-- like sign backgrounds
  //  gPad->SetLogy();
  for (int ipT = 0; ipT < background_dist_plus.size(); ++ipT) {
    // First draw both BG on same plot panel
    background_dist_plus[ipT]->SetTitle(";p_{T} [GeV];Like-Sign Background Distributions");
    background_dist_plus[ipT]->SetLineColor(kBlue+2);
    background_dist_plus[ipT]->Draw("hist");
    background_dist_plus[ipT]->SetLineColor(kRed + 2);
    background_dist_minus[ipT]->Draw("hist same");
    
    drawText(Form("#bf{STAR} %s",dataset_label_formal[use_data]), 0.15, 0.83);
    drawText("#it{pp} #sqrt{s} = 200 GeV", 0.15, 0.775);
    pT_bin_tex[ipT]->Draw();
    
    TLegend* leg = new TLegend(0.7, 0.7, 0.95, 0.85);
    leg->SetLineWidth(0);
    leg->AddEntry(background_dist_plus[ipT], Form("#pi^{+}K^{+} (#times%.2f)",1.0), "l");
    leg->AddEntry(background_dist_minus[ipT], Form("#pi^{-}K^{-} (#times%.2f)",1.0), "l");
    leg->Draw();
    
    c->SaveAs(Form("../plots/reconstruction/diff_background/%s_bgplot_ptbin%i.pdf",
                   dataset_label[use_data], ipT));
    
    
    
    // Draw ratio (++)/(--)
    TH1F* ratio = RatioTH1F(background_dist_plus[ipT], background_dist_minus[ipT]);
    ratio->SetTitle(";p_{T} [GeV];Like-Sign Ratio (#pi^{+}K^{+})/(#pi^{-}K^{-})");
    ratio->GetYaxis()->SetRangeUser(0.25, 1.75);
    ratio->Draw("hist");
    
    TLine* unity_line = new TLine();
    unity_line->SetLineStyle(7);
    unity_line->SetLineColor(kGray+1);
    unity_line->DrawLine(ratio->GetBinLowEdge(1), 1,
                         ratio->GetBinLowEdge(ratio->GetNbinsX()+1), 1);
    
    drawText(Form("#bf{STAR} %s",dataset_label_formal[use_data]), 0.15, 0.83);
    drawText("#it{pp} #sqrt{s} = 200 GeV", 0.15, 0.775);
    pT_bin_tex[ipT]->Draw();
    
    c->SaveAs(Form("../plots/reconstruction/diff_background/%s_ratiobg_ptbin%i.pdf",
                   dataset_label[use_data], ipT));
    
    
    
    // Draw difference (++) - (--)
    background_dist_plus[ipT]->Add(background_dist_plus[ipT], background_dist_minus[ipT], 1, -1);
    
    background_dist_plus[ipT]->SetTitle(";p_{T} [GeV];Like-Sign Difference (#pi^{+}K^{+}) - (#pi^{-}K^{-})");
    
    background_dist_plus[ipT]->Draw("hist");
    drawText(Form("#bf{STAR} %s",dataset_label_formal[use_data]), 0.15, 0.83);
    drawText("#it{pp} #sqrt{s} = 200 GeV", 0.15, 0.775);
    pT_bin_tex[ipT]->Draw();
    
    c->SaveAs(Form("../plots/reconstruction/diff_background/%s_diffbg_ptbin%i.pdf",
                   dataset_label[use_data], ipT));
  }// End of like sign comparison
#endif
}// End of extract_signal::main




//========================================================================== Sub-Macros



void plotSpectra() {
  
  return;
}// End of extract_signal::plotSpectra

