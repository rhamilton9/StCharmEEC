/// A ROOT macro for producing a contour of statistical likelihood
/// to successfully perform a dead cone extraction, dependent on
/// signal/noise ratio and dataset luminosity as the x/y axes.
///
/// R. Hamilton 6/10/2026

#include "pythia_config.hpp"

#include "../utils/root_draw_tools.h"

// Settings and config
bool use_charged_jets = false;

// Use a restricted pThat range for this feasibility study which should be reasonable!
double pt_binedge[2] = {15, 35};     // pThat range to generate hard scatterings in PYTHIA
double EEC_pt_binedge[2] = {20, 30}; // pT range of final clustered jets to store in the EEC histograms
double cross_section = 3.98029761e-4; // Estimated cross section for PYTHIA generation

// Available options to compare the histograms at each iteration of the MC generation:
//   0 - Kalmogorov-Smirnov test
//   1 - Chi2 test
// Add more as desired.
int switch_comparison_mode = 1;

// Label strings
char comparison_mode_string[2][20] = {"KS dist.", "#chi^{2}/NDF"};
char chargejetstring[2][10] = {"full","char"};
char chargejetstringformal[2][10] = {"Full","Charged"};

// Control how many iterations the MC variation undertakes
const int n_iter = 1000;


// Global pointers to key histograms
TH1D* gEEC[3];
TH1D* gConstituentMultiplicity[3];
TCanvas* gCanvas;

// Forward declare helper methods
double getFeasibilitySingleBin(double signal_to_noise, double Njet);





// ------------------------------------------------------Main
void make_feasibility_contour() {
  // File I/O
  double snn_tev = ((double) sqrt_s)/1000;
  char gen_info_string[50];
  snprintf(gen_info_string, 50, "%.2fTeV_r%.1fjet", snn_tev, jetRadius);
  TFile *datafile = new TFile(Form("PYTHIA/out_data/pythpp_merged_%s.root",gen_info_string));
  if (datafile->IsZombie()) {
    std::cout << "\033[31mError\033[39m :: File could not be found!" << std::endl;
    return;
  }
  
  // Get the necessary histograms from the parent file
  // EEC histograms
  std::cout << "Attempting to find {" << Form("EEC_%sjet_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]) << "}..." << std::endl;
  gEEC[0] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gEEC[1] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_D0tagged",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gEEC[2] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_nocharm",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  // Jet constituent multiplicity histograms
  gConstituentMultiplicity[0] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gConstituentMultiplicity[1] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_D0tagged",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gConstituentMultiplicity[2] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_nocharm",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  int njet[3] = {
    static_cast<int>(gConstituentMultiplicity[0]->GetEntries()),
    static_cast<int>(gConstituentMultiplicity[1]->GetEntries()),
    static_cast<int>(gConstituentMultiplicity[2]->GetEntries()),
  };
  
  
  // Canvas setup
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  gCanvas = new TCanvas();
  gCanvas->SetCanvasSize(800, 500);
  gCanvas->SetMargin(0.09, 0.03, 0.08, 0.12);
  gPad->SetTicks(1,1);
  
  // Label strings
  char hist_same_switch[2][10] = {"hist","hist same"};
  char EEC_label_formal[3][15] = {"Inclusive", "D^{0}-Tagged", "D^{0}-Free"};
  
  // Draw EEC histograms in a clean manner first
  gPad->SetLogx();
  TLegend* leg_EEC = new TLegend(0.16, 0.50, 0.33, 0.70);
  leg_EEC->SetLineWidth(0);
  for (int i_hist = 0; i_hist < 3; ++i_hist) {
    // Adjust axes
    gEEC[i_hist]->GetYaxis()->SetTitle("1/N^{jet} d#Sigma/dR_{L}");
    gEEC[i_hist]->GetYaxis()->SetTitleOffset(1.25);
    gEEC[i_hist]->GetXaxis()->SetTitleOffset(0.8);
    
    // Rescale to 1/njet and divide by bin width
    gEEC[i_hist]->Scale(1./njet[i_hist]);
//    gEEC[i_hist]->Scale(1./njet[i_hist], "width");
    
    // Set aesthetic settings
    gEEC[i_hist]->SetLineColor(plot_colors[i_hist]);
    
    // Draw, add to legend
    gEEC[i_hist]->Draw(hist_same_switch[i_hist != 0]);
    leg_EEC->AddEntry(gEEC[i_hist], EEC_label_formal[i_hist],"l");
  }
  // Add some labels and legend
  leg_EEC->Draw();
  drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", snn_tev), gPad->GetLeftMargin(), 0.95);
  drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
  drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
  drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
  
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/%s_EEC.pdf", gen_info_string));
  
  
  
  // Analyze feasibility in a grid of signal:background, Njet.
  getFeasibilitySingleBin(50.0, 1000);
  
  
  return;
}// End of make_feasibility_contour::main






// ------------------------------------------------------Helper Methods

// Find the feasibility of the measurement
double getFeasibilitySingleBin(double signal_to_noise, double Njet) {
  int njet[3] = {
    static_cast<int>(Njet),
    static_cast<int>((signal_to_noise / (1 + signal_to_noise)) * Njet),
    static_cast<int>((1               / (1 + signal_to_noise)) * Njet)
  };
  std::cout << "Njet {tot,sig,bkg} :: {" << njet[0] << ',' << njet[1] << ',' << njet[2] << "}." << std::endl;
  
  // Histogram of deviation
  double comparison_range[2] = {50, -50};
  TH1D* hist_comparison_metric = new TH1D("hist_comparison_metric",
                                          Form(";%s;Simulated Experiment Counts",comparison_mode_string[switch_comparison_mode]),
                                          1000, -50, 50);
  
  
  //==========================================Simulated Experiment Generation
  
  
  // Begin event loop
  double integral_scale = gEEC[1]->Integral("width");
  TH1D* hist_reco_EEC = static_cast<TH1D*>(gEEC[0]->Clone());
  TH1D* hist_sideband = static_cast<TH1D*>(gEEC[0]->Clone());
  for (int i_iter = 0; i_iter < n_iter; ++i_iter) {
    // Reset the EEC reconstructed histograms for this generation
    hist_reco_EEC->Reset();
    hist_sideband->Reset();
    
    // Fill signal jets for this histogram
    for (int i_sig = 0; i_sig <= njet[1]; ++i_sig) {
      // Sample constituent histogram for the number of splittings
      int n_constit = std::floor(gConstituentMultiplicity[1]->GetRandom());
      
      // Sample EEC fills 1/2(n)(n-1) combinations among n constituents
      for (int i_pair = 0; i_pair <= n_constit*(n_constit-1)/2; ++i_pair) {
        hist_reco_EEC->Fill(gEEC[1]->GetRandom());
      }
    }// End of signal fill
    
    
    // Fill background jets for signal hist and sideband hist
    for (int i_bkg = 0; i_bkg <= njet[2]; ++i_bkg) {
      // Sample constituent histogram for the number of splittings
      int n_constit = std::floor(gConstituentMultiplicity[2]->GetRandom());
      
      // Sample EEC fills 1/2(n)(n-1) combinations among n constituents
      for (int i_pair = 0; i_pair <= n_constit*(n_constit-1)/2; ++i_pair) {
        hist_reco_EEC->Fill(gEEC[2]->GetRandom());
      }
      
      // Do the same for the sideband histogram
      n_constit = std::floor(gConstituentMultiplicity[2]->GetRandom());
      for (int i_pair = 0; i_pair <= n_constit*(n_constit-1)/2; ++i_pair) {
        hist_reco_EEC->Fill(gEEC[2]->GetRandom());
      }
    }// End of background fill
    
    
    // Perform the subtraction of the sideband and renormalize
    TH1D* diff = static_cast<TH1D*>(gEEC[0]->Clone());
    diff->Reset();
    diff->Add(hist_reco_EEC, hist_sideband, 1, -1);
    diff->Scale(integral_scale/(diff->Integral("width")));
    
    // Analyze the difference in the subtracted histograms
    double comparison_metric = 0;
    switch (switch_comparison_mode) {
      case 0: default: // Kolmogorov-Smirnov Test
        comparison_metric = diff->KolmogorovTest(gEEC[1], "M");
        break;
      case 1: // Chi2 test
        comparison_metric = diff->Chi2Test(gEEC[1], "CHI2/NDF");
        break;
    }
    hist_comparison_metric->Fill(comparison_metric);
    if (comparison_metric < comparison_range[0]) comparison_range[0] = comparison_metric;
    if (comparison_metric > comparison_range[1]) comparison_range[1] = comparison_metric;
    // End of comparison evaluation
    
    
    // Draw the iterations occasionally to compare
    if (i_iter % (n_iter / 10) == 0) {
      std::cout << "Finished with simul. experiment generation #\033[32m" << i_iter << "\033[39m" << std::endl;
      gEEC[1]->Draw("hist");
      diff->Draw("hist same");
      
      // Legend
      TLegend* leg_iter = new TLegend(0.16, 0.50, 0.33, 0.65);
      leg_iter->SetLineWidth(0);
      leg_iter->AddEntry(diff, "Sampled", "l");
      leg_iter->AddEntry(gEEC[1], "D0-Tagged", "l");
      leg_iter->Draw();
      
      // Generation text
      drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", ((double) sqrt_s)/1000), gPad->GetLeftMargin(), 0.95);
      drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
      drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
      drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
      
      // Text for the comparison
      drawText(Form("N^{jet} = %i, N^{signal} = %i", njet[0], njet[1]), 0.05 + gPad->GetLeftMargin(), 0.925 - gPad->GetTopMargin());
      drawText(Form("%s :: #color[2]{%.4f}", comparison_mode_string[switch_comparison_mode], comparison_metric), 0.05 + gPad->GetLeftMargin(), 0.88 - gPad->GetTopMargin());
      drawText(Form("Iteration #bf{%i}", i_iter), 0.05 + gPad->GetLeftMargin(), 0.83 - gPad->GetTopMargin());
      
      
      gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/iteration_check/iter_%i.pdf",i_iter / (n_iter / 10) ) );
    }// End of iteration drawing
    
  }// End of MC loop
  
  //==========================================Evaluation and final feasibility computation
  
  // Evaluate the p value of agreement from the monte carlo
  
  
  
  // Plot the histogram of the comparison metric results
  gPad->SetLogx(0);
  hist_comparison_metric->GetXaxis()->SetRangeUser(comparison_range[0] - 1, comparison_range[1] + 1);
  hist_comparison_metric->Draw("hist");
  
  // Generation text
  drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", ((double) sqrt_s)/1000), gPad->GetLeftMargin(), 0.95);
  drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
  drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
  drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
  
  // Simulation text
  drawText(Form("N^{jet} = %i, N^{signal} = %i", njet[0], njet[1]), 0.05 + gPad->GetLeftMargin(), 0.925 - gPad->GetTopMargin());
  drawText(Form("Total iterations: %i", n_iter), 0.05 + gPad->GetLeftMargin(), 0.88 - gPad->GetTopMargin());
  
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/metric_plots/feasibility_SN%.2f_Njet%i.pdf",signal_to_noise, static_cast<int>(Njet)));
  
  
  
  
  return 0;
}// End of make_feasibility_contour::getFeasibilitySingleBin
