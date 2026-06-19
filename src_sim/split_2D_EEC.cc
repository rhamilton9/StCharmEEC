/// A ROOT macro for dividing out a 2D histogram of EEC/pT
/// into several histograms in different pT bins
///
/// R. Hamilton 6/17/2026

#include "pythia_config.hpp"

#include "../utils/root_draw_tools.h"

// Settings and config
bool use_charged_jets = true;

// Use a restricted pThat range for this feasibility study which should be reasonable!
double pt_binedge[2] = {5, -1};     // pThat range to generate hard scatterings in PYTHIA
double EEC_pt_binedge[2] = {5, 50}; // pT range of final clustered jets to store in the EEC histograms
double cross_section = 3.98029761e-4; // Estimated cross section for PYTHIA generation


// =============== * Global Variables

// Label strings
char comparison_mode_string[4][30] = {"KS dist.", "#chi^{2}/NDF", "KS #it{p}-value", "#chi^{2} Test #it{p}-value"};
char chargejetstring[2][10] = {"full","char"};
char chargejetstringformal[2][10] = {"Full","Charged"};

// ------------------------------------------------------Main
void split_2D_EEC() {
  // File I/O
  double snn_tev = ((double) sqrt_s)/1000;
  char gen_info_string[50];
  snprintf(gen_info_string, 50, "%.2fTeV_r%.1fjet", snn_tev, jetRadius);
  TFile *datafile = new TFile(Form("PYTHIA/out_data/pythia_%s_merged_lowbias.root",gen_info_string), "update");
  if (datafile->IsZombie()) {
    std::cout << "\033[31mError\033[39m :: File could not be found!" << std::endl;
    return;
  }
  
  
  // Get the necessary histograms from the parent file
  // EEC histograms
  TH2D* EEC_2D[3];
  TH1D* gConstituentMultiplicity[3];
  std::cout << "Attempting to find {" << Form("EEC_%sjet_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]) << "}..." << std::endl;
  EEC_2D[0] = static_cast<TH2D*>(datafile->Get<TH2D>(Form("EEC_%sjet_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]))->Clone());
  EEC_2D[1] = static_cast<TH2D*>(datafile->Get<TH2D>(Form("EEC_%sjet_ptbin%.f-%.f_D0tagged",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]))->Clone());
  EEC_2D[2] = static_cast<TH2D*>(datafile->Get<TH2D>(Form("EEC_%sjet_ptbin%.f-%.f_nocharm",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]))->Clone());
  // Jet constituent multiplicity histograms
  gConstituentMultiplicity[0] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gConstituentMultiplicity[1] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_D0tagged",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  gConstituentMultiplicity[2] = datafile->Get<TH1D>(Form("%sjet_constituents_ptbin%.f-%.f_nocharm",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  int njet[3] = {
    static_cast<int>(gConstituentMultiplicity[0]->GetEntries()),
    static_cast<int>(gConstituentMultiplicity[1]->GetEntries()),
    static_cast<int>(gConstituentMultiplicity[2]->GetEntries()),
  };
  
  
//  TFile* outfile_hist = new TFile(Form("PYTHIA/out_data/pythia_%s_merged_lowbias_outhist.root",gen_info_string), "recreate");
  TH1D* EEC_1D[3];
  
  
  // Canvas setup
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  TCanvas* gCanvas = new TCanvas();
  gPad->SetTicks(1,1);
  gPad->SetLogx();
  gCanvas->SetCanvasSize(800, 500);
  gCanvas->SetMargin(0.10, 0.03, 0.08, 0.12);
  
  TCanvas* gCanvas_2D = new TCanvas();
  gPad->SetTicks(1,1);
  gPad->SetLogx();
  gPad->SetLogz();
  gCanvas_2D->SetCanvasSize(800, 450);
  gCanvas_2D->SetMargin(0.09, 0.13, 0.08, 0.12);
  
  
  
  
  // Label strings
  char hist_same_switch[2][10] = {"hist","hist same"};
  char EEC_label_short[3][15] = {"inclusive", "D0tagged", "nocharm"};
  char EEC_label_formal[3][15] = {"Inclusive", "D^{0}-Tagged", "D^{0}-Free"};
  
  // Draw EEC histograms in a clean manner first
  TLegend* leg_EEC = new TLegend(0.16, 0.50, 0.33, 0.70);
  leg_EEC->SetLineWidth(0);
  for (int i_hist = 0; i_hist < 3; ++i_hist) {
    // Adjust axes
    gCanvas_2D->cd();
    EEC_2D[i_hist]->GetZaxis()->SetTitle("1/N^{jet} d^{2}#Sigma/dlogR_{L}dp_{T}");
    EEC_2D[i_hist]->GetYaxis()->SetTitleOffset(1.25);
    EEC_2D[i_hist]->GetXaxis()->SetTitleOffset(0.8);
    
    // Rescale to 1/njet and plot
    EEC_2D[i_hist]->Scale(1./njet[i_hist]);
    EEC_2D[i_hist]->Draw("colz");
    drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", snn_tev), gPad->GetLeftMargin(), 0.95);
    drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
    drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
    drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
    gCanvas_2D->SaveAs(Form("PYTHIA/out_plots/feasibility/EEC/%s%s_EEC_pT_%s.pdf",chargejetstring[use_charged_jets],gen_info_string, EEC_label_short[i_hist]));
    
    
    // *------ 1D histograms
    
    // Set up 1d histograms using same binning as is done for the
    gCanvas->cd();
    std::vector<double> bins_EEC;
    int nbin_EEC = EEC_2D[i_hist]->GetXaxis()->GetNbins();
    for (int i_bin_EEC = 1; i_bin_EEC <= nbin_EEC + 1; ++i_bin_EEC) {
      bins_EEC.push_back(EEC_2D[i_hist]->GetXaxis()->GetBinLowEdge(i_bin_EEC));
    }
    EEC_1D[i_hist] = new TH1D(Form("EEC_%sjet_ptbin%.f-%.f_%s_proj",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1],EEC_label_short[i_hist]),
                              ";R_{L};1/N^{jet} d#Sigma/dlog(R_{L})",
                              nbin_EEC, bins_EEC.data());
    EEC_1D[i_hist]->SetLineColor(kBlack);
    
    // Perform pT projection
    int nbin_pt = EEC_2D[i_hist]->GetYaxis()->GetNbins();
    for (int i_bin_EEC = 1; i_bin_EEC <= nbin_EEC; ++i_bin_EEC) {
      for (int i_bin_pt = 0; i_bin_pt <= nbin_pt; ++i_bin_pt) {
        EEC_1D[i_hist]->SetBinContent(i_bin_EEC, EEC_1D[i_hist]->GetBinContent(i_bin_EEC) + EEC_2D[i_hist]->GetBinContent(i_bin_EEC, i_bin_pt));
      }
    }
    
    // Set aesthetic settings
    EEC_1D[i_hist]->SetLineColor(plot_colors[i_hist]);
    
    // Draw, add to legend
    EEC_1D[i_hist]->Draw(hist_same_switch[i_hist != 0]);
    leg_EEC->AddEntry(EEC_1D[i_hist], EEC_label_formal[i_hist],"l");
  }
  // Add some labels and legend
  leg_EEC->Draw();
  drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", snn_tev), gPad->GetLeftMargin(), 0.95);
  drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
  drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
  drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
  
  
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/EEC/%s%s_EEC.pdf", chargejetstringformal[use_charged_jets], gen_info_string));
  
  // Write projected histograms to file
//  outfile_hist->cd();
  for (int i_hist = 0; i_hist < 3; ++i_hist) {
    std::cout << "debug" << std::endl;
    EEC_1D[i_hist]->Write(EEC_1D[i_hist]->GetName(), TObject::kOverwrite);
  }
  
  return;
}// End of split_2D_EEC::main
