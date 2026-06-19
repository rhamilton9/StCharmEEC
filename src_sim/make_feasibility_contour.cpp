/// A ROOT macro for producing a contour of statistical likelihood
/// to successfully perform a dead cone extraction, dependent on
/// signal/noise ratio and dataset luminosity as the x/y axes.
///
/// R. Hamilton 6/10/2026

#include "pythia_config.hpp"

#include "../utils/root_draw_tools.h"

// Default--use 0 inf
//#define use2030
//#define use5inf

#define rejection_test

// Settings and config
bool use_charged_jets = true;
bool _debug_global = false;

// Use a restricted pThat range for this feasibility study which should be reasonable!
#ifdef use2030
double pt_binedge[2] = {15, 35};     // pThat range to generate hard scatterings in PYTHIA
double EEC_pt_binedge[2] = {20, 30}; // pT range of final clustered jets to store in the EEC histograms
#else
#ifdef use5inf
double pt_binedge[2] = {5, -1};     // pThat range to generate hard scatterings in PYTHIA
#else
double pt_binedge[2] = {0, -1};     // pThat range to generate hard scatterings in PYTHIA
#endif
double EEC_pt_binedge[2] = {5, 50}; // pT range of final clustered jets to store in the EEC histograms
#endif
double cross_section = 3.98029761e-4; // Estimated cross section for PYTHIA generation


// =============== * Settings for the experimental signal/noise input D meson yields

// {D0, D*+}
double Dmeson_yields_raw[2] = {4085, 364};
double Dmeson_signalnoise[2] = {0.00625, 0.06066};
double Dmeson_ptrange[2][2] = {
  {0.6, 2.2},
  {2.0, 6.0}
};


// =============== * Experiment simulation settings

// Available options to compare the histograms at each iteration of the MC generation:
//   0 - Kalmogorov-Smirnov test
//   1 - Chi2 test
//   2 - Kalmogorov-Smirnov p-value
//   3 - Chi2 test p-value
// Add more as desired.
int switch_comparison_mode = 3;

// Control how many iterations the MC variation undertakes
const int n_iter = 500;

// =============== * Feasibility Histogram controls

int nbin_signalnoise = 8;
double range_signalnoise[2] = {1e-3, 1e1};
int nbin_jetcount = 4;
double range_jetcount[2] = {50, 5e4};

// =============== * Global Variables

// Label strings
char comparison_mode_string[4][30] = {"KS dist.", "#chi^{2}/NDF", "KS #it{p}-value", "#chi^{2} Test #it{p}-value"};
char chargejetstring[2][10] = {"full","char"};
char chargejetstringformal[2][10] = {"Full","Charged"};

// Global pointers to key histograms
TH1D* gEEC[3];
TH1D* gConstituentMultiplicity[3];
TCanvas* gCanvas;
bool thumbnail_mode = false;


// =============== * Struct for storing data

struct feasibility_data_node {
  double signal_noise;
  double njet;
  double p_value;
  feasibility_data_node* next;
};

// Global pointers
feasibility_data_node* head = NULL;
feasibility_data_node* tail = NULL;

void appendNode(double sn, double nj, double pv);
void printNodeList();

// =============== * Forward Declarations

// General helper methods
std::vector<double> getUniformBinningLogScale(int nbin, double low, double high);

// File I/O to avoid repeated computation
void getPreviouslyMeasuredData(TH2D* hist_to_fill, double bin_matching_resolution);
void writeDataToFile();

// Computation helper
double getFeasibilitySingleBin(double signal_to_noise, double Njet);




// ------------------------------------------------------Main
void make_feasibility_contour() {
  // File I/O
  double snn_tev = ((double) sqrt_s)/1000;
  char gen_info_string[50];
  snprintf(gen_info_string, 50, "%.2fTeV_r%.1fjet", snn_tev, jetRadius);
  
  char datafile_name[100];
#ifdef use2030
  snprintf(datafile_name, 100, "PYTHIA/out_data/pythia_%s_merged_smallsample.root", gen_info_string);
#else
#ifdef use5inf
  snprintf(datafile_name, 100, "PYTHIA/out_data/pythia_%s_merged_lowbias.root", gen_info_string);
#else
  snprintf(datafile_name, 100, "PYTHIA/out_data/pythia_%s_merged_lowbias_0.root", gen_info_string);
#endif
#endif
  TFile *datafile = new TFile(datafile_name);
  if (datafile->IsZombie()) {
    std::cout << "\033[31mError\033[39m :: File could not be found!" << std::endl;
    return;
  }
  
  if (_debug_global) datafile->ls();
  
  TFile* outfile_hist = new TFile("PYTHIA/out_data/feasibility_hist.root", "recreate");
  
  // Get the necessary histograms from the parent file
  // EEC histograms
  std::cout << "Attempting to find {" << Form("EEC_%sjet_ptbin%.f-%.f_inclusive_proj",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]) << "}..." << std::endl;
  std::string check2D(datafile_name);
  if (check2D.find("lowbias") != std::string::npos) {
    gEEC[0] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_inclusive_proj",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
    gEEC[1] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_D0tagged_proj",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
    gEEC[2] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_nocharm_proj",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
    if (_debug_global) gEEC[1]->Print();
  } else {
    gEEC[0] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_inclusive",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
    gEEC[1] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_D0tagged",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
    gEEC[2] = datafile->Get<TH1D>(Form("EEC_%sjet_ptbin%.f-%.f_nocharm",chargejetstring[use_charged_jets],EEC_pt_binedge[0],EEC_pt_binedge[1]));
  }// Jet constituent multiplicity histograms
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
  if (!_debug_global) gErrorIgnoreLevel = kWarning;
  gCanvas = new TCanvas();
  gPad->SetTicks(1,1);
  if (thumbnail_mode) {
    gCanvas->SetCanvasSize(500, 400);
    gCanvas->SetMargin(0.15, 0.03, 0.10, 0.03);
  } else {
    gCanvas->SetCanvasSize(800, 500);
    gCanvas->SetMargin(0.09, 0.03, 0.08, 0.12);
  }
  
  // Label strings
  char ptbin_inf_string[2][20] = {"", "#infty)"};
  snprintf(ptbin_inf_string[0], 20, "%.f]",pt_binedge[1]);
  char hist_same_switch[2][10] = {"hist","hist same"};
  char EEC_label_formal[3][15] = {"Inclusive", "D-Tagged", "D-Free"};
  
  // Draw EEC histograms in a clean manner first
  gPad->SetLogx();
  TLegend* leg_EEC;
  if (thumbnail_mode) leg_EEC = new TLegend(0.16, 0.35, 0.55, 0.85);
  else                leg_EEC = new TLegend(0.16, 0.50, 0.33, 0.70);
  leg_EEC->SetLineWidth(0);
  for (int i_hist = 0; i_hist < 3; ++i_hist) {
    // Adjust axes
    gEEC[i_hist]->GetYaxis()->SetTitle("1/N^{jet} d#Sigma/dlog(R_{L})");
    gEEC[i_hist]->GetYaxis()->SetTitleOffset(1.25);
    gEEC[i_hist]->GetXaxis()->SetTitleOffset(0.8);
    if (thumbnail_mode) {
      gEEC[i_hist]->SetLineWidth(3);
      gEEC[i_hist]->GetXaxis()->SetLabelSize(0.06);
      gEEC[i_hist]->GetXaxis()->SetTitleSize(0.065);
      gEEC[i_hist]->GetXaxis()->SetTitleOffset(0.6);
      gEEC[i_hist]->GetYaxis()->SetLabelSize(0.06);
      gEEC[i_hist]->GetYaxis()->SetTitleSize(0.06);
      gEEC[i_hist]->GetYaxis()->SetTitleOffset(0.9);
    }
    
    // Rescale to 1/njet and divide by bin width
    gEEC[i_hist]->Scale(1./gEEC[i_hist]->Integral("width"));
//    gEEC[i_hist]->Scale(1./njet[i_hist], "width");
    
    // Set aesthetic settings
    gEEC[i_hist]->SetLineColor(plot_colors[i_hist]);
    
    // Draw, add to legend
    gEEC[i_hist]->Draw(hist_same_switch[i_hist != 0]);
    leg_EEC->AddEntry(gEEC[i_hist], Form("%s",EEC_label_formal[i_hist]),"l");
  }
  // Add some labels and legend
  leg_EEC->Draw();
  if (!thumbnail_mode) {
    drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", snn_tev), gPad->GetLeftMargin(), 0.95);
    drawText(Form("#hat{#it{p}}_{T} #in  [%.0f,%s, #it{p}_{T}^{jet} #in  [%.0f,%0.f]",
                  pt_binedge[0],ptbin_inf_string[pt_binedge[1]==-1],EEC_pt_binedge[0],EEC_pt_binedge[1]),
             1.0 - gPad->GetRightMargin(), 0.95, true);
    drawText(Form("%s, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
    drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
  } else {
    drawText(Form("#hat{#it{p}}_{T} #in  [%.0f,%s, #it{p}_{T}^{jet} #in  [%.0f,%0.f]",
                  pt_binedge[0],ptbin_inf_string[pt_binedge[1]==-1],EEC_pt_binedge[0],EEC_pt_binedge[1]),
             0.05 + gPad->GetLeftMargin(), 0.87, false, kBlack, 0.065);
  }
  
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/%s_EEC.pdf", gen_info_string));
  
  
  // Determine the estimated jet yield from the hadron yield
  double Dmeson_yields_jetcorrected[2];
#ifdef use2030
  // Did not have these histograms in the 20-30 GeV run
  Dmeson_yields_jetcorrected[0] = Dmeson_yields_raw[0];
  Dmeson_yields_jetcorrected[1] = Dmeson_yields_raw[1];
#else
  // Estiamte the fraction in the pT range using the D0 jet partitioned histogram
  TH2D* hist_Dmeson_jet_separated = datafile->Get<TH2D>("histJetDMesonMomentum");
  int nbin_pt = hist_Dmeson_jet_separated->GetXaxis()->GetNbins();
  std::cout << "\033[36mDebug!\033[39m" << std::endl;
  
  // Gather the projections for each meson
  for (int i_Dmeson = 0; i_Dmeson < 2; ++i_Dmeson) {
    
    // Find which bins to loop over
    int lowbin_Dpt = hist_Dmeson_jet_separated->GetYaxis()->FindBin(Dmeson_ptrange[i_Dmeson][0]);
    int topbin_Dpt = hist_Dmeson_jet_separated->GetYaxis()->FindBin(Dmeson_ptrange[i_Dmeson][1]);
    
    
    double count_Dmeson_injet = 0;
    double count_Dmeson_all = 0;
    for (int i_bin_pt = 1; i_bin_pt <= nbin_pt; ++i_bin_pt) {
      for (int i_bin_Dpt = lowbin_Dpt; i_bin_Dpt <= topbin_Dpt; ++i_bin_Dpt) {
        // Count all yields for total
        count_Dmeson_all += hist_Dmeson_jet_separated->GetBinContent(i_bin_pt, i_bin_Dpt);
        
        // Count Dmesons in jet
        if (hist_Dmeson_jet_separated->GetXaxis()->GetBinCenter(i_bin_pt) < 0) continue;
        
        count_Dmeson_injet += hist_Dmeson_jet_separated->GetBinContent(i_bin_pt, i_bin_Dpt);
      }// End of D pT bin loop
      
    }// End of pT bin loop
    
    
    // Set the Dmeson yield as the ratio multiplied into the raw yield
    Dmeson_yields_jetcorrected[i_Dmeson] = Dmeson_yields_raw[i_Dmeson] * count_Dmeson_injet / count_Dmeson_all;
    if (i_Dmeson) std::cout << "\033[36mD*\033[39m";
    else          std::cout << "\033[32mD0\033[39m";
    std::cout << " meson jet yield fraction :: " << count_Dmeson_injet / count_Dmeson_all << std::endl;
  }// End of D0/D* loop
  
#endif
  
  
  
  
  
  
  
  // Analyze feasibility in a grid of signal:background, Njet.
  
  // Set up histogram with uniform binning in log and check for existing data
  std::vector<double> binedge_signalnoise = getUniformBinningLogScale(nbin_signalnoise,
                                                                      range_signalnoise[0],
                                                                      range_signalnoise[1]);
  std::vector<double> binedge_jetcount = getUniformBinningLogScale(nbin_jetcount,
                                                                   range_jetcount[0],
                                                                   range_jetcount[1]);
#ifdef rejection_test
  TH2D* hist_feasibility = new TH2D("hist_feasibility", ";Signal/Noise Ratio;D-Tagged Jet Yield;Average Inclusive Rejection #it{p}-Value",
                                    nbin_signalnoise, binedge_signalnoise.data(),
                                    nbin_jetcount, binedge_jetcount.data());
#else
  TH2D* hist_feasibility = new TH2D("hist_feasibility", ";Signal/Noise Ratio;D-Tagged Jet Yield;Average D-Tagged Observation #it{p}-Value",
                                    nbin_signalnoise, binedge_signalnoise.data(),
                                    nbin_jetcount, binedge_jetcount.data());
#endif
  getPreviouslyMeasuredData(hist_feasibility, 0.05);
  hist_feasibility->GetZaxis()->SetRangeUser(-.001,1.001);
  
  // Loop over any data not yet collected
  for (int i_sn = 1; i_sn <= nbin_signalnoise; ++i_sn) {
    double bincenter_x = hist_feasibility->GetXaxis()->GetBinCenter(i_sn);
    double last_result = 0;
    for (int i_jc= 1; i_jc <= nbin_jetcount; ++i_jc) {
      double bincenter_y = hist_feasibility->GetYaxis()->GetBinCenter(i_jc);
      
      // Check if the data was stored in the file
      if (hist_feasibility->GetBinContent(i_sn, i_jc) != 0) continue;
      
      // If the last jet count result was near 1.0, this one will be too--skip it
      if (switch_comparison_mode >= 2 &&
          last_result == 1) {
        hist_feasibility->SetBinContent(i_sn, i_jc, 10);
        continue;
      }
      
      last_result = getFeasibilitySingleBin(bincenter_x, bincenter_y);
      hist_feasibility->SetBinContent(i_sn, i_jc, last_result);
      
      // Append result to the linked list of nodes
      appendNode(bincenter_x, bincenter_y, last_result);
    }
  }// End of feasibility collection
  
  gPad->SetMargin(0.07, 0.13, 0.075, 0.12);
  gPad->SetLogx();
  gPad->SetLogy();
//  gPad->SetLogz();
  hist_feasibility->GetYaxis()->SetTitleOffset(0.9);
  hist_feasibility->Draw("colz");
  
  // Draw lines for the D0 signal/noise ratio and jet count in previous data
  TLine* D0line[2];
  D0line[0] = new TLine();
  D0line[0]->SetLineColor(plot_colors[1]);
  D0line[1] = new TLine();
  D0line[1]->SetLineColor(plot_colors[1]);
  D0line[1]->SetLineStyle(7);
  for (int i_D0_line = 0; i_D0_line < 2; ++i_D0_line) {
    D0line[i_D0_line]->DrawLine(Dmeson_signalnoise[i_D0_line], range_jetcount[0],
                                Dmeson_signalnoise[i_D0_line], range_jetcount[1]);
    D0line[i_D0_line]->DrawLine(range_signalnoise[0], Dmeson_yields_jetcorrected[i_D0_line],
                                range_signalnoise[1], Dmeson_yields_jetcorrected[i_D0_line]);
  }
  
  
  
  // Draw some text about the setup
  drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", snn_tev), gPad->GetLeftMargin(), 0.95);
  drawText(Form("#hat{#it{p}}_{T} #in  [%.0f,%s, #it{p}_{T}^{jet} #in  [%.0f,%0.f]",
                pt_binedge[0],ptbin_inf_string[pt_binedge[1]==-1],EEC_pt_binedge[0],EEC_pt_binedge[1]),
           1.0 - gPad->GetRightMargin(), 0.95, true);
  drawText(Form("%s, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
  drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
  
#ifdef rejection_test
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/hist_feasibility_iter%i_%.fGeV_rejection.pdf", n_iter, pt_binedge[0]));
#else
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/hist_feasibility_iter%i_%.fGeV.pdf", n_iter, pt_binedge[0]));
#endif
  
  
  
  hist_feasibility->Write(hist_feasibility->GetName(), TObject::kOverwrite);
  outfile_hist->Close();
  
  
  // Write file of data
  if (_debug_global) printNodeList();
  writeDataToFile();
  
  
//  getFeasibilitySingleBin(1, 1000);
//
//  std::cout << "Comparison of base EEC histograms :: (KS " << gEEC[1]->KolmogorovTest(gEEC[0]) << ")\n";
//  gEEC[1]->Chi2Test(gEEC[0], "WWP");
  return;
}// End of make_feasibility_contour::main







// ------------------------------------------------------Struct helpers

// Append a new node to the list
void appendNode(double sn, double nj, double pv) {
  if (head == NULL) {
    head = new feasibility_data_node;
    head->signal_noise = sn;
    head->njet = nj;
    head->p_value = pv;
    tail = head;
    return;
  }
  
  tail->next = new feasibility_data_node;
  tail = tail->next;
  tail->signal_noise = sn;
  tail->njet = nj;
  tail->p_value = pv;
  return;
}// End of make_feasibility_contour::appendNode


// Print the node list
void printNodeList() {
  feasibility_data_node* print = head;
  while (print != NULL) {
    std::cout << "-------------------------" << std::endl;
    std::cout << "| S/N = " << print->signal_noise << "     \t|" << std::endl;
    std::cout << "| Njet = " << print->njet << "      \t|" << std::endl;
    std::cout << "| Pval = " << print->p_value << "      \t|" << std::endl;
    std::cout << "-------------------------" << std::endl;
    std::cout << "            ||           " << std::endl;
    std::cout << "            ||           " << std::endl;
    std::cout << "            \\/           " << std::endl;
    print = print->next;
  }return;
}

// ------------------------------------------------------Helper Methods

// Produce a uniform binning as would be viewed on a logarithmic scale
std::vector<double> getUniformBinningLogScale(int nbin, double low, double high) {
  std::vector<double> logbinedge = std::vector<double>();
  double logbinrange[2] = {std::log(low), std::log(high)};
  double delta_log_binedge = (logbinrange[1] - logbinrange[0]) / nbin;
  for (int i = 0; i <= nbin; ++i) logbinedge.push_back(TMath::Exp(logbinrange[0] + i * delta_log_binedge));
  return logbinedge;
}// End of make_feasibility_contour::getUniformBinningLogScale()


// Find any relevant values from previously gathered data
void getPreviouslyMeasuredData(TH2D* hist_to_fill, double bin_matching_resolution) {
  
  bool _debug = true;
  
  char filename[100];
#ifdef rejection_test
  snprintf(filename, 100, "PYTHIA/out_data/feasibility_data_iter%i_%.fGeV_rejection.txt", n_iter, pt_binedge[0]);
#else
  snprintf(filename, 100, "PYTHIA/out_data/feasibility_data_iter%i_%.fGeV.txt", n_iter, pt_binedge[0]);
#endif
  
  // Check the file exists
  std::ifstream prevdata(filename);
  if (prevdata.fail()) {
    std::cout << "No pre-existing data file! Creating a new data file for N_iter = \033[32m" << n_iter << "\033[39m" << std::endl;
    return;
  }
  
  // Gather existing data
  std::string line;
  while (getline(prevdata, line)) {
    if (line[0] == '#') continue; //Comment line
    if (_debug) std::cout << "Line to read :: " << line << std::endl;
    
    // Get formatted data from table
    std::stringstream linestream(line);
    std::string entry;
    getline(linestream, entry, '\t');
    if (_debug) std::cout << "s/n :: " << entry;
    double signal_noise = std::stod(entry);
    getline(linestream, entry, '\t');
    if (_debug) std::cout << ", jc :: " << entry;
    double njet = std::stod(entry);
    getline(linestream, entry);
    if (_debug) std::cout << ", p_val :: " << entry << std::endl;
    double p_value = std::stod(entry);
    
    // Check if this data is relevant to the table
    int bin_sn = hist_to_fill->GetXaxis()->FindBin(signal_noise);
    int bin_nj = hist_to_fill->GetYaxis()->FindBin(njet);
    if (std::fabs(std::log(signal_noise) - 
                  std::log(hist_to_fill->GetXaxis()->GetBinCenter(bin_sn)))
          <= std::fabs(bin_matching_resolution) &&
        std::fabs(std::log(njet) -
                  std::log(hist_to_fill->GetYaxis()->GetBinCenter(bin_nj)))
          <= std::fabs(bin_matching_resolution)) {
      std::cout << "File bin (s/n,Njet)[log] = (" << std::log(signal_noise) << ',' << std::log(njet) << ")[log] matched to bin with center (";
      std::cout << std::log(hist_to_fill->GetXaxis()->GetBinCenter(bin_sn)) << ',' << std::log(hist_to_fill->GetYaxis()->GetBinCenter(bin_nj)) << ")[log], ";
      std::cout << "p_value = " << p_value << std::endl;
      
      hist_to_fill->SetBinContent(bin_sn, bin_nj, p_value);
    } // End of successful bin match
    
    // Append data to linked list
    appendNode(signal_noise, njet, p_value);
  }// End of bin match/data gathering
  
  return;
}// End of make_feasibility_contour::getPreviouslyMeasuredData



// Write linked list to file
void writeDataToFile() {
  
  // Overwrite file
  char filename[100];
#ifdef rejection_test
  snprintf(filename, 100, "PYTHIA/out_data/feasibility_data_iter%i_%.fGeV_rejection.txt", n_iter, pt_binedge[0]);
#else
  snprintf(filename, 100, "PYTHIA/out_data/feasibility_data_iter%i_%.fGeV.txt", n_iter, pt_binedge[0]);
#endif
  std::ofstream outfile(filename);
  outfile << "# Feasibility test for " << n_iter << " iterations\n#\n# S/N\tNjet\tmean p-value" << std::endl;
  
  // Write all available data to file
  feasibility_data_node* write = head;
  while (write != NULL) {
    outfile << write->signal_noise << '\t';
    outfile << write->njet << '\t';
    outfile << write->p_value << std::endl;
    
    write = write->next;
  }// End of data write loop
  
  outfile.close();
  return;
}



// Find the feasibility of the measurement
double getFeasibilitySingleBin(double signal_to_noise, double Njet) {
  // If input Njet is the total number of jets
//  int njet[3] = {
//    static_cast<int>(Njet),
//    static_cast<int>((signal_to_noise / (1. + signal_to_noise)) * Njet),
//    static_cast<int>((1.               / (1. + signal_to_noise)) * Njet)
//  };
  // If the input Njet is the total number of D0 jets
  int njet[3] = {
    static_cast<int>(((1. + signal_to_noise) / signal_to_noise) * Njet),
    static_cast<int>(Njet),
    static_cast<int>(Njet / signal_to_noise)
  };
  
  
  char ptbin_inf_string[2][20] = {"", "#infty)"};
  snprintf(ptbin_inf_string[0], 20, "%.f]",pt_binedge[1]);
  std::cout << "Njet {tot,sig,bkg} :: {" << njet[0] << ',' << njet[1] << ',' << njet[2] << "}";
  std::cout << " from S/N = \033[31m" << signal_to_noise*1e3 << "*10^{-3}\033[39m, Njet = \033[31m" << Njet << "\033[39m." << std::endl;
  
 
  
  // Histogram of deviation
  // Source the hist when the range is known better
  TH1D* hist_comparison_metric;
  const bool run_logscale = false; /* (switch_comparison_mode >= 2) */
  double comparison_range[2] = {50, 0};
  
  
  
  //==========================================Simulated Experiment Generation
  
  
  // Begin event loop
  double integral_scale = gEEC[1]->Integral();
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
        hist_sideband->Fill(gEEC[2]->GetRandom());
      }
    }// End of background fill
    
    
    // Perform the subtraction of the sideband and renormalize
    TH1D* diff = static_cast<TH1D*>(gEEC[0]->Clone());
    diff->Reset();
    diff->Add(hist_reco_EEC, hist_sideband, 1, -1);
    diff->Scale(integral_scale/(diff->Integral()));
    
    // Analyze the difference in the subtracted histograms
    double comparison_metric = 0;
    char printmode[2] = {' ','P'};
    switch (switch_comparison_mode) {
      case 0: default: // Kolmogorov-Smirnov Test
        comparison_metric = diff->KolmogorovTest(gEEC[1], "WW M");
        break;
      case 1: // Chi2 test
        comparison_metric = diff->Chi2Test(gEEC[1], "CHI2/NDF");
        break;
      case 2:// Kolmogorov-Smirnov Test (p-value)
#ifdef rejection_test
        comparison_metric = diff->KolmogorovTest(gEEC[0]);
#else
        comparison_metric = diff->KolmogorovTest(gEEC[1]);
#endif
        break;
      case 3: // Chi2 test (p-value)
#ifdef rejection_test
        comparison_metric = diff->Chi2Test(gEEC[0], Form("WW%c",printmode[_debug_global]));
#else
        comparison_metric = diff->Chi2Test(gEEC[1], Form("WW%c",printmode[_debug_global]));
#endif
        break;
    }// End of comparison evaluation
    
    
    // Check for max if not using a logscale method
    if (!run_logscale && i_iter == 0) {
      hist_comparison_metric = new TH1D(Form("hist_comparison_metric_sn%.2f_Njet%.f", signal_to_noise, Njet),
                                        Form(";%s;Simulated Experiment Counts",comparison_mode_string[switch_comparison_mode]),
                                        15, 0, 1.00001);
    } else if (i_iter == 0) {
      // If using log scale relevant methods, compute the binedges
      // of the p-value plots using the first value as a rough estimate
      int nbin_logrange = 100;
      std::vector<double> logbinedge = getUniformBinningLogScale(nbin_logrange, comparison_metric / 20, TMath::Min(1., comparison_metric * 20));
      hist_comparison_metric = new TH1D(Form("hist_comparison_metric_sn%.2f_Njet%.f", signal_to_noise, Njet),
                                        Form(";%s;Simulated Experiment Counts",comparison_mode_string[switch_comparison_mode]),
                                        nbin_logrange, logbinedge.data());
      
    }// End of hist setup
    hist_comparison_metric->Fill(comparison_metric);
    
    
    
    
    // Draw the iterations occasionally to compare
    if (i_iter % (n_iter / 10) == 0) {
      gPad->SetLogx();
      hist_reco_EEC->SetLineColor(plot_colors[2]);
      hist_sideband->SetLineColor(plot_colors[3]);
      hist_reco_EEC->Draw("hist");
      hist_sideband->Draw("hist same");
      gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/iteration_check/basehist_iter_%i.pdf",i_iter / (n_iter / 10) ) );
      
      
      std::cout << "Finished with simul. experiment generation #\033[32m" << i_iter << "\033[39m" << std::endl;
 
      
#ifdef rejection_test
      diff->GetYaxis()->SetRangeUser(0, 1.2*TMath::Max(diff->GetMaximum(), gEEC[0]->GetMaximum()));
      diff->Draw("hist");
      diff->SetLineColor(kViolet);
      gEEC[0]->Draw("hist same");
#else
      diff->GetYaxis()->SetRangeUser(0, 1.2*TMath::Max(diff->GetMaximum(), gEEC[1]->GetMaximum()));
      diff->Draw("hist");
      gEEC[1]->Draw("hist same");
#endif
      
      // Legend
      TLegend* leg_iter;
      if (thumbnail_mode) leg_iter = new TLegend(0.20, 0.63, 0.55, 0.85);
      else                leg_iter = new TLegend(0.16, 0.50, 0.33, 0.70);
      leg_iter->SetLineWidth(0);
      leg_iter->AddEntry(diff, "Sampled", "l");
      leg_iter->AddEntry(gEEC[1], "D0-Tagged", "l");
      leg_iter->Draw();
      
      // text if not plotting thumbnail/small plot mode
      if (!thumbnail_mode) {
        // Generation text
        drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", ((double) sqrt_s)/1000), gPad->GetLeftMargin(), 0.95);
        drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
        drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
        drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
        
        // Text for the comparison
        drawText(Form("N^{jet} = %i, N^{signal} = %i", njet[0], njet[1]), 0.05 + gPad->GetLeftMargin(), 0.925 - gPad->GetTopMargin());
        drawText(Form("%s :: #color[2]{%.4f}", comparison_mode_string[switch_comparison_mode], comparison_metric), 0.05 + gPad->GetLeftMargin(), 0.88 - gPad->GetTopMargin());
        drawText(Form("Iteration #bf{%i}", i_iter), 0.05 + gPad->GetLeftMargin(), 0.83 - gPad->GetTopMargin());
      } else {
        gCanvas->SetBottomMargin(0.10);
        drawText(Form("#hat{#it{p}}_{T} #in  [%.0f,%s, #it{p}_{T}^{jet} #in  [%.0f,%0.f]",
                      pt_binedge[0],ptbin_inf_string[pt_binedge[1]==-1],EEC_pt_binedge[0],EEC_pt_binedge[1]),
                 0.05 + gPad->GetLeftMargin(), 0.87, false, kBlack, 0.065);
        drawText(Form("#chi^{2} #it{p}-value: #color[2]{%.4f}", comparison_metric),
                 0.20, 0.55, false, kBlack, 0.065);
      }
      
      
      gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/iteration_check/iter_%i.pdf",i_iter / (n_iter / 10) ) );
    }// End of iteration drawing
    
  }// End of MC loop
  
  //==========================================Evaluation and final feasibility computation
  
  // Get average of the distribution
  double mean_p_value = hist_comparison_metric->GetMean();
  
  
  // Plot the histogram of the comparison metric results
  if (run_logscale) gPad->SetLogx();
  else              gPad->SetLogx(0);
  hist_comparison_metric->GetYaxis()->SetRangeUser(0, hist_comparison_metric->GetMaximum() * 1.2);
  // Enlarge axes if drawing in thumbnail mode
  if (thumbnail_mode) {
    gCanvas->SetBottomMargin(0.13);
    hist_comparison_metric->GetXaxis()->SetLabelSize(0.05);
    hist_comparison_metric->GetYaxis()->SetLabelSize(0.05);
    hist_comparison_metric->GetXaxis()->SetTitleSize(0.06);
    hist_comparison_metric->GetYaxis()->SetTitleSize(0.06);
  }
  
  hist_comparison_metric->GetYaxis()->SetTitleOffset(1.15);
  hist_comparison_metric->SetLineColor(kBlack);
  hist_comparison_metric->SetFillColorAlpha(kBlack, 0.15);
  hist_comparison_metric->Draw("hist");
  
  // text if not plotting thumbnail/small plot mode
  if (!thumbnail_mode) {
    // Generation text
    
    drawText(Form("#it{#bf{PYTHIA}} #sqrt{s_{NN}} = %.2f TeV", ((double) sqrt_s)/1000), gPad->GetLeftMargin(), 0.95);
    drawText(Form("#it{p}_{T}^{full jet} in range [%.0f,%0.f]",EEC_pt_binedge[0],EEC_pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
    drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jetRadius), 1.0 - gPad->GetRightMargin(), 0.90, true);
    drawText(Form("EEC Formed with %s Constituents",chargejetstringformal[use_charged_jets]), gPad->GetLeftMargin(), 0.90);
    
    // Simulation text
    drawText(Form("N^{jet} = %i, N^{signal} = %i", njet[0], njet[1]), 0.05 + gPad->GetLeftMargin(), 0.925 - gPad->GetTopMargin());
    drawText(Form("Total iterations: %i", n_iter), 0.05 + gPad->GetLeftMargin(), 0.88 - gPad->GetTopMargin());
  } else {
    drawText(Form("#hat{#it{p}}_{T} #in  [%.0f,%s, #it{p}_{T}^{jet} #in  [%.0f,%0.f]",
                  pt_binedge[0],ptbin_inf_string[pt_binedge[1]==-1],EEC_pt_binedge[0],EEC_pt_binedge[1]),
             0.05 + gPad->GetLeftMargin(), 0.87, false, kBlack, 0.065);
    
    // Simulation text
    drawText(Form("N^{jet} = %i, N^{signal} = %i", njet[0], njet[1]), 0.05 + gPad->GetLeftMargin(), 0.79 - gPad->GetTopMargin(), false, kBlack, 0.060);
    drawText(Form("Total iterations: %i", n_iter), 0.05 + gPad->GetLeftMargin(), 0.70 - gPad->GetTopMargin(), false, kBlack, 0.060);
  }
  
  
  // Draw line for average p-valie
  TLine* pval_line = new TLine();
  pval_line->SetLineWidth(2);
  pval_line->SetLineColor(plot_colors[1]);
  pval_line->DrawLine(mean_p_value, 0, mean_p_value, hist_comparison_metric->GetMaximum());
  drawText(Form("Avg. #it{p}-value = %.4f", mean_p_value),
           gPad->GetLeftMargin() + (1.0 - gPad->GetLeftMargin() - gPad->GetRightMargin())*(mean_p_value - 0.02 + 0.06*(mean_p_value <= 0.3)),
           0.25 + gPad->GetBottomMargin(), false, plot_colors[1])->SetTextAngle(90);
  
  
  gCanvas->SaveAs(Form("PYTHIA/out_plots/feasibility/metric_plots/feasibility_SN%.2f_Njet%i.pdf",signal_to_noise, static_cast<int>(Njet)));
  
  std::cout << "Result " << chargejetstringformal[use_charged_jets] << " :: " << std::scientific << hist_comparison_metric->GetMean() << std::fixed << std::endl;
  
  // Return the average performance
  return mean_p_value;
}// End of make_feasibility_contour::getFeasibilitySingleBin
