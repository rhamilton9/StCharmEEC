// ROOT macro. Should be soon incorperated into pythia_summary
// Just getting an early sense of the spectra shape

#include "../../utils/root_draw_tools.h"
#include "../../utils/hist_tools.h"
#include "config_local.h"
#include "ptbin.h"

const bool debug = false;

// Color and aesthetic settings for root plots
Int_t palette_full = kCherry;
Int_t palette_ch = kLake;
Int_t markerStyle_full = 20;
Int_t markerStyle_ch = 33;
Float_t markerSize_full = 0.75;
Float_t markerSize_ch = 1.0;

// bool local_flag_xsec_scale = false;

TPad* pad_main[2];
TPad* pads_left[5];

void spectra_summary() {
  // pTHat bin setup
  const int total_threads = getNbins(sqrt_s);
  float collision_energy_TeV = (float) sqrt_s / 1000.;
  char gen_info_string[100];
  snprintf(gen_info_string, 100, "%.2fTeV_r%.1fjet", collision_energy_TeV, jet_radius);
  
  // Jet spectra setup
  TH1D* hist_fulljet[total_threads];
  TH1D* hist_charjet[total_threads];
  TH1D* hist_alljet[2];
  int max_pt;
  if (sqrt_s < 900) max_pt = 60;
  else              max_pt = 360;
  
  hist_alljet[0] = new TH1D("hist_pT_fulljet",";p_{T};1/N^{event} d^{2}N^{jet}/dp_{T}d#eta",20, 0, max_pt);
  hist_alljet[1] = new TH1D("hist_pT_charjet",";p_{T};1/N^{event} d^{2}N^{charged jet}/dp_{T}d#eta",20, 0, max_pt);
  
  // Counters for the jets counted in each case
  // Each variable is indexed [fulljet, chargedjet]
  // Updated 6/7/25 to work for files with all processes included
  double final_total_jets[2] = {0,0};
  double D0_jets_gluonsourced[2] = {0,0};
  double D0_jets_lightsourced[2] = {0,0};
  double D0_jets_charmsourced[2] = {0,0};
  double D0_jets_beautysourced[2] = {0,0};
  double final_D0_fraction[2] = {0,0};
  double final_missed_charm_fraction[2] = {0,0};
  
  // Cross-section setup
  TH1D* hist_xsec_given;     // From Isaac, Laura
  TH1D* hist_xsec_pythia[5]; // {total, light, gluon, charm, beauty}
  TH1D* hist_charmfrac[3];   // Charm source fractions
  TH1D* hist_beautyfrac[3];  // Beauty source fractions
  
  char xsec_type_string[5][10] = {"total","light","gluon","charm","beauty"};
  char xsec_type_string_short[5][5] = {"tot","uds","g","c","b"};
  char xsec_hf_source[3][7] = {"ccbar","cg","cq"};
  int markerStyle[5] = {53, 33, 47, 20, 21};
  int flavorColor[5] = {TColor::GetColor("#3B5A84"),TColor::GetColor("#C6D8A0"), // Lighter yellow "#e7de9f" darker "#ccc489"
    TColor::GetColor("#97BE9D"),TColor::GetColor("#6FA0A3"),TColor::GetColor("#537993")};
  
  // set up histograms for cross sections
  hist_xsec_given = new TH1D("hist_xsec_ptHatBins_given",";p_{T};#sigma(p_{T}^{Hard}) [mb]",total_threads,0,total_threads);
  hist_xsec_given->SetLineColor(TColor::GetColor("#1A1E3E")); // Ernst 0
  int color_index[3] = {3, 2, 1};
  int bar_index[5] = {2, 0, 1, 3, 4};
  for (int f = 0; f < 5; ++f) {
    // {total, light, gluon, charm, beauty}
    hist_xsec_pythia[f] = new TH1D(Form("hist_xsec_ptHatBins_%s",xsec_type_string[f]),
                                   Form(";p_{T};d^{2}#sigma^{%s}/dp_{T}^{Hard}d#eta [mb]",xsec_type_string_short[f]),
                                   total_threads,0,total_threads);
    hist_xsec_pythia[f]->SetLineColor(flavorColor[f]);
    hist_xsec_pythia[f]->SetMarkerColor(flavorColor[f]);
    hist_xsec_pythia[f]->SetMarkerStyle(markerStyle[f]);
    hist_xsec_pythia[f]->SetMarkerStyle(markerStyle[f]);
    hist_xsec_pythia[f]->SetBarWidth(1./5);
    hist_xsec_pythia[f]->SetBarOffset(bar_index[f]/5.);
    if (f >= 3) continue;
    hist_charmfrac[f] = new TH1D(Form("hist_cfrac_source:%s",xsec_hf_source[f]),
                                 ";p_{T};Charm source process",
                                 total_threads,0,total_threads);
    hist_beautyfrac[f] = new TH1D(Form("hist_bfrac_source:%s",xsec_hf_source[f]),
                                  ";p_{T};Beauty source process",
                                  total_threads,0,total_threads);
    hist_charmfrac[f]->SetFillColor(flavorColor[color_index[f]]);
    hist_beautyfrac[f]->SetFillColor(flavorColor[color_index[f]]);
  }
  hist_beautyfrac[0]->SetFillColor(flavorColor[4]);
  
  // TFile I/O
  const char process_label[5][5] = {"uds","cc","bb","gg","all"};
  const char process_legend[4][5] = {"uds","g","c","b"};
  TFile* write_file = new TFile(Form("out_data/%s/pythpp_%s_%s.root",gen_info_string,process_label[flavor_flag],gen_info_string),"recreate");
  
  // Merged TTree for storing jet information in each event
  int                     pt_hard_bin;        // Which bin this event belongs to
  float                   pt_bin_xsec;        // Corresponding xsec of that bin
  int                     process_merge;      // Process at PYTHIA truth level
  int                     pcharm_merge;       // process charm multiplicity
//  int                     pbeauty_merge;      // process beauty multiplicity
  int                     mult_merge;         // Total event track multiplicity
  int                     cmult_merge;        // Event charged track multiplicity
  float                   delta_E_merge;      // Deviation from energy conservation
  int                     ntot_jet_merge;     // Number of full jets in an event
  std::vector<int>        ntot_cst_merge;     // Array of full jet constituent multiplicties
  std::vector<float>      jet_pT_merge;       // Array of full jet p_T
  std::vector<float>      jet_y_merge;        // Array of full jet rapidity
  std::vector<float>      jet_eta_merge;      // Array of full jet pseudorapidity
  std::vector<float>      jet_phi_merge;      // Array of full jet azimuth
  std::vector<float>      jet_area_merge;     // Array of full jet area
  int                     ntot_cjet_merge;    // Number of charged jets in an event
  std::vector<int>        ntot_ccst_merge;    // Array of charged jet constituent multiplicities
  std::vector<float>      chjet_pT_merge;     // Array of charged jet p_T
  std::vector<float>      chjet_y_merge;      // Array of charged jet rapidity
  std::vector<float>      chjet_eta_merge;    // Array of charged jet pseudorapidity
  std::vector<float>      chjet_phi_merge;    // Array of charged jet azimuth
  std::vector<float>      chjet_area_merge;   // Array of charged jet area
  
  // Which jets contain D0 mesons?
  std::vector<int>        jet_charm_status_merge;
  std::vector<int>        chjet_charm_status_merge;
  // Note:
  //    0: No charm in jet history
  //    1: Contains reconstructed D0
  //    2: Contains another charmed hadron
  
  TTree* pythia_tree_merged = new TTree("pythia_tree","pythia_tree");
  pythia_tree_merged->Branch("pt_bin",     &pt_hard_bin);
  pythia_tree_merged->Branch("bin_xsec",   &pt_bin_xsec);
  pythia_tree_merged->Branch("process",    &process_merge);
  pythia_tree_merged->Branch("pcharm",     &pcharm_merge);
//  pythia_tree_merged->Branch("pbeauty",    &pbeauty_merge);
  pythia_tree_merged->Branch("ntotal",     &mult_merge);
  pythia_tree_merged->Branch("ncharge",    &cmult_merge);
  pythia_tree_merged->Branch("delta_E",    &delta_E_merge);
  pythia_tree_merged->Branch("njet",       &ntot_jet_merge);
  pythia_tree_merged->Branch("nchjet",     &ntot_cjet_merge);
  pythia_tree_merged->Branch("jet_n",      &ntot_cst_merge);
  pythia_tree_merged->Branch("chjet_n",    &ntot_ccst_merge);
  pythia_tree_merged->Branch("jet_pT",     &jet_pT_merge);
  pythia_tree_merged->Branch("chjet_pT",   &chjet_pT_merge);
  pythia_tree_merged->Branch("jet_y",      &jet_y_merge);
  pythia_tree_merged->Branch("chjet_y",    &chjet_y_merge);
  pythia_tree_merged->Branch("jet_eta",    &jet_eta_merge);
  pythia_tree_merged->Branch("chjet_eta",  &chjet_eta_merge);
  pythia_tree_merged->Branch("jet_phi",    &jet_phi_merge);
  pythia_tree_merged->Branch("chjet_phi",  &chjet_phi_merge);
  pythia_tree_merged->Branch("jet_area",   &jet_area_merge);
  pythia_tree_merged->Branch("chjet_area", &chjet_area_merge);
  pythia_tree_merged->Branch("jet_charm_status",   &jet_charm_status_merge);
  pythia_tree_merged->Branch("chjet_charm_status", &chjet_charm_status_merge);
  
  // Merged cross section tree
  double xsec_total_given;    // Isaac+Laura's cross section values
  double xsec_total_final;    // Cross section for any hard scattering at this pTHat
  double xsec_light_final;    // Total light flavor jet cross section (uds)
  double xsec_gluon_final;    // Total gluon jet cross section
  double xsec_charm_final;    // Total charm jet cross section
  double xsec_beauty_final;   // Total beauty jet cross section
  double xsec_total_error;    // Error on full cross section for any hard scattering
  double xsec_light_error;    // Error on light flavor jet cross section
  double xsec_gluon_error;    // Error on gluon jet cross section
  double xsec_charm_error;    // Error on charm jet cross section
  double xsec_beauty_error;   // Error on beauty jet cross section
  
  double frac_charm_fromcc;   // Fraction of charm jets which come from ccbar production
  double frac_charm_fromgq;   // Fraction of charm jets which come from cg scattering
  double frac_charm_fromqQ;   // Fraction of charm jets which come from cq scattering
  double frac_beauty_frombb;  // Fraction of beauty jets which come from bbbar production
  double frac_beauty_fromgq;  // Fraction of beauty jets which come from bg scattering
  double frac_beauty_fromqQ;  // Fraction of beauty jets which come from bq scattering
  
  
  TTree* pythia_xsec_tree = new TTree("xsec_tree","xsec_tree");
  pythia_xsec_tree->Branch("given",     &xsec_total_given);
  pythia_xsec_tree->Branch("xtotal",    &xsec_total_final);
  pythia_xsec_tree->Branch("xlight",    &xsec_light_final);
  pythia_xsec_tree->Branch("xgluon",    &xsec_gluon_final);
  pythia_xsec_tree->Branch("xcharm",    &xsec_charm_final);
  pythia_xsec_tree->Branch("xbeauty",   &xsec_beauty_final);
  pythia_xsec_tree->Branch("etotal",    &xsec_total_error);
  pythia_xsec_tree->Branch("elight",    &xsec_light_error);
  pythia_xsec_tree->Branch("egluon",    &xsec_gluon_error);
  pythia_xsec_tree->Branch("echarm",    &xsec_charm_error);
  pythia_xsec_tree->Branch("ebeauty",   &xsec_beauty_error);
  pythia_xsec_tree->Branch("fcharmT",   &frac_charm_fromcc);
  pythia_xsec_tree->Branch("fcharmQ",   &frac_charm_fromqQ);
  pythia_xsec_tree->Branch("fcharmG",   &frac_charm_fromgq);
  pythia_xsec_tree->Branch("fbeautyT",  &frac_beauty_frombb);
  pythia_xsec_tree->Branch("fbeautyQ",  &frac_beauty_fromqQ);
  pythia_xsec_tree->Branch("fbeautyG",  &frac_beauty_fromgq);
  
  
  // Print some information about this dataset
  std::cout << Form("PYTHIA pp sqrt_s_NN = %.2f TeV, pT_hadron >= %.2f GeV, |eta_hadron| < %.2f",
                    collision_energy_TeV, pTmin_hadron, max_eta) << std::endl;
  std::cout << Form("FastJet3 R = %.1f %s Jets, |eta_jet| < %.2f, jet_A >= %.2f*pi*R^2",
                    jet_radius,algo_string,max_eta - jet_radius,min_area/(3.141592*jet_radius*jet_radius)) << std::endl;
  std::cout << Form("pT_jet >= %.2f, pT_charged_jet >= %.2f, pT_leading >= %.2f",
                    jetcut_minpT_full,jetcut_minpT_chrg, pTmin_jetcore) << std::endl;
  std::cout << "----------------------------------------------------------------------------" << std::endl;
  
  // TCanvas Setup
  gStyle->SetOptStat(0);
  TCanvas* canvas_spectra = new TCanvas();
  canvas_spectra->SetCanvasSize(450,500);
  std::vector<std::vector<TPad*>> pads = divideFlush(canvas_spectra, 1, 2);
  TLegend* leg_pTHard = new TLegend(0.7, 0.35, 0.85, 0.95);
  leg_pTHard->SetLineWidth(0);
  
  TCanvas* canvas_xsec = new TCanvas();
  canvas_xsec->SetCanvasSize(350, 500);
  TPad* pads_xsec[3];
  pads_xsec[0] = buildPad("xsec_pad0", 0, 0.40, 1, 1.0,   0.1, 0.05, 0.01, 0.05);
  pads_xsec[1] = buildPad("xsec_pad0", 0, 0.20, 1, 0.4,   0.1, 0.05, 0.025, 0.25);
  pads_xsec[2] = buildPad("xsec_pad0", 0, 0.0, 1, 0.20,   0.1, 0.05, 0.25, 0.025);
  
  // Colors from root spectra
  gStyle->SetPalette(kBird);
  Int_t color_sample_fulljet[total_threads];
  for (Int_t i = 0; i < total_threads; ++i) color_sample_fulljet[i] = TColor::GetPalette().At(TColor::GetPalette().GetSize() * i/(1.2*total_threads));
  gStyle->SetPalette(kBird);
  Int_t color_sample_charjet[total_threads];
  for (Int_t i = 0; i < total_threads; ++i) color_sample_charjet[i] = TColor::GetPalette().At(TColor::GetPalette().GetSize() * i/(1.2*total_threads));
  
  double xsec_range[2] = {10, 1e-10};
  int nevent_min = 0;
  for (int this_thread = total_threads-1; this_thread >= 0; --this_thread) {
    
    //-----------------------------------------------Thread-level setup
    
    // Get pTHat Bin info for this generation
    pt_hard_bin = this_thread;
    pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
    const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
    pt_bin_xsec = settings.cross_section;
    
    // Set up histograms for jet spectra
    hist_fulljet[this_thread] = new TH1D(Form("hist_pT_fulljet_thread%i",this_thread),
                                         ";p_{T};1/N^{event} dN^{jet}/dp_{T}",
                                         20, 0, max_pt);
    hist_fulljet[this_thread]->SetLineColor(color_sample_fulljet[this_thread]);
    hist_fulljet[this_thread]->SetLineWidth(2);
    hist_fulljet[this_thread]->SetFillColorAlpha(color_sample_fulljet[this_thread], 0.2);
    hist_charjet[this_thread] = new TH1D(Form("hist_pT_charjet_thread%i",this_thread),
                                         ";p_{T};1/N^{event} dN^{charged jet}/dp_{T}",
                                         20, 0, max_pt);
    hist_charjet[this_thread]->SetLineColor(color_sample_charjet[this_thread]);
    hist_charjet[this_thread]->SetLineWidth(2);
    hist_charjet[this_thread]->SetFillColorAlpha(color_sample_charjet[this_thread], 0.2);

    
    // File for this thread
    TFile* current_file = new TFile(Form("out_data/%s/pythpp_%s_%s_thread%i.root",gen_info_string,process_label[flavor_flag],gen_info_string,this_thread), "read");
    if (current_file->IsZombie() || current_file->GetNkeys() == 0) {
      std::cout << "\033[31mWarning\033[39m in spectra_summary.c: Unable to read deta from thread " << this_thread << ". Continuing..." << std::endl;
      hist_fulljet[this_thread]->SetLineColor(kRed);
      hist_fulljet[this_thread]->SetFillColorAlpha(kRed, 0.2);
      leg_pTHard->AddEntry(hist_fulljet[this_thread], Form("[%.f,%.f]",pt_binedge[0],pt_binedge[1]), "f");
      continue;
    }
    
    
    //-----------------------------------------------Gather event-by-event data + spectra
    
    // Find event data tree in this pT bin
    char tree_name[50];
    snprintf(tree_name, 50, "pythia_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
    const int nevent = current_file->Get<TTree>(tree_name)->GetEntries();
    if (nevent_min == 0 || nevent < nevent_min) nevent_min = nevent;
    TTreeReader *reader = new TTreeReader(tree_name, current_file);
    if (!reader) {std::cout << "Warning: Tree not found in bin " << this_thread << std::endl; continue;}
    TTreeReaderValue<Int_t>                     process(*reader, "process");        // PYTHIA label for the hard scattering process
    TTreeReaderValue<Int_t>                     pcharm(*reader, "pcharm");          // Count of produced charm quarks in hard scattering
//    TTreeReaderValue<Int_t>                     pbeauty(*reader, "pbeauty");        // Count of produced beauty quarks in hard scattering
    TTreeReaderValue<Int_t>                     ntotal(*reader, "ntotal");          // Total event track multiplicity
    TTreeReaderValue<Int_t>                     ncharge(*reader, "ncharge");        // Event charged track multiplicity
    TTreeReaderValue<Float_t>                   deltaE(*reader, "delta_E");         // Deviation from energy conservation
    TTreeReaderValue<Int_t>                     njet(*reader, "njet");              // Number of full jets in an event
    TTreeReaderValue<Int_t>                     nchjet(*reader, "nchjet");          // Number of charged jets in an event
    TTreeReaderValue<std::vector<Int_t>>        jet_n(*reader, "jet_n");            // Array of full jet constituent multiplicties
    TTreeReaderValue<std::vector<Int_t>>        chjet_n(*reader, "chjet_n");        // Array of charged jet constituent multiplicities
    TTreeReaderValue<std::vector<Float_t>>      jet_pT(*reader, "jet_pT");          // Array of full jet p_T
    TTreeReaderValue<std::vector<Float_t>>      chjet_pT(*reader, "chjet_pT");      // Array of charged jet p_T
    TTreeReaderValue<std::vector<Float_t>>      jet_y(*reader, "jet_y");            // Array of full jet rapidity
    TTreeReaderValue<std::vector<Float_t>>      chjet_y(*reader, "chjet_y");        // Array of charged jet rapidity
    TTreeReaderValue<std::vector<Float_t>>      jet_eta(*reader, "jet_eta");        // Array of full jet pseudorapidity
    TTreeReaderValue<std::vector<Float_t>>      chjet_eta(*reader, "chjet_eta");    // Array of charged jet pseudorapidity
    TTreeReaderValue<std::vector<Float_t>>      jet_phi(*reader, "jet_phi");        // Array of full jet azimuth
    TTreeReaderValue<std::vector<Float_t>>      chjet_phi(*reader, "chjet_phi");    // Array of charged jet azimuth
    TTreeReaderValue<std::vector<Float_t>>      jet_area(*reader, "jet_area");      // Array of full jet areas
    TTreeReaderValue<std::vector<Float_t>>      chjet_area(*reader, "chjet_area");  // Array of charged jet areas
    
    TTreeReaderValue<std::vector<Int_t>>        jet_charm_status(*reader, "jet_charm_status");      // Array of full jet charm status
    TTreeReaderValue<std::vector<Int_t>>        chjet_charm_status(*reader, "chjet_charm_status");    // Array of charged jet charm status
    
    // D0/charm status tallies
    int total_reconstructable_charmfree[2] = {0,0};
    int total_missed_charm[2] = {0,0};
    int total_reconstructable_D0[2] = {0,0};
    double total_D0_fromgluon[2] = {0,0};
    double total_D0_fromlight[2] = {0,0};
    double total_D0_fromcharm[2] = {0,0};
    double total_D0_frombeauty[2] = {0,0};
    
    // Gather data
    int current_event = -1;
    int ntot_jet = 0;
    int ntot_cjet = 0;
    while (reader->Next()) {
      if (debug) std::cout << "Event " << ++current_event << " Njet = " << *njet << ", NCjet = " << *nchjet << std::endl;
      ntot_jet += *njet;
      ntot_cjet += *nchjet;
      
      // Assign data from current branh into the merged tree
      process_merge =     *process;
      pcharm_merge =      *pcharm;
//      pbeauty_merge =     *pbeauty;
      mult_merge =        *ntotal;
      cmult_merge =       *ncharge;
      delta_E_merge =     *deltaE;
      ntot_jet_merge =    *njet;
      ntot_cst_merge =    *jet_n;
      jet_pT_merge =      *jet_pT;
      jet_y_merge =       *jet_y;
      jet_eta_merge =     *jet_eta;
      jet_phi_merge =     *jet_phi;
      jet_area_merge =    *jet_area;
      jet_charm_status_merge = *jet_charm_status;
      ntot_cjet_merge =   *nchjet;
      ntot_ccst_merge =   *chjet_n;
      chjet_pT_merge =    *chjet_pT;
      chjet_y_merge =     *chjet_y;
      chjet_eta_merge =   *chjet_eta;
      chjet_phi_merge =   *chjet_phi;
      chjet_charm_status_merge = *chjet_charm_status;
      pythia_tree_merged->Fill();
      
      // Fill histograms
      for (int i_jet = 0; i_jet < *njet; ++i_jet) {
        hist_fulljet[this_thread]->Fill(jet_pT->at(i_jet));
      } for (int i_jet = 0; i_jet < *nchjet; ++i_jet) {
        hist_charjet[this_thread]->Fill(chjet_pT->at(i_jet));
      }
      
      // Tally charm data, for debugging/comparing to event displays
      if (debug) std::cout << "fulljet:";
      for (int i_jet = 0; i_jet < *njet; ++i_jet) {
        if (debug) std::cout << jet_charm_status->at(i_jet) << ' ';
        switch (jet_charm_status->at(i_jet)) {
          case 0: // no charm
            ++total_reconstructable_charmfree[0];
            break;
          case 1: // Reconstructable D0 tagged jet
            ++total_reconstructable_D0[0];
            switch (*process) {// Temp until PYTHIA output is jet level rather than event level.
              case 111: case 115: // final state gg
                total_D0_fromgluon[0] += 1;
                break;
              case 112: case 116: // final state ll
                total_D0_fromlight[0] += 1;
                break;
              case 113: // final state qg
//                total_D0_fromgluon[0] += 0.5;
//                total_D0_fromcharm[0] += 0.25;
//                total_D0_frombeauty[0] += 0.15;
                break;
              case 114: // final state qq
//                total_D0_fromlight[0] += 0.75;
//                total_D0_fromcharm[0] += 0.25;
                break;
              case 121: case 122:
                total_D0_fromcharm[0] += 1;
                break;
              case 123: case 124:
                total_D0_frombeauty[0] += 1;
                break;
            }
            
            break;
          default: // missed charm
            ++total_missed_charm[0];
        }
      }
      if (debug) std::cout << std::endl << "charjet:";
      for (int i_jet = 0; i_jet < *nchjet; ++i_jet) {
        if (debug) std::cout << chjet_charm_status->at(i_jet) << ' ';
        switch (chjet_charm_status->at(i_jet)) {
          case 0: // no charm
            ++total_reconstructable_charmfree[1];
            break;
          case 1: // Reconstructable D0 tagged jet
            ++total_reconstructable_D0[1];
            switch (*process) {// Temp until PYTHIA output is jet level rather than event level.
              case 111: case 115: // final state gg
                total_D0_fromgluon[1] += 1;
                break;
              case 112: case 116: // final state ll
                total_D0_fromlight[1] += 1;
                break;
              case 113: // final state qg
//                total_D0_fromgluon[1] += 0.5;
//                total_D0_fromcharm[1] += 0.25;
//                total_D0_frombeauty[1] += 0.15;
                break;
              case 114: // final state qq
//                total_D0_fromlight[1] += 0.75;
//                total_D0_fromcharm[1] += 0.25;
                break;
              case 121: case 122:
                total_D0_fromcharm[1] += 1;
                break;
              case 123: case 124:
                total_D0_frombeauty[1] += 1;
                break;
            }
            break;
          default: // missed charm
            ++total_missed_charm[1];
        }
      } if (debug) std::cout << std::endl;
    }// end of tree reader loop
    
    // Jet count info for this bin
    std::cout << Form("pTHat bin [%.f,%.f]    \tn_event = %i\tn_jet = %i\tn_charged_jet = %i", pt_binedge[0],pt_binedge[1],nevent,ntot_jet,ntot_cjet) << std::endl;
    
    // D0/charm information for this bin
    std::cout << "Jet Conservation (Thread" << this_thread << ") :: ";
    if (total_reconstructable_charmfree[0] + total_reconstructable_D0[0] + total_missed_charm[0] == ntot_jet)
      std::cout << "\033[32mO\033[39m" << std::endl;
    else
      std::cout << "\033[31mX\033[39m" << std::endl;
    std::cout << "Charged Jet Conservation (Thread" << this_thread << ") :: ";
    if (total_reconstructable_charmfree[1] + total_reconstructable_D0[1] + total_missed_charm[1] == ntot_cjet)
      std::cout << "\033[32mO\033[39m" << std::endl;
    else
      std::cout << "\033[31mX\033[39m" << std::endl;
    
    std::cout << "Tot\tD0\tXc" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << ntot_jet << "\t" << total_reconstructable_D0[0] << "\t" << total_missed_charm[0] << std::endl;
    std::cout << ntot_cjet << "\t" << total_reconstructable_D0[1] << "\t" << total_missed_charm[1] << std::endl;
  
    // xsec-weight the source production fractions
    final_total_jets[0] += ntot_jet * pt_bin_xsec / nevent;
    final_total_jets[1] += ntot_cjet * pt_bin_xsec / nevent;
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      final_D0_fraction[iCharge] += total_reconstructable_D0[iCharge] * pt_bin_xsec / nevent;
      final_missed_charm_fraction[iCharge] += total_missed_charm[iCharge] * pt_bin_xsec / nevent;
      D0_jets_gluonsourced[iCharge] += total_D0_fromgluon[iCharge] * pt_bin_xsec / nevent;
      D0_jets_lightsourced[iCharge] += total_D0_fromlight[iCharge] * pt_bin_xsec / nevent;
      D0_jets_charmsourced[iCharge] += total_D0_fromcharm[iCharge] * pt_bin_xsec / nevent;
      D0_jets_beautysourced[iCharge] += total_D0_frombeauty[iCharge] * pt_bin_xsec / nevent;
    }
    
    
    //-----------------------------------------------Gather overall cross section information
    
    // Find cross section data tree in this pT bin
    char xsec_tree_name[50];
    snprintf(xsec_tree_name, 50, "xsec_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
    TTreeReader *reader_xsec = new TTreeReader(xsec_tree_name, current_file);
    if (!reader_xsec) std::cout << "Warning: Xsec tree not found in bin " << this_thread << std::endl;
    TTreeReaderValue<Double_t>    read_xgiven(*reader_xsec, "given");
    TTreeReaderValue<Double_t>    read_xtotal(*reader_xsec, "xtotal");
    TTreeReaderValue<Double_t>    read_xlight(*reader_xsec, "xlight");
    TTreeReaderValue<Double_t>    read_xgluon(*reader_xsec, "xgluon");
    TTreeReaderValue<Double_t>    read_xcharm(*reader_xsec, "xcharm");
    TTreeReaderValue<Double_t>    read_xbeauty(*reader_xsec, "xbeauty");
    TTreeReaderValue<Double_t>    read_etotal(*reader_xsec, "etotal");
    TTreeReaderValue<Double_t>    read_elight(*reader_xsec, "elight");
    TTreeReaderValue<Double_t>    read_egluon(*reader_xsec, "egluon");
    TTreeReaderValue<Double_t>    read_echarm(*reader_xsec, "echarm");
    TTreeReaderValue<Double_t>    read_ebeauty(*reader_xsec, "ebeauty");
    TTreeReaderValue<Double_t>    read_fcharmT(*reader_xsec, "fcharmT");
    TTreeReaderValue<Double_t>    read_fcharmQ(*reader_xsec, "fcharmQ");
    TTreeReaderValue<Double_t>    read_fcharmG(*reader_xsec, "fcharmG");
    TTreeReaderValue<Double_t>    read_fbeautyT(*reader_xsec, "fbeautyT");
    TTreeReaderValue<Double_t>    read_fbeautyQ(*reader_xsec, "fbeautyQ");
    TTreeReaderValue<Double_t>    read_fbeautyG(*reader_xsec, "fbeautyG");
    
    reader_xsec->Next(); // each tree should have just one entry
    xsec_total_given   = *read_xgiven;
    xsec_total_final   = *read_xtotal;
    xsec_light_final   = *read_xlight;
    xsec_gluon_final   = *read_xgluon;
    xsec_charm_final   = *read_xcharm;
    xsec_beauty_final  = *read_xbeauty;
    xsec_total_error   = *read_etotal;
    xsec_light_error   = *read_elight;
    xsec_gluon_error   = *read_egluon;
    xsec_charm_error   = *read_echarm;
    xsec_beauty_error  = *read_ebeauty;
    if (*read_fcharmT <= 1)  frac_charm_fromcc  = *read_fcharmT;
    else frac_charm_fromcc = 0;
    frac_charm_fromgq  = *read_fcharmG;
    frac_charm_fromqQ  = *read_fcharmQ;
    frac_beauty_frombb = *read_fbeautyT;
    if (*read_fbeautyG <= 1) frac_beauty_fromgq = *read_fbeautyG;
    else frac_beauty_fromgq = 0;
    if (*read_fbeautyQ <= 1) frac_beauty_fromqQ = *read_fbeautyQ;
    else frac_beauty_fromqQ = 0;
    pythia_xsec_tree->Fill();
    
    hist_xsec_given->SetBinContent(this_thread+1, xsec_total_given);
    hist_xsec_given->SetBinError(this_thread+1, 0);
    if (this_thread == total_threads - 1)
      hist_xsec_given->GetXaxis()->SetBinLabel(this_thread+1, Form("%.f-#infty", pt_binedge[0]));
    else
      hist_xsec_given->GetXaxis()->SetBinLabel(this_thread+1, Form("%.f-%.f", pt_binedge[0],pt_binedge[1]));
    hist_xsec_pythia[0]->SetBinContent(this_thread+1, xsec_total_final);
    hist_xsec_pythia[1]->SetBinContent(this_thread+1, xsec_light_final);
    hist_xsec_pythia[2]->SetBinContent(this_thread+1, xsec_gluon_final);
    hist_xsec_pythia[3]->SetBinContent(this_thread+1, xsec_charm_final);
    hist_xsec_pythia[4]->SetBinContent(this_thread+1, xsec_beauty_final);
    hist_xsec_pythia[0]->SetBinError(this_thread+1, xsec_total_error);
    hist_xsec_pythia[1]->SetBinError(this_thread+1, xsec_light_error);
    hist_xsec_pythia[2]->SetBinError(this_thread+1, xsec_gluon_error);
    hist_xsec_pythia[3]->SetBinError(this_thread+1, xsec_charm_error);
    hist_xsec_pythia[4]->SetBinError(this_thread+1, xsec_beauty_error);
    
    hist_charmfrac[0]->SetBinContent(this_thread+1, frac_charm_fromcc);
    hist_charmfrac[1]->SetBinContent(this_thread+1, frac_charm_fromgq);
    hist_charmfrac[2]->SetBinContent(this_thread+1, frac_charm_fromqQ);
    hist_beautyfrac[0]->SetBinContent(this_thread+1, frac_beauty_frombb);
    hist_beautyfrac[1]->SetBinContent(this_thread+1, frac_beauty_fromgq);
    hist_beautyfrac[2]->SetBinContent(this_thread+1, frac_beauty_fromqQ);
    
    // Candidate max cross sections
    if (xsec_range[1] < xsec_total_final) xsec_range[1] = xsec_total_final;
    if (xsec_range[1] < xsec_total_given) xsec_range[1] = xsec_total_given;
    
    // candidate min cross sections
    if (xsec_range[0] > xsec_charm_final && xsec_charm_final != 0) xsec_range[0] = xsec_charm_final;
    if (xsec_range[0] > xsec_beauty_final && xsec_beauty_final != 0) xsec_range[0] = xsec_beauty_final;
    
    //------------------------------------------ Plotting the spectra in this bin with xsec weighting
    
    //Normalize by n_event
    float max_range;
    if (sqrt_s < 900) max_range = 0.15;
    else max_range = 0.025;
    hist_fulljet[this_thread]->Scale(1./((float)nevent), "width");
    hist_fulljet[this_thread]->GetYaxis()->SetRangeUser(0, max_range);
    hist_charjet[this_thread]->Scale(1./((float)nevent), "width");
    hist_charjet[this_thread]->GetYaxis()->SetRangeUser(0, max_range);
    
    // Add to net spectrum
    hist_alljet[0]->Add(hist_alljet[0], hist_fulljet[this_thread], 1, pt_bin_xsec);
    hist_alljet[1]->Add(hist_alljet[1], hist_charjet[this_thread], 1, pt_bin_xsec);
    
    // Plot this thread's xsec weighted spectra
    pads.at(0).at(0)->cd();
    if (this_thread == total_threads-1) {
      gPad->SetTicks(1,1);
      hist_fulljet[this_thread]->Draw("hist");
    } else
      hist_fulljet[this_thread]->Draw("hist same");
    
    
    pads.at(1).at(0)->cd();
    if (this_thread == total_threads-1) {
      gPad->SetTicks(1,1);
      hist_charjet[this_thread]->Draw("hist");
    } else
      hist_charjet[this_thread]->Draw("hist same");
    
    leg_pTHard->AddEntry(hist_fulljet[this_thread], Form("[%.f,%.f]",pt_binedge[0],pt_binedge[1]), "f");
    
    std::cout << "Finished with thread #" << this_thread << std::endl;
    
    // caused problems, not fully sure why
//    current_file->Close();
//    delete current_file;
  }// End of thread data loop
  
  //------------------------------------------Final state plotting (xsec and spectra legend)
  
  
  // Spectra plot labels
  pads.at(0).at(0)->cd();
  drawText(Form("Full Jet p_{T} #geq %.1f",jetcut_minpT_full), 0.15, 0.82, false);
  drawText(Form("#it{PYTHIA} #bf{pp} #sqrt{s_{NN}} = %.2f", (float) sqrt_s / 1000.), 0.9, 0.82, true);
  drawText(Form("#it{FastJet} %s R = %.1f", algo_string, jet_radius), 0.9, 0.76, true);
  drawText(Form("#left|#eta_{hadron}#right| #leq %.1f", max_eta), 0.9, 0.70, true);
  drawText(Form("p_{T}^{leading} #geq %.1f", pTmin_jetcore), 0.9, 0.62, true);
  
  pads.at(1).at(0)->cd();
  drawText(Form("Charged Jet p_{T} #geq %.1f",jetcut_minpT_chrg), 0.15, 0.92, false);
  leg_pTHard->SetHeader("p_{T}^{Hard} Bin\n");
  leg_pTHard->Draw();
  canvas_spectra->SaveAs(Form("out_plots/%s/spectra_normalized_%s.pdf",gen_info_string,gen_info_string));
  
  // Xsec plots
  
  // Make sure axes labels are the same
  for (int this_thread = 0; this_thread < total_threads; ++this_thread) {
    hist_charmfrac[0]->GetXaxis()->SetBinLabel(this_thread+1, hist_xsec_given->GetXaxis()->GetBinLabel(this_thread+1));
    hist_beautyfrac[0]->GetXaxis()->SetBinLabel(this_thread+1, hist_xsec_given->GetXaxis()->GetBinLabel(this_thread+1));
  }
  
  pads_xsec[0]->cd();
  gPad->SetLogy();
  gPad->SetTicks(1,1);
  TLegend* leg_xsec = new TLegend(0.73, 0.60, 0.90, 0.90);
  leg_xsec->SetLineWidth(0);
  hist_xsec_given->GetXaxis()->SetNdivisions(total_threads);
  hist_xsec_given->GetXaxis()->SetLabelOffset(10);
  hist_xsec_given->GetYaxis()->SetRangeUser(0.2*xsec_range[0], 5*xsec_range[1]);
  hist_xsec_given->GetYaxis()->SetNdivisions(5);
  hist_xsec_given->GetYaxis()->SetTitleOffset(1.3);
  hist_xsec_given->SetLineStyle(7);
  leg_xsec->AddEntry(hist_xsec_given, "Ref. Total", "l");
  hist_xsec_given->Draw("b");
  hist_xsec_pythia[0]->SetMarkerStyle(markerStyle[0]);
  hist_xsec_pythia[0]->SetMarkerColor(flavorColor[0]);
  leg_xsec->AddEntry(hist_xsec_pythia[0], "Total #sigma", "pl");
  hist_xsec_pythia[0]->Draw("hist p e0 same");
  for (int f = 1; f < 5; ++f) {
    leg_xsec->AddEntry(hist_xsec_pythia[f], Form("#it{pp} #rightarrow %s",process_legend[f-1]), "p");
    hist_xsec_pythia[f]->Draw("b p e x0 same");
  }
  
  drawText(Form("#it{PYTHIA #bf{pp}} #sqrt{s} = %.2f TeV", (float) sqrt_s / 1000.), 0.7, 0.87, true, kBlack);
  drawText(Form("N^{event} #geq %i per Bin", nevent_min), 0.7, 0.81, true, kBlack);
  leg_xsec->Draw();
  
  pads_xsec[1]->cd();
  gPad->SetTicks(1,1);
  THStack* stack_charm = new THStack("stack_charm", "");
  for (int s = 0; s < 3; ++s) stack_charm->Add(hist_charmfrac[s]);
  stack_charm->Draw("b");
  stack_charm->SetTitle(";p_{T}^{Hard} Bin;Charm Source Fraction");
  stack_charm->GetXaxis()->SetNdivisions(total_threads);
  stack_charm->GetXaxis()->SetLabelOffset(10);
  stack_charm->GetYaxis()->SetNdivisions(6);
  stack_charm->GetYaxis()->SetLabelSize(0.07);
  stack_charm->GetYaxis()->SetTitleSize(0.065);
  stack_charm->GetYaxis()->SetTitleOffset(0.6);
  stack_charm->SetMinimum(0);
  stack_charm->SetMaximum(1);
  stack_charm->Draw("b");
  
  // Add division legend
  TMarker* fraction_marker = new TMarker();
  fraction_marker->SetMarkerStyle(21);
  fraction_marker->SetMarkerSize(1.5);
  
  int nbin = stack_charm->GetXaxis()->GetNbins();
  float upperbound = stack_charm->GetXaxis()->GetBinUpEdge(nbin);
  
  
  fraction_marker->SetMarkerColor(flavorColor[1]);
  fraction_marker->DrawMarker(upperbound/11., 1.17);
  drawText("qg #rightarrow qg", 0.205, 0.85, false, kBlack, 0.07);
  
  fraction_marker->SetMarkerColor(flavorColor[2]);
  fraction_marker->DrawMarker(upperbound*3.5/11., 1.17);
  drawText("qq' #rightarrow qq'", 0.4, 0.85, false, kBlack, 0.07);
  
  fraction_marker->SetMarkerColor(flavorColor[3]);
  fraction_marker->DrawMarker(upperbound*6./11., 1.17);
  drawText("#splitline{gg #rightarrow c#bar{c}}{q#bar{q} #rightarrow c#bar{c}}", 0.595, 0.85, false, kBlack, 0.06);
  
  fraction_marker->SetMarkerColor(flavorColor[4]);
  fraction_marker->DrawMarker(upperbound*8.5/11., 1.17);
  drawText("#splitline{gg #rightarrow b#bar{b}}{q#bar{q} #rightarrow b#bar{b}}", 0.79, 0.85, false, kBlack, 0.06);
  
  pads_xsec[2]->cd();
  gPad->SetTicks(1,1);
  THStack* stack_beauty = new THStack("stack_beauty", hist_beautyfrac[0]->GetTitle());
  for (int s = 0; s < 3; ++s) stack_beauty->Add(hist_beautyfrac[s]);
  stack_beauty->Draw("b");
  stack_beauty->SetTitle(";p_{T}^{Hard} Bin;Beauty Source Fraction");
  stack_beauty->GetXaxis()->SetNdivisions(total_threads);
  stack_beauty->GetXaxis()->SetLabelSize(0.11);
  if (nbin > 15) stack_beauty->GetXaxis()->SetLabelOffset(stack_beauty->GetXaxis()->GetLabelOffset()*2.2);
  stack_beauty->GetXaxis()->SetTitleSize(0.08);
  stack_beauty->GetXaxis()->SetTitleOffset(1.35);
  stack_beauty->GetYaxis()->SetNdivisions(6);
  stack_beauty->GetYaxis()->SetLabelSize(0.07);
  stack_beauty->GetYaxis()->SetTitleSize(0.065);
  stack_beauty->GetYaxis()->SetTitleOffset(0.6);
  stack_beauty->SetMinimum(0);
  stack_beauty->SetMaximum(1);
  stack_beauty->Draw("b");
  
  
  canvas_xsec->SaveAs(Form("out_plots/%s/process_xsec_%s.pdf",gen_info_string,gen_info_string));
  
  
  // Draw the net sigma scaled spectrum (in black)
  for (int this_thread = 0; this_thread < total_threads; ++this_thread) {
    pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
    const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
    pt_bin_xsec = settings.cross_section;
    const int nevent = (int) (settings.nevent * nevent_scale_factor); // scaled event number for testing
    float collision_energy_TeV = (float) sqrt_s / 1000.;
    
    hist_fulljet[this_thread]->SetTitle(";p_{T};#sigma/N^{event} dN^{jet}/dp_{T}");
    hist_charjet[this_thread]->SetTitle(";p_{T};#sigma/N^{event} dN^{charged jet}/dp_{T}");
    
    // scale by cross section
    hist_fulljet[this_thread]->Scale(pt_bin_xsec);
    hist_fulljet[this_thread]->GetYaxis()->SetRangeUser(1e-9, 1.5e-2);
    hist_charjet[this_thread]->Scale(pt_bin_xsec);
    hist_charjet[this_thread]->GetYaxis()->SetRangeUser(1e-9, 1.5e-2);
    
    // Plot spectra
    pads.at(0).at(0)->cd();
    if (this_thread == 0) {
      gPad->SetLogy();
      hist_fulljet[this_thread]->Draw("hist");
    } else
      hist_fulljet[this_thread]->Draw("hist same");
    
    
    pads.at(1).at(0)->cd();
    if (this_thread == 0) {
      gPad->SetLogy();
      hist_charjet[this_thread]->Draw("hist");
    } else
      hist_charjet[this_thread]->Draw("hist same");
  }
  
  
  pads.at(0).at(0)->cd();
  hist_alljet[0]->SetLineColor(kBlack);
  hist_alljet[0]->Draw("hist same");
  drawText(Form("Full Jet p_{T} #geq %.1f",jetcut_minpT_full), 0.15, 0.82, false);
  drawText(Form("#it{PYTHIA} #bf{pp} #sqrt{s_{NN}} = %.2f", (float) sqrt_s / 1000.), 0.9, 0.82, true);
  drawText(Form("#it{FastJet} %s R = %.1f", algo_string, jet_radius), 0.9, 0.76, true);
  drawText(Form("#left|#eta_{hadron}#right| #leq %.1f", max_eta), 0.9, 0.70, true);
  drawText(Form("p_{T}^{leading} #geq %.1f", pTmin_jetcore), 0.9, 0.62, true);
  
  pads.at(1).at(0)->cd();
  hist_alljet[1]->SetLineColor(kBlack);
  hist_alljet[1]->Draw("hist same");
  leg_pTHard->AddEntry(hist_alljet[1], "Net Spectrum", "l");
  drawText(Form("Charged Jet p_{T} #geq %.1f",jetcut_minpT_chrg), 0.15, 0.92, false);
  leg_pTHard->SetHeader("p_{T}^{Hard} Bin\n");
  leg_pTHard->Draw();
  canvas_spectra->SaveAs(Form("out_plots/%s/spectra_sigmascaled_%s.pdf",gen_info_string,gen_info_string));
  
  
  //------------------------------------------Save/write/report to terminal
  
  // Save the completed spectra and full tree to a composite output file
  write_file->cd();
  pythia_tree_merged->Print();
  pythia_tree_merged->Write(pythia_tree_merged->GetName(), TObject::kOverwrite);
  
  pythia_xsec_tree->Print();
  pythia_xsec_tree->Write(pythia_xsec_tree->GetName(), TObject::kOverwrite);
  
  hist_alljet[0]->Write(hist_alljet[0]->GetName(), TObject::kOverwrite);
  for (int i_thread = 0; i_thread < total_threads; ++i_thread)
    hist_fulljet[i_thread]->Write(hist_fulljet[i_thread]->GetName(), TObject::kOverwrite);
  hist_alljet[1]->Write(hist_alljet[1]->GetName(), TObject::kOverwrite);
  for (int i_thread = 0; i_thread < total_threads; ++i_thread)
    hist_charjet[i_thread]->Write(hist_charjet[i_thread]->GetName(), TObject::kOverwrite);
  
  // Charm reco information
  std::cout << "\n\nFinal jet/charm reconstruction predictions for this system (generation_mode = " << flavor_flag << "):" << std::endl;
  std::cout << "Type\tFulljet\t\tCharged Jet" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "Tot\t" << final_total_jets[0] << "  \t" << final_total_jets[1] << std::endl;
  std::cout << "D0\t" << final_D0_fraction[0] << "  \t" << final_D0_fraction[1] << std::endl;
  std::cout << "Xc\t" << final_missed_charm_fraction[0] << "  \t" << final_missed_charm_fraction[1] << std::endl;
  std::cout << "--------------------------------------------Relevant ratios:" << std::endl;
  std::cout << "D0:\t\t" << final_D0_fraction[0]/final_total_jets[0] << "  \t" << final_D0_fraction[1]/final_total_jets[1] << std::endl;
  std::cout << "Missed Charm:\t" << final_missed_charm_fraction[0]/final_total_jets[0] << "  \t" << final_missed_charm_fraction[1]/final_total_jets[1] << std::endl;
  
  
  double frac_baseline_fulljet = D0_jets_gluonsourced[0] + D0_jets_lightsourced[0] + D0_jets_charmsourced[0] + D0_jets_beautysourced[0];
  double frac_baseline_charjet = D0_jets_gluonsourced[1] + D0_jets_lightsourced[1] + D0_jets_charmsourced[1] + D0_jets_beautysourced[1];
  
  std::cout << "\n\nD0 jet source information with updated xsec and B feeddown (generation_mode = " << flavor_flag << "):" << std::endl;
  std::cout << "Type\tFulljet\t\tCharged Jet" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  std::cout << "g\t" <<D0_jets_gluonsourced[0]/frac_baseline_fulljet << "  \t" <<  D0_jets_gluonsourced[1]/frac_baseline_charjet << std::endl;
  std::cout << "uds\t" <<D0_jets_lightsourced[0]/frac_baseline_fulljet << "  \t" <<  D0_jets_lightsourced[1]/frac_baseline_charjet << std::endl;
  std::cout << "c\t" <<D0_jets_charmsourced[0]/frac_baseline_fulljet << "  \t" <<  D0_jets_charmsourced[1]/frac_baseline_charjet << std::endl;
  std::cout << "b\t" <<D0_jets_beautysourced[0]/frac_baseline_fulljet << "  \t" <<  D0_jets_beautysourced[1]/frac_baseline_charjet << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  
  delete write_file;
  return;
}
