/// A ROOT macro for producing EEC plots from output TTrees gathering PYTHIA data.
/// Written by R J Hamilton on Oct. 3, 2024
///
/// Modified 6/6/2025 to include more sensible PYTHIA inputs.

#include "../config.hpp"

#include "../utils/root_draw_tools.h"

void pythia_EEC() {
  // File I/O
  double snn_tev = ((double) sqrt_s)/1000;
  char gen_info_string[50];
  snprintf(gen_info_string, 50, "%.2fTeV_r%.1fjet", snn_tev, jetRadius);
  TFile *datafile = new TFile(Form("PYTHIA/out_data/%s/pythpp_all_%s.root", gen_info_string,gen_info_string));
  
  
  // Set up EEC histograms
  char cjet_string[2][3] = {"", "c"};
  char cjet_longstring[2][10] = {"full", "charged"};
  Double_t binedge_array[11] = {
    2e-3, 5e-3, 7e-3, 1e-2, 2e-2, 4e-2,
    7e-2, 1e-1, 2e-1, 4e-1, 7e-1
  };
  Int_t jetpt_binedge[4] = {20, 40, 60, 80};
  Int_t color_flavor[4] = {kGreen+2, kRed+2, kBlue+2, kGray+2};
  char char_flavor[4] = {'0', 'c', 'b', 'g'};
  Int_t marker_flavor[4] = {24, 25, 28, 30};
  TH1D* EEC_hist[2][4][3];
  Int_t count_jets[2][4][3];
  for (int iCharge = 0; iCharge < 2; ++iCharge) {
    for (int iFlavor = 0; iFlavor < 4; ++iFlavor) {
      for (int ipt = 0; ipt < 3; ++ipt) {
        EEC_hist[iCharge][iFlavor][ipt] = new TH1D(Form("hist_%sjet_%cEEC_pT%i-%i",cjet_string,char_flavor[iFlavor],jetpt_binedge[ipt],jetpt_binedge[ipt+1]),
                                                   ";R_{L};1/N_{jet} d#Sigma/dR_{L}", 10, binedge_array);
        EEC_hist[iCharge][iFlavor][ipt]->SetLineColor(color_flavor[iFlavor]);
        EEC_hist[iCharge][iFlavor][ipt]->SetMarkerColor(color_flavor[iFlavor]);
        EEC_hist[iCharge][iFlavor][ipt]->SetMarkerStyle(marker_flavor[iFlavor]);
        count_jets[iCharge][iFlavor][ipt] = 0;
      }
    }
  }
  
  // Canvas setup
  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas();
  canvas->SetCanvasSize(400, 800);
  std::vector<std::vector<TPad*>> pads = divideFlush(gPad, 1, 3, 0.1, 0.05, 0.05, 0.02);
  
  TCanvas* canvas_solo = new TCanvas();
  canvas_solo->SetCanvasSize(400, 300);
  canvas_solo->SetRightMargin(0.05);
  canvas_solo->SetTopMargin(0.05);
  
  const int nevent = datafile->Get<TTree>("pythia_tree")->GetBranch("ntotal")->GetEntries();
  TTreeReader* reader_all = new TTreeReader("pythia_tree", datafile);
  TTreeReaderValue<int>                 pTHat_bin(*reader, "pt_bin");     // pT hard bin
  TTreeReaderValue<float>               pTHat_xsec(*reader, "bin_xsec");  // pT hard xsec
  TTreeReaderValue<int>                 process(*reader, "process");      // PYTHIA truth process
  TTreeReaderValue<int>                 pcharm(*reader, "pcharm");        // count of charm quarks in the hard scattering
//  TTreeReaderValue<int>                 pbeauty(*reader, "pbeauty");      // count of charm quarks in the hard scattering
  TTreeReaderValue<int>                 mult(*reader, "ntotal");          // final state hadron multiplicity
  TTreeReaderValue<int>                 cmult(*reader, "ncharge");        // final state charged hadron multiplicity
  TTreeReaderValue<int>                 njet(*reader, "njet");            // jet multiplicity in this event
  TTreeReaderValue<int>                 nchjet(*reader, "nchjet");        // charged jet multiplicity in this event
  TTreeReaderValue<std::vector<int>>    jet_n(*reader, "jet_n");          // count of constituents in each jet
  TTreeReaderValue<std::vector<int>>    chjet_n(*reader, "chjet_n");      // count of charged constituents in each jet
  TTreeReaderValue<std::vector<float>>  jet_pt(*reader, "jet_pt");        // transverse momentum of each jet constituent
  TTreeReaderValue<std::vector<float>>  chjet_pt(*reader, "chjet_pt");    // transverse momentum of each charged jet constitiuent
//  TTreeReaderValue<std::vector<float>>  jet_y(*reader, "jet_y");          // rapidity of each jet constituent
//  TTreeReaderValue<std::vector<float>>  chjet_y(*reader, "chjet_y");      // rapidity of each charged jet constitiuent
  TTreeReaderValue<std::vector<float>>  jet_eta(*reader, "jet_eta");      // pseudorapidity of each jet constituent
  TTreeReaderValue<std::vector<float>>  chjet_eta(*reader, "chjet_eta");  // pseudorapidity of each charged jet constitiuent
  TTreeReaderValue<std::vector<float>>  jet_phi(*reader, "jet_phi");      // azimuth of each jet constituent
  TTreeReaderValue<std::vector<float>>  chjet_phi(*reader, "chjet_phi");  // azimuth of each charged jet constitiuent
  
  // Charm status (at the level of D0 tagging)
  //    0: No charm in jet history
  //    1: Contains reconstructed D0
  //    2: Contains another charmed hadron
  TTreeReaderValue<std::vector<int>>    jet_charm_status(*reader, "jet_charm_status");      // charm status of
  TTreeReaderValue<std::vector<int>>    chjet_charm_status(*reader, "chjet_charm_status");
  
  
  
  pythia_event_tree->Branch("jet_charm_status",   &jet_charm_status);
  pythia_event_tree->Branch("chjet_charm_status", &chjet_charm_status);
  
  
  
  int                     mult;         // Total event track multiplicity
  int                     cmult;        // Event charged track multiplicity
  float                   delta_E;      // Deviation from energy conservation
  int                     ntot_jet;     // Number of full jets in an event
  std::vector<int>        ntot_cst;     // Array of full jet constituent multiplicties
  std::vector<float>      jet_pT;       // Array of full jet p_T
  std::vector<float>      jet_y;        // Array of full jet rapidity
  std::vector<float>      jet_eta;      // Array of full jet pseudorapidity
  std::vector<float>      jet_phi;      // Array of full jet azimuth
  std::vector<float>      jet_area;     // Array of full jet areas
  int                     ntot_cjet;    // Number of charged jets in an event
  std::vector<int>        ntot_ccst;    // Array of charged jet constituent multiplicities
  std::vector<float>      chjet_pT;     // Array of charged jet p_T
  std::vector<float>      chjet_y;      // Array of charged jet rapidity
  std::vector<float>      chjet_eta;    // Array of charged jet pseudorapidity
  std::vector<float>      chjet_phi;    // Array of charged jet azimuth
  std::vector<float>      chjet_area;   // Array of charged jet areas
  
  
  for (int iCharge = 0; iCharge < 2; ++iCharge) {
    
    
    
    
  }// End of loop over charge state
  
  
  //====================================================================================Light Jets (u,d,s quark initiated)
  
  // Gather tree level info for light jets
  const Int_t ntotal_jet_light = datafile->Get<TTree>(Form("%sjetTree_light",cjet_string[charged_jets]))->GetBranch("jet_mult")->GetEntries();
  TTreeReader* reader_light = new TTreeReader(Form("%sjetTree_light",cjet_string[charged_jets]), datafile);
  TTreeReaderValue<Int_t> jet_mult_light(*reader_light, "jet_mult");
  TTreeReaderValue<Float_t> jet_pT_light(*reader_light, "jet_pT");
  TTreeReaderValue<Float_t> jet_eta_light(*reader_light, "jet_eta");
  TTreeReaderArray<Float_t> constituent_pT_light(*reader_light, "cont_pT");
  TTreeReaderArray<Float_t> constituent_eta_light(*reader_light, "cont_eta");
  TTreeReaderArray<Float_t> constituent_phi_light(*reader_light, "cont_phi");
  
  // Compute EEC for these particles
  while (reader_light->Next()) {
    if (TMath::Abs(*jet_eta_light) > etajet_latecut) continue;
    Int_t cjet_nconstituents = *jet_mult_light;
    for (Int_t i = 0; i < cjet_nconstituents; ++i) {
      for (Int_t j = i+1; j < cjet_nconstituents; ++j) {
        Double_t RL = TMath::Sqrt(TMath::Sq(constituent_eta_light[i] - constituent_eta_light[j]) +
                                  TMath::Sq(constituent_phi_light[i] - constituent_phi_light[j]) );
        if (*jet_pT_light > jetpt_binedge[0] && *jet_pT_light < jetpt_binedge[1]) {
          EEC_hist[0][0]->Fill(RL, constituent_pT_light[i]*constituent_pT_light[j] / TMath::Sq(*jet_pT_light));
        } else if (*jet_pT_light > jetpt_binedge[1] && *jet_pT_light < jetpt_binedge[2]) {
          EEC_hist[0][1]->Fill(RL, constituent_pT_light[i]*constituent_pT_light[j] / TMath::Sq(*jet_pT_light));
        } else if (*jet_pT_light > jetpt_binedge[2] && *jet_pT_light < jetpt_binedge[3]) {
          EEC_hist[0][2]->Fill(RL, constituent_pT_light[i]*constituent_pT_light[j] / TMath::Sq(*jet_pT_light));
        }
      }
    }
    
    // Tally jets in each pT bin
    for (Int_t ipt = 0; ipt < 3; ++ipt)
      if (*jet_pT_light > jetpt_binedge[ipt] && *jet_pT_light < jetpt_binedge[ipt+1])
        ++count_jets[0][ipt];
    
  } // End of light jet tree loop
//  
//  //====================================================================================Charm Jets (c quark initiated)
//  
//  // Gather tree level info for light jets
//  const Int_t ntotal_jet_charm = datafile->Get<TTree>(Form("%sjetTree_charm",cjet_string[charged_jets]))->GetBranch("jet_mult")->GetEntries();
//  TTreeReader* reader_charm = new TTreeReader(Form("%sjetTree_charm",cjet_string[charged_jets]), datafile);
//  TTreeReaderValue<Int_t> jet_mult_charm(*reader_charm, "jet_mult");
//  TTreeReaderValue<Float_t> jet_pT_charm(*reader_charm, "jet_pT");
//  TTreeReaderValue<Float_t> jet_eta_charm(*reader_charm, "jet_eta");
//  TTreeReaderArray<Float_t> constituent_pT_charm(*reader_charm, "cont_pT");
//  TTreeReaderArray<Float_t> constituent_eta_charm(*reader_charm, "cont_eta");
//  TTreeReaderArray<Float_t> constituent_phi_charm(*reader_charm, "cont_phi");
//  
//  // Compute EEC for these particles
//  while (reader_charm->Next()) {
//    if (TMath::Abs(*jet_eta_charm) > etajet_latecut) continue;
//    Int_t cjet_nconstituents = *jet_mult_charm;
//    for (Int_t i = 0; i < cjet_nconstituents; ++i) {
//      for (Int_t j = i+1; j < cjet_nconstituents; ++j) {
//        Double_t RL = TMath::Sqrt(TMath::Sq(constituent_eta_charm[i] - constituent_eta_charm[j]) +
//                                  TMath::Sq(constituent_phi_charm[i] - constituent_phi_charm[j]) );
//        if (*jet_pT_charm > jetpt_binedge[0] && *jet_pT_charm < jetpt_binedge[1]) {
//          EEC_hist[1][0]->Fill(RL, constituent_pT_charm[i]*constituent_pT_charm[j] / TMath::Sq(*jet_pT_charm));
//        } else if (*jet_pT_charm > jetpt_binedge[1] && *jet_pT_charm < jetpt_binedge[2]) {
//          EEC_hist[1][1]->Fill(RL, constituent_pT_charm[i]*constituent_pT_charm[j] / TMath::Sq(*jet_pT_charm));
//        } else if (*jet_pT_charm > jetpt_binedge[2] && *jet_pT_charm < jetpt_binedge[3]) {
//          EEC_hist[1][2]->Fill(RL, constituent_pT_charm[i]*constituent_pT_charm[j] / TMath::Sq(*jet_pT_charm));
//        }
//      }
//    }
//    
//    // Tally jets in each pT bin
//    for (Int_t ipt = 0; ipt < 3; ++ipt)
//      if (*jet_pT_charm > jetpt_binedge[ipt] && *jet_pT_charm < jetpt_binedge[ipt+1])
//        ++count_jets[1][ipt];
//    
//  } // End of charm jet tree loop
//  
//  //====================================================================================Beauty Jets (b quark initiated)
//  
//  // Gather tree level info for light jets
//  const Int_t ntotal_jet_beauty = datafile->Get<TTree>(Form("%sjetTree_beauty",cjet_string[charged_jets]))->GetBranch("jet_mult")->GetEntries();
//  TTreeReader* reader_beauty = new TTreeReader(Form("%sjetTree_beauty",cjet_string[charged_jets]), datafile);
//  TTreeReaderValue<Int_t> jet_mult_beauty(*reader_beauty, "jet_mult");
//  TTreeReaderValue<Float_t> jet_pT_beauty(*reader_beauty, "jet_pT");
//  TTreeReaderValue<Float_t> jet_eta_beauty(*reader_beauty, "jet_eta");
//  TTreeReaderArray<Float_t> constituent_pT_beauty(*reader_beauty, "cont_pT");
//  TTreeReaderArray<Float_t> constituent_eta_beauty(*reader_beauty, "cont_eta");
//  TTreeReaderArray<Float_t> constituent_phi_beauty(*reader_beauty, "cont_phi");
//  
//  // Compute EEC for these particles
//  while (reader_beauty->Next()) {
//    if (TMath::Abs(*jet_eta_beauty) > etajet_latecut) continue;
//    Int_t cjet_nconstituents = *jet_mult_beauty;
//    for (Int_t i = 0; i < cjet_nconstituents; ++i) {
//      for (Int_t j = i+1; j < cjet_nconstituents; ++j) {
//        Double_t RL = TMath::Sqrt(TMath::Sq(constituent_eta_beauty[i] - constituent_eta_beauty[j]) +
//                                  TMath::Sq(constituent_phi_beauty[i] - constituent_phi_beauty[j]) );
//        if (*jet_pT_beauty > jetpt_binedge[0] && *jet_pT_beauty < jetpt_binedge[1]) {
//          EEC_hist[2][0]->Fill(RL, constituent_pT_beauty[i]*constituent_pT_beauty[j] / TMath::Sq(*jet_pT_beauty));
//        } else if (*jet_pT_beauty > jetpt_binedge[1] && *jet_pT_beauty < jetpt_binedge[2]) {
//          EEC_hist[2][1]->Fill(RL, constituent_pT_beauty[i]*constituent_pT_beauty[j] / TMath::Sq(*jet_pT_beauty));
//        } else if (*jet_pT_beauty > jetpt_binedge[2] && *jet_pT_beauty < jetpt_binedge[3]) {
//          EEC_hist[2][2]->Fill(RL, constituent_pT_beauty[i]*constituent_pT_beauty[j] / TMath::Sq(*jet_pT_beauty));
//        }
//      }
//    }
//    
//    // Tally jets in each pT bin
//    for (Int_t ipt = 0; ipt < 3; ++ipt)
//      if (*jet_pT_beauty > jetpt_binedge[ipt] && *jet_pT_beauty < jetpt_binedge[ipt+1])
//        ++count_jets[2][ipt];
//    
//  } // End of beauty jet tree loop
//  
//  //====================================================================================Gluon Jets (gluon initiated)
//  
//  // Gather tree level info for light jets
//  const Int_t ntotal_jet_gluon = datafile->Get<TTree>(Form("%sjetTree_gluon",cjet_string[charged_jets]))->GetBranch("jet_mult")->GetEntries();
//  TTreeReader* reader_gluon = new TTreeReader(Form("%sjetTree_gluon",cjet_string[charged_jets]), datafile);
//  TTreeReaderValue<Int_t> jet_mult_gluon(*reader_gluon, "jet_mult");
//  TTreeReaderValue<Float_t> jet_pT_gluon(*reader_gluon, "jet_pT");
//  TTreeReaderValue<Float_t> jet_eta_gluon(*reader_gluon, "jet_eta");
//  TTreeReaderArray<Float_t> constituent_pT_gluon(*reader_gluon, "cont_pT");
//  TTreeReaderArray<Float_t> constituent_eta_gluon(*reader_gluon, "cont_eta");
//  TTreeReaderArray<Float_t> constituent_phi_gluon(*reader_gluon, "cont_phi");
//  
//  // Compute EEC for these particles
//  while (reader_gluon->Next()) {
//    if (TMath::Abs(*jet_eta_gluon) > etajet_latecut) continue;
//    Int_t cjet_nconstituents = *jet_mult_gluon;
//    for (Int_t i = 0; i < cjet_nconstituents; ++i) {
//      for (Int_t j = i+1; j < cjet_nconstituents; ++j) {
//        Double_t RL = TMath::Sqrt(TMath::Sq(constituent_eta_gluon[i] - constituent_eta_gluon[j]) +
//                                  TMath::Sq(constituent_phi_gluon[i] - constituent_phi_gluon[j]) );
//        if (*jet_pT_gluon > jetpt_binedge[0] && *jet_pT_gluon < jetpt_binedge[1]) {
//          EEC_hist[3][0]->Fill(RL, constituent_pT_gluon[i]*constituent_pT_gluon[j] / TMath::Sq(*jet_pT_gluon));
//        } else if (*jet_pT_gluon > jetpt_binedge[1] && *jet_pT_gluon < jetpt_binedge[2]) {
//          EEC_hist[3][1]->Fill(RL, constituent_pT_gluon[i]*constituent_pT_gluon[j] / TMath::Sq(*jet_pT_gluon));
//        } else if (*jet_pT_gluon > jetpt_binedge[2] && *jet_pT_gluon < jetpt_binedge[3]) {
//          EEC_hist[3][2]->Fill(RL, constituent_pT_gluon[i]*constituent_pT_gluon[j] / TMath::Sq(*jet_pT_gluon));
//        }
//      }
//    }
//    
//    // Tally jets in each pT bin
//    for (Int_t ipt = 0; ipt < 3; ++ipt)
//      if (*jet_pT_gluon > jetpt_binedge[ipt] && *jet_pT_gluon < jetpt_binedge[ipt+1])
//        ++count_jets[3][ipt];
//    
//  } // End of gluon jet tree loop
//  
//  //====================================================================================Cleanup and Plotting
//  
//  Int_t count_jets_total[3] = {0,0,0};
//  
//  for (Int_t iFlavor = 0; iFlavor < 4; ++iFlavor) {
//    for (int ipt = 0; ipt < 3; ++ipt) {
//      std::cout << Form("%c-Jets in bin %i-%i: %i", char_flavor[iFlavor], jetpt_binedge[ipt], jetpt_binedge[ipt+1], count_jets[iFlavor][ipt]) << std::endl;
//      count_jets_total[ipt] += count_jets[iFlavor][ipt];
//    }
//  }
//  for (int ipt = 0; ipt < 3; ++ipt) {
//    for (int iFlavor = 0; iFlavor < 4; ++iFlavor) {
//      if (flavor_dependent_norm) EEC_hist[iFlavor][ipt]->Scale(1./count_jets[iFlavor][ipt], "width");
//      else                       EEC_hist[iFlavor][ipt]->Scale((crosssec_process[iFlavor]/crosssec_process[0])/count_jets[iFlavor][ipt], "width");
//    }
//  }
//  
//  TLegend* legend_flavor = new TLegend(0.7, 0.7, 0.9, 0.95);
//  legend_flavor->AddEntry(EEC_hist[0][0], "Light Jet", "lp");
//  legend_flavor->AddEntry(EEC_hist[1][0], "Charm Jet", "lp");
//  legend_flavor->AddEntry(EEC_hist[2][0], "Beauty Jet", "lp");
//  legend_flavor->AddEntry(EEC_hist[3][0], "Gluon Jet", "lp");
//  legend_flavor->SetLineWidth(0);
//  
//  for (int ipt = 0; ipt < 3; ++ipt) {
//    pads[ipt][0]->cd();
//    gPad->SetLogx();
//    
//    EEC_hist[0][ipt]->GetYaxis()->SetRangeUser(0, 3.3);
//    EEC_hist[0][ipt]->Draw("hist p l");
//    EEC_hist[1][ipt]->Draw("hist p l same");
//    EEC_hist[2][ipt]->Draw("hist p l same");
//    EEC_hist[3][ipt]->Draw("hist p l same");
//    
//    drawText(Form("%i < p_{T}^{jet} < %i GeV", jetpt_binedge[ipt], jetpt_binedge[ipt+1]), 0.15, 0.9 - gPad->GetTopMargin(), false, kBlack, 0.05);
//    if (ipt == 2) legend_flavor->Draw();
//    if (ipt == 0) { //Kinematic labels
//      drawText(Form("#bf{PYTHIA} #it{pp} #sqrt{s} = %.2f TeV", sqrt_s/1000), 0.9, 0.92 - gPad->GetTopMargin(), true);
//      drawText(Form("#bf{FastJet3} anti-k_{T} %s jet R = %.1f", cjet_longstring[charged_jets], jetRadius), 0.9, 0.85 - gPad->GetTopMargin(), true);
//      drawText(Form("#left|#eta^{jet}#right| < %.1f, p_{T, min}^{hadron} = %.1f GeV", etajet_latecut, pTmin_hadron), 0.9, 0.77 - gPad->GetTopMargin(), true);
////      drawText(Form("Minimum total invariant #hat{p}_{T, min} = %.1f GeV", pTmin_jet),0.9, 0.7 - gPad->GetTopMargin(), true);
//    }
//  }
//  canvas->SaveAs("EEC.pdf");
//  
//  canvas_solo->cd();
//  gPad->SetLogx();
//  
//  
//  TLegend* legend_flavor_solo = new TLegend(0.15, 0.45, 0.4, 0.75);
//  legend_flavor_solo->AddEntry(EEC_hist[0][0], "Light Jet", "l p");
//  legend_flavor_solo->AddEntry(EEC_hist[1][0], "Charm Jet", "l p");
//  legend_flavor_solo->AddEntry(EEC_hist[2][0], "Beauty Jet", "l p");
//  legend_flavor_solo->AddEntry(EEC_hist[3][0], "Gluon Jet", "l p");
//  legend_flavor_solo->SetLineWidth(0);
//  
//  EEC_hist[0][0]->Draw("hist p c");
//  EEC_hist[1][0]->Draw("hist p c same");
//  EEC_hist[2][0]->Draw("hist p c same");
//  EEC_hist[3][0]->Draw("hist p c same");
//  
//  drawText(Form("%i < p_{T}^{jet} < %i GeV", jetpt_binedge[0], jetpt_binedge[1]), 0.15, 0.9 - gPad->GetTopMargin(), false, kBlack, 0.05);
//  legend_flavor_solo->Draw();
//  
//  drawText(Form("#bf{PYTHIA} #it{pp} #sqrt{s} = %.2f TeV", sqrt_s/1000), 0.9, 0.92 - gPad->GetTopMargin(), true);
//  drawText(Form("#bf{FastJet3} anti-k_{T} %s jet R = %.1f", cjet_longstring[charged_jets], jetRadius), 0.9, 0.85 - gPad->GetTopMargin(), true);
//  drawText(Form("#left|#eta^{jet}#right| < %.1f, p_{T, min}^{hadron} = %.1f GeV", etajet_latecut, pTmin_hadron), 0.9, 0.77 - gPad->GetTopMargin(), true);
//  canvas_solo->SaveAs("EEC_solo.pdf");
  
  return;
}
