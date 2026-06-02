/// Code for the generation of proton-proton collision events using PYTHIA
/// The data of the events produced in the Monte Carlo will be processed and
/// compressed into some essential information like spectrum histograms.
/// These histograms can be used for data analysis and MC comparisons
/// against data.
///
/// Requires cross section information from ptbin.h
/// Designed for running on the Grace high-performance computing cluster
/// events are divided into many pTHatMin (pTHard) bins that are generated
/// concurrently (multithreaded)
///
/// Important kinematic settings like sqrt{s} and detector rapidity window can
/// be modified in the local config.h. Settings for PYTHIA can be changed in the
/// eventgen\_jets.cmnd file.
///
/// Written by Ryan Hamilton with cross section data from Isaac Mooney and Laura Havener
///
/// This code utilizes the PYTHIA event generator.
/// Copyright (C) 2023 Torbjorn Sjostrand.
/// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
/// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <string>
#include <chrono>
#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TVectorT.h"
#include "TString.h"
#include "TH2D.h"
#include "TMath.h"
#include "TPave.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h" 
#include "TRandom3.h" 
#include "TStyle.h"  
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "../../utils/root_draw_tools.h"
#include "config_local.h"
#include "ptbin.h"
 
using namespace Pythia8;

const bool debug = false;

// Forward declare relevant methods. These will be defined at the end of the file
fastjet::JetAlgorithm getJetAlgorithm(const char*);

double getSquareEtaPhiDistance(double eta1, double eta2, double phi1, double phi2);

bool isKnownProcess(int process);

TMarker* drawPythiaParticleMarker(const Pythia8::Particle&, int, int, double);


//================================================================== Begin main


// This method is designed for multithreading
// by dividing pTHard bins among many cores.
// The bash call expects an argument which
// represents the thread index 0-(tot_threads - 1).
// It should be called as ./eventgen_jets [THREAD_INDEX]
//
// The max thread index depends on the number of pTHatMin bins to run.
// More info on this setting is contained in the ptbin.h file.
int main(int argc, char *argv[]) {
  auto starttime = std::chrono::high_resolution_clock::now(); // to time the run
  
  // Check that inputs are provided correctly
  if (argc == 1) {
    std::cerr << "Error in <eventgen_Dzero_jets>: Not enough arguments!" << std::endl;
    std::cout << "Requires one argument for the pT_Hard bin to run over." << std::endl;
    return 0;
  }
  
  // Thread setup, throw error if thread input is out of range
  int this_thread = atoi(argv[1]);
  const int total_threads = getNbins(sqrt_s);
  if (this_thread >= total_threads) {
    std::cerr << "Error in <eventgen_jets.cc>: requested pTHatBin " << this_thread << " out of bin range " << total_threads << std::endl;
    return 0;
  }
  
  // Gather bin settings for this thread
  pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
  const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
  const float cross_section = settings.cross_section;
  const int nevent = (int) (settings.nevent * nevent_scale_factor); // scaled event number for testing
  float collision_energy_TeV = (float) sqrt_s / 1000.;
  
  // Set longitudinal parameters based on rapidity or pseudorapidity, whichever is requested
  const double longitudinal_acceptance[2] = {max_eta, max_y};
  const char longitudinal_label[2][5] = {"#eta", "y"};
  const char longitudinal_longlabel[2][25] = {"Pseudorapidity #eta", "Rapidity y"};
  
  // Create output file
  char gen_info_string[100];
  const char process_label[5][5] = {"uds","cc","bb","gg","all"};
  snprintf(gen_info_string, 100, "%.2fTeV_r%.1fjet",collision_energy_TeV,jet_radius);
  TFile *outfile = new TFile(Form("out_data/%s/pythpp_%s_%s_thread%i.root",gen_info_string,
                                  process_label[flavor_flag],gen_info_string,this_thread),"recreate");
  
  
  // Initialize Pythia, adjust necessary settings
  Pythia pythia;
  Event& cEvent    = pythia.event;            // Actual event with particles
  Event& cProcess  = pythia.process;          // Info from just the hard scattering
  const Pythia8::Info* gInfo = &pythia.info;  // For getting info about the PYTHA process.
  // Kinematic settings
  pythia.readString(Form("PhaseSpace:pTHatMax = %.f", pt_binedge[1]));
  pythia.readString(Form("PhaseSpace:pTHatMin = %.f", pt_binedge[0]));
  pythia.readString(Form("Beams:eCM = %.f", sqrt_s));
  pythia.readString(Form("Main:numberOfEvents = %i", nevent));
  // Generator settings
  pythia.readFile("eventgen_jets.cmnd");
  pythia.readString(Form("Main:timesAllowErrors = %i", times_allow_errors));
  pythia.readString(Form("Tune:pp = %i", pythia_tune));
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %i", pythia_seed));
  // Allowed QCD Processes
  switch (flavor_flag) {
    case 0: // Light quark (u, d, s)
      pythia.readString("HardQCD:nQuarkNew = 3");
      pythia.readString("HardQCD:gg2qqbar = on");
      pythia.readString("HardQCD:qqbar2qqbarNew = on");
      break;
    case 1: // Charm
      pythia.readString("HardQCD:hardccbar = on");
      break;
    case 2: // Beauty
      pythia.readString("HardQCD:hardbbbar = on");
      break;
    case 3: // Gluon only
      pythia.readString("HardQCD:gg2gg = on");
      pythia.readString("HardQCD:qqbar2gg = on");
      break;
    default: // All processes
      pythia.readString("HardQCD:all = on"); // Good to run when testing cross sections
  }
  
  // initialize
  pythia.init();
   
  std::cout << "Beginning generation for thread " << this_thread << " of " << total_threads - 1 << std::endl;
  
  
  // TTree for storing jet information in each event
  int                     process;      // PYTHIA truth level process
  int                     pcharm;       // process charm multiplicity
  int                     pbeauty;      // process beauty multiplicity
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
  
  // Which jets contain D0 mesons?
  std::vector<int>        jet_charm_status;
  std::vector<int>        chjet_charm_status;
  // Note:
  //    0: No charm in jet history
  //    1: Contains reconstructed D0
  //    2: Contains another charmed hadron
  
  char tree_name[50];
  snprintf(tree_name, 50, "pythia_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
  TTree* pythia_event_tree = new TTree(tree_name,tree_name);
  pythia_event_tree->Branch("process",    &process);
  pythia_event_tree->Branch("pcharm",    &process);
  pythia_event_tree->Branch("ntotal",     &mult);
  pythia_event_tree->Branch("ncharge",    &cmult);
  pythia_event_tree->Branch("delta_E",    &delta_E);
  pythia_event_tree->Branch("njet",       &ntot_jet);
  pythia_event_tree->Branch("nchjet",     &ntot_cjet);
  pythia_event_tree->Branch("jet_n",      &ntot_cst);
  pythia_event_tree->Branch("chjet_n",    &ntot_ccst);
  pythia_event_tree->Branch("jet_pT",     &jet_pT);
  pythia_event_tree->Branch("chjet_pT",   &chjet_pT);
  pythia_event_tree->Branch("jet_y",      &jet_y);
  pythia_event_tree->Branch("chjet_y",    &chjet_y);
  pythia_event_tree->Branch("jet_eta",    &jet_eta);
  pythia_event_tree->Branch("chjet_eta",  &chjet_eta);
  pythia_event_tree->Branch("jet_phi",    &jet_phi);
  pythia_event_tree->Branch("chjet_phi",  &chjet_phi);
  pythia_event_tree->Branch("jet_area",   &jet_area);
  pythia_event_tree->Branch("chjet_area", &chjet_area);
  
  pythia_event_tree->Branch("jet_charm_status",   &jet_charm_status);
  pythia_event_tree->Branch("chjet_charm_status", &chjet_charm_status);
  
  // Construct histograms for storing information on the single-particle level
  const int n_hist_type = 3;
  const int n_particle_type = 4;
  char histtype[n_hist_type][15] = {"hist_pT", "hist_y_phi", "hist_eta_phi"};
  char particletype[n_particle_type][10] = {"hadron","chadron","const","chconst"};
  char plottitle[2*n_hist_type][100] = {
    ";#it{p}_{T} [GeV];1/N^{event} dN^{hadron}/d#it{p}_{T}",
    ";y;#phi;1/N^{event} d^{2}N^{hadron}/dyd#phi",
    ";#eta;#phi;1/N^{event} d^{2}N^{hadron}/d#etad#phi",
    ";#it{p}_{T} [GeV];1/N^{event} dN^{constituent}/d#it{p}_{T}",
    ";y;#phi;1/N^{event} d^{2}N^{constituent}/dyd#phi",
    ";#eta;#phi;1/N^{event} d^{2}N^{constituent}/d#etad#phi"
  };
  
  TH1D* hist_pT[n_particle_type];
  TH2D* hist_2D[2][n_particle_type];
  for (int i_particle_type = 0; i_particle_type < n_particle_type; ++i_particle_type) {
    hist_pT[i_particle_type]    = new TH1D(Form("%s_%s_ptbin%.f-%.f",histtype[0],particletype[i_particle_type],pt_binedge[0],pt_binedge[1]),
                                           plottitle[3*(i_particle_type/2)], nbins_pT, 0, pTmax_hadron_eventdisplay);
    hist_2D[0][i_particle_type] = new TH2D(Form("%s_%s_ptbin%.f-%.f",histtype[1],particletype[i_particle_type],pt_binedge[0],pt_binedge[1]),
                                           plottitle[3*(i_particle_type/2)+1],nbins_rap,-max_y,max_y,nbins_phi,-TMath::Pi(),TMath::Pi());
    hist_2D[1][i_particle_type] = new TH2D(Form("%s_%s_ptbin%.f-%.f",histtype[2],particletype[i_particle_type],pt_binedge[0],pt_binedge[1]),
                                           plottitle[3*(i_particle_type/2)+2],nbins_rap,-max_eta,max_eta,nbins_phi,-TMath::Pi(),TMath::Pi());
  }
  
  // Setup jet definiton for fastjet
  // Using algorithm and jet radius according to config file
  fastjet::JetAlgorithm  jet_algo         = getJetAlgorithm(algo_string_short);
  fastjet::JetDefinition jet_definition   = fastjet::JetDefinition(jet_algo, jet_radius);
  fastjet::JetDefinition chjet_definition = fastjet::JetDefinition(jet_algo, jet_radius);
  
  // 2D Histograms that will display solid regions for reconstructed jets
  // Purely for illustrative purposes to show that the jet finder is working properly
  const double dA = 4*TMath::Pi()*longitudinal_acceptance[cylinder_uses_rapidity]/((double)(nbins_rap_display*nbins_phi_display));
  TH2D* jet_diagram2D[2];
  char jettype[2][10] = {"jet","chjet"};
  for (int i_charged = 0; i_charged < 2; ++i_charged) {
    jet_diagram2D[i_charged] = new TH2D(Form("%s_display",jettype[i_charged]),
                                        Form(";%s;#phi [rad]; #it{p}_{T}^{jet} [GeV]",longitudinal_longlabel[cylinder_uses_rapidity]),
                                        nbins_rap_display,-max_eta,max_eta,nbins_phi_display,-TMath::Pi(),TMath::Pi());
    jet_diagram2D[i_charged]->GetYaxis()->SetTitleOffset(0.8);
    jet_diagram2D[i_charged]->GetZaxis()->SetTitleOffset(1.2);
  }
   
  // Setup for root plotting
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainBow);
  TCanvas* canvas = new TCanvas();
  canvas->SetCanvasSize(1750, 600);
  canvas->SetMargin(0, 0, 0, 0);
  TPad* pads[2];
  pads[0] = buildPad("pad_1", 0, 0, 0.5, 1, 0.08, 0.15, 0.1, 0.075);
  pads[1] = buildPad("pad_2", 0.5, 0, 1, 1, 0.06, 0.15, 0.1, 0.075);
  pads[0]->SetLogz();
  jet_diagram2D[0]->GetZaxis()->SetRangeUser(jetcut_minpT_full, 1000);
  pads[1]->SetLogz();
  jet_diagram2D[1]->GetZaxis()->SetRangeUser(jetcut_minpT_chrg, 1000);
  
  // Open pdf file for event display
  TString plotfile = Form("out_plots/%s/%s_eventdisplay_ptbin%.f-%.f.pdf",
                          gen_info_string,gen_info_string,pt_binedge[0],pt_binedge[1]);
  canvas->Print(plotfile + "[");
  
  
  //========================================================================== Start of event loop
  
  int total_reco_D0[2] = {0,0};
  int total_missed_charm[2] = {0,0};
  
  bool flag_draw_event;
  int ntotal_plotted[4] = {0,0,0,0};
  int ntotal_skew_plotted[2] = {0,0};
  int iError = 0; 
  int grandtotal_jet[2] = {0,0};
  int count_charm_in_gq[3] = {0, 0, 0}; // for process ids 113, 114, 115
  int count_beauty_in_gq[3] = {0, 0, 0}; // for process ids 113, 114, 115
  std::map<int, int> process_count;
  for (int iEvent = 0; iEvent < nevent; ++iEvent) {
    
    // Test for errors, throw if necessary
    if (!pythia.next()) {
      if (++iError < times_allow_errors) continue;
      // else
      cerr << "Event generation aborted due to error.\n";
      break;
    }
    
    if (iEvent == 0) {
      // Print initial cross sections
      std::cout << "Initial cross sections:" << std::endl;
      for (int i = 111; i <= 124; ++i) {
        if (i > 116 && i < 121) continue;
        std::cout << i << " :: " << std::scientific << gInfo->sigmaGen(i) << " +/- " << gInfo->sigmaErr(i) << std::endl;
      }
    }
    
    // PYTHIA event-level process info
    process = gInfo->code();
    ++process_count[process];
    if (!isKnownProcess(process)) {
      std::cout << "Unrecognized Process: ID " << process << " corresponds to " << gInfo->name() << std::endl;
      std::cout << "Process listing:" << std::endl;
      for (Particle p : cProcess) std::cout << p.id() << " with status " << p.status();
    } if (gInfo->weight() != 1) {
      std::cout << "Warning in PYTHIA Eventgen :: nonzero weight = " << gInfo->weight() << " found!" << std::endl;
    }
    
    // Count charm/beauty quarks in the process
    pcharm = 0; pbeauty = 0;
    for (Particle p : cProcess) {
      if (std::abs(p.id()) == 4) ++pcharm;
      if (std::abs(p.id()) == 5) ++pbeauty;
    }
    if (pcharm > 0) {
      std::cout << pcharm << " charm quarks in process " << process << "!" << std::endl;
      if (process >= 113 && process <= 115) count_charm_in_gq[process-113] += pcharm/2; // due to flavor conservation
    }
    if (pbeauty > 0) {
      std::cout << pbeauty << " beauty quarks in process " << process << "!" << std::endl;
      if (process >= 113 && process <= 115) count_beauty_in_gq[process-113] += pbeauty/2; // due to flavor conservation
    }
    

    //========================================================================== Gather PYTHIA particles into FastJet vectors
    
    // Extract particles from current pythia event
    std::vector<Particle> particles;
    std::vector<Particle> reco_D0;
    std::vector<Particle> missed_charm;
    std::vector<fastjet::PseudoJet> stable_particles, stable_charged_particles;
    Vec4 pSum;
    for (int i = 0; i < cEvent.size(); ++i) { // Get info from final state particles
      Particle &p = cEvent[i];
      
      // PDG MC ID for D0 is 421.
      if (std::abs(p.id()) == 421 && p.pT() >= D0_minpT
          && std::abs(p.eta()) < max_eta) {
        // Add the D0 regardless of status (decay daughter or other)
        // From what I can tell investigating this in PYTHIA, there isn't
        // a status setting for a particle directly connected to an original
        // hard scattering... I could be missing something though...
        // Still though, maybe this is better, to be more MC-agnostic?
        
        reco_D0.push_back(p);
      } else if ((std::abs(p.id()) == 421 && p.pT() < D0_minpT) || // low pT D0, missed.
                 (std::abs(p.id()) == 421 && std::abs(p.eta()) > max_eta) || // D0 outside detector acceptance
                 (p.id() >= 400 && p.id() < 500) || // other charmed meson
                 (p.id() >= 4000 && p.id() < 5000)) { // a charmed baryon
        missed_charm.push_back(p);
      }
      
      if (!cEvent[i].isFinal()) continue; // Only add final state hadrons to clustering
      
      particles.push_back(p);
      pSum += Vec4(p.px(), p.py(), p.pz(), p.e());
      if (cEvent[i].pT() < pTmin_hadron) continue; // ignore hadrons with pT is lower than detector efficiency threshold
      
      stable_particles.push_back(fastjet::PseudoJet(p.px(),p.py(),p.pz(),p.e()));
      if (p.isCharged()) stable_charged_particles.push_back(fastjet::PseudoJet(p.px(),p.py(),p.pz(),p.e()));
      
      // Fill histograms at the track level
      hist_pT[0]->Fill(p.pT());
      hist_2D[0][0]->Fill(p.y(), p.phi());
      hist_2D[1][0]->Fill(p.eta(),p.phi());
      if (p.isCharged()) {
        hist_pT[1]->Fill(p.pT());
        hist_2D[0][1]->Fill(p.y(), p.phi());
        hist_2D[1][1]->Fill(p.eta(),p.phi());
      }
    }
    
    mult = stable_particles.size();
    cmult = stable_charged_particles.size();
    pSum /= cEvent[0].e(); // total event energy
    delta_E = std::abs(pSum.e() - 1);  // normalized energy loss fraction DeltaE / E
    
    // Make a grid of ghost particles, add them to the stack
    fastjet::PseudoJet ghost;
    double ghost_pt = 1e-100;
    for (int ir = 1; ir <= nbins_rap_display; ++ir) {
      // For massless particles like ghosts, rapidity = pseudorapidity
      double longitude = jet_diagram2D[0]->GetXaxis()->GetBinCenter(ir);
      for (int iphi = 1; iphi <= nbins_phi_display; ++iphi) {
        double phi = jet_diagram2D[0]->GetYaxis()->GetBinCenter(iphi);
        ghost.reset_momentum_PtYPhiM(ghost_pt, longitude, phi, 0);
        stable_particles.push_back(ghost);
        stable_charged_particles.push_back(ghost);
      }
    }
    
    //========================================================================== Full Jet Clustering
    
    // debug
    if (debug) {
      std::cout << "Event " << iEvent << ", Missed charm = " << missed_charm.size() << std::endl;
      std::cout << "Event " << iEvent << ", Reco D0 = " << reco_D0.size() << std::endl;
    }

    // Run jet clustering with fastjet -- All Jets
    // Note that the min pT cut is handled by fastjet::inclusive_jets)
    fastjet::ClusterSequence clustering(stable_particles, jet_definition);
    std::vector<fastjet::PseudoJet> final_jets = sorted_by_pt(clustering.inclusive_jets(jetcut_minpT_full));
    
    ntot_jet = 0;
    for (fastjet::PseudoJet jet : final_jets) { // Full jets that pass the pT cut
      
      // Cut this jet if its radius extends outside the detector acceptance
      double jet_long[2] = {jet.eta(), jet.rap()};
      if (abs(jet_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]-jet_radius) continue;
      
      // Cut this jet if its core (highest energy particle) is too soft
      std::vector<fastjet::PseudoJet> constituents_ptsort = sorted_by_pt(jet.constituents());
      if ((constituents_ptsort.at(0)).pt() < pTmin_jetcore) {
        cout << "Full jet cut with core pT :: " << (constituents_ptsort.at(0)).pt() << endl;
        continue;
      }
      
      // Process jet data
      int ctot_constituents = 0;
      double ctot_area = 0;
      for (fastjet::PseudoJet c : jet.constituents()) {
        if (c.pt() < 1e-50) ctot_area += dA; // grid based jet area calculation
        else { // Fill constituent jet histograms with non-ghost particles
          hist_pT[2]->Fill(c.pt());
          hist_2D[0][2]->Fill(c.rap(), c.phi_std());
          hist_2D[1][2]->Fill(c.eta(), c.phi_std());
          ++ctot_constituents;
        }
      }
       
      // Cut jet if it fails to meet the area cut
      if (ctot_area < min_area) {
        cout << "Full jet cut low area :: " << ctot_area/(3.14159*jet_radius*jet_radius) << "*pi*R^2" << endl;
        continue;
      }
      
      // Check for D0/charm within the jet
      int charm_status = 0;
      for (int mc = 0; mc < missed_charm.size(); ++mc) {
        if (jet_radius*jet_radius > getSquareEtaPhiDistance(missed_charm[mc].eta(), jet.eta(),
                                                            missed_charm[mc].phi(), jet.phi())) {
          if (debug) std::cout << "Missed Charm (" << missed_charm[mc].eta() << "," << missed_charm[mc].phi() << ") grouped with jet (" << jet.eta() << "," << jet.phi() << ")." << std::endl;
          charm_status = 2; break;
        }
      } for (int d0 = 0; d0 < reco_D0.size(); ++d0) {
        if (jet_radius*jet_radius > getSquareEtaPhiDistance(reco_D0[d0].eta(), jet.eta(),
                                                            reco_D0[d0].phi(), jet.phi())) {
          if (debug) std::cout << "Reco D0 (" << reco_D0[d0].eta() << "," << reco_D0[d0].phi() << ") grouped with jet (" << jet.eta() << "," << jet.phi() << ")." << std::endl;
          charm_status = 1; break;
        }
      }// not separated since a jet may have a D0 and another charmed hadron
      
      
      // debug
      switch (charm_status) {
        case 1:
          ++total_reco_D0[0];
          if (debug) std::cout << "Full r D0 in event " << iEvent << std::endl;
          break;
        case 2:
          ++total_missed_charm[0];
          if (debug) std::cout << "Full Xc in event " << iEvent << std::endl;
        default:
	  break;
      }
      
      // Update jet counters for pythia_event_tree
      ntot_cst.push_back(ctot_constituents);
      jet_pT.push_back(jet.pt());
      jet_y.push_back(jet.rap());
      jet_eta.push_back(jet.eta());
      jet_phi.push_back(jet.phi_std());
      jet_area.push_back(ctot_area);
      jet_charm_status.push_back(charm_status);
      ++ntot_jet;
      
    }// end of jet loop
    
    //========================================================================== Charged Jet Clustering
    
    // Run jet clustering with fastjet -- Charged Jets only
    fastjet::ClusterSequence clustering_charged(stable_charged_particles, chjet_definition);
    std::vector<fastjet::PseudoJet> final_chjets = sorted_by_pt(clustering_charged.inclusive_jets(jetcut_minpT_chrg));
    
    ntot_cjet = 0;
    for (fastjet::PseudoJet jet:final_chjets) { // Charged jets that pass the pT cut
      
      // Cut this jet if its radius extends outside the detector acceptance
      double jet_long[2] = {jet.eta(), jet.rap()};
      if (abs(jet_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]-jet_radius) continue;
      
      // Cut this jet if its core (highest energy particle) is too soft
      std::vector<fastjet::PseudoJet> constituents_ptsort = sorted_by_pt(jet.constituents());
      if ((constituents_ptsort.at(0)).pt() < pTmin_jetcore) {
        cout << "Ch jet cut with core pT :: " << (constituents_ptsort.at(0)).pt() << endl;
        continue;
      }
      
      // Process jet data
      int ctot_constituents = 0;
      double ctot_area = 0;
      for (fastjet::PseudoJet c : jet.constituents()) {
        if (c.pt() < 1e-50) ctot_area += dA;
        else { // Fill constituent jet histograms with non-ghost particles
          hist_pT[3]->Fill(c.pt());
          hist_2D[0][3]->Fill(c.rap(), c.phi_std());
          hist_2D[1][3]->Fill(c.eta(), c.phi_std());
          ++ctot_constituents;
        }
      }
      
      // Cut jet if it fails to meet the area cut
      if (ctot_area < min_area) {
        cout << "Charged jet cut low area :: " << ctot_area/(3.14159*jet_radius*jet_radius) << "*pi*R^2" << endl;
        continue;
      }
      
      // Check for D0/charm if the jet pT is above the threshold
      int charm_status = 0;
      for (int mc = 0; mc < missed_charm.size(); ++mc) {
        if (jet_radius*jet_radius > getSquareEtaPhiDistance(missed_charm[mc].eta(), jet.eta(),
                                                            missed_charm[mc].phi(), jet.phi())) {
          if (debug) std::cout << "Missed Charm (" << missed_charm[mc].eta() << "," << missed_charm[mc].phi() << ") grouped with cjet (" << jet.eta() << "," << jet.phi() << ")." << std::endl;
          charm_status = 2; break;
        }
      } for (int d0 = 0; d0 < reco_D0.size(); ++d0) {
        if (jet_radius*jet_radius > getSquareEtaPhiDistance(reco_D0[d0].eta(), jet.eta(),
                                                            reco_D0[d0].phi(), jet.phi())) {
          if (debug) std::cout << "Reco D0 (" << reco_D0[d0].eta() << "," << reco_D0[d0].phi() << ") grouped with cjet (" << jet.eta() << "," << jet.phi() << ")." << std::endl;
          charm_status = 1; break;
        }
      }
      
      
      // debug
      switch (charm_status) {
        case 1:
          ++total_reco_D0[1];
          if (debug) std::cout << "Charged r D0 in event " << iEvent << std::endl;
          break;
        case 2:
          ++total_missed_charm[1];
          if (debug) std::cout << "Charged Xc in event " << iEvent << std::endl;
          break;
        default:
          break;
      }
      
      // Update jet counters for pythia_event_tree and gendata
      ntot_ccst.push_back(ctot_constituents);
      chjet_pT.push_back(jet.pt());
      chjet_y.push_back(jet.rap());
      chjet_eta.push_back(jet.eta());
      chjet_phi.push_back(jet.phi_std());
      chjet_area.push_back(ctot_area);
      chjet_charm_status.push_back(charm_status);
      ++ntot_cjet;
    }// end of charged jet loop
    
    //========================================================================== Event Display
    
    // Draw at least 10 events regardless of their jet content
    flag_draw_event = iEvent % (nevent/10) == 0;
    
    // Plot some events with various jet/cjet content, no more than 3 per type
    for (int i_jet = 0; i_jet < 4; ++i_jet) if (ntot_jet == i_jet + 2 && ntotal_plotted[i_jet] < 3) {
      ++ntotal_plotted[i_jet];
      flag_draw_event = true;
    } if (ntot_jet == 1 && ntot_cjet == 0 && ntotal_skew_plotted[0] < 3) {
      ++ntotal_skew_plotted[0];
      flag_draw_event = true;
    } else if (ntot_jet == 0 && ntot_cjet == 1 && ntotal_skew_plotted[1] < 3) {
      ++ntotal_skew_plotted[1];
      flag_draw_event = true;
    }
    
    // plot eta-phi event display of the current event, if flagged to do so
    if (debug || flag_draw_event) {
      jet_diagram2D[0]->Reset();
      // loop back over jets and fill the display histograms
      for (fastjet::PseudoJet jet : final_jets) {
        // Do not draw jets that fail cuts. May change this in the future (see todo at top of file)
        
        double jet_long[2] = {jet.eta(), jet.rap()};
        if (abs(jet_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]-jet_radius) continue;
        
        std::vector<fastjet::PseudoJet> constituents_ptsort = sorted_by_pt(jet.constituents());
        if ((constituents_ptsort.at(0)).pt() < pTmin_jetcore) continue;
        
        if (jet.pt() < jetcut_minpT_full) continue;
        
        for (fastjet::PseudoJet c : jet.constituents()) if (c.pt() < 1e-50) {
          double constituent_long[2] = {c.eta(), c.rap()};
          jet_diagram2D[0]->Fill(constituent_long[cylinder_uses_rapidity], c.phi_std(), jet.pt());
        }
      }
      
      pads[0]->cd();
      jet_diagram2D[0]->Draw("colz");
      drawText(Form("#it{PYTHIA} Event %i, #sqrt{s_{NN}} = %.2f TeV", iEvent, collision_energy_TeV), 0.08, 0.95);
      drawText(Form("%s R = %.1f, #it{p}_{T}^{Hard} #in [%.0f,%0.f]",algo_string,jet_radius,pt_binedge[0],pt_binedge[1]), 0.9, 0.95, true);
      
      
      // plot eta-phi jet diagram of the current event -- Charged Jets Only
      jet_diagram2D[1]->Reset();
      // loop back over jets and fill the display histograms
      for (fastjet::PseudoJet jet : final_chjets) {
        // Do not draw jets that fail cuts. May change this in the future (see todo at top of file)
        
        double jet_long[2] = {jet.eta(), jet.rap()};
        if (abs(jet_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]-jet_radius) continue;
        
        std::vector<fastjet::PseudoJet> constituents_ptsort = sorted_by_pt(jet.constituents());
        if ((constituents_ptsort.at(0)).pt() < pTmin_jetcore) continue;
        
        if (jet.pt() < jetcut_minpT_chrg) continue;
        
        for (fastjet::PseudoJet c : jet.constituents()) if (c.pt() < 1e-50) {
          double constituent_long[2] = {c.eta(), c.rap()};
          jet_diagram2D[1]->Fill(constituent_long[cylinder_uses_rapidity], c.phi_std(), jet.pt());
        }
      }
      
      pads[1]->cd();
      jet_diagram2D[1]->Draw("colz");
      drawText("#it{FastJet} ver. 3.4.1", 0.06, 0.95);
      drawText(Form("charged jet %s R = %.1f, #it{p}_{T}^{Hard} #in [%.0f,%0.f]",algo_string,jet_radius,pt_binedge[0],pt_binedge[1]), 0.9, 0.95, true);
      
      // Draw markers for the event hadrons
      for (int i_diagram = 0; i_diagram < 2; ++i_diagram) {
        pads[i_diagram]->cd();
        for (Particle &p:particles) {
          if (i_diagram == 1 && p.charge() == 0) continue;
          double particle_long[2] = {p.eta(), p.y()};
          if (abs(particle_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]) continue;
          
          if (!p.isFinal()) drawPythiaParticleMarker(p, kBlack, -1, 1); // Should never happen, this is a way of checking if nonfinal particles are leaking
          else if (p.pT() < pTmin_hadron) drawPythiaParticleMarker(p, kGray, -1, 1);
          else drawPythiaParticleMarker(p, -1, -1, 1);
        }
        
        for (Particle &p : reco_D0) drawPythiaParticleMarker(p, kBlue-10, 42, 2);
        for (Particle &p : missed_charm) drawPythiaParticleMarker(p, kRed-7, 5, 2);
      }
      
      // Save the canvas to the display file
      canvas->Print(plotfile);
    }
    
    //========================================================================== Write/Store Data
    
    // Fill tree and clear jet vectors
    pythia_event_tree->Fill();
    ntot_cst.clear();
    jet_pT.clear();
    jet_y.clear();
    jet_eta.clear();
    jet_phi.clear();
    jet_area.clear(); 
    jet_charm_status.clear();
    ntot_ccst.clear();
    chjet_pT.clear();
    chjet_y.clear();
    chjet_eta.clear();
    chjet_phi.clear();
    chjet_area.clear();
    chjet_charm_status.clear();
    grandtotal_jet[0] += ntot_jet;
    grandtotal_jet[1] += ntot_cjet;
  }//========================================================================== End of event loop
  canvas->Print(plotfile + "]"); // close event display file
  
  // Display tree info and write to file
  pythia_event_tree->Print();
  pythia_event_tree->Write(pythia_event_tree->GetName(), TObject::kOverwrite);
  
  // Print post-generation cross section information
  std::cout << "Final cross sections (process level):" << std::endl;
  for (std::map<int,int>::iterator it = process_count.begin(); it != process_count.end(); ++it) {
    std::cout << it->first << " :: " << std::scientific << gInfo->sigmaGen(it->first) << " +/- " << gInfo->sigmaErr(it->first) << std::endl;
  }
  
  // Construct and write a tree with cross section information
  double xsec_total_given = static_cast<double>(cross_section);
  double xsec_total_final = 0;
  double xsec_light_final = 0;
  double xsec_gluon_final = 0;
  double xsec_charm_final = 0;
  double xsec_beauty_final = 0;
  double xsec_total_error = 0;
  double xsec_light_error = 0;
  double xsec_gluon_error = 0;
  double xsec_charm_error = 0;
  double xsec_beauty_error = 0;
  for (std::map<int,int>::iterator it = process_count.begin(); it != process_count.end(); ++it) {
    xsec_total_final += gInfo->sigmaGen(it->first);
    // Quadrature sum of errors
    xsec_total_error += std::pow(gInfo->sigmaErr(it->first), 2);
    
    switch (it->first) {
      case 111: // gg->gg
        xsec_gluon_final += gInfo->sigmaGen(it->first);
        xsec_gluon_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
      case 112: // gg->llbar
        xsec_light_final += gInfo->sigmaGen(it->first);
        xsec_light_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
      case 113: // gq->gq
        // half are gluons, by definition
        xsec_gluon_final += 0.5*gInfo->sigmaGen(it->first);
        xsec_gluon_error += std::pow( 0.5*gInfo->sigmaErr(it->first), 2);
        
        // charm fraction
        xsec_charm_final += ( ((float) count_charm_in_gq[0])/(2.*it->second)) * gInfo->sigmaGen(it->first);
        xsec_charm_error += std::pow( ((float) count_charm_in_gq[0])/(2.*it->second) * gInfo->sigmaErr(it->first), 2);
        
        // beauty fraction
        xsec_beauty_final += ( ((float) count_beauty_in_gq[0])/(2.*it->second)) * gInfo->sigmaGen(it->first);
        xsec_beauty_error += std::pow( ((float) count_beauty_in_gq[0])/(2.*it->second) * gInfo->sigmaErr(it->first), 2);
        
        // light quarks are not charm or beauty
        xsec_light_final += ( ((float) (it->second
                                        - count_charm_in_gq[0]
                                        - count_beauty_in_gq[0]))
                             / (2.*it->second) ) * gInfo->sigmaGen(it->first);
        xsec_light_error += std::pow( ( ((float) (it->second
                                                  - count_charm_in_gq[0]
                                                  - count_beauty_in_gq[0]))
                                      / (2.*it->second) ) * gInfo->sigmaErr(it->first), 2);
        // doesn't have factor of two on it->second in numerator since other half is g.
        break;
      case 114: // qq'->qq' includes all quark/antiquark interactions
        // charm fraction
        xsec_charm_final += ( ((float) count_charm_in_gq[1])/(2.*it->second)) * gInfo->sigmaGen(it->first);
        xsec_charm_error += std::pow( ((float) count_charm_in_gq[1])/(2.*it->second) * gInfo->sigmaErr(it->first), 2);
        
        // beauty fraction
        xsec_beauty_final += ( ((float) count_beauty_in_gq[1])/(2.*it->second)) * gInfo->sigmaGen(it->first);
        xsec_beauty_error += std::pow( ((float) count_beauty_in_gq[1])/(2.*it->second) * gInfo->sigmaErr(it->first), 2);
        
        // light quarks are not charm or beauty
        xsec_light_final += ( ((float) (2.*it->second
                                        - count_charm_in_gq[1]
                                        - count_beauty_in_gq[1]))
                             / (2.*it->second) ) * gInfo->sigmaGen(it->first);
        xsec_light_error += std::pow( ( ((float) (2.*it->second
                                                  - count_charm_in_gq[1]
                                                  - count_beauty_in_gq[1]))
                                       / (2.*it->second) ) * gInfo->sigmaErr(it->first), 2);
        break;
      case 115: // qqbar->gg
        // result is always gluons
        xsec_gluon_final += gInfo->sigmaGen(it->first);
        xsec_gluon_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
      case 116: // llbar->llbar
        xsec_light_final += gInfo->sigmaGen(it->first);
        xsec_light_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
      case 121: case 122: // charm pair production
        xsec_charm_final += gInfo->sigmaGen(it->first);
        xsec_charm_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
      case 123: case 124: // beauty pair production
        xsec_beauty_final += gInfo->sigmaGen(it->first);
        xsec_beauty_error += std::pow(gInfo->sigmaErr(it->first), 2);
        break;
    }
  }
  
  // Complete the quadrature
  xsec_total_error = std::sqrt(xsec_total_error);
  xsec_light_error = std::sqrt(xsec_light_error);
  xsec_gluon_error = std::sqrt(xsec_gluon_error);
  xsec_charm_error = std::sqrt(xsec_charm_error);
  xsec_beauty_error = std::sqrt(xsec_beauty_error);
  
  int total_charm = 2 * (process_count[121] + process_count[122]) + count_charm_in_gq[0] + count_charm_in_gq[1];
  double frac_charm_fromgq;
  if (process_count[113] != 0)
    frac_charm_fromgq = ((float) count_charm_in_gq[0])/total_charm;
  else frac_charm_fromgq = 0;
  double frac_charm_fromqQ;
  if (process_count[114] != 0)
    frac_charm_fromqQ = ((float) count_charm_in_gq[1])/total_charm;
  else frac_charm_fromqQ = 0;
  double frac_charm_fromcc = 0;
  if (process_count[121] != 0)
    frac_charm_fromcc += ((float) 2*process_count[121])/total_charm;
  if (process_count[122] != 0)
    frac_charm_fromcc += ((float) 2*process_count[122])/total_charm;
  
  int total_beauty = 2 * (process_count[123] + process_count[124]) + count_beauty_in_gq[0] + count_beauty_in_gq[1];
  double frac_beauty_fromgq;
  if (process_count[113] != 0)
    frac_beauty_fromgq = ((float) count_beauty_in_gq[0])/total_beauty;
  else frac_beauty_fromgq = 0;
  double frac_beauty_fromqQ;
  if (process_count[114] != 0)
    frac_beauty_fromqQ = ((float) count_beauty_in_gq[1])/total_beauty;
  else frac_beauty_fromqQ = 0;
  double frac_beauty_frombb = 0;
  if (process_count[123] != 0)
    frac_beauty_frombb += ((float) 2*process_count[123])/total_beauty;
  if (process_count[124] != 0)
    frac_beauty_frombb += ((float) 2*process_count[124])/total_beauty;
  
  
  char xsec_tree_name[50];
  snprintf(xsec_tree_name, 50, "xsec_tree_pTHard%.f-%.f", pt_binedge[0], pt_binedge[1]);
  TTree* pythia_xsec_tree = new TTree(xsec_tree_name,xsec_tree_name);
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
  
  pythia_xsec_tree->Fill();
  
  std::cout << "\nPythia parton-level cross sections in p_T^hat bin [" << pt_binedge[0] << "," << pt_binedge[1] << "]:" << std::endl;
  pythia_xsec_tree->Show(0);
  std::cout << "against charm fractions 113 :: " << count_charm_in_gq[0] << ";\t114 :: " << count_charm_in_gq[1] << std::endl;
  std::cout << "...and beauty fractions 113 :: " << count_beauty_in_gq[0] << ";\t114 :: " << count_beauty_in_gq[1] << std::endl;
  pythia_xsec_tree->Write();
  
  
  // Write histograms to file
  for (int i_particle_type = 0; i_particle_type < n_particle_type; ++i_particle_type) {
    hist_pT[i_particle_type]->Write(hist_pT[i_particle_type]->GetName(), TObject::kOverwrite);
    hist_2D[0][i_particle_type]->Write(hist_2D[0][i_particle_type]->GetName(), TObject::kOverwrite);
    hist_2D[1][i_particle_type]->Write(hist_2D[1][i_particle_type]->GetName(), TObject::kOverwrite);
  }
  
  delete outfile;
  
  
  
  // Print some information about the PYTHIA processes generated.
  std::cout << "\nPYTHIA Hard Process generation:" << std::endl;
  std::cout << "ID\t# Events" << std::endl;
  std::cout << "------------------------" << std::endl;
  for (std::map<int,int>::iterator it = process_count.begin(); it != process_count.end(); ++it)
    std::cout << it->first << "\t" << it->second << std::endl;
  std::cout << "========================" << std::endl;
  std::cout << "all\t" << nevent << std::endl;
  std::cout << "------------------------" << std::endl << std::endl;
  
  // Print informaton about D0/tagging containment inside the jets.
  std::cout << "Tot\tD0\tXc" << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << grandtotal_jet[0] << "\t" << total_reco_D0[0] << "\t" << total_missed_charm[0] << std::endl;
  std::cout << grandtotal_jet[1] << "\t" << total_reco_D0[1] << "\t" << total_missed_charm[1] << std::endl;
  std::cout << "------------------------" << std::endl << std::endl;
  
  
  // Print some valuable information about the jet generation, with an extrapolation if a scale factor is used.
  std::cout << "Total jets in this generation    :: " << grandtotal_jet[0] << ",\t" << grandtotal_jet[1] << std::endl;
  std::cout << "Extrapolated jets in full sample :: " << static_cast<int>(grandtotal_jet[0]*1./nevent_scale_factor) << ",\t" << static_cast<int>(grandtotal_jet[1]*1./nevent_scale_factor) << std::endl;
  
  // Print information about the runtime, with an extrapolation if a scale factor is used.
  auto endtime = std::chrono::high_resolution_clock::now();
  int runtime = std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count();
  char timeprint[100];
  snprintf(timeprint, 100, "%i Events generated in runtime\t:: %02d:%02d:%02d.%03d", nevent, runtime/3600000, runtime/60000 - 60*(runtime/3600000), runtime/1000 - 60*(runtime/60000), runtime%1000);
  std::cout << timeprint << std::endl;
  int extrap_runtime = 200*runtime;
  char extrap_timeprint[100];
  snprintf(extrap_timeprint, 100, "Extrapolated Runtime\t\t\t:: %02d:%02d:%02d", extrap_runtime/3600000, extrap_runtime/60000 - 60*(extrap_runtime/3600000), extrap_runtime/1000 - 60*(extrap_runtime/60000));
  std::cout << extrap_timeprint << std::endl;
  return 0;
}

//========================================================================== Helper methods

// Get fastjet algorithm based on the short string provided in config_local
fastjet::JetAlgorithm getJetAlgorithm(const char* algo_string_short) {
  std::string algo_str = static_cast<std::string>(algo_string_short);
  
  if (algo_str.compare("a-kT") == 0) return fastjet::antikt_algorithm;
  else if (algo_str.compare("kT") == 0) return fastjet::kt_algorithm;
  else if (algo_str.compare("cam") == 0) return fastjet::cambridge_algorithm;
  else if (algo_str.compare("g-kTn") == 0) return fastjet::genkt_algorithm;
  
  std::cerr << "Error in <getJetAlgorithm>: could not identify specified jet algorithm. Check config_local.h" << std::endl;
  return fastjet::antikt_algorithm;
}


// Find the distance between two objects in the eta-phi plane
// Some care must be taken to handle certain azimuthal angle cases
// to make sure the distance on the cylinder is accurately measured.
double getSquareEtaPhiDistance(double eta1, double eta2,
                               double phi1, double phi2) {
  // Project [-pi, pi) onto [0, 2pi).
  phi1 = std::fmod(phi1+TMath::TwoPi(),TMath::TwoPi());
  phi2 = std::fmod(phi2+TMath::TwoPi(),TMath::TwoPi());
  
  // Cylinder wrapping cases: distance should be the smallest wrapping
  double angular_separation = std::abs(phi1 - phi2);
  if (std::abs(phi1 - TMath::TwoPi() - phi2) < angular_separation)
      angular_separation = std::abs(phi1 - TMath::TwoPi() - phi2);
  if (std::abs(phi2 - TMath::TwoPi() - phi1) < angular_separation)
    angular_separation = std::abs(phi2 - TMath::TwoPi() - phi1);
  
  return angular_separation*angular_separation + TMath::Sq(eta1 - eta2);
}

bool isKnownProcess(int process) {
  switch (process) {
    case 111:
    case 112:
    case 113:
    case 114:
    case 115:
    case 116:
    case 121:
    case 122:
    case 123:
    case 124:
      return true;
    default:
      return false;
  }
}

// Draws a root TMarker object at the position of an input PYTHIA particle in the y-phi plane
// Default inputs for color and style are -1, will use particle information to format.
TMarker* drawPythiaParticleMarker(const Pythia8::Particle& p,
                                  int color = -1,
                                  int style = -1,
                                  double size = 1) {
  static TMarker* marker = new TMarker();
  if (style == -1) {
    if (p.charge() > 0)       marker->SetMarkerStyle(2);
    else if (p.charge() < 0)  marker->SetMarkerStyle(5);
    else                      marker->SetMarkerStyle(24);
  } else marker->SetMarkerStyle(style);
  if (color == -1) {
    if (p.charge() > 0)       marker->SetMarkerColor(kAzure);
    else if (p.charge() < 0)  marker->SetMarkerColor(kRed+1);
    else                      marker->SetMarkerColor(kGreen+2);
  } else marker->SetMarkerColor(color);
  if (p.charge() == 0) marker->SetMarkerSize(0.8*size);
  else marker->SetMarkerSize(size);
  double particle_longitude[2] = {p.eta(), p.y()};
  marker->DrawMarker(particle_longitude[cylinder_uses_rapidity], p.phi());
  return marker;
}

