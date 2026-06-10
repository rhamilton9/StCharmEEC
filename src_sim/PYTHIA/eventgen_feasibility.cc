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
/// Changelog ::
///   - This file broken away from eventjen\_Dzero\_jets.cc on 6/3/2026
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
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "../../utils/root_draw_tools.h"
#include "config_local.h"
#include "ptbin.h"
 
using namespace Pythia8;

const bool debug = false;

const double pion_mass = 0.13957039;

// Use a restricted pThat range for this feasibility study which should be reasonable!
double pt_binedge[2] = {15, 35};     // pThat range to generate hard scatterings in PYTHIA
double EEC_pt_binedge[2] = {20, 30}; // pT range of final clustered jets to store in the EEC histograms
double cross_section = 3.98029761e-4;
int nevent = 1000;


// Energy correlator settings
int nbin_EEC = 100;
double EEC_range[2] = {2e-4, 7e-1};



// Define a user class for info associated with pseudojets
class LocalInfo: public fastjet::PseudoJet::UserInfoBase {
  int _pdg_id;
  bool _is_charged;
public:
  // Constructor
  LocalInfo(int id) {_pdg_id = id; _is_charged = false;}
  LocalInfo(int id, bool charge) {_pdg_id = id; _is_charged = charge;}
  // PDG ID getter
  int pdg_id() const {return _pdg_id;}
  bool isCharged() const {return _is_charged;}
};



// Forward declare relevant methods. These will be defined at the end of the file
fastjet::JetAlgorithm getJetAlgorithm(const char*);

double getSquareEtaPhiDistance(double eta1, double eta2, double phi1, double phi2);

bool isKnownProcess(int process);

std::string getProcessString(Event& process);

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
    std::cerr << "Error in <eventgen_feasibility>: Not enough arguments!" << std::endl;
    std::cout << "arg1: Thread to seed randomizer (for parallelization)." << std::endl;
    return 0;
  }
  
  // Thread setup, throw error if thread input is out of range
  int this_thread = std::atoi(argv[1]);
//  const int total_threads = getNbins(sqrt_s);
//  if (this_thread >= total_threads) {
//    std::cerr << "Error in <eventgen_jets.cc>: requested pTHatBin " << this_thread << " out of bin range " << total_threads << std::endl;
//    return 0;
//  }
  
  // Gather bin settings for this thread
  pTBinSettings settings = getBinSettings(sqrt_s, this_thread, local_flag_xsec_scale);
//  const float pt_binedge[2] = {settings.binedge_low, settings.binedge_high};
//  const float cross_section = settings.cross_section;
//  const int nevent = (int) (settings.nevent * nevent_scale_factor); // scaled event number for testing
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
  pythia.readString(Form("Random:seed = %i", this_thread + 1000));
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
  
  // Don't allow D0 mesons to decay, so that we can cluster 
  // them into our jets as we would in data
  pythia.readString("421:mayDecay = off");
  
  // Only allow decays that occur within the primary vertex pointing resolution
  // This simulates the removal/reconstruction of weak decays
//  pythia.readString("ParticleDecays:limitRadius = on");
//  pythia.readString("ParticleDecays:rMax = 10"); // 100 micron resolution
  // Maybe has a wierd bug?? Increases D0-Tagged jet counts a lot when making the radius smaller which seems odd
  
  // For now just remove weak decays manually
  pythia.readString("130:mayDecay = off");  // K0L
  pythia.readString("310:mayDecay = off");  // K0S
  pythia.readString("3122:mayDecay = off"); // Lambda
  pythia.readString("3222:mayDecay = off"); // Sigma+
  pythia.readString("3212:mayDecay = off"); // Sigma0
  pythia.readString("3112:mayDecay = off"); // Sigma-
  pythia.readString("3322:mayDecay = off"); // Xi0
  pythia.readString("3312:mayDecay = off"); // Xi-
  pythia.readString("3334:mayDecay = off"); // Omega
  
  // initialize
  pythia.init();
   
  std::cout << "Beginning generation for thread " << this_thread << std::endl;
  
  
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
  
  
  
  // Generate uniform binning in log so it looks reasonable when plotted on logscale
  // This is important to get a fine resolution on the features of the EEC
  std::vector<double> EEC_binedge = std::vector<double>();
  double logbinedge[2] = {TMath::Log(EEC_range[0]), TMath::Log(EEC_range[1])};
  double delta_log_binedge = (logbinedge[1] - logbinedge[0]) / nbin_EEC;
  for (int i = 0; i <= nbin_EEC; ++i) EEC_binedge.push_back(TMath::Exp(logbinedge[0] + i * delta_log_binedge));
  
  // Record constituent, EEC information as a histogram
  // {inclusive, D0-tagged, non-D0-tagged}
  TH1D* hist_event_jetmult[3];
  TH1D* hist_fulljet_constituents[3];
  TH1D* hist_charjet_constituents[3];
  TH1D* EEC_fulljet[3];
  TH1D* EEC_charjet[3];
  char D0tag_string[3][20] = {"inclusive","D0tagged","nocharm"};
  for (int i = 0; i < 3; ++i) {
    // Jet multiplicity for each event, important to know how many jets there are
    // Note that since neutral energy is included in the charged jets, the multiplicity of full/char jets is equal
    hist_event_jetmult[i] = new TH1D(Form("eventlevel_jetmult_ptbin%.f-%.f_%s",EEC_pt_binedge[0],EEC_pt_binedge[1],D0tag_string[i]),
                                     ";Count of Jets in Event;Event Count",10, 0, 10);
    
    // Jet constituent multiplicities, important to know how many times to sample the EEC
    // These differ strongly for neutral/charged jets and so are stored separately
    hist_fulljet_constituents[i] = new TH1D(Form("fulljet_constituents_ptbin%.f-%.f_%s",EEC_pt_binedge[0],EEC_pt_binedge[1],D0tag_string[i]),
                                            ";Count of Constituents in Full Jet;Jet Count",50, 0, 50);
    hist_charjet_constituents[i] = new TH1D(Form("charjet_constituents_ptbin%.f-%.f_%s",EEC_pt_binedge[0],EEC_pt_binedge[1],D0tag_string[i]),
                                            ";Count of Constituents in Charged Jet;Jet Count",50, 0, 50);
    
    
    // EEC--use log uniform binning
    EEC_fulljet[i] = new TH1D(Form("EEC_fulljet_ptbin%.f-%.f_%s",EEC_pt_binedge[0],EEC_pt_binedge[1],D0tag_string[i]),
                              ";R_{L};d#Sigma^{Full Incl.}/dR_{L}",
                              nbin_EEC, EEC_binedge.data());
    EEC_charjet[i] = new TH1D(Form("EEC_charjet_ptbin%.f-%.f_%s",EEC_pt_binedge[0],EEC_pt_binedge[1],D0tag_string[i]),
                              ";R_{L};d#Sigma^{Charged Incl.}/dR_{L}",
                              nbin_EEC, EEC_binedge.data());
  }// End of D0-tagged hist construction
  
  
  
  // Setup jet definiton for fastjet
  // Using algorithm and jet radius according to config file
  // Recombination Scheme options: See
  //   https://fastjet.fr/repo/doxygen-3.5.1/namespacefastjet.html#a46fcc48dcb00a10557d773e328153bcb
  fastjet::JetAlgorithm  jet_algo         = getJetAlgorithm(algo_string_short);
  fastjet::JetDefinition jet_definition   = fastjet::JetDefinition(jet_algo, jet_radius, fastjet::E_scheme);
  
  // 2D Histograms that will display solid regions for reconstructed jets
  // Purely for illustrative purposes to show that the jet finder is working properly
  const double dA = 4*TMath::Pi()*longitudinal_acceptance[cylinder_uses_rapidity]/((double)(nbins_rap_display*nbins_phi_display));
  TH2D* jet_diagram2D;
  jet_diagram2D = new TH2D(Form("jet_display"),
                                      Form(";%s;Azimuth #phi [rad]; #it{p}_{T}^{jet} [GeV]",longitudinal_longlabel[cylinder_uses_rapidity]),
                                      nbins_rap_display,-max_eta,max_eta,nbins_phi_display,-TMath::Pi(),TMath::Pi());
  jet_diagram2D->GetYaxis()->SetTitleOffset(0.8);
  jet_diagram2D->GetZaxis()->SetTitleOffset(0.4);
  
   
  // Setup for root plotting
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBlueRedYellow);
  TCanvas* canvas = new TCanvas();
  canvas->SetCanvasSize(1200, 600);
  canvas->SetMargin(0.07, 0.35, 0.1, 0.12);
  gPad->SetLogz();
  gPad->SetTicks(1,1);
  jet_diagram2D->GetZaxis()->SetRangeUser(jetcut_minpT_full, 300);
  
  // Open pdf file for event display
  TString plotfile = Form("out_plots/%s/%s_eventdisplay_ptbin%.f-%.f-thread%i.pdf",
                          gen_info_string,gen_info_string,pt_binedge[0],pt_binedge[1],this_thread);
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
    int count_D0 = 0;
    std::vector<Particle> particles;
    std::vector<Particle> reco_D0;
    std::vector<Particle> missed_charm;
    std::vector<fastjet::PseudoJet> stable_particles;
    Vec4 pSum;
    for (int i = 0; i < cEvent.size(); ++i) { // Get info from final state particles
      Particle &p = cEvent[i];
      
      // Only add final state hadrons to clustering
      if (!cEvent[i].isFinal()) continue;
      
      // Record the total final state 4-momentum to check against momentum conservation
      pSum += Vec4(p.px(), p.py(), p.pz(), p.e());
      
      // Check that the particle is within the eta range
      if (std::fabs(p.eta()) > max_eta) continue;
      
      // Ignore hadrons with pT is lower than detector efficiency threshold
      // D0 are excluded from this since they are reconstructed from daughters that will
      // generally have pT that pass this threshold (though admittedly not always)
      //
      // PDG MC ID for D0 is 421.
      // Note that D0 decay is currently turned off, so these D0 will
      // be regarded as final state hadrons from PYTHIA's point of view
      if (std::abs(p.id()) == 421) {
        if (p.pT() < D0_minpT) continue;
        ++count_D0;
        std::cout << "Push back valid final state D0!" << std::endl;
      } else if (cEvent[i].pT() < pTmin_hadron) continue; // Other hadrons
      
      // Append to particle vector for jet clustering
      particles.push_back(p);
      fastjet::PseudoJet fj_particle(p.px(),p.py(),p.pz(),p.e());
      fj_particle.set_user_info(new LocalInfo(p.id(), p.isCharged())); // Store PDG ID in local class
      
      stable_particles.push_back(fj_particle);
      
      
      // Fill histograms at the track/particle level
      hist_pT[0]->Fill(p.pT());
      hist_2D[0][0]->Fill(p.y(), p.phi());
      hist_2D[1][0]->Fill(p.eta(),p.phi());
    }// End of PYTHIA particle loop
    
    // Record event multiplicity
    mult = stable_particles.size();
    pSum /= cEvent[0].e(); // total event energy
    delta_E = std::abs(pSum.e() - 1);  // normalized energy loss fraction DeltaE / E
    
    // Make a grid of ghost particles, add them to the stack
    // This is necessary to find the jet area, defined as the grid
    // area of ghost particles that are clustered into a jet.
    fastjet::PseudoJet ghost;
    double ghost_pt = 1e-100;
    for (int ir = 1; ir <= nbins_rap_display; ++ir) {
      // For massless particles like ghosts, rapidity = pseudorapidity
      double longitude = jet_diagram2D->GetXaxis()->GetBinCenter(ir);
      for (int iphi = 1; iphi <= nbins_phi_display; ++iphi) {
        double phi = jet_diagram2D->GetYaxis()->GetBinCenter(iphi);
        ghost.reset_momentum_PtYPhiM(ghost_pt, longitude, phi, 0);
        stable_particles.push_back(ghost);
      }
    }// End of ghost loop
    
    // Default to no drawing this event
    flag_draw_event = false;
    
    //========================================================================== Jet Clustering
    
    
    // Run jet clustering with fastjet -- All Jets
    // Note that the min pT cut is handled by fastjet::inclusive_jets)
    fastjet::ClusterSequence clustering(stable_particles, jet_definition);
    std::vector<fastjet::PseudoJet> final_jets = sorted_by_pt(clustering.inclusive_jets(jetcut_minpT_full));
    
    ntot_jet = 0;
    int njet_D0tagged = 0;
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
      int ctot_chargedconstituents = 0;
      double ctot_area = 0;
      int charm_status = 0;
      for (fastjet::PseudoJet c : jet.constituents()) {
        // Check for D0
        if (c.has_user_info() && std::abs(c.user_info<LocalInfo>().pdg_id()) == 421) {
          std::cout << "D0 associated with jet; pT = \033[31m" << jet.pt() << "\033[39m, (eta,phi) = (";
          std::cout << jet.eta() << ',' << jet.phi() << ")." << std::endl;
          charm_status = 1;
          flag_draw_event = true;
        }// End of D0 check
        
        // Tabulate jet area using ghost grid, fill constituent histograms
        if (c.pt() < 1e-50) ctot_area += dA; // grid based jet area calculation
        else { // Fill constituent jet histograms with non-ghost particles
          hist_pT[2]->Fill(c.pt());
          hist_2D[0][2]->Fill(c.rap(), c.phi_std());
          hist_2D[1][2]->Fill(c.eta(), c.phi_std());
          ++ctot_constituents;
          if (!c.has_user_info() || !c.user_info<LocalInfo>().isCharged()) continue;
          ++ctot_chargedconstituents;
        }// End of area/constituent fill
      }// End of consituent loop
      
      // Cut jet if it fails to meet the area cut
      if (ctot_area < min_area) {
        cout << "Full jet cut low area :: (" << ctot_area/(3.14159*jet_radius*jet_radius) << ")*pi*R_{jet}^2" << endl;
        continue;
      }
      
      // record the jet count/information for this event
      ++ntot_jet;
      hist_fulljet_constituents[0]->Fill(ctot_constituents);
      hist_charjet_constituents[0]->Fill(ctot_chargedconstituents);
      if (charm_status == 1) {
        ++njet_D0tagged;
        hist_fulljet_constituents[1]->Fill(ctot_constituents);
        hist_charjet_constituents[1]->Fill(ctot_chargedconstituents);
      } else {
        hist_fulljet_constituents[2]->Fill(ctot_constituents);
        hist_charjet_constituents[2]->Fill(ctot_chargedconstituents);
      }
      
      
      
      //========================================================================== Compute Energy-Energy Correlators
      
      // Add to energy correlator histogram--full jets
      if (jet.pt() > EEC_pt_binedge[0] && jet.pt() < EEC_pt_binedge[1]) {
        double sq_jet_pt = TMath::Sq(jet.pt()) / 2; // 2 for double counting
        for (int i = 0; i < ctot_constituents; ++i) {
          
          if (constituents_ptsort[i].pt() == 1e-100) std::cout << "Warning :: Ghost included in EEC!" << std::endl;
          
          // Mass conditions: assume pion if charged, photon if neutral
          Double_t mass_assumption_i = pion_mass*(constituents_ptsort[i].has_user_info() && constituents_ptsort[i].user_info<LocalInfo>().isCharged());
          
          for (int j = i + 1; j < ctot_constituents; ++j) {
            
            // Compute eta-phi difference carefully with 2pi modulus for cylinder wrapping
            Double_t R_L = getSquareEtaPhiDistance(constituents_ptsort[i].eta(), constituents_ptsort[j].eta(),
                                                   constituents_ptsort[i].phi(), constituents_ptsort[j].phi());
            
            
            // Mass conditions: assume pion if charged, photon if neutral
            Double_t mass_assumption_j = pion_mass*(constituents_ptsort[j].has_user_info() && constituents_ptsort[j].user_info<LocalInfo>().isCharged());
            
            // Get EEC weight for this pair
            Double_t EEC_weight = (TMath::Hypot(mass_assumption_i, constituents_ptsort[i].pt()) *
                                   TMath::Hypot(mass_assumption_j, constituents_ptsort[j].pt()) ) / sq_jet_pt;
            
            // Fill inclusive
            EEC_fulljet[0]->Fill(R_L, EEC_weight);
            
            // Fill D0-tagged hist if there is a D0 in the jet, otherwise fill non-D0
            if (charm_status == 1) {
              EEC_fulljet[1]->Fill(R_L, EEC_weight);
            } else {
              EEC_fulljet[2]->Fill(R_L, EEC_weight);
            }
            
            // Add only charged pairs to Charged EEC, but still use full jet energy as pT weight
            if (!constituents_ptsort[i].has_user_info() || !constituents_ptsort[i].user_info<LocalInfo>().isCharged() ||
                !constituents_ptsort[j].has_user_info() || !constituents_ptsort[j].user_info<LocalInfo>().isCharged() ) continue;
            
            // Fill inclusive
            EEC_charjet[0]->Fill(R_L, EEC_weight);
            
            // Fill D0-tagged hist if there is a D0 in the jet, otherwise fill non-D0
            if (charm_status == 1) {
              EEC_charjet[1]->Fill(R_L, EEC_weight);
            } else {
              EEC_charjet[2]->Fill(R_L, EEC_weight);
            }
            
          }// End of j consitutent loop
        }// End of i constituent loop
      }// End of EEC Fill for this jet
      
    }// End of final state full jet loop
    
    
    //========================================================================== Event Display
    
    // Choose whether to draw the event based on its characteristics
    
    // Draw at least 10 events regardless of their jet content
    flag_draw_event = flag_draw_event || (iEvent % (nevent/10) == 0) || iEvent == 44;
    
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
      jet_diagram2D->Reset();
      Int_t jet_fill_color = kWhite;
      TPaletteAxis *palette = static_cast<TPaletteAxis*>(jet_diagram2D->GetListOfFunctions()->FindObject("palette"));
      
      // loop back over jets and fill the display histograms
      for (fastjet::PseudoJet jet : final_jets) {
        // Do not draw jets that fail cuts.
        
        double jet_long[2] = {jet.eta(), jet.rap()};
        if (abs(jet_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]-jet_radius) continue;
        
        std::vector<fastjet::PseudoJet> constituents_ptsort = sorted_by_pt(jet.constituents());
        if ((constituents_ptsort.at(0)).pt() < pTmin_jetcore) continue;
        
        if (jet.pt() < jetcut_minpT_full) continue;
        
        for (fastjet::PseudoJet c : jet.constituents()) if (c.pt() < 1e-50) {
          double constituent_long[2] = {c.eta(), c.rap()};
          jet_diagram2D->Fill(constituent_long[cylinder_uses_rapidity], c.phi_std(), jet.pt());
        }
        
        // Doesn't work -- method not implemented in TPaletteAxis yet???
//        jet_fill_color = palette->GetBinColor(jet_diagram2D->GetXaxis()->FindBin(jet_long[cylinder_uses_rapidity]),
//                                              jet_diagram2D->GetYaxis()->FindBin(jet.phi_std()));
      }// End of final jet loop
      
      // Draw jets as area grids with their ghosts
      jet_diagram2D->Draw("colz");
      
      // Add description text
      drawText(Form("#it{#bf{PYTHIA}} Event #color[2]{%i}, #sqrt{s_{NN}} = %.2f TeV", iEvent, collision_energy_TeV), gPad->GetLeftMargin(), 0.95);
      drawText(Form("#it{p}_{T}^{Hard} in range [%.0f,%0.f]",pt_binedge[0],pt_binedge[1]), 1.0 - gPad->GetRightMargin(), 0.95, true);
      drawText(Form("Clustering : #it{%s}, R = %.1f, E#minusScheme",algo_string,jet_radius), 1.0 - gPad->GetRightMargin(), 0.90, true);
      drawText("#it{#bf{FastJet}} ver. 3.4.1", gPad->GetLeftMargin(), 0.90);
      drawText(Form("Event Process Code #it{%i}", gInfo->code()), 0.75 , 0.95, false);
      drawText(Form("(Parton Scattering %s)", getProcessString(cProcess).c_str()), 0.75, 0.90, false);
      
      // Make a legend of particle markers
      TLegend* leg_event_display = new TLegend(0.75, 0.25, 0.95, 0.75);
      leg_event_display->SetLineWidth(0);
      jet_diagram2D->SetFillColor(jet_fill_color);
      leg_event_display->AddEntry(jet_diagram2D, "Clustered Jets", "f");
      
      // Style is {positive, negative, neutral, D0}
      int marker_style[4] = {2, 5, 24, 42};
      int marker_color[4] = {kAzure, kRed+1, kGreen+2, kBlue-10};
      char leg_description[4][20] = {"Positive","Negative","Neutral","D^{0} Meson"};
      TMarker* leg_markers[4];
      for (int i_marker = 0; i_marker < 4; ++i_marker) {
        leg_markers[i_marker] = new TMarker();
        leg_markers[i_marker]->SetMarkerSize(2);
        leg_markers[i_marker]->SetMarkerColor(marker_color[i_marker]);
        leg_markers[i_marker]->SetMarkerStyle(marker_style[i_marker]);
        leg_event_display->AddEntry(leg_markers[i_marker], leg_description[i_marker], "p");
      }leg_event_display->Draw();
      
      
      // Draw markers for the event hadrons
      for (Particle &p:particles) {
        
        // Check if the particle is within the accptance
        double particle_long[2] = {p.eta(), p.y()};
        if (abs(particle_long[cylinder_uses_rapidity]) > longitudinal_acceptance[cylinder_uses_rapidity]) continue;
        
        // Draw D0 always
        if (std::abs(p.id()) == 421) {drawPythiaParticleMarker(p, kBlue-10, 42, 2); continue;}
        
        // Draw other hadrons
        if (!p.isFinal()) drawPythiaParticleMarker(p, kBlack, -1, 1); // Should never happen, this is a way of checking if nonfinal particles are leaking
        else if (p.pT() < pTmin_hadron) drawPythiaParticleMarker(p, kGray, -1, 1);
        else drawPythiaParticleMarker(p, -1, -1, 1);
      }// End of particle loop
        
      
      // Save the canvas to the display file
      canvas->Print(plotfile);
    }// End of event display draw
    
    
    //========================================================================== Write/Store Event-level Data
    
    // Store event level information to histograms
    hist_event_jetmult[0]->Fill(ntot_jet);
    hist_event_jetmult[1]->Fill(njet_D0tagged);
    hist_event_jetmult[2]->Fill(ntot_jet - njet_D0tagged);
    grandtotal_jet[0] += ntot_jet;
    grandtotal_jet[1] += njet_D0tagged;
    
  }//========================================================================== End of event loop
  canvas->Print(plotfile + "]"); // close event display file
  
  
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
  }// End of process loop
  
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
  
  for (int i_D0_tagged = 0; i_D0_tagged < 3; ++i_D0_tagged) {
    hist_event_jetmult[i_D0_tagged]->Write(hist_event_jetmult[i_D0_tagged]->GetName(), TObject::kOverwrite);
    hist_fulljet_constituents[i_D0_tagged]->Write(hist_fulljet_constituents[i_D0_tagged]->GetName(), TObject::kOverwrite);
    hist_charjet_constituents[i_D0_tagged]->Write(hist_charjet_constituents[i_D0_tagged]->GetName(), TObject::kOverwrite);
    EEC_fulljet[i_D0_tagged]->Write(EEC_fulljet[i_D0_tagged]->GetName(), TObject::kOverwrite);
    EEC_charjet[i_D0_tagged]->Write(EEC_charjet[i_D0_tagged]->GetName(), TObject::kOverwrite);
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
  std::cout << "Jet D0 Tagging:" << std::endl;
  std::cout << "Incl.\tD0-Tagged" << std::endl;
  std::cout << "------------------------" << std::endl;
  std::cout << grandtotal_jet[0] << "\t" << grandtotal_jet[1] << std::endl;
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
}// End of eventgen_feasibility::main











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

// Check if the PYTHIA hard process is a known one
bool isKnownProcess(int process) {
  switch (process) {
    case 111: // gg -> gg                       Gluon scattering
    case 112: // gg -> qqbar (uds)              Gluon fusion to light quark/antiquark
    case 113: // qg -> qg (any flavor)          Quark-gluon scattering
    case 114: // qq' -> qq'                     Quark scattering
    case 115: // qqbar -> gg                    Quark-Antiquark fusion to gluon pair
    case 116: // qqbar -> q'q'bar (q' uds)      Quark-Antiquark fusion to light quark/antiquark
    case 121: // gg -> ccbar                    Gluon fusion to charm/anticharm pair
    case 122: // qqbar -> ccbar                 Quark-Antiquark fusion to charm/anticharm pair
    case 123: // gg -> bbbar                    Gluon fusion to bottom/antibottom pair
    case 124: // qqbar -> bbbar                 Quark-Antiquark fusion to bottom/antibottom pair
      return true;
    default:
      return false;
  }
}

// Convert the PYTHIA hard process to a legible string
// e.g. cdbar -> dcbar, gg->uubar
// Format is TLatex, for overline
// Note that the method assumes 2->2 processes (of the type listed in isKnownProcess()
std::string getProcessString(Event& process) {
  char process_string[100];
  for (int i_process_candidate = 3; i_process_candidate <= 6; ++i_process_candidate) {
    if (i_process_candidate == 5) snprintf(process_string, 100, "%s#rightarrow", process_string);
    switch (process[i_process_candidate].id()) {
      case 1: // down quark
        snprintf(process_string, 100, "%sd", process_string); break;
      case -1: // anti-down quark
        snprintf(process_string, 100, "%s#bar{d}", process_string); break;
      case 2: // up quark
        snprintf(process_string, 100, "%su", process_string); break;
      case -2: // anti-up quark
        snprintf(process_string, 100, "%s#bar{u}", process_string); break;
      case 3: // strange quark
        snprintf(process_string, 100, "%ss", process_string); break;
      case -3: // anti-strange quark
        snprintf(process_string, 100, "%s#bar{s}", process_string); break;
      case 4: // charm quark
        snprintf(process_string, 100, "%sc", process_string); break;
      case -4: // anti-charm quark
        snprintf(process_string, 100, "%s#bar{c}", process_string); break;
      case 5: // bottom quark
        snprintf(process_string, 100, "%sb", process_string); break;
      case -5: // anti-bottom quark
        snprintf(process_string, 100, "%s#bar{b}", process_string); break;
      case 6: // top quark
        snprintf(process_string, 100, "%st", process_string); break;
      case -6: // anti-top quark
        snprintf(process_string, 100, "%s#bar{t}", process_string); break;
      case 21: // gluon
        snprintf(process_string, 100, "%sg", process_string); break;
      default:
        snprintf(process_string, 100, "%s?", process_string); break;
    }
  }// End of process loop
  std::string to_return(process_string);
  return to_return;
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

