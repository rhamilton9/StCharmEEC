// Generate some figures summarizing the vertex cuts and vertex quality from analyzed data.
//
// ----------------- Changelog -----------------
//  - 4/30/2026  ::  Created (with directory structure of extract_signal.cc)

// Inclusions
#include "config.h"

// Compiler flags
//#define include_MC


// Summarize the primary vertex reconstruction and cut efficiencies
//
// The input to the method is the modifier on the file wich should be
//         ../data.nosync/test_pp2012_[modifier].root
void analyze_primary_vertex(const char *modifier = "VertexCut") {
  
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
  const int n_directories = 1;
  char list_directories[n_directories][50] = {
    "primary_vertex"
  };
  const int n_subdirectories = 6;
  char list_subdirectories[n_directories][n_subdirectories][50] = {
    {"vertex_position", "vertex_cut_efficiency", "vertex_tracks", "", "", ""}
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
  
  TCanvas* canvas_solo = new TCanvas();
  gPad->SetTicks(1,1);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.12);
  
  TCanvas* canvas_duos = new TCanvas();
  canvas_duos->SetCanvasSize(1600, 600);
  canvas_duos->Divide(2,1);
  for (int i = 0; i < 2; ++i) {
    canvas_duos->cd(i+1);
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.9);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.11);
  }// End duo canvas setup
  
  TCanvas* canvas_quad = new TCanvas();
  canvas_quad->SetCanvasSize(800, 600);
  canvas_quad->Divide(2,2);
  for (int i = 0; i < 4; ++i) {
    canvas_quad->cd(i+1);
    gPad->SetTicks(1,1);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.105);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.11);
  }// End quad canvas setup
  
  // Colors from input hex codes
  const int n_colors = sizeof(plot_colors) / sizeof(plot_colors[0]);
  Int_t plot_TColors[n_colors];
  for (int c = 0; c < n_colors; ++c) plot_TColors[c] = getTColorFromHex(plot_colors[c]);
  double plot_fill_alpha = 0.15;
  
  // Other helpful plotting settings and variables
  char coord_label[4] = {'x','y','z','r'};
  char coord_label_2D[2][2] = {{'x','y'},{'z','r'}};
  
  // Good settings for pp2012
//  double axis_range_draw_limits[4][2] = {
//    {-0.12, 0.25}, // x
//    {-0.12, 0.25}, // y
//    {-83, 83}, // z
//    {0, 0.275}, // r
//  };
  
  // Good settings for dAu2016
  double axis_range_draw_limits[4][2] = {
    {-0.50, 0.50}, // x
    {-0.50, 0.50}, // y
    {-83, 83}, // z
    {-0.50, 0.50}, // r
  };
  
  
  
  // Get histograms from file
  // Vertex position histograms
  TH1F* hist_1D_vertex_coordinate[4];
  hist_1D_vertex_coordinate[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("TPC_Vx");
  hist_1D_vertex_coordinate[1] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("TPC_Vy");
  hist_1D_vertex_coordinate[2] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("TPC_Vz");
  hist_1D_vertex_coordinate[3] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("TPC_Vr");
  
  TH2F* hist_2D_vertex_coordinate[2];
  hist_2D_vertex_coordinate[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH2F>("TPC_PV_XY");
  hist_2D_vertex_coordinate[1] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH2F>("TPC_PV_ZR");
  
  // Vertex cut efficiency/information histograms
  TH1F* hist_vertex_cut_counts = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("VertexCutLevels");
  TH2F* hist_TPC_VPD_corr[2];
  hist_TPC_VPD_corr[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH2F>("TPCVz_VPDVz_picoDST");
  hist_TPC_VPD_corr[1] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH2F>("TPCVz_VPDVz_cut");
  
  // Vertex track information histograms
  TH1F* hist_vertex_tracks[3];
  hist_vertex_tracks[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("NTracks");
  hist_vertex_tracks[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("NTracksMain");
  hist_vertex_tracks[0] = datafile->Get<TDirectoryFile>("PrimaryVertices")->Get<TH1F>("NTracksPileup");
  
  
  // *--------------------------- Plotting :: Vertex position histograms
  
  // Plot 1D vertex position information histograms
  float z_cut[2];
  TLine* zero_line = new TLine();
  zero_line->SetLineStyle(7);
  zero_line->SetLineColor(plot_TColors[3]);
  for (int i_solo = 0; i_solo < 2; ++i_solo) {
    for (int i_hist = 0; i_hist < 4; ++i_hist) {
      if (i_solo) canvas_solo->cd();
      else        canvas_quad->cd(i_hist + 1);
      
      // Add plot settings and constrain to the relevant regions
      hist_1D_vertex_coordinate[i_hist]->GetXaxis()->SetRangeUser(axis_range_draw_limits[i_hist][0],
                                                                  axis_range_draw_limits[i_hist][1]);
      switch (i_hist) {
        case 2: // Longutidinal coordinate
          hist_1D_vertex_coordinate[i_hist]->SetLineColor(plot_TColors[1]);
          hist_1D_vertex_coordinate[i_hist]->SetFillColorAlpha(plot_TColors[1], plot_fill_alpha);
          break;
        default: // Transverse coordinates
          hist_1D_vertex_coordinate[i_hist]->SetLineColor(plot_TColors[0]);
          hist_1D_vertex_coordinate[i_hist]->SetFillColorAlpha(plot_TColors[0], plot_fill_alpha);
      }// End of plot region adjustment
      
      // Modify axes and draw text to the plot panels
      
      hist_1D_vertex_coordinate[i_hist]->SetTitle("");
      hist_1D_vertex_coordinate[i_hist]->GetXaxis()->SetTitle(Form("TPC V_{%c} [cm]",coord_label[i_hist]));
      hist_1D_vertex_coordinate[i_hist]->GetYaxis()->SetTitle("Final PV Counts");
      hist_1D_vertex_coordinate[i_hist]->GetYaxis()->SetTitleOffset(1.42);
      
      // Plot to populate the panel
      hist_1D_vertex_coordinate[i_hist]->Draw("hist");
      
      // Draw label text on the dataset
      double delta_v = 0.055;
      drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.3*(i_hist == 4) + 0.15, 0.83);
      drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*1);
      drawText(Form("#it{N}^{PV}= %.f",hist_1D_vertex_coordinate[i_hist]->Integral()), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*2);
      drawText(Form("#mu = %.4f",hist_1D_vertex_coordinate[i_hist]->GetMean(1)), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*3.5);
      drawText(Form("#sigma = %.4f",hist_1D_vertex_coordinate[i_hist]->GetStdDev(1)), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*4.5);
      
      drawText(Form("PV #bf{#it{%c}}",coord_label[i_hist]), 0.93, 0.8, true, kBlack, 0.06);
      
      drawText(Form("Modifier #it{%s}",modifier),
               1.0 - gPad->GetRightMargin(),
               1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
      
      // Draw line for the zero axis
      if (i_hist != 3) zero_line->DrawLine(0, 0, 0,
                                           hist_1D_vertex_coordinate[i_hist]->GetMaximum()*1.05);
      
      
      // Save solo
      if (i_solo) canvas_solo->SaveAs(Form("%s/%s/%s/%s_%s_V%c.pdf",
                                           base_directory,                    // ../data/[modifier]
                                           list_directories[0],               // /primary_vertex
                                           list_subdirectories[0][0],         // /vertex_position
                                           dataset_label,
                                           list_subdirectories[0][0],
                                           coord_label[i_hist] ) );
    }// End of hist loop
  }// End of solo/composite plotting loop
  
  canvas_quad->SaveAs(Form("%s/%s/%s/%s_%s_all.pdf",
                           base_directory,                    // ../data/[modifier]
                           list_directories[0],               // /primary_vertex
                           list_subdirectories[0][0],         // /vertex_position
                           dataset_label,
                           list_subdirectories[0][0]) );
  
  
  
  // 2D vertex plots (x,y) and (z,r)
  canvas_solo->SetRightMargin(0.12);
  canvas_solo->SetLeftMargin(0.10);
  for (int i_solo = 0; i_solo < 2; ++i_solo) {
    for (int i_hist = 0; i_hist < 2; ++i_hist) {
      if (i_solo) canvas_solo->cd();
      else        canvas_duos->cd(i_hist + 1);
      
      // Add plot settings and constrain to the relevant regions
      if (!i_hist) { // 0 :: (x, y)
        hist_2D_vertex_coordinate[i_hist]->GetXaxis()->SetRangeUser(axis_range_draw_limits[0][0],
                                                                    axis_range_draw_limits[0][1]);
        hist_2D_vertex_coordinate[i_hist]->GetYaxis()->SetRangeUser(axis_range_draw_limits[1][0],
                                                                    axis_range_draw_limits[1][1]);
      } else { // 1 :: (z, r)
        hist_2D_vertex_coordinate[i_hist]->GetXaxis()->SetRangeUser(axis_range_draw_limits[2][0],
                                                                    axis_range_draw_limits[2][1]);
        hist_2D_vertex_coordinate[i_hist]->GetYaxis()->SetRangeUser(axis_range_draw_limits[3][0],
                                                                    axis_range_draw_limits[3][1]);
      }// End of plot region adjustment
      
      // Modify axes and draw text to the plot panels
      hist_2D_vertex_coordinate[i_hist]->SetTitle("");
      hist_2D_vertex_coordinate[i_hist]->GetXaxis()->SetTitle(Form("TPC V_{%c} [cm]",coord_label_2D[i_hist][0]));
      hist_2D_vertex_coordinate[i_hist]->GetYaxis()->SetTitle(Form("TPC V_{%c} [cm]",coord_label_2D[i_hist][1]));
      hist_2D_vertex_coordinate[i_hist]->GetZaxis()->SetTitle("Final PV Counts");
      hist_2D_vertex_coordinate[i_hist]->GetZaxis()->SetTitleOffset(1.42);
      
      // Plot to populate the panel
      hist_2D_vertex_coordinate[i_hist]->Draw("colz");
      
      // Draw label text on the dataset
      double delta_v = 0.055;
      drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.3*(i_hist == 4) + 0.15, 0.83);
      drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*1);
      drawText(Form("#it{N}^{PV}= %.f",hist_1D_vertex_coordinate[i_hist]->Integral()), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*2);
      drawText(Form("#mu = %.4f",hist_1D_vertex_coordinate[i_hist]->GetMean(1)), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*3.5);
      drawText(Form("#sigma = %.4f",hist_1D_vertex_coordinate[i_hist]->GetStdDev(1)), 0.3*(i_hist == 4) + 0.15, 0.83 - delta_v*4.5);
      
      drawText(Form("PV (#bf{#it{%c}}, #bf{#it{%c}})",
                    coord_label_2D[i_hist][0],
                    coord_label_2D[i_hist][1]),
               0.95 - gPad->GetRightMargin() - i_hist*0.24,
               0.8 - gPad->GetTopMargin(), true, kBlack, 0.06);
      
      drawText(Form("Modifier #it{%s}",modifier),
               1.0 - gPad->GetRightMargin(),
               1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
      
      // Draw lines for the axes
      zero_line->DrawLine(0, axis_range_draw_limits[2*i_hist+1][0],
                          0, axis_range_draw_limits[2*i_hist+1][1]);
      if (!i_hist)
        zero_line->DrawLine(axis_range_draw_limits[0][0], 0, axis_range_draw_limits[0][1], 0);
      
      
      // Save solo
      if (i_solo) canvas_solo->SaveAs(Form("%s/%s/%s/%s_%s2D_V%cV%c.pdf",
                                           base_directory,                    // ../data/[modifier]
                                           list_directories[0],               // /primary_vertex
                                           list_subdirectories[0][0],         // /vertex_position
                                           dataset_label,
                                           list_subdirectories[0][0],
                                           coord_label_2D[i_hist][0],
                                           coord_label_2D[i_hist][1] ) );
    }// End of hist loop
  }// End of solo/composite plotting loop
  
  canvas_duos->SaveAs(Form("%s/%s/%s/%s_%s2D_all.pdf",
                           base_directory,                    // ../data/[modifier]
                           list_directories[0],               // /primary_vertex
                           list_subdirectories[0][0],         // /vertex_position
                           dataset_label,
                           list_subdirectories[0][0]) );
  
  
  // *--------------------------- Plotting :: Vertex cut efficiency/information histograms
  
  
  // 2D VPD / TPC Vertex plots
  char label_vertexcut[2][10] = {"raw","cut"};
  for (int i_solo = 0; i_solo < 2; ++i_solo) {
    for (int i_hist = 0; i_hist < 2; ++i_hist) {
      if (i_solo) canvas_solo->cd();
      else        canvas_duos->cd(i_hist + 1);
      gPad->SetLogz();
      
      // Add plot settings and constrain to the relevant regions
      int text_TColor;
      double base = 0.83;
      if (!i_hist) { // 0 :: uncut/raw/all PV
        text_TColor = plot_TColors[4];
        base = 0.80;
      } else { // 1 :: After applying V_z matching cuts
        text_TColor = kBlack;
      }// End of plot region adjustment
      
      // Modify axes and draw text to the plot panels
      hist_TPC_VPD_corr[i_hist]->SetTitle("");
      hist_TPC_VPD_corr[i_hist]->GetXaxis()->SetTitle("TPC V_{z} [cm]");
      hist_TPC_VPD_corr[i_hist]->GetYaxis()->SetTitle("VPD V_{z} [cm]");
      hist_TPC_VPD_corr[i_hist]->GetZaxis()->SetTitle("Final PV Counts");
      hist_TPC_VPD_corr[i_hist]->GetZaxis()->SetTitleOffset(1.42);
      
      // Plot to populate the panel
      hist_TPC_VPD_corr[i_hist]->Draw("colz");
      
      // Draw label text on the dataset
      double delta_v = 0.055;
      drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.3*(i_hist == 4) + 0.2, base, false, text_TColor);
      drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.3*(i_hist == 4) + 0.2, base - delta_v*1, false, text_TColor);
      drawText(Form("#it{N}^{PV}= %.f",hist_1D_vertex_coordinate[i_hist]->Integral()), 0.3*(i_hist == 4) + 0.2, base - delta_v*2, false, text_TColor);
      drawText(Form("#mu = %.4f",hist_1D_vertex_coordinate[i_hist]->GetMean(1)), 0.3*(i_hist == 4) + 0.2, base - delta_v*3.5, false, text_TColor);
      drawText(Form("#sigma = %.4f",hist_1D_vertex_coordinate[i_hist]->GetStdDev(1)), 0.3*(i_hist == 4) + 0.2, base - delta_v*4.5, false, text_TColor);
      
      drawText("PV (#bf{#it{z}}^{TPC}, #bf{#it{z}}^{VPD})",
               0.95 - gPad->GetRightMargin(),
               0.9 - gPad->GetTopMargin(), true, kBlack, 0.06);
      
      drawText(Form("Modifier #it{%s}",modifier),
               1.0 - gPad->GetRightMargin(),
               1.03 - gPad->GetTopMargin(), true, kBlack, 0.05);
      
      // Draw lines for the axes
//      zero_line->DrawLine(0, axis_range_draw_limits[2*i_hist+1][0],
//                          0, axis_range_draw_limits[2*i_hist+1][1]);
//      if (!i_hist)
//        zero_line->DrawLine(axis_range_draw_limits[0][0], 0, axis_range_draw_limits[0][1], 0);
      
      
      // Save solo
      if (i_solo) canvas_solo->SaveAs(Form("%s/%s/%s/%s_TPC_VPD_%s.pdf",
                                           base_directory,                    // ../data/[modifier]
                                           list_directories[0],               // /primary_vertex
                                           list_subdirectories[0][1],         // /vertex_efficiency
                                           dataset_label,
                                           label_vertexcut[i_hist] ) );
    }// End of hist loop
  }// End of solo/composite plotting loop
  
  canvas_duos->SaveAs(Form("%s/%s/%s/%s_TPC_VPD_all.pdf",
                           base_directory,                    // ../data/[modifier]
                           list_directories[0],               // /primary_vertex
                           list_subdirectories[0][1],         // /vertex_efficiency
                           dataset_label) );
  
  
  
  // Vertex Cut Levels count information
  bool draw_logscale_vertexcuts = true;
  gPad->SetBottomMargin(0.15);
  const ULong64_t known_total_events =  1231716967;
  char bin_cutlevel_labels_5[5][20] = {"Processed Events","VPDMB","VPD/TPC Matched","V_{z} Rank > 0","V_z #in [-50,50]"};
  char bin_cutlevel_labels_6[6][20] = {"Processed Events","Good Runs","VPDMB","VPD/TPC Matched","V_{z} Rank > 0","V_z #in [-50,50]"};
  for (int i_bin = 1; i_bin <= hist_vertex_cut_counts->GetXaxis()->GetNbins(); ++i_bin) {
    
    // Setup hist for bar chart drawing
    hist_vertex_cut_counts->SetBinContent(i_bin, hist_vertex_cut_counts->GetBinContent(i_bin));
    switch (hist_vertex_cut_counts->GetXaxis()->GetNbins()) {
      case 5:
        hist_vertex_cut_counts->GetXaxis()->SetBinLabel(i_bin, bin_cutlevel_labels_5[i_bin-1]);
        break;
      case 6:
        hist_vertex_cut_counts->GetXaxis()->SetBinLabel(i_bin, bin_cutlevel_labels_6[i_bin-1]);
        break;
    }// End of bin label switch
  }
  
  canvas_solo->cd();
  if (draw_logscale_vertexcuts) {
    gPad->SetLogy();
    hist_vertex_cut_counts->GetYaxis()->SetRangeUser(1e6, 6*known_total_events);
  } else {
    hist_vertex_cut_counts->GetYaxis()->SetRangeUser(0, 1.25*known_total_events);
  }
  
  hist_vertex_cut_counts->GetXaxis()->SetTitleOffset(1.7);
  hist_vertex_cut_counts->SetFillColorAlpha(kBlack, 0.05);
  hist_vertex_cut_counts->Draw("hist");
  
  zero_line->DrawLine(0, known_total_events, hist_vertex_cut_counts->GetXaxis()->GetNbins(), known_total_events);
  drawText("All Events",
           0.95 - gPad->GetRightMargin(),
           0.87 - gPad->GetTopMargin(), true, kGray+2, 0.05);
  
  drawText(Form("#bf{STAR} %s",dataset_label_plotstring), 0.15, 0.92-gPad->GetTopMargin(), false, kBlack);
  drawText(Form("#it{%s} #sqrt{s} = %.f GeV",collision_system, sqrt_s), 0.15, 0.87-gPad->GetTopMargin(), false, kBlack);
  
  // Draw event counts in each bin
//  if (draw_logscale_vertexcuts) {
//    for (int i_bin = 0; i_bin < hist_vertex_cut_counts->GetXaxis()->GetNbins(); ++i_bin)
//  } else {
//    
//    
//  }
  
  
  canvas_solo->SaveAs(Form("%s/%s/%s/%s_PV_cut_summary.pdf",
                           base_directory,                    // ../data/[modifier]
                           list_directories[0],               // /primary_vertex
                           list_subdirectories[0][1],         // /vertex_efficiency
                           dataset_label) );
  
  // *--------------------------- Plotting :: Vertex track information histograms
  
  
  // TODO :: Vertex Track plots
  
  return;
}

