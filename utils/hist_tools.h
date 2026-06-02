// A collection of small macros for histogram analysis
// Macros are organized roughly by function

// Variables for storing data that may be of interest in the daughter macro

//#include "interpolation_tools.h"

double lostdata_translate;
double lostdata_rotate;

//========================================================================== Search/Retrieval/Optimization

// Returns an integer which signals the type of TH1 inheriting class.
// Used to determine qualities of the histogram when not specified.
int getHistType(TH1* hist) {
  // Test most common types first
  if (TH1D* hist_TH1D = dynamic_cast<TH1D*>(hist))  return 11;   // TH1D
  if (TH1F* hist_TH1F = dynamic_cast<TH1F*>(hist))  return 12;   // TH1F
  
  if (TH2D* hist_TH2D = dynamic_cast<TH2D*>(hist))  return 21;   // TH2D
  if (TH2F* hist_TH2F = dynamic_cast<TH2F*>(hist))  return 22;   // TH2F
  
  if (TH3D* hist_TH3D = dynamic_cast<TH3D*>(hist))  return 31;   // TH3D
  if (TH3F* hist_TH3F = dynamic_cast<TH3F*>(hist))  return 32;   // TH3F
  
  // Test for more obscure types
  if (TH1I* hist_TH1I = dynamic_cast<TH1I*>(hist))  return 13;   // TH1I
  if (TH1C* hist_TH1C = dynamic_cast<TH1C*>(hist))  return 14;   // TH1C
  if (TH1S* hist_TH1S = dynamic_cast<TH1S*>(hist))  return 15;   // TH1S
  if (TH1K* hist_TH1K = dynamic_cast<TH1K*>(hist))  return 16;   // TH1K
  
  if (TH2I* hist_TH2I = dynamic_cast<TH2I*>(hist))  return 23;   // TH2I
  if (TH2C* hist_TH2C = dynamic_cast<TH2C*>(hist))  return 24;   // TH2C
  if (TH2S* hist_TH2S = dynamic_cast<TH2S*>(hist))  return 25;   // TH2S
  
  if (TH3I* hist_TH3I = dynamic_cast<TH3I*>(hist))  return 33;   // TH3I
  if (TH3C* hist_TH3C = dynamic_cast<TH3C*>(hist))  return 34;   // TH3C
  if (TH3S* hist_TH3S = dynamic_cast<TH3S*>(hist))  return 35;   // TH3S
  
  // Not recognized. (could be TProfile, TH2Poly, TGLTH3Composition)
  return 0;
}

// Returns the largest value from two input histograms
double getMaxFromHists(TH1* hist1, TH1* hist2) {
  double max1 = hist1->GetMaximum();
  double max2 = hist2->GetMaximum();
  if (max2 > max1) return max2;
  return max1;
}

// Returns the usual calculus definite integral of a hist within a region
// Treats bin contents as uniformly distributed within the bin
double definiteIntegral(TH1* hist, double lowbound = INT_MIN, double topbound = INT_MAX) {
  TAxis* axis = hist->GetXaxis();
  int nbin = axis->GetNbins();
  
  // Edge cases & signed integral considerations
  int int_sign = 1;
  if (axis->FindBin(lowbound) == axis->FindBin(topbound)) return hist->GetBinContent(axis->FindBin(lowbound)) * (topbound - lowbound);
  if (lowbound > topbound) {double temp = topbound; topbound = lowbound; lowbound = temp; int_sign = -1;}
  if (lowbound == INT_MIN || lowbound < axis->GetBinLowEdge(1))       lowbound = axis->GetBinLowEdge(1);
  if (topbound == INT_MAX || topbound > axis->GetBinLowEdge(nbin+1))  topbound = axis->GetBinLowEdge(nbin+1);
  
  // Compute integral
  int startbin = axis->FindBin(lowbound);
  double integral = hist->GetBinContent(startbin)*(axis->GetBinLowEdge(startbin+1) - lowbound);
  for (int ibin = startbin+1; ibin <= nbin; ++ibin) {
    if (axis->GetBinLowEdge(ibin+1) >= topbound) {
      integral += hist->GetBinContent(ibin) * (topbound - axis->GetBinLowEdge(ibin));
      break;
    }integral += hist->GetBinContent(ibin)*axis->GetBinWidth(ibin);
  }return integral*int_sign;
}

// Returns the coordinates for the center-of-mass of an input hist.
// Output is a std::vector<double> of coordinates.
// Length of output varies with histogram dimension.
// Flagging hasUniformBinning increases computation speed.
std::vector<double> getCM(TH1* hist,
                          int histflag = -1,
                          bool hasUniformBinning = false) {
  // Variable setup
  std::vector<double> cm;
  TAxis* ref_axis[3];
  double binVol = 1;
  
  if (histflag == -1) histflag = getHistType(hist);
  int hist_dim = histflag/10;
  switch (histflag) {
    case 31: case 32: case 33:
      ref_axis[2] = hist->GetZaxis();
      if (hasUniformBinning) binVol *= ref_axis[2]->GetBinWidth(1);
    case 21: case 22: case 23:
      ref_axis[1] = hist->GetYaxis();
      if (hasUniformBinning) binVol *= ref_axis[1]->GetBinWidth(1);
    case 11: case 12: case 13:
      ref_axis[0] = hist->GetXaxis();
      if (hasUniformBinning) binVol *= ref_axis[0]->GetBinWidth(1);
      break;
    default:
      cout << "Error in getCM: Input Histogram type not recognized." << endl;
      cout << "Input must be of type THnD, THnF or THnI for n = 1,2,3." << endl;
      return cm;
  }
  
  // Compute CM
  double sum;
  double temp;
  double weighted_sum[3] = {0,0,0};
  for (int iBinX = 1; iBinX <= ref_axis[0]->GetNbins(); ++iBinX) {
    if (hist_dim > 1) {
      for (int iBinY = 1; iBinY <= ref_axis[1]->GetNbins(); ++iBinY) {
        if (hist_dim > 2) {
          for (int iBinZ = 1; iBinZ <= ref_axis[2]->GetNbins(); ++iBinZ) {
            // TH3
            if (hasUniformBinning)
              temp = (hist->GetBinContent(iBinX, iBinY, iBinZ)
                      * binVol);
            else
              temp = (hist->GetBinContent(iBinX, iBinY, iBinZ)
                      * ref_axis[0]->GetBinWidth(iBinX)
                      * ref_axis[1]->GetBinWidth(iBinY)
                      * ref_axis[2]->GetBinWidth(iBinZ));
            sum += temp;
            weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
            weighted_sum[1] += temp * ref_axis[1]->GetBinCenter(iBinY);
            weighted_sum[2] += temp * ref_axis[2]->GetBinCenter(iBinZ);
          }
        } else {
          // TH2
          if (hasUniformBinning)
            temp = (hist->GetBinContent(iBinX, iBinY)
                    * binVol);
          else
            temp = (hist->GetBinContent(iBinX, iBinY)
                    * ref_axis[0]->GetBinWidth(iBinX)
                    * ref_axis[1]->GetBinWidth(iBinY));
          sum += temp;
          weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
          weighted_sum[1] += temp * ref_axis[1]->GetBinCenter(iBinY);
        }
      }
    } else {
      // TH1
      if (hasUniformBinning)
        temp = (hist->GetBinContent(iBinX)
                * binVol);
      else
        temp = (hist->GetBinContent(iBinX)
                * ref_axis[0]->GetBinWidth(iBinX));
      sum += temp;
      weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
    }
  }
  
  // Push to CM vector and return
  for (int iDim = 0; iDim < hist_dim; ++iDim) {
    cm.push_back(weighted_sum[iDim] / sum);
  }return cm;
}

std::vector<double> getCM_fast(TH2D* hist) {
  // Variable setup
  TAxis* ref_axis[2];
  ref_axis[0] = hist->GetXaxis();
  ref_axis[1] = hist->GetYaxis();
  double binVol = ref_axis[0]->GetBinWidth(1)*ref_axis[1]->GetBinWidth(1);
  
  // Compute CM
  double sum;
  double temp;
  double weighted_sum[2] = {0,0};
  for (int iBinX = 1; iBinX <= ref_axis[0]->GetNbins(); ++iBinX)
    for (int iBinY = 1; iBinY <= ref_axis[1]->GetNbins(); ++iBinY) {
      temp = (hist->GetBinContent(iBinX, iBinY) * binVol);
      
      sum += temp;
      weighted_sum[0] += temp * ref_axis[0]->GetBinCenter(iBinX);
      weighted_sum[1] += temp * ref_axis[1]->GetBinCenter(iBinY);
    }
  
  std::vector<double> cm = {weighted_sum[0]/sum, weighted_sum[1]/sum};
  return cm;
}

//========================================================================== Placeholder Methods

// Simple translate method for TH2D, assumes uniform binning.
// Only translates histogram within the resolution of bins,
// i.e. moving by the nearest bin division to input
TH2D* translateHist_simple(TH2D* hist,
                           double delta_x = 0,
                           double delta_y = 0,
                           bool suppress_warnings = false) {
  if (delta_x == 0 && delta_y == 0) {
    std::vector<double> cm = getCM(hist, 21, true);
    delta_x = -cm.at(0);
    delta_y = -cm.at(1);
    if (!suppress_warnings) cout << "Defaulting to CM at x = " << cm.at(0) << ", y = " << cm.at(1) << endl;
  }
  
  double width_x = hist->GetXaxis()->GetBinWidth(1);
  double width_y = hist->GetYaxis()->GetBinWidth(1);
  int nbins_x = hist->GetXaxis()->GetNbins();
  int nbins_y = hist->GetYaxis()->GetNbins();
  double histbounds[2][2] = {
    {hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinLowEdge(nbins_x+1)},
    {hist->GetYaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(nbins_y+1)} };
  if (TMath::Abs(delta_x) <= width_x/2 && TMath::Abs(delta_y) <= width_y/2) {
    if (!suppress_warnings) cout << "Shift too small to affect hist. Returned original hist." << endl;
    return hist;
  }
  
  TH2D* out_hist = static_cast<TH2D*>(hist->Clone());
  out_hist->Reset();
  
  double cx, cy, cfill;
  double lostdata = 0;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    for (int iy = 1; iy <= nbins_y; ++iy) {
      cfill = hist->GetBinContent(ix, iy);
      cx = hist->GetXaxis()->GetBinCenter(ix) + delta_x;
      cy = hist->GetYaxis()->GetBinCenter(iy) + delta_y;
      
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        lostdata += cfill*width_x * width_y;
        continue;
      }out_hist->Fill(cx, cy, cfill);
    }
  }
  
  // Report on lost data, store in public variable
  if (!suppress_warnings && lostdata > 0)
    cout << Form("Warning in translateHist_simple: %.5f%% of data lost to translation over edge.",
                 100.*lostdata / hist->Integral("width")) << endl;
  lostdata_translate = lostdata;
  
  return out_hist;
}


TH2D* rotateHist2D_simple(TH2D* hist, double theta, bool suppress_warnings = false) {
  double width_x = hist->GetXaxis()->GetBinWidth(1);
  double width_y = hist->GetYaxis()->GetBinWidth(1);
  int nbins_x = hist->GetXaxis()->GetNbins();
  int nbins_y = hist->GetYaxis()->GetNbins();
  double histbounds[2][2] = {
    {hist->GetXaxis()->GetBinLowEdge(1), hist->GetXaxis()->GetBinLowEdge(nbins_x+1)},
    {hist->GetYaxis()->GetBinLowEdge(1), hist->GetYaxis()->GetBinLowEdge(nbins_y+1)} };
  if (TMath::Abs(theta) <= TMath::Abs(TMath::ATan(width_y / (width_x * nbins_x)))) {
    if (!suppress_warnings) cout << "Warning in rotatehist2D_simple: Angle too small to affect hist; returned original hist." << endl;
    return hist;
  }
  
  TH2D* out_hist = static_cast<TH2D*>(hist->Clone());
  out_hist->Reset();
  
  double cx, cy, cr, ctheta, cfill;
  double lostdata = 0;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    for (int iy = 1; iy <= nbins_y; ++iy) {
      // Tally lost data by manually translating each item
      cfill = hist->GetBinContent(ix, iy);
      cx = hist->GetXaxis()->GetBinCenter(ix);
      cy = hist->GetYaxis()->GetBinCenter(iy);
      cr = TMath::Sqrt(cx*cx + cy*cy);
      ctheta = TMath::ATan2(cy, cx);
      cx = cr*TMath::Cos(ctheta + theta);
      cy = cr*TMath::Sin(ctheta + theta);
      
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        lostdata += cfill * width_x * width_y;
      }
      
      // Fill hist by forward rotating to get bin content
      // This is necessary to avoid 2 bins filling the same rotated bin
      cx = cr*TMath::Cos(ctheta - theta);
      cy = cr*TMath::Sin(ctheta - theta);
      if (cx <= histbounds[0][0] || cx >= histbounds[0][1] ||
          cy <= histbounds[1][0] || cy >= histbounds[1][1]) {
        continue;
      }
      cfill = hist->GetBinContent(hist->GetXaxis()->FindBin(cx),
                                  hist->GetYaxis()->FindBin(cy));
      out_hist->SetBinContent(ix, iy, cfill);
    }
  }
  
  // Report on lost data, store in public variable
  if (!suppress_warnings && lostdata > 0)
    cout << Form("Warning in rotateHist2D_simple: %.5f%% of data lost to rotation over edge.",
                 100.*lostdata / hist->Integral("width")) << endl;
  lostdata_rotate = lostdata;
  
  return out_hist;
}

// Returns gradient-norm map of input 2D histogram.
// Gradient is computed using finite difference approximation.
//
// CONSIDERATIONS :: Could be implemented more elegantly, or include
//    higher order corrections to finite difference upon request.
TH2D* gradientNorm(TH2D* hist) {
  TH2D* out_hist = (TH2D*) hist->Clone();
  out_hist->Reset();
  
  const int nbins_x = hist->GetXaxis()->GetNbins();
  const int nbins_y = hist->GetYaxis()->GetNbins();
  
  double ux, uy;
  double cx, cy;
  double lx, ly;
  double hx, hy;
  double dx, dy;
  for (int ix = 1; ix <= nbins_x; ++ix) {
    if (ix == 1) {
      lx = ix;
      ux = ix + 1;
      cx = ix;
      hx = hist->GetXaxis()->GetBinWidth(ix);
    } else if (ix == nbins_x) {
      lx = ix-1;
      cx = ix;
      ux = ix;
      hx = hist->GetXaxis()->GetBinWidth(ix);
    } else {
      lx = ix - 1;
      cx = ix;
      ux = ix + 1;
      hx = hist->GetXaxis()->GetBinWidth(ix-1)
      +hist->GetXaxis()->GetBinWidth(ix);
    }
    
    for (int iy = 1; iy <= nbins_y; ++iy) {
      if (iy == 1) {
        ly = iy;
        uy = iy + 1;
        cy = iy;
        hy = hist->GetYaxis()->GetBinWidth(iy);
      } else if (iy == nbins_y) {
        ly = iy-1;
        cy = iy;
        uy = iy;
        hy = hist->GetYaxis()->GetBinWidth(iy);
      } else {
        ly = iy - 1;
        cy = iy;
        uy = iy + 1;
        hy = hist->GetYaxis()->GetBinWidth(iy-1)
        +hist->GetYaxis()->GetBinWidth(iy);
      }
      
      // compute norm gradient with finite difference
      out_hist->Fill(hist->GetXaxis()->GetBinCenter(ix),
                     hist->GetXaxis()->GetBinCenter(iy),
                     TMath::Sqrt(TMath::Power((hist->GetBinContent(ux, cy) - hist->GetBinContent(lx, cy)) / hx, 2) +
                                 TMath::Power((hist->GetBinContent(cx, uy) - hist->GetBinContent(cx, ly)) / hy, 2)));
    }
  }return out_hist;
}
