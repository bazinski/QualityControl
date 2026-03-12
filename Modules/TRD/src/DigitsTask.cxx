// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>

#include "TRD/DigitsTask.h"
#include "TRD/TRDHelpers.h"
#include "DataFormatsTRD/Constants.h"
#include "DataFormatsTRD/Digit.h"
#include "DataFormatsTRD/HelperMethods.h"
#include "DataFormatsTRD/TriggerRecord.h"
#include "QualityControl/ObjectsManager.h"
#include "QualityControl/TaskInterface.h"
#include "QualityControl/QcInfoLogger.h"
#include "Common/Utils.h"
#include <Framework/InputRecord.h>
#include <Framework/InputRecordWalker.h>
#include <gsl/span>
#include <map>
#include <tuple>

using namespace o2::quality_control_modules::common;
using namespace o2::trd::constants;
using Helper = o2::trd::HelperMethods;

namespace o2::quality_control_modules::trd
{

void DigitsTask::initialize(o2::framework::InitContext& /*ctx*/)
{
  ILOG(Debug, Devel) << "initialize TRDDigitQcTask" << ENDM; // QcInfoLogger is used. FairMQ logs will go to there as well.

  // this is how to get access to custom parameters defined in the config file at qc.tasks.<task_name>.taskParameters
  mDoClusterize = getFromConfig<bool>(mCustomParameters, "doClusterize", false);
  mPulseHeightPeakRegion.first = getFromConfig<float>(mCustomParameters, "peakRegionStart", 0.f);
  mPulseHeightPeakRegion.second = getFromConfig<float>(mCustomParameters, "peakRegionEnd", 5.f);
  mPulseHeightThreshold = getFromConfig<unsigned int>(mCustomParameters, "phThreshold", 400u);
  mChambersToIgnore = getFromConfig<std::string>(mCustomParameters, "ignoreChambers", "16_3_0");
  mClsCutoff = getFromConfig<int>(mCustomParameters, "clsCutoff", 1000);
  mAdcBaseline = getFromConfig<int>(mCustomParameters, "adcBaseline", 10);

  buildChamberIgnoreBP();
  buildHistograms();
}

void DigitsTask::buildHistograms()
{
  constexpr int nLogBins = 100;
  float xBins[nLogBins + 1];
  float xBinLogMin = 0.f;
  float xBinLogMax = 8.f;
  float logBinWidth = (xBinLogMax - xBinLogMin) / nLogBins;
  for (int iBin = 0; iBin <= nLogBins; ++iBin) {
    xBins[iBin] = TMath::Power(10, xBinLogMin + iBin * logBinWidth);
  }
  mDigitsPerEvent.reset(new TH1F("digitsperevent", "Number of digits per Event", nLogBins, xBins));
  getObjectsManager()->startPublishing(mDigitsPerEvent.get());

  mDigitHCID.reset(new TH1F("digithcid", "Digit distribution over Halfchambers", 1080, -0.5, 1079.5));
  getObjectsManager()->startPublishing(mDigitHCID.get());

  mDigitsSizevsTrackletSize.reset(new TH1F("digitsvstracklets", "Digit count vs. Tracklet count; N_Digits / N_Tracklets;counts", 100, 0, 100));
  getObjectsManager()->startPublishing(mDigitsSizevsTrackletSize.get());

  mADCvalue.reset(new TH1F("ADCvalue", "ADC values for all digits;ADC value;Counts", 1024, -0.5, 1023.5));
  getObjectsManager()->startPublishing(mADCvalue.get());

  for (Int_t iSM = 0; iSM < NSECTOR; ++iSM) {
    mHCMCM[iSM].reset(new TH2F(Form("DigitsPerMCM_%i", iSM), Form("Digits per MCM in sector %i;pad row; MCM column", iSM), 76, -0.5, 75.5, 8, -0.5, 7.5));
    getObjectsManager()->startPublishing(mHCMCM[iSM].get());
    getObjectsManager()->setDefaultDrawOptions(mHCMCM[iSM]->GetName(), "COLZ");
  }

  if (mDoClusterize) {
    mClsNTb.reset(new TH1F("Cluster/ClsNTb", "ClsNTb", 30, -0.5, 29.5));
    mClsNTb->SetTitle("Clusters per time bin;Number of clusters;Timebin");
    getObjectsManager()->startPublishing(mClsNTb.get());
    mClsChargeTb.reset(new TH1F("Cluster/ClsChargeTb", "ClsChargeTb;", 30, -0.5, 29.5));
    getObjectsManager()->startPublishing(mClsChargeTb.get());
    mClsAmp.reset(new TH1F("Cluster/ClsAmp", "Amplitude of clusters;Amplitdu(ADC);Counts", 200, -0.5, 1999.5));
    getObjectsManager()->startPublishing(mClsAmp.get());
    mNCls.reset(new TH1F("Cluster/NCls", "Total number of clusters per sector", 18, -0.5, 17.5));
    mNCls->SetTitle("Total number of clusters per sector;Sector;Counts");
    getObjectsManager()->startPublishing(mNCls.get());
    mClsTb.reset(new TH2F("Cluster/ChargeTB", "Cluster charge;time bin;cluster charge", 30, -0.5, 29.5, 200, 0, 2000));
    getObjectsManager()->startPublishing(mClsTb.get());
  }

  mPulseHeight.reset(new TH1F("PulseHeight1D", "PH spectrum 1D;time bin;ADC sum", 30, -0.5, 29.5));
  drawLinesOnPulseHeight(mPulseHeight.get());
  getObjectsManager()->startPublishing(mPulseHeight.get());
  mPulseHeight.get()->GetYaxis()->SetTickSize(0.01);

  mPulseHeightHD.reset(new TH1F("PulseHeightHD1D", "PH spectrum 1D (with phase);time bin;ADC sum", 120, -0.5, 29.5));
  drawLinesOnPulseHeight(mPulseHeightHD.get());
  getObjectsManager()->startPublishing(mPulseHeightHD.get());
  mPulseHeightHD.get()->GetYaxis()->SetTickSize(0.01);
  
  mPulseHeightHDcorrected.reset(new TH1F("PulseHeightHD1Dcorrected", "PH spectrum 1D (with phase corrected);time bin;ADC sum", 120, -0.5, 29.5));
  drawLinesOnPulseHeight(mPulseHeightHDcorrected.get());
  getObjectsManager()->startPublishing(mPulseHeightHDcorrected.get());
  mPulseHeightHDcorrected.get()->GetYaxis()->SetTickSize(0.01);

  mTotalPulseHeight2D.reset(new TH2F("PulseHeight2D", "PH spectrum 2D;time bin;ADC sum", 30, 0., 30., 200, 0., 200.));
  getObjectsManager()->startPublishing(mTotalPulseHeight2D.get());
  getObjectsManager()->setDefaultDrawOptions(mTotalPulseHeight2D->GetName(), "COLZ");

  mTotalPulseHeightHD2D.reset(new TH2F("PulseHeightHD2D", "PH spectrum 2D (with phase);time bin;ADC sum", 120, 0., 30., 200, 0., 200.));
  getObjectsManager()->startPublishing(mTotalPulseHeightHD2D.get());
  getObjectsManager()->setDefaultDrawOptions(mTotalPulseHeightHD2D->GetName(), "COLZ");

  mTotalPulseHeightHD2Dcorrected.reset(new TH2F("PulseHeightHD2Dcorrected", "PH spectrum 2D (with phase corrected);time bin;ADC sum", 120, 0., 30., 200, 0., 200.));
  getObjectsManager()->startPublishing(mTotalPulseHeightHD2Dcorrected.get());
  getObjectsManager()->setDefaultDrawOptions(mTotalPulseHeightHD2Dcorrected->GetName(), "COLZ");

  mPulseHeightpro.reset(new TProfile("PulseHeightProfile", "PH spectrum for all chambers combined;time bin;ADC sum", 30, -0.5, 29.5));
  mPulseHeightpro.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightpro.get());
  
  mPulseHeightHDpro.reset(new TProfile("PulseHeightHDProfile", "PH spectrum for all chambers combined (with phase);time bin;ADC sum", 120, -0.5, 29.5));
  mPulseHeightHDpro.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightpro.get());
 
  mPulseHeightHDprocorrected.reset(new TProfile("PulseHeightHDProfilecorrected", "PH spectrum for all chambers combined (with phase corrected);time bin;ADC sum", 120, -0.5, 29.5));
  mPulseHeightHDprocorrected.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightprocorrected.get());

  mPulseHeightperchamber.reset(new TProfile2D("PulseHeightPerChamber", "PH spectrum for all chambers;time bin;chamber", 30, -0.5, 29.5, 540, -0.5, 539.5));
  mPulseHeightperchamber.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightperchamber.get());
  getObjectsManager()->setDefaultDrawOptions(mPulseHeightperchamber.get()->GetName(), "colz");
 
  mPulseHeightHDperchamber.reset(new TProfile2D("PulseHeightHDPerChamber", "PH spectrum for all chambers (with phase);time bin;chamber", 120, -0.5, 29.5, 540, -0.5, 539.5));
  mPulseHeightHDperchamber.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightHDperchamber.get());
  getObjectsManager()->setDefaultDrawOptions(mPulseHeightHDperchamber.get()->GetName(), "colz");
  
  mPulseHeightHDperchambercorrected.reset(new TProfile2D("PulseHeightHDPerChambercorrected", "PH spectrum for all chambers (with phase corrected);time bin;chamber", 120, -0.5, 29.5, 540, -0.5, 539.5));
  mPulseHeightHDperchambercorrected.get()->Sumw2();
  getObjectsManager()->startPublishing(mPulseHeightHDperchambercorrected.get());
  getObjectsManager()->setDefaultDrawOptions(mPulseHeightHDperchambercorrected.get()->GetName(), "colz");

  for (int iSec = 0; iSec < NSECTOR; ++iSec) {
    mPulseHeight2DperSM[iSec].reset(new TH1F(Form("PulseHeight_%i", iSec), Form("PH spectrum for sector %i;time bin;ADC sum count", iSec), 30, -0.5, 29.5));
    getObjectsManager()->startPublishing(mPulseHeight2DperSM[iSec].get());
    mPulseHeightHD2DperSM[iSec].reset(new TH1F(Form("PulseHeightHD_%i", iSec), Form("PH spectrum for sector %i (with phase);time bin;ADC sum count", iSec), 120, -0.5, 29.5));
    getObjectsManager()->startPublishing(mPulseHeightHD2DperSM[iSec].get());
    mPulseHeightHD2DperSMcorrected[iSec].reset(new TH1F(Form("PulseHeightHD_%i", iSec), Form("PH spectrum for sector %i (with phase);time bin;ADC sum count", iSec), 120, -0.5, 29.5));
    getObjectsManager()->startPublishing(mPulseHeightHD2DperSMcorrected[iSec].get());
  }

  // Build digits layers
  int unitsPerSection = NCOLUMN;
  for (int iLayer = 0; iLayer < NLAYER; ++iLayer) {
    mLayers[iLayer].reset(new TH2F(Form("DigitsPerLayer_%i", iLayer), Form("Digit count per pad in layer %i;glb pad row;glb pad col", iLayer),
                                   76, -0.5, 75.5, unitsPerSection * NSECTOR, -0.5, unitsPerSection * NSECTOR - 0.5));
    mLayers[iLayer]->SetStats(0);
    TRDHelpers::addChamberGridToHistogram(mLayers[iLayer], unitsPerSection);
    getObjectsManager()->startPublishing(mLayers[iLayer].get());
    getObjectsManager()->setDefaultDrawOptions(mLayers[iLayer]->GetName(), "COLZ");
    getObjectsManager()->setDisplayHint(mLayers[iLayer].get(), "logz");
  }
}

void DigitsTask::monitorData(o2::framework::ProcessingContext& ctx)
{
  // Load CCDB objects (needs to be done only once)
  if (!mNoiseMap) {
    auto ptr = ctx.inputs().get<o2::trd::NoiseStatusMCM*>("noiseMap");
    mNoiseMap = ptr.get();
  }

  auto ptr = ctx.inputs().get<std::array<int, MAXCHAMBER>*>("fedChamberStatus");

  // fill histograms
  auto digits = ctx.inputs().get<gsl::span<o2::trd::Digit>>("digits");
  auto triggerrecords = ctx.inputs().get<gsl::span<o2::trd::TriggerRecord>>("triggers");
  for (auto& trigger : triggerrecords) {
    auto nDigits = trigger.getNumberOfDigits();
    auto iDigitFirst = trigger.getFirstDigit();
    int iDigitLast = iDigitFirst + nDigits;
    if (nDigits == 0) {
      continue;
    }
    mDigitsPerEvent->Fill(nDigits);
    if (trigger.getNumberOfTracklets() > 0) {
      mDigitsSizevsTrackletSize->Fill(nDigits / trigger.getNumberOfTracklets());
    }

    for (int iDigit = iDigitFirst; iDigit < iDigitLast; ++iDigit) {
      const auto& digit = digits[iDigit];
      if (digit.isSharedDigit()) {
        continue;
      }

      int detector = digit.getDetector();
      if (detector >= MAXCHAMBER) {
        // stupid fix to avoid crashes below on corrupted digits recorded in 2022
        // observed e.g. in run 523308
        continue;
      }
      int sector = Helper::getSector(detector);
      int layer = Helper::getLayer(detector);
      int stack = Helper::getStack(detector);
      int row = digit.getPadRow();
      int col = digit.getPadCol();
      int rowGlb = FIRSTROW[stack] + row;
      int colGlb = col + sector * NCOLUMN;
      mLayers[layer]->Fill(rowGlb, colGlb);
      mHCMCM[sector]->Fill(rowGlb, digit.getMCMCol());
      mDigitHCID->Fill(digit.getHCId());
      for (auto adc : digit.getADC()) {
        mADCvalue->Fill(adc);
      }
      if (iDigit == iDigitFirst || iDigit == iDigitLast - 1) {
        // above histograms are filled for all digits, now we are looking for clusters
        continue;
      }

      const auto* digitLeft = &digits[iDigit + 1];
      const auto* digitRight = &digits[iDigit - 1];
      // due to shared digits the neighbouring element might still be in the same pad column
      if (digitLeft->getPadCol() == col && iDigit < iDigitLast - 2) {
        digitLeft = &(digits[iDigit + 2]);
      }
      if (digitRight->getPadCol() == col && iDigit > iDigitFirst + 2) {
        digitRight = &(digits[iDigit - 2]);
      }
      if (!digitLeft->isNeighbour(digit) || !digitRight->isNeighbour(digit)) {
        continue;
      }

      if (mDoClusterize) {
        for (int tb = 1; tb < TIMEBINS - 1; ++tb) {
          auto adc = digit.getADC()[tb];
          if (adc > digitLeft->getADC()[tb] && adc > digitRight->getADC()[tb]) {
            // local maximum
            int valueLU = digitLeft->getADC()[tb - 1] > mAdcBaseline ? digitLeft->getADC()[tb - 1] - mAdcBaseline : 0;
            int valueRU = digitRight->getADC()[tb - 1] > mAdcBaseline ? digitRight->getADC()[tb - 1] - mAdcBaseline : 0;
            int valueU = digit.getADC()[tb - 1] - mAdcBaseline; // why not protect against negative values here?

            int valueLD = digitLeft->getADC()[tb + 1] > mAdcBaseline ? digitLeft->getADC()[tb + 1] - mAdcBaseline : 0;
            int valueRD = digitRight->getADC()[tb + 1] > mAdcBaseline ? digitRight->getADC()[tb + 1] - mAdcBaseline : 0;
            int valueD = digit.getADC()[tb + 1] - mAdcBaseline; // same question as above

            int valueL = digitLeft->getADC()[tb] > mAdcBaseline ? digitLeft->getADC()[tb] - mAdcBaseline : 0;
            int valueR = digitRight->getADC()[tb] > mAdcBaseline ? digitRight->getADC()[tb] - mAdcBaseline : 0;

            int sum = adc + valueL + valueR;
            int sumU = valueU + valueLU + valueRU;
            int sumD = valueD + valueLD + valueRD;

            if (sumU < mAdcBaseline || sumD < mAdcBaseline) {
              continue;
            }
            if (TMath::Abs(1. * sum / sumU - 1) < 0.01) {
              continue;
            }
            if (TMath::Abs(1. * sum / sumD - 1) < 0.01) {
              continue;
            }
            if (TMath::Abs(1. * sumU / sumD - 1) < 0.01) {
              continue;
            }
            mNCls->Fill(sector);
            mClsTb->Fill(tb, sum);
            mClsAmp->Fill(sum);
            if (sum > mAdcBaseline && sum < mClsCutoff) {
              mClsChargeTb->Fill(tb, sum);
              mClsNTb->Fill(tb);
            }
          } // local maximum
        }   // loop over time bins
      }

      auto sumLeft = digitLeft->getADCsum();
      auto sumCenter = digit.getADCsum();
      auto sumRight = digitRight->getADCsum();
      if (!mChambersToIgnoreBP.test(detector)) {
        if (sumCenter > sumLeft && sumCenter > sumRight) {
          auto lowestSum = (sumLeft > sumRight) ? sumRight : sumLeft;
          if (lowestSum > mPulseHeightThreshold) {
            for (int tb = 0; tb < TIMEBINS; tb++) {
              int phVal = (digit.getADC()[tb] + digitLeft->getADC()[tb] + digitRight->getADC()[tb]);
              mPulseHeight->Fill(tb, phVal);
              mTotalPulseHeight2D->Fill(tb, phVal);
              mPulseHeight2DperSM[sector]->Fill(tb, phVal);
              mPulseHeightpro->Fill(tb, phVal);
              mPulseHeightperchamber->Fill(tb, detector, phVal);
            } // loop over time bins
            for (int tb = 0; tb < TIMEBINS; tb++) {
              int phValHD = digit.getADC()[tb] + digitLeft->getADC()[tb] + digitRight->getADC()[tb];
              float tbphase=(tb<<2)+ digit.getPreTrigPhase();
              //LOGP(info,"Pretrig phase is : {} tb {} tbphase {}", digit.getPreTrigPhase(), tb, tbphase);
              tbphase= tbphase/4; // get the tb back into 0-30, its currently in 0-120;
              mPulseHeightHD->Fill(tbphase, phValHD);
              mTotalPulseHeightHD2D->Fill(tbphase, phValHD);
              mPulseHeightHD2DperSM[sector]->Fill(tbphase, phValHD);
              mPulseHeightHDpro->Fill(tbphase, phValHD);
              mPulseHeightHDperchamber->Fill(tbphase, detector, phValHD);
            } // loop over time bins with phase
          }   // lower ADC sum above threshold
        }     // local ADC maximum
      }       // chamber OK
    }         // for loop over digits
  }           // loop over triggers
}

void DigitsTask::drawLinesOnPulseHeight(TH1F* h)
{
  TLine* lmin = new TLine(mPulseHeightPeakRegion.first, 0, mPulseHeightPeakRegion.first, 1e9);
  TLine* lmax = new TLine(mPulseHeightPeakRegion.second, 0, mPulseHeightPeakRegion.second, 1e9);
  lmin->SetLineStyle(2);
  lmax->SetLineStyle(2);
  lmin->SetLineColor(kRed);
  lmax->SetLineColor(kRed);
  h->GetListOfFunctions()->Add(lmin);
  h->GetListOfFunctions()->Add(lmax);
}

void DigitsTask::drawLinesOnPulseHeightHD(TH1F* h)
{
  //the same as above but hopefully lines extending from top to bottom of y-axis
  double ymin = h->GetMinimum();
  double ymax = h->GetMaximum();
  TLine* lstart = new TLine(mPulseHeightPeakRegion.first, ymin, mPulseHeightPeakRegion.first, ymax);
  TLine* lend = new TLine(mPulseHeightPeakRegion.second, ymin, mPulseHeightPeakRegion.second, ymax);
  lstart->SetLineStyle(2);
  lend->SetLineStyle(2);
  lstart->SetLineColor(kRed);
  lend->SetLineColor(kRed);
  h->GetListOfFunctions()->Add(lstart);
  h->GetListOfFunctions()->Add(lend);
}

void DigitsTask::buildChamberIgnoreBP()
{
  mChambersToIgnoreBP.reset();
  // Vector of string to save tokens
  std::vector<std::string> tokens;
  // stringstream class check1
  std::stringstream check1(mChambersToIgnore);
  std::string intermediate;
  // Tokenizing w.r.t. space ','
  while (getline(check1, intermediate, ',')) {
    tokens.push_back(intermediate);
  }
  // Printing the token vector
  for (auto& token : tokens) {
    // token now holds something like 16_3_0
    std::vector<std::string> parts;
    std::stringstream indexcheck(token);
    std::string tokenpart;
    // Tokenizing w.r.t. space ','
    while (getline(indexcheck, tokenpart, '_')) {
      parts.push_back(tokenpart);
    }
    // now flip the bit related to the sector:stack:layer stored in parts.
    int sector = std::stoi(parts[0]);
    int stack = std::stoi(parts[1]);
    int layer = std::stoi(parts[2]);
    mChambersToIgnoreBP.set(Helper::getDetector(sector, stack, layer));
  }
}

void DigitsTask::startOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "startOfActivity" << ENDM;
} // set stats/stacs

void DigitsTask::startOfCycle()
{
  ILOG(Debug, Devel) << "startOfCycle" << ENDM;
}

void DigitsTask::endOfCycle()
{
  ILOG(Debug, Devel) << "endOfCycle" << ENDM;
  /*
  std::shared_ptr<TH1F> mPulseHeightHD = nullptr; 
  std::shared_ptr<TH1F> mPulseHeightHDcorrected = nullptr; **
  std::shared_ptr<TH2F> mTotalPulseHeightHD2Dcorrected = nullptr; **
  std::array<std::shared_ptr<TH1F>, o2::trd::constants::NSECTOR> mPulseHeightHD2DperSM;
  std::shared_ptr<TProfile> mPulseHeightHDpro = nullptr; 
  std::shared_ptr<TProfile> mPulseHeightHDprocorrected = nullptr; **
 
 std::shared_ptr<TProfile2D> mPulseHeightHDperchambercorrected = nullptr;**
*/
  //fix the corrected HD spectra for the prevelance (efficiency) of the relative phases.
  std::array<float,4> phasetotals{};
  std::array<float,4> phasescaling{};
  for(int i=0;i<120;++i){
    phasetotals[i/4]+= mPulseHeightHD->GetBinContent(i);
  }
  const auto maxval =std::max_element(phasetotals.begin(),phasetotals.end());
  if(*maxval>0){
    for(int i=0;i<4;++i){
      phasescaling[i]=phasetotals[i]/ (*maxval);
    }
  } 

  for(int i=0;i<120;++i){
    if(phasetotals[i/4]>0){
      mPulseHeightHDcorrected->SetBinContent(i, mPulseHeightHD->GetBinContent(i)/phasescaling[i/4]);
      for(int j=0;j<??;++j){
        mTotalPulseHeightHD2Dcorrected->SetBinContent(i, j, mTotalPulseHeightHD2Dcorrected->GetBinContent(i,j)/phasescaling[i/4]);
      }
      //weighted graphs
      mPulseHeightHDprocorrected->SetBinContent(i, mPulseHeightHDpro->GetBinContent(i)/phasescaling[i/4]);
      mPulseHeightHDprocorrected->SetBinError(i, mPulseHeightHDpro->GetBinError(i)/phasescaling[i/4]);
    }
  }
  // now for the per chamber version:
  for(int chamber=0;chamber<540;++chamber){
    std::array<float,4> phasetotals{};
    std::array<float,4> phasescaling{};
    for(int i=0;i<120;++i){
      phasetotals[i/4]+= mPulseHeightHDperchamber->GetBinContent(i,chamber);
    }
    const auto maxval =std::max_element(phasetotals.begin(),phasetotals.end());
    if(*maxval>0){
      for(int i=0;i<4;++i){
        phasescaling[i]=phasetotals[i]/ (*maxval);
      }
      
    } 
    for(int i=0;i<120;++i){
      if(phasetotals[i/4]>0){
        mPulseHeightHDperchambercorrected->SetBinContent(i, chamber, mPulseHeightHDperchamber->GetBinContent(i)/phasescaling[i/4]);
        mPulseHeightHDperchambercorrected->SetBinError(i, chamber,mPulseHeightHDperchamber->GetBinError(i)/phasescaling[i/4]);
      }
    }
  }


}

void DigitsTask::endOfActivity(const Activity& /*activity*/)
{
  ILOG(Debug, Devel) << "endOfActivity" << ENDM;
}

void DigitsTask::finaliseCCDB(o2::framework::ConcreteDataMatcher& matcher, void* obj)
{
  if (matcher == o2::framework::ConcreteDataMatcher("TRD", "FCHSTATUS", 0)) {
    // LB: no half chamber distribution map for Digits, pass it as null pointer
    TRDHelpers::drawChamberStatusOnHistograms(static_cast<std::array<int, MAXCHAMBER>*>(obj), nullptr, mLayers, NCOLUMN);
  }
}

void DigitsTask::reset()
{
  // clean all the monitor objects here
  ILOG(Debug, Devel) << "Resetting the histogram" << ENDM;
  mDigitsPerEvent->Reset();
  mDigitHCID->Reset();
  mADCvalue->Reset();
  mDigitsSizevsTrackletSize->Reset();
  mTotalPulseHeight2D->Reset();
  mPulseHeight->Reset();
  mPulseHeightpro->Reset();
  mPulseHeightperchamber->Reset();
  mTotalPulseHeightHD2D->Reset();
  mPulseHeightHD->Reset();
  mPulseHeightHDpro->Reset();
  mPulseHeightHDperchamber->Reset();
  mPulseHeightHDcorrected->Reset();
  mPulseHeightHDprocorrected->Reset();
  mPulseHeightHDperchambercorrected->Reset();
  for (auto& h : mPulseHeight2DperSM) {
    h->Reset();
  }
  for (auto& h : mPulseHeightHD2DperSM) {
    h->Reset();
  }
  for (auto& h : mPulseHeightHD2DperSMcorrected) {
    h->Reset();
  }
  for (auto& h : mHCMCM) {
    h->Reset();
  }
  for (auto& h : mLayers) {
    h->Reset();
  }
  if (mDoClusterize) {
    mNCls->Reset();
    mClsTb->Reset();
    mClsChargeTb->Reset();
    mClsNTb->Reset();
    mClsAmp->Reset();
  }
}
} // namespace o2::quality_control_modules::trd
