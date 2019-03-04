/*
 *  JHistos.cxx
 *  
 *
 *  Created by Sami Rasanen on 2/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "JHistos.h"
#include "TMath.h"

JHistos::JHistos(){

    // If you change this, remember to change the numbers in JHistos.h also as that was left alone.
    int const static nRuns = 3; //0: detector, 1: detector without corrections, 2: detector with corrections
    TString sType[nRuns] = {"True", "DetWoCorr","DetWCorr"};
	
    // pT -histograms, taken from JCORRAN version on 11th February 2013
    int NBINS=150;
    double LogBinsX[NBINS+1], LimL=0.1, LimH=100;
    double logBW = (TMath::Log(LimH)-TMath::Log(LimL))/NBINS;
    for(int ij=0;ij<=NBINS;ij++) LogBinsX[ij]=LimL*exp(ij*logBW);
    
    // xT -histograms, use also here logarithmic bins
    int binsQ = 400, binsQ4 = 800, binsQ6 = 1000;
    double LogQBinsX[binsQ+1], LogQ2BinsX[binsQ+1], LogQ4BinsX[binsQ4+1], LogQ6BinsX[binsQ6+1];
    double LimQL=0.1, LimQH=2.e1, LimQ2H=2.e2, LimQ4H=2.e4, LimQ6H=2.e6;
    double logQBW  = (TMath::Log(LimQH)-TMath::Log(LimQL))/binsQ;
    double logQ2BW = (TMath::Log(LimQ2H)-TMath::Log(LimQL))/binsQ;
    double logQ4BW = (TMath::Log(LimQ4H)-TMath::Log(LimQL))/binsQ4;
    double logQ6BW = (TMath::Log(LimQ6H)-TMath::Log(LimQL))/binsQ6;
    for(int ij=0;ij<=binsQ;ij++){
        LogQBinsX[ij] =LimQL*exp(ij*logQBW);
        LogQ2BinsX[ij]=LimQL*exp(ij*logQ2BW);
    }
    for(int ij=0;ij<=binsQ4;ij++) LogQ4BinsX[ij]=LimQL*exp(ij*logQ4BW);
    for(int ij=0;ij<=binsQ6;ij++) LogQ6BinsX[ij]=LimQL*exp(ij*logQ6BW);
    
    for(int iType=0; iType<nRuns; iType++) {
        hPt[iType] = new TH1D(Form("hPtT%02i",iType),Form("pT - %s",sType[iType].Data()), NBINS, LogBinsX); hPt[iType]->Sumw2();
        hPhi[iType] = new TH1D(Form("hPhiT%02i",iType),Form("phi - %s",sType[iType].Data()), 129,-3.2,3.2); hPhi[iType]->Sumw2();
        hEta[iType] = new TH1D(Form("hEtaT%02i",iType),Form("eta - %s",sType[iType].Data()), 81,-0.8,0.8); hEta[iType]->Sumw2();
        hPhiEta[iType] = new TH2D(Form("hPhiEtaT%02i",iType),Form("phi-eta - %s",sType[iType].Data()),129,-3.2,3.2,81,-0.8,0.8); hPhiEta[iType]->Sumw2();
        hMultiplicity[iType] = new TH1D(Form("hMultiplicityT%02i",iType),Form("Multiplicity - %s",sType[iType].Data()), 125, 0.0, 2500.); hMultiplicity[iType]->Sumw2();
        hSqrtSumWeights[iType] = new TH1D(Form("hSqrtSumWeightsT%02i",iType),Form("sqrt of sum of weights squares - %s",sType[iType].Data()), 240, 0.0, 60.); hSqrtSumWeights[iType]->Sumw2();
        hSqrtSumWeightsA[iType] = new TH1D(Form("hSqrtSumWeightsAT%02i",iType),Form("sqrt of sum of weights squares in subevent A - %s",sType[iType].Data()), 240, 0.0, 60.); hSqrtSumWeightsA[iType]->Sumw2();
        hSqrtSumWeightsB[iType] = new TH1D(Form("hSqrtSumWeightsBT%02i",iType),Form("sqrt of sum of weights squares in subevent B - %s",sType[iType].Data()), 240, 0.0, 60.); hSqrtSumWeightsB[iType]->Sumw2();
        // integrated
        for(int i = 0; i < 5; i++){
            hQx[iType][i] = new TH1D(Form("hQxT%02iH%02i",iType,i+1),Form("hQx%s%02i",sType[iType].Data(),i+1),401,-25.0,25.0); hQx[iType][i]->Sumw2();
            hQy[iType][i] = new TH1D(Form("hQyT%02iH%02i",iType,i+1),Form("hQy%s%02i",sType[iType].Data(),i+1),401,-25.0,25.0); hQy[iType][i]->Sumw2();
            hQ[iType][i]  = new TH1D(Form("hQT%02iH%02i",iType,i+1),Form("hQ%s%02i",sType[iType].Data(),i+1),binsQ, LogQBinsX); hQ[iType][i]->Sumw2();
            hQ2[iType][i] = new TH1D(Form("hQ2T%02iH%02i",iType,i+1),Form("hQ2%s%02i",sType[iType].Data(),i+1),binsQ, LogQ2BinsX); hQ2[iType][i]->Sumw2();
            hQ4[iType][i] = new TH1D(Form("hQ4T%02iH%02i",iType,i+1),Form("hQ4%s%02i",sType[iType].Data(),i+1),binsQ4, LogQ4BinsX); hQ4[iType][i]->Sumw2();
            hQ6[iType][i] = new TH1D(Form("hQ6T%02iH%02i",iType,i+1),Form("hQ6%s%02i",sType[iType].Data(),i+1),binsQ6, LogQ6BinsX); hQ6[iType][i]->Sumw2();
            hvObs[iType][i] = new TH1D(Form("hvObsT%02iH%02i",iType,i+1),Form("hvObs%s%02i",sType[iType].Data(),i+1),100, -1.0, 1.0); hvObs[iType][i]->Sumw2();
            hrSub[iType][i] = new TH1D(Form("hrSubT%02iH%02i",iType,i+1),Form("hrSub%s%02i",sType[iType].Data(),i+1),100, -1.0, 1.0); hrSub[iType][i]->Sumw2();
            for(int j=0;j<2;j++) {
                hPsiAB[iType][i][j] = new TH2D(Form("hPsiABT%02iH%02iE%02i",iType,i+1,j),Form("hPsiAB%s%02i%02i",sType[iType].Data(),i+1,j),100, -1.1*TMath::Pi(), 1.1*TMath::Pi(), 100, -1.1*TMath::Pi(), 1.1*TMath::Pi()); hPsiAB[iType][i][j]->Sumw2();
                hPsiDiff[iType][i][j] = new TH1D(Form("hPsiDiffT%02iH%02iE%02i",iType,i+1,j),Form("hPsiDiff%s%02i%02i",sType[iType].Data(),i+1,j),200, -4.0*TMath::Pi(), 4.0*TMath::Pi()); hPsiDiff[iType][i][j]->Sumw2();
                hPsiDiffN[iType][i][j] = new TH1D(Form("hPsiDiffNT%02iH%02iE%02i",iType,i+1,j),Form("hPsiDiffN%s%02i%02i",sType[iType].Data(),i+1,j),200, -4.0*TMath::Pi(), 4.0*TMath::Pi()); hPsiDiffN[iType][i][j]->Sumw2();
            }
            hSPnominator[iType][i] = new TH1D(Form("hSPnominatorVT%02iH%02i",iType,i+1),Form("hSPnominatorV%s%02i",sType[iType].Data(),i+1),binsQ, LogQ2BinsX); hSPnominator[iType][i]->Sumw2();
            hSPdenominator[iType][i] = new TH1D(Form("hSPdenominatorVT%02iH%02i",iType,i+1),Form("hSPdenominatorV%s%02i",sType[iType].Data(),i+1),binsQ, LogQ2BinsX); hSPdenominator[iType][i]->Sumw2();
            hEPnominator[iType][i] = new TH1D(Form("hEPnominatorVT%02iH%02i",iType,i+1),Form("hEPnominatorV%s%02i",sType[iType].Data(),i+1),400,0.0,20.0); hEPnominator[iType][i]->Sumw2();
            hEPdenominator[iType][i] = new TH1D(Form("hEPdenominatorVT%02iH%02i",iType,i+1),Form("hEPdenominatorV%s%02i",sType[iType].Data(),i+1),100,0.0,1.0); hEPdenominator[iType][i]->Sumw2();
        }
    }
}
