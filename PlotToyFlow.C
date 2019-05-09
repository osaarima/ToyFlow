//
//  LoptToyFlow.c
//  
//
//  Created by Räsänen Sami on 13.4.2017.
//
//

#include <stdio.h>

double CalculateV2Using4Cumul( double Q4, double Q2, double weight );
double CalculateV2Using6Cumul( double Q6, double Q4, double Q2, double weight );
double CalculateV2UsingEPorSP( double up, double down, double weight );
double CalculateV2ErrorUsing4Cumul( double Q4, double Q2, double eQ4, double eQ2, double weight, double weightError );
double CalculateV2ErrorUsing6Cumul( double Q6, double Q4, double Q2, double eQ6, double eQ4, double eQ2, double weight, double weightError );
double CalculateV2ErrorUsingEPorSP( double up, double down, double eUp, double eDown, double weight, double weightError );
double CalculateResolutionKOne( double Xi );
double fFitMeanQ( double *x, double *p );
void checkUnderOverFlow( TH1 *h );
double AcceptanceFunc(double *x, double *p);
bool isInAcc(double phi, double detMax, double detMin);
double v2PtDependence(double *x, double *p);
double v2PtDependenceFun(double pT, double pTmax, double alpha);

void PlotToyFlow(int iType=2, bool bDrawNegRHisto = false, bool bUseWeightning = false, bool bCheckAllHistos = false, bool bPrintGraphs = false ){
    int nType = 3;
    
    int iDrawNegRHisto=0;
    if(bDrawNegRHisto) iDrawNegRHisto=1;
    //TString sFileName = "toyFlow_weight_randomPsi_ptDep_dNdeta-1000_nEvents-1000-TestWPt.root";
    //TString sFileName = "toyFlow_weight_randomPsi_ptDep_dNdeta-1000_nEvents-1000-Test.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-1000-Test.root";
    //TString sFileName = "toyFlow_weight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-TestWPt.root";
    //TString sFileName = "toyFlow_weight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-Test.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-Test.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-ptbinkokeilu.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-300-overflowTest.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-100-newPtDepTest.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-1000-PtDep_noWeight_1100.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-100-ttestt.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-EPtest.root";
    //TString sFileName = "toyFlow_weight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-asdasd.root";
    TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000-asdasd.root";
    //run1
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-100000-run1.root";
    //TString sFileName = "toyFlow_weight_randomPsi_noptDep_dNdeta-1000_nEvents-100000-run1.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-5000-run1.root";
    //TString sFileName = "toyFlow_weight_randomPsi_ptDep_dNdeta-1000_nEvents-5000-run1.root";
    //puck run
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000000-overflowTest.root";
    //TString sFileName = "toyFlow_weight_randomPsi_noptDep_dNdeta-1000_nEvents-1000000-overflowTest.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-100000-overflowTest.root";
    //TString sFileName = "toyFlow_weight_randomPsi_ptDep_dNdeta-1000_nEvents-100000-overflowTest.root";
    //puck run 2
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-1000.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-100000.root";
    //puck run 3
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-1000_nEvents-1000000-EP.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_noptDep_dNdeta-10000_nEvents-1000000-EP-HiMulti.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-1000_nEvents-1000000-EP.root";
    //TString sFileName = "toyFlow_noWeight_randomPsi_ptDep_dNdeta-10000_nEvents-1000000-EP-HiMulti.root";

    TFile *fIn = TFile::Open(sFileName,"read");
    if(fIn==0x0) {
        cout << "File not found, terminating." << endl;
        return -1;
    }

    
    const int nCoef = 5;
    const int nPtBins = 10;
    double inputFlow[nCoef] = {0.0};
    
    TH1D *hInputNumbers = (TH1D*)fIn->Get("hInputNumbers");
    double haddScaling = hInputNumbers->GetBinContent(22);
    double nEvents = hInputNumbers->GetBinContent(1);
    double dNdeta = hInputNumbers->GetBinContent(2)/haddScaling;
    double etaRange = hInputNumbers->GetBinContent(3)/haddScaling;
    double nMult = hInputNumbers->GetBinContent(4)/haddScaling;
    for(int i = 0; i < nCoef; i++) inputFlow[i] = hInputNumbers->GetBinContent(5+i)/haddScaling;
    double Tdec = hInputNumbers->GetBinContent(10)/haddScaling;
    double vr = hInputNumbers->GetBinContent(11)/haddScaling;
    double Teff = hInputNumbers->GetBinContent(12)/haddScaling;
    double slope = hInputNumbers->GetBinContent(13)/haddScaling;
    double const detAMax = hInputNumbers->GetBinContent(14)/haddScaling;
    double const detAMin = hInputNumbers->GetBinContent(15)/haddScaling;
    double const detBMax = hInputNumbers->GetBinContent(16)/haddScaling;
    double const detBMin = hInputNumbers->GetBinContent(17)/haddScaling;
    double const detAEff = hInputNumbers->GetBinContent(18)/haddScaling;
    double const detBEff = hInputNumbers->GetBinContent(19)/haddScaling;
    double const alpha = hInputNumbers->GetBinContent(20)/haddScaling;
    double const pTMax = hInputNumbers->GetBinContent(21)/haddScaling;
    cout << "alpha: " << alpha << ", pTMax: " << pTMax << endl;

    
    TH1D *hFlowIn = new TH1D("hFlowIn", "hFlowIn", nCoef, 0.5, double(nCoef)+0.5);
    hFlowIn->SetLineStyle(1);
    hFlowIn->SetLineColor(1);
    hFlowIn->SetLineWidth(2);
    for(int i = 0; i < nCoef; i++) hFlowIn->SetBinContent(i+1,inputFlow[i]);
    
    TH1D *hPhi[nType];
    hPhi[iType] = (TH1D*)fIn->Get(Form("hPhiT%02i",iType));
    hPhi[iType]->Scale(1./nEvents,"width");
    hPhi[iType]->SetMarkerStyle(24);
    hPhi[iType]->SetMarkerColor(1);
    hPhi[iType]->SetMarkerSize(0.8);
    
    TH1D *hPt[nType];
    hPt[iType] = (TH1D*)fIn->Get(Form("hPtT%02i",iType));
    hPt[iType]->Scale(1./nEvents,"width");
    hPt[iType]->SetMarkerStyle(24);
    hPt[iType]->SetMarkerColor(1);
    hPt[iType]->SetMarkerSize(0.8);
    
    TF1 *fInputPt = new TF1("fInputPt", "[0]*exp(-[1]*x)", 0.0, 10.0);
    fInputPt->SetParameters(10.0, slope);
    fInputPt->SetLineStyle(1);
    fInputPt->SetLineColor(1);
    fInputPt->SetLineWidth(2);
    hPt[iType]->Fit("fInputPt","NRO");
    
    TH1D *hEta[nType];
    hEta[iType] = (TH1D*)fIn->Get(Form("hEtaT%02i",iType));
    hEta[iType]->Scale(1./nEvents,"width");
    hEta[iType]->SetMarkerStyle(24);
    hEta[iType]->SetMarkerColor(1);
    hEta[iType]->SetMarkerSize(0.8);
    
    TH1D *hMultiplicity[nType];
    hMultiplicity[0] = (TH1D*)fIn->Get(Form("hMultiplicityT%02i",0));
    hMultiplicity[0]->Scale(1./nEvents,"width");
    hMultiplicity[0]->SetMarkerStyle(24);
    hMultiplicity[0]->SetMarkerColor(1);
    hMultiplicity[0]->SetMarkerSize(0.8);
    if(iType!=0) {
        hMultiplicity[iType] = (TH1D*)fIn->Get(Form("hMultiplicityT%02i",iType));
        hMultiplicity[iType]->Scale(1./nEvents,"width");
        hMultiplicity[iType]->SetMarkerStyle(24);
        hMultiplicity[iType]->SetMarkerColor(1);
        hMultiplicity[iType]->SetMarkerSize(0.8);
    }
    
    TH1D *hSqrtSumWeights[nType];
    hSqrtSumWeights[iType] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsT%02i",iType));
    hSqrtSumWeights[iType]->Scale(1./nEvents,"width");
    hSqrtSumWeights[iType]->SetMarkerStyle(24);
    hSqrtSumWeights[iType]->SetMarkerColor(1);
    hSqrtSumWeights[iType]->SetMarkerSize(0.8);
    
    TH1D *hSqrtSumWeightsA[nType];
    hSqrtSumWeightsA[iType] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsAT%02i",iType));
    hSqrtSumWeightsA[iType]->Scale(1./nEvents,"width");
    hSqrtSumWeightsA[iType]->SetMarkerStyle(24);
    hSqrtSumWeightsA[iType]->SetMarkerColor(1);
    hSqrtSumWeightsA[iType]->SetMarkerSize(0.8);
    
    TH1D *hSqrtSumWeightsB[nType];
    hSqrtSumWeightsB[iType] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsBT%02i",iType));
    hSqrtSumWeightsB[iType]->Scale(1./nEvents,"width");
    hSqrtSumWeightsB[iType]->SetMarkerStyle(24);
    hSqrtSumWeightsB[iType]->SetMarkerColor(1);
    hSqrtSumWeightsB[iType]->SetMarkerSize(0.8);

    TH1D *hSqrtSumWeightsPtBins[nType][nPtBins];
    for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        hSqrtSumWeightsPtBins[iType][iPtBin] = (TH1D*)fIn->Get(Form("hSqrtSumWeightsPtBinsT%02iPtB%02i",iType,iPtBin));
        hSqrtSumWeightsPtBins[iType][iPtBin]->Scale(1./nEvents,"width");
        hSqrtSumWeightsPtBins[iType][iPtBin]->SetMarkerStyle(24);
        hSqrtSumWeightsPtBins[iType][iPtBin]->SetMarkerColor(1);
        hSqrtSumWeightsPtBins[iType][iPtBin]->SetMarkerSize(0.8);
    }
    
    double meanMultiplicity = hMultiplicity[iType]->GetMean();
    double meanMultiplicityTrue = hMultiplicity[0]->GetMean();
    double meanMultiplicityError = hMultiplicity[iType]->GetMeanError();
    
    double weight = hSqrtSumWeights[iType]->GetMean();
    double weightError = hSqrtSumWeights[iType]->GetMeanError();
    double weightPtBins[nPtBins];
    double weightErrorPtBins[nPtBins];
    for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        weightPtBins[iPtBin] = hSqrtSumWeightsPtBins[iType][iPtBin]->GetMean();
        weightErrorPtBins[iPtBin] = hSqrtSumWeightsPtBins[iType][iPtBin]->GetMeanError();
    }

    if( bUseWeightning ){
        weight *= TMath::Sqrt( meanMultiplicity );
        weightError *= TMath::Sqrt( meanMultiplicity ); // Need to update error from mean multiplicity itself
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
            weightPtBins[iPtBin] *= TMath::Sqrt( meanMultiplicity ); // Need to update error from mean multiplicity itself
            weightErrorPtBins[iPtBin] *= TMath::Sqrt( meanMultiplicity ); // Need to update error from mean multiplicity itself
        }
    }

    TF1 *fAcceptanceFunc = new TF1("fAcceptanceFunc", AcceptanceFunc, -TMath::Pi(), TMath::Pi(), 10);
    fAcceptanceFunc->SetParameter(0,1.0);
    fAcceptanceFunc->SetParameter(1,detAEff);
    fAcceptanceFunc->SetParameter(2,detBEff);
    fAcceptanceFunc->SetParameter(3,detAMax);
    fAcceptanceFunc->SetParameter(4,detAMin);
    fAcceptanceFunc->SetParameter(5,detBMax);
    fAcceptanceFunc->SetParameter(6,detBMin);
    //const double accIntegral = fAcceptanceFunc->Integral(-pi,pi);
    fAcceptanceFunc->SetParameter(0,meanMultiplicityTrue);
    fAcceptanceFunc->SetLineStyle(1);
    fAcceptanceFunc->SetLineColor(1);
    fAcceptanceFunc->SetLineWidth(2);
    
    
    // integrated flow
    TH1D *hQx[nType][nCoef];
    TH1D *hQy[nType][nCoef];
    TH1D *hQ[nType][nCoef];
    TH1D *hQ2[nType][nCoef];
    TH1D *hQ4[nType][nCoef];
    TH1D *hQ6[nType][nCoef];
    TH1D *hvObs[nType][nCoef];
    TH1D *hvObsPtBins[nType][nCoef][nPtBins];
    TH1D *hEPnominatorPtBins[nType][nCoef][nPtBins];
    TH1D *hrSub[nType][nCoef];
    TH1D *hTrueReso[nType][nCoef];
    TH1D *hTrueResoPtBins[nType][nCoef][nPtBins];
    TH2D *hPsiAB[nType][nCoef][2];
    TH1D *hPsiDiff[nType][nCoef][2];
    TH1D *hPsiDiffN[nType][nCoef][2];
    TH1D *hSPnominator[nType][nCoef];
    TH1D *hSPdenominator[nType][nCoef];
    TH1D *hEPnominator[nType][nCoef];
    TH1D *hEPdenominator[nType][nCoef];
    
    for(int i = 0; i < nCoef; i++){
        hQx[iType][i] = (TH1D*)fIn->Get( Form("hQxT%02iH%02i",iType,i+1) );
        hQx[iType][i]->Scale(1./hQx[iType][i]->GetEntries(),"width");
        hQx[iType][i]->SetMarkerStyle(20);
        hQx[iType][i]->SetMarkerColor(2);
        hQx[iType][i]->SetMarkerSize(0.8);
        
        hQy[iType][i] = (TH1D*)fIn->Get( Form("hQyT%02iH%02i",iType,i+1) );
        hQy[iType][i]->Scale(1./hQy[iType][i]->GetEntries(),"width");
        hQy[iType][i]->SetMarkerStyle(20);
        hQy[iType][i]->SetMarkerColor(2);
        hQy[iType][i]->SetMarkerSize(0.8);
        
        hQ[iType][i]  = (TH1D*)fIn->Get( Form("hQT%02iH%02i",iType,i+1) );
        hQ[iType][i]->Scale(1./hQ[iType][i]->GetEntries(),"width");
        hQ[iType][i]->SetMarkerStyle(20);
        hQ[iType][i]->SetMarkerColor(2);
        hQ[iType][i]->SetMarkerSize(0.8);
        
        hQ2[iType][i] = (TH1D*)fIn->Get( Form("hQ2T%02iH%02i",iType,i+1) );
        hQ2[iType][i]->Scale(1./hQ2[iType][i]->GetEntries(),"width");
        hQ2[iType][i]->SetMarkerStyle(20);
        hQ2[iType][i]->SetMarkerColor(2);
        hQ2[iType][i]->SetMarkerSize(0.8);
        
        hQ4[iType][i] = (TH1D*)fIn->Get( Form("hQ4T%02iH%02i",iType,i+1) );
        hQ4[iType][i]->Scale(1./hQ4[iType][i]->GetEntries(),"width");
        hQ4[iType][i]->SetMarkerStyle(20);
        hQ4[iType][i]->SetMarkerColor(2);
        hQ4[iType][i]->SetMarkerSize(0.8);
        
        hQ6[iType][i] = (TH1D*)fIn->Get( Form("hQ6T%02iH%02i",iType,i+1) );
        hQ6[iType][i]->Scale(1./hQ6[iType][i]->GetEntries(),"width");
        hQ6[iType][i]->SetMarkerStyle(20);
        hQ6[iType][i]->SetMarkerColor(2);
        hQ6[iType][i]->SetMarkerSize(0.8);
        
        hvObs[iType][i] = (TH1D*)fIn->Get( Form("hvObsT%02iH%02i",iType,i+1) );
        hvObs[iType][i]->Scale(1./hvObs[iType][i]->GetEntries(),"width");
        hvObs[iType][i]->SetMarkerStyle(20);
        hvObs[iType][i]->SetMarkerColor(2);
        hvObs[iType][i]->SetMarkerSize(0.8);
        checkUnderOverFlow(hvObs[iType][i]);
        
        hrSub[iType][i] = (TH1D*)fIn->Get( Form("hrSubT%02iH%02i",iType,i+1) );
        hrSub[iType][i]->Scale(1./hrSub[iType][i]->GetEntries(),"width");
        hrSub[iType][i]->SetMarkerStyle(20);
        hrSub[iType][i]->SetMarkerColor(2);
        hrSub[iType][i]->SetMarkerSize(0.8);
        checkUnderOverFlow(hrSub[iType][i]);
        
        hTrueReso[iType][i] = (TH1D*)fIn->Get( Form("hTrueResoT%02iH%02i",iType,i+1) );
        hTrueReso[iType][i]->Scale(1./hTrueReso[iType][i]->GetEntries(),"width");
        hTrueReso[iType][i]->SetMarkerStyle(20);
        hTrueReso[iType][i]->SetMarkerColor(1);
        hTrueReso[iType][i]->SetMarkerSize(0.8);
        checkUnderOverFlow(hTrueReso[iType][i]);
        
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
            hvObsPtBins[iType][i][iPtBin] = (TH1D*)fIn->Get( Form("hvObsPtBinsT%02iH%02iPtB%02i",iType,i+1,iPtBin) );
            hvObsPtBins[iType][i][iPtBin]->Scale(1./hvObsPtBins[iType][i][iPtBin]->GetEntries(),"width");
            hvObsPtBins[iType][i][iPtBin]->SetMarkerStyle(20);
            hvObsPtBins[iType][i][iPtBin]->SetMarkerColor(1+i);
            hvObsPtBins[iType][i][iPtBin]->SetMarkerSize(0.8);
            checkUnderOverFlow(hvObsPtBins[iType][i][iPtBin]);

            hEPnominatorPtBins[iType][i][iPtBin] = (TH1D*)fIn->Get( Form("hEPnominatorPtBinsT%02iH%02iPtB%02i",iType,i+1,iPtBin) );
            hEPnominatorPtBins[iType][i][iPtBin]->Scale(1./hEPnominatorPtBins[iType][i][iPtBin]->GetEntries(),"width");
            hEPnominatorPtBins[iType][i][iPtBin]->SetMarkerStyle(20);
            hEPnominatorPtBins[iType][i][iPtBin]->SetMarkerColor(1+i);
            hEPnominatorPtBins[iType][i][iPtBin]->SetMarkerSize(0.8);
            checkUnderOverFlow(hEPnominatorPtBins[iType][i][iPtBin]);

            hTrueResoPtBins[iType][i][iPtBin] = (TH1D*)fIn->Get( Form("hTrueResoPtBinsT%02iH%02iPtB%02i",iType,i+1,iPtBin) );
            hTrueResoPtBins[iType][i][iPtBin]->Scale(1./hTrueResoPtBins[iType][i][iPtBin]->GetEntries(),"width");
            hTrueResoPtBins[iType][i][iPtBin]->SetMarkerStyle(20);
            hTrueResoPtBins[iType][i][iPtBin]->SetMarkerColor(1+i);
            hTrueResoPtBins[iType][i][iPtBin]->SetMarkerSize(0.8);
            checkUnderOverFlow(hTrueResoPtBins[iType][i][iPtBin]);
        }
        
        for(int j=0; j<2; j++) {
            hPsiAB[iType][i][j] = (TH2D*)fIn->Get( Form("hPsiABT%02iH%02iE%02i",iType,i+1,j) );
            //hPsiAB[iType][i][j]->Scale(1./hPsiAB[iType][i][j]->GetEntries(),"width");
            //hPsiAB[iType][i][j]->SetMarkerStyle(20+4*j);
            hPsiAB[iType][i][j]->SetMarkerColor(2+2*j);
            //hPsiAB[iType][i][j]->SetMarkerSize(0.8);

            hPsiDiff[iType][i][j] = (TH1D*)fIn->Get( Form("hPsiDiffT%02iH%02iE%02i",iType,i+1,j) );
            hPsiDiff[iType][i][j]->Scale(1./hPsiDiff[iType][i][0]->GetEntries(),"width");
            hPsiDiff[iType][i][j]->SetMarkerStyle(20+4*j);
            hPsiDiff[iType][i][j]->SetMarkerColor(2+2*j);
            hPsiDiff[iType][i][j]->SetMarkerSize(0.8);
            checkUnderOverFlow(hPsiDiff[iType][i][j]);

            hPsiDiffN[iType][i][j] = (TH1D*)fIn->Get( Form("hPsiDiffNT%02iH%02iE%02i",iType,i+1,j) );
            hPsiDiffN[iType][i][j]->Scale(1./hPsiDiffN[iType][i][0]->GetEntries(),"width");
            hPsiDiffN[iType][i][j]->SetMarkerStyle(20+4*j);
            hPsiDiffN[iType][i][j]->SetMarkerColor(2+2*j);
            hPsiDiffN[iType][i][j]->SetMarkerSize(0.8);
            checkUnderOverFlow(hPsiDiffN[iType][i][j]);
        }

        hSPnominator[iType][i] = (TH1D*)fIn->Get( Form("hSPnominatorVT%02iH%02i",iType,i+1) );
        hSPnominator[iType][i]->Scale(1./hSPnominator[iType][i]->GetEntries(),"width");
        hSPnominator[iType][i]->SetMarkerStyle(20);
        hSPnominator[iType][i]->SetMarkerColor(2);
        hSPnominator[iType][i]->SetMarkerSize(0.8);
        
        hSPdenominator[iType][i] = (TH1D*)fIn->Get( Form("hSPdenominatorVT%02iH%02i",iType,i+1) );
        hSPdenominator[iType][i]->Scale(1./hSPdenominator[iType][i]->GetEntries(),"width");
        hSPdenominator[iType][i]->SetMarkerStyle(20);
        hSPdenominator[iType][i]->SetMarkerColor(2);
        hSPdenominator[iType][i]->SetMarkerSize(0.8);
        
        hEPnominator[iType][i] = (TH1D*)fIn->Get( Form("hEPnominatorVT%02iH%02i",iType,i+1) );
        hEPnominator[iType][i]->Scale(1./hEPnominator[iType][i]->GetEntries(),"width");
        hEPnominator[iType][i]->SetMarkerStyle(20);
        hEPnominator[iType][i]->SetMarkerColor(2);
        hEPnominator[iType][i]->SetMarkerSize(0.8);
        
        hEPdenominator[iType][i] = (TH1D*)fIn->Get( Form("hEPdenominatorVT%02iH%02i",iType,i+1) );
        hEPdenominator[iType][i]->Scale(1./hEPdenominator[iType][i]->GetEntries(),"width");
        hEPdenominator[iType][i]->SetMarkerStyle(20);
        hEPdenominator[iType][i]->SetMarkerColor(2);
        hEPdenominator[iType][i]->SetMarkerSize(0.8);
    }
    
    // means and mean errors from simulation
    double meanQx[nCoef] = {0}, meanErrorQx[nCoef] = {0};
    double meanQy[nCoef] = {0}, meanErrorQy[nCoef] = {0};
    double meanQ[nCoef]  = {0}, meanErrorQ[nCoef]  = {0};
    double meanQ2[nCoef] = {0}, meanErrorQ2[nCoef] = {0};
    double meanQ4[nCoef] = {0}, meanErrorQ4[nCoef] = {0};
    double meanQ6[nCoef] = {0}, meanErrorQ6[nCoef] = {0};
    double meanvObs[nCoef] = {0}, meanErrorvObs[nCoef] = {0};
    double meanvObsPtBins[nCoef][nPtBins] = {{0}}, meanErrorvObsPtBins[nCoef][nPtBins] = {{0}};
    double meanEPnomPtBins[nCoef][nPtBins] = {{0}}, meanErrorEPnomPtBins[nCoef][nPtBins] = {{0}};
    double meanvRealEP[nCoef] = {0}, meanErrorvRealEP[nCoef] = {0};
    double meanvRealEPTrue[nCoef] = {0}, meanErrorvRealEPTrue[nCoef] = {0};
    double meanvRealEPTruePtBins[nCoef][nPtBins] = {0}, meanErrorvRealEPTruePtBins[nCoef][nPtBins] = {0};
    double meanrSub[nCoef] = {0}, meanErrorrSub[nCoef] = {0};
    double meanrTrue[nCoef] = {0}, meanErrorrTrue[nCoef] = {0};
    double meanrTruePtBins[nCoef][nPtBins] = {{0}}, meanErrorrTruePtBins[nCoef][nPtBins] = {{0}};
    double meanrSubFull[nCoef] = {0}, meanErrorrSubFull[nCoef] = {0};
    double meanSPnom[nCoef] = {0}, meanErrorSPnom[nCoef] = {0};
    double meanSPdenom[nCoef] = {0}, meanErrorSPdenom[nCoef] = {0};
    double meanEPnom[nCoef] = {0}, meanErrorEPnom[nCoef] = {0};
    double meanEPdenom[nCoef] = {0}, meanErrorEPdenom[nCoef] = {0};
    // various estimates for vn's calculated from the the means above
    double vFitQdist[nCoef], vFitQdistError[nCoef]; // mean of the Q-vector
    double v1cumul[nCoef], v1cumulError[nCoef]; // mean of the Q-vector
    double v4cumul[nCoef], v4cumulError[nCoef]; // 4th order cumulants
    double v6cumul[nCoef], v6cumulError[nCoef]; // 6th order cumulants
    double vEP[nCoef], vEPError[nCoef];         // Event Plane -method
    double vSP[nCoef], vSPError[nCoef];         // Scalar Product -method
    double vEPPtBins[nCoef][nPtBins], vEPErrorPtBins[nCoef][nPtBins];         // Event Plane -method
    
    // Extract flow from the fit to the Q
    TF1 *fitQ = new TF1("fitQ", fFitMeanQ, 0.0, 15.0, 3);
    fitQ->SetLineColor(1);
    fitQ->SetLineStyle(1);
    fitQ->SetLineWidth(2);
    int iFig = 1;
    for(int i = 0; i < nCoef; i++){
        
        fitQ->SetParameters(1.0, hQ[iType][i]->GetBinCenter( hQ[iType][i]->GetMaximumBin() ), 1.0);
        //fitQ->FixParameter(0,1.0);
        hQ[iType][i]->Fit("fitQ","RNO");
        meanQ[i]  = TMath::Max(0.0,fitQ->GetParameter(1));  meanErrorQ[i]  = fitQ->GetParError(1);

        mc(iFig++);
        
        gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
        gPad->SetLogx(0); gPad->SetLogy(0);
        //mpad->SetGridx(0); mpad->SetGridy(0);
        
        TH2F *hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, 0.0, 12.0, 1, 0.0, 0.9);
        myhset( *hfr, "Q", "1/N_{entries} dN/dQ",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfr->Draw();
        
        hQ[iType][i]->Draw("same,p");
        fitQ->DrawCopy("same,l");
        
        TLegend *leg = new TLegend(0.50, 0.68, 0.80, 0.93,"","brNDC");
        leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
        leg->AddEntry(hQ[iType][i],"toy simulation","p");
        leg->AddEntry(fitQ,"fit, resulting parameters:","l");
        leg->AddEntry(fitQ,Form("<Q> = %4.3f +- %4.3f",fitQ->GetParameter(1), TMath::Max(0.001,fitQ->GetParError(1))),"");
        leg->AddEntry(fitQ,Form("#sigma = %4.3f +- %4.3f",fitQ->GetParameter(2), TMath::Max(0.001,fitQ->GetParError(2))),"");
        leg->AddEntry(fitQ,Form("n = %d",i+1),"");
        leg->Draw();
    }

    for(int i = 0; i < nCoef; i++){
        mc(iFig++);
        gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
        float mleft = 0.15, mbott = 0.15, mtop= 0.06, mright = 0.15;
        TPad* mpad = new TPad("mpad", "mpad", 0.01, 0.01, 0.99, 0.99, 0,0,0);
        mpad->SetLeftMargin(mleft);
        mpad->SetBottomMargin(mbott);
        mpad->SetTopMargin(mtop);
        mpad->SetRightMargin(mright);
        mpad->SetLogx(0);
        mpad->SetLogy(0);
        mpad->SetLogz(0);
        mpad->Draw();
        mpad->cd();

        //TH2F *hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -3.2, 3.2, 1, -3.2, 3.2);
        //
        //hfr->Draw();

        hPsiAB[iType][i][0]->Draw("same");
        hPsiAB[iType][i][1]->Draw("same");
        myhset( *hPsiAB[iType][i][iDrawNegRHisto], Form("%d * #Psi_{A}",i+1), Form("%d * #Psi_{B}",i+1),0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        TLegend *leg = new TLegend(0.50, 0.68, 0.80, 0.93,"","brNDC");
        leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
        leg->AddEntry(hPsiAB[iType][i][iDrawNegRHisto],Form("v_%d",i+1),"");
        leg->Draw("same");
    }

    for(int i = 0; i < nCoef; i++){
        mc(iFig++);
        gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
        gPad->SetLogx(0); gPad->SetLogy(0);

        TH2F *hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -6.4, 6.4, 1, 0.0, 3.5);
        myhset( *hfr, "#Psi_{A} - #Psi_{B}", "Probability",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfr->Draw();
        
        hPsiDiff[iType][i][0]->Draw("same");
        hPsiDiff[iType][i][1]->Draw("same");

        TLegend *leg = new TLegend(0.50, 0.68, 0.80, 0.93,"","brNDC");
        leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
        leg->AddEntry(hPsiAB[iType][i][iDrawNegRHisto],Form("v_%d",i+1),"");
        leg->Draw("same");
    }
    
    for(int i = 0; i < nCoef; i++){
        mc(iFig++);
        gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
        gPad->SetLogx(0); gPad->SetLogy(0);

        TH2F *hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -6.4, 6.4, 1, 0.0, 3.5);
        myhset( *hfr, Form("%d * (#Psi_{A} - #Psi_{B})",i+1), "Probability",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfr->Draw();
        
        hPsiDiffN[iType][i][0]->Draw("same");
        hPsiDiffN[iType][i][1]->Draw("same");

        TLegend *leg = new TLegend(0.50, 0.68, 0.80, 0.93,"","brNDC");
        leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
        leg->AddEntry(hPsiAB[iType][i][iDrawNegRHisto],Form("v_%d",i+1),"");
        leg->Draw("same");
    }
    double arg;
    double iterDiff, iterRSub, tryXi, prevIterDiff;
    for(int i = 0; i < nCoef; i++){
        // Get means and mean errors
        meanQx[i] = hQx[iType][i]->GetMean(); meanErrorQx[i] = hQx[iType][i]->GetMeanError();
        meanQy[i] = hQy[iType][i]->GetMean(); meanErrorQy[i] = hQy[iType][i]->GetMeanError();
        meanQ2[i] = hQ2[iType][i]->GetMean(); meanErrorQ2[i] = hQ2[iType][i]->GetMeanError();
        meanQ4[i] = hQ4[iType][i]->GetMean(); meanErrorQ4[i] = hQ4[iType][i]->GetMeanError();
        meanQ6[i] = hQ6[iType][i]->GetMean(); meanErrorQ6[i] = hQ6[iType][i]->GetMeanError();
        meanvObs[i] = hvObs[iType][i]->GetMean(); meanErrorvObs[i] = hvObs[iType][i]->GetMeanError();
        meanrSub[i] = sqrt(hrSub[iType][i]->GetMean()); meanErrorrSub[i] = hrSub[iType][i]->GetMeanError()/(2*sqrt(meanrSub[i]));
        meanrTrue[i] = hTrueReso[iType][i]->GetMean(); meanErrorrTrue[i] = hTrueReso[iType][i]->GetMeanError();
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
            meanvObsPtBins[i][iPtBin] = hvObsPtBins[iType][i][iPtBin]->GetMean(); meanErrorvObsPtBins[i][iPtBin] = hvObsPtBins[iType][i][iPtBin]->GetMeanError();
            meanrTruePtBins[i][iPtBin] = hTrueResoPtBins[iType][i][iPtBin]->GetMean(); meanErrorrTruePtBins[i][iPtBin] = hTrueResoPtBins[iType][i][iPtBin]->GetMeanError();
            meanEPnomPtBins[i][iPtBin] = hEPnominatorPtBins[iType][i][iPtBin]->GetMean(); meanErrorEPnomPtBins[i][iPtBin] = hEPnominatorPtBins[iType][i][iPtBin]->GetMeanError();
        }
        tryXi=0.001;
        iterDiff=1.0;
        iterRSub = 0.0;
        while(iterDiff>0.001) {
            iterRSub = CalculateResolutionKOne(tryXi);
            if(iterRSub < meanrSub[i]) tryXi += 0.01*iterDiff;
            else tryXi -= 0.01*iterDiff;
            iterDiff = TMath::Abs(iterRSub - meanrSub[i]);
            prevIterDiff = iterDiff;
        }
        //meanrSubFull[i] = meanrSubSquared[i]<0.0 ? 0.01 : sqrt(2)*sqrt(meanrSubSquared[i]); meanErrorrSubFull[i] = sqrt(2)*sqrt(meanErrorrSubSquared[i]);
        meanrSubFull[i] = CalculateResolutionKOne(sqrt(2)*tryXi); meanErrorrSubFull[i] = sqrt(2)*sqrt(meanErrorrSub[i]);
        cout << "meanrSubFull[" << i << "]: " << meanrSubFull[i] << endl;
        meanvRealEP[i] = meanvObs[i]/meanrSubFull[i];
        meanErrorvRealEP[i] = meanvRealEP[i] * TMath::Sqrt( meanErrorvObs[i]*meanErrorvObs[i]/meanvObs[i]/meanvObs[i] +
                                                            meanErrorvRealEP[i]*meanErrorvRealEP[i]/meanvRealEP[i]/meanvRealEP[i] ); 
        meanvRealEPTrue[i] = meanvObs[i]/meanrTrue[i];
        meanErrorvRealEPTrue[i] = meanvRealEPTrue[i] * TMath::Sqrt( meanErrorvObs[i]*meanErrorvObs[i]/meanvObs[i]/meanvObs[i] +
                                                                    meanErrorvRealEPTrue[i]*meanErrorvRealEPTrue[i]/meanvRealEPTrue[i]/meanvRealEPTrue[i] ); 
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
            meanvRealEPTruePtBins[i][iPtBin] = meanvObsPtBins[i][iPtBin]/meanrTruePtBins[i][iPtBin];
            meanErrorvRealEPTruePtBins[i][iPtBin] = meanvRealEPTruePtBins[i][iPtBin] *
                                                    TMath::Sqrt( meanErrorvObsPtBins[i][iPtBin]*meanErrorvObsPtBins[i][iPtBin]/meanvObsPtBins[i][iPtBin]/meanvObsPtBins[i][iPtBin] +
                                                                 meanErrorvRealEPTruePtBins[i][iPtBin]*meanErrorvRealEPTruePtBins[i][iPtBin]/meanvRealEPTruePtBins[i][iPtBin]/meanvRealEPTruePtBins[i][iPtBin] ); 
        }
        mc(iFig++);
        TH2F *hfrR = new TH2F(Form("hfr%02d",iFig)," ", 1, -1.1, 1.1, 1, 0.01, hvObs[iType][i]->GetMaximum()+1);
        myhset( *hfrR, Form("v^{obs}_{%i}",i+1), "Probability",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfrR->Draw();
        hvObs[iType][i]->Draw("same");
        if(hrSub[iType][i]->GetBinContent(0)>0) cout << "rSub underflow bin not empty!" << endl;

        mc(iFig++);
        hfrR = new TH2F(Form("hfr%02d",iFig)," ", 1, -1.1, 1.1, 1, 0.01, hrSub[iType][i]->GetMaximum()+1);
        myhset( *hfrR, Form("R_{%i}",i+1), "Probability",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
        hfrR->Draw();
        hrSub[iType][i]->Draw("same");
        hTrueReso[iType][i]->Draw("same");

        if(hrSub[iType][i]->GetBinContent(hrSub[iType][i]->GetXaxis()->GetNbins()+1)>0) cout << "rSub overflow bin not empty!" << endl;
        meanSPnom[i] = hSPnominator[iType][i]->GetMean(); meanErrorSPnom[i] = hSPnominator[iType][i]->GetMeanError();
        meanSPdenom[i] = hSPdenominator[iType][i]->GetMean(); meanErrorSPdenom[i] = hSPdenominator[iType][i]->GetMeanError();
        meanEPnom[i] = hEPnominator[iType][i]->GetMean(); meanErrorEPnom[i] = hEPnominator[iType][i]->GetMeanError();
        meanEPdenom[i] = hEPdenominator[iType][i]->GetMean(); meanErrorEPdenom[i] = hEPdenominator[iType][i]->GetMeanError();
        // calculate various estimates for vn's
        v1cumul[i] = meanQ[i] / weight;
        v4cumul[i] = CalculateV2Using4Cumul( meanQ4[i], meanQ2[i], weight );
        v6cumul[i] = CalculateV2Using6Cumul( meanQ6[i], meanQ4[i], meanQ2[i], weight );
        vEP[i] = CalculateV2UsingEPorSP(  meanEPnom[i], meanEPdenom[i], weight );
        vSP[i] = CalculateV2UsingEPorSP(  meanSPnom[i], meanSPdenom[i], weight );
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) vEPPtBins[i][iPtBin] = CalculateV2UsingEPorSP(  meanEPnomPtBins[i][iPtBin], meanEPdenom[i], weightPtBins[iPtBin] );
        // errors
        v1cumulError[i] = meanErrorQ[i]/weight + meanQ[i]/weight/weight*weightError;
        v4cumulError[i] = CalculateV2ErrorUsing4Cumul( meanQ4[i], meanQ2[i], meanErrorQ4[i], meanErrorQ2[i], weight, weightError );
        v6cumulError[i] = CalculateV2ErrorUsing6Cumul( meanQ6[i], meanQ4[i], meanQ2[i], meanErrorQ6[i], meanErrorQ4[i], meanErrorQ2[i], weight, weightError );
        vEPError[i] = CalculateV2ErrorUsingEPorSP(  meanEPnom[i], meanEPdenom[i], meanErrorEPnom[i], meanErrorEPdenom[i], weight, weightError );
        vSPError[i] = CalculateV2ErrorUsingEPorSP(  meanSPnom[i], meanSPdenom[i], meanErrorSPnom[i], meanErrorSPdenom[i], weight, weightError );
        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) vEPErrorPtBins[i][iPtBin] = CalculateV2ErrorUsingEPorSP(  meanEPnomPtBins[i][iPtBin], meanEPdenom[i], meanErrorEPnomPtBins[i][iPtBin], meanErrorEPdenom[i], weightPtBins[iPtBin], weightErrorPtBins[iPtBin] ); //TODO: implement pt-binned weight
        cout << i << "  " << v1cumul[i] << "  " << v4cumul[i] << "  " << v6cumul[i] << "  " << vEP[i] << "  " << vSP[i] << endl;
    }
    
    TF1 *fPhiInput = new TF1("fPhiInput","[0] * (1. + 2.*[1]*cos(x) + 2.*[2]*cos(2.*x) + 2.*[3]*cos(3.*x) + 2.*[4]*cos(4.*x) + 2.*[5]*cos(5.*x) )",-3.2,3.2);
    fPhiInput->SetParameters(1.0, inputFlow[0], inputFlow[1], inputFlow[2], inputFlow[3], inputFlow[4]);
    double idealNormForPhi = meanMultiplicity / 2. / pi;
    fPhiInput->SetParameter(0, idealNormForPhi);
    fPhiInput->FixParameter(0, idealNormForPhi);
    fPhiInput->FixParameter(1, inputFlow[0]);
    fPhiInput->FixParameter(2, inputFlow[1]);
    fPhiInput->FixParameter(3, inputFlow[2]);
    fPhiInput->FixParameter(4, inputFlow[3]);
    fPhiInput->FixParameter(5, inputFlow[4]);
    fPhiInput->SetLineStyle(2);
    fPhiInput->SetLineColor(4);
    fPhiInput->SetLineWidth(2);
    
    TF1 *fitPhi = new TF1("fitPhi","[0] * (1. + 2.*[1]*cos(x) + 2.*[2]*cos(2.*x) + 2.*[3]*cos(3.*x) + 2.*[4]*cos(4.*x) + 2.*[5]*cos(5.*x) )",-3.1,3.1);
    fPhiInput->SetParameters(idealNormForPhi, inputFlow[0], inputFlow[1], inputFlow[2], inputFlow[3], inputFlow[4]);
    fitPhi->SetLineStyle(1);
    fitPhi->SetLineColor(2);
    fitPhi->SetLineWidth(2);
    
    hPhi[iType]->Fit("fitPhi","NO");
    
    cout << Form("FIT: v1 = %5.4f +- %5.4f, v2 = %5.4f +- %5.4f, v3 = %5.4f +- %5.4f, v4 = %5.4f +- %5.4f, v5 = %5.4f +- %5.4f", fitPhi->GetParameter(1), fitPhi->GetParError(1), fitPhi->GetParameter(2), fitPhi->GetParError(2), fitPhi->GetParameter(3), fitPhi->GetParError(3), fitPhi->GetParameter(4), fitPhi->GetParError(4), fitPhi->GetParameter(5), fitPhi->GetParError(5) ) << endl;
    
    // Final results
    TGraphErrors *gv2with1cumul = new TGraphErrors(nCoef);
    gv2with1cumul->SetMarkerStyle(20);
    gv2with1cumul->SetMarkerColor(2);
    gv2with1cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2with4cumul = new TGraphErrors(nCoef);
    gv2with4cumul->SetMarkerStyle(21);
    gv2with4cumul->SetMarkerColor(4);
    gv2with4cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2with6cumul = new TGraphErrors(nCoef);
    gv2with6cumul->SetMarkerStyle(33);
    gv2with6cumul->SetMarkerColor(1);
    gv2with6cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2withEP = new TGraphErrors(nCoef);
    gv2withEP->SetMarkerStyle(34);
    gv2withEP->SetMarkerColor(9);
    gv2withEP->SetMarkerSize(0.8);
    TGraphErrors *gv2withSP = new TGraphErrors(nCoef);
    gv2withSP->SetMarkerStyle(29);
    gv2withSP->SetMarkerColor(6);
    gv2withSP->SetMarkerSize(0.8);
    TGraphErrors *gv2Obs = new TGraphErrors(nCoef);
    gv2Obs->SetMarkerStyle(28);
    gv2Obs->SetMarkerColor(1);
    gv2Obs->SetMarkerSize(0.8);
    TGraphErrors *gv2withRealEP = new TGraphErrors(nCoef);
    gv2withRealEP->SetMarkerStyle(28);
    gv2withRealEP->SetMarkerColor(7);
    gv2withRealEP->SetMarkerSize(0.8);
    TGraphErrors *gv2withRealEPTrue = new TGraphErrors(nCoef);
    gv2withRealEPTrue->SetMarkerStyle(28);
    gv2withRealEPTrue->SetMarkerColor(8);
    gv2withRealEPTrue->SetMarkerSize(0.8);
    for(int i = 0; i < nCoef; i++){
        gv2with1cumul->SetPoint(i, double(i+1)-0.35, v1cumul[i]);
        gv2with4cumul->SetPoint(i,double(i+1)-0.25, v4cumul[i]);
        gv2with6cumul->SetPoint(i,double(i+1)-0.15, v6cumul[i]);
        gv2withEP->SetPoint(i,double(i+1)-0.05, vEP[i]);
        gv2withSP->SetPoint(i,double(i+1)+0.05, vSP[i]);
        gv2Obs->SetPoint(i,double(i+1)+0.15, meanvObs[i]);
        gv2withRealEP->SetPoint(i,double(i+1)+0.25, meanvRealEP[i]);
        gv2withRealEPTrue->SetPoint(i,double(i+1)+0.35, meanvRealEPTrue[i]);
        gv2with1cumul->SetPointError(i, 0.0, v1cumulError[i]);
        gv2with4cumul->SetPointError(i, 0.0, v4cumulError[i]);
        gv2with6cumul->SetPointError(i, 0.0, v6cumulError[i]);
        gv2withEP->SetPointError(i, 0.0, vEPError[i]);
        gv2withSP->SetPointError(i, 0.0, vSPError[i]);
        gv2Obs->SetPointError(i, 0.0, meanErrorvObs[i]);
        gv2withRealEP->SetPointError(i, 0.0, meanErrorvRealEP[i]);
        gv2withRealEPTrue->SetPointError(i, 0.0, meanErrorvRealEPTrue[i]);
    }
    
    if(bPrintGraphs){
        gv2with1cumul->Print();
        gv2with4cumul->Print();
        gv2with6cumul->Print();
        gv2withEP->Print();
        gv2withSP->Print();
        gv2Obs->Print();
        gv2withRealEP->Print();
        gv2withRealEPTrue->Print();
    }
    
    TGraphErrors *gv2ratioWith1cumul = new TGraphErrors(nCoef);
    gv2ratioWith1cumul->SetMarkerStyle(20);
    gv2ratioWith1cumul->SetMarkerColor(2);
    gv2ratioWith1cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWith4cumul = new TGraphErrors(nCoef);
    gv2ratioWith4cumul->SetMarkerStyle(21);
    gv2ratioWith4cumul->SetMarkerColor(4);
    gv2ratioWith4cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWith6cumul = new TGraphErrors(nCoef);
    gv2ratioWith6cumul->SetMarkerStyle(33);
    gv2ratioWith6cumul->SetMarkerColor(1);
    gv2ratioWith6cumul->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWithEP = new TGraphErrors(nCoef);
    gv2ratioWithEP->SetMarkerStyle(34);
    gv2ratioWithEP->SetMarkerColor(9);
    gv2ratioWithEP->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWithSP = new TGraphErrors(nCoef);
    gv2ratioWithSP->SetMarkerStyle(29);
    gv2ratioWithSP->SetMarkerColor(6);
    gv2ratioWithSP->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioObs = new TGraphErrors(nCoef);
    gv2ratioObs->SetMarkerStyle(28);
    gv2ratioObs->SetMarkerColor(1);
    gv2ratioObs->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWithRealEP = new TGraphErrors(nCoef);
    gv2ratioWithRealEP->SetMarkerStyle(28);
    gv2ratioWithRealEP->SetMarkerColor(7);
    gv2ratioWithRealEP->SetMarkerSize(0.8);
    TGraphErrors *gv2ratioWithRealEPTrue = new TGraphErrors(nCoef);
    gv2ratioWithRealEPTrue->SetMarkerStyle(28);
    gv2ratioWithRealEPTrue->SetMarkerColor(8);
    gv2ratioWithRealEPTrue->SetMarkerSize(0.8);
    for(int i = 0; i < nCoef; i++){
        double x, y, ey, ratio;
        double in = inputFlow[i];
        
        gv2with1cumul->GetPoint(i, x, y );
        ey = gv2with1cumul->GetErrorY(i);
        gv2ratioWith1cumul->SetPoint(i, x, y - in );
        gv2ratioWith1cumul->SetPointError(i, 0.0, ey );
        
        gv2with4cumul->GetPoint(i, x, y );
        ey = gv2with4cumul->GetErrorY(i);
        gv2ratioWith4cumul->SetPoint(i, x, y - in );
        gv2ratioWith4cumul->SetPointError(i, 0.0, ey );
        
        gv2with6cumul->GetPoint(i, x, y );
        ey = gv2with6cumul->GetErrorY(i);
        gv2ratioWith6cumul->SetPoint(i, x, y - in );
        gv2ratioWith6cumul->SetPointError(i, 0.0, ey );
        
        gv2withEP->GetPoint(i, x, y );
        ey = gv2withEP->GetErrorY(i);
        gv2ratioWithEP->SetPoint(i, x, y - in );
        gv2ratioWithEP->SetPointError(i, 0.0, ey );
        
        gv2withSP->GetPoint(i, x, y );
        ey = gv2withSP->GetErrorY(i);
        gv2ratioWithSP->SetPoint(i, x, y - in );
        gv2ratioWithSP->SetPointError(i, 0.0, ey );

        gv2Obs->GetPoint(i, x, y );
        ey = gv2Obs->GetErrorY(i);
        gv2ratioObs->SetPoint(i, x, y - in );
        gv2ratioObs->SetPointError(i, 0.0, ey );

        gv2withRealEP->GetPoint(i, x, y );
        ey = gv2withRealEP->GetErrorY(i);
        gv2ratioWithRealEP->SetPoint(i, x, y - in );
        gv2ratioWithRealEP->SetPointError(i, 0.0, ey );

        gv2withRealEPTrue->GetPoint(i, x, y );
        ey = gv2withRealEPTrue->GetErrorY(i);
        gv2ratioWithRealEPTrue->SetPoint(i, x, y - in );
        gv2ratioWithRealEPTrue->SetPointError(i, 0.0, ey );
    }

    double const ptBins[11] = {0.0,0.2,0.6,1.0,1.5,1.8,2.2,2.8,3.5,4.5,6.0};
    TF1 *fPtDep = new TF1("fPhiDistribution", v2PtDependence, 0.0, 6.0, 3);
    fPtDep->SetParameters(alpha,pTMax,inputFlow[1]);
    TGraphErrors *gv2RealEPPtDep = new TGraphErrors(nPtBins);
    TGraphErrors *gv2EPPtDep = new TGraphErrors(nPtBins);
    gv2RealEPPtDep->SetMarkerStyle(28);
    gv2RealEPPtDep->SetMarkerColor(1);
    gv2RealEPPtDep->SetMarkerSize(0.8);
    gv2EPPtDep->SetMarkerStyle(28);
    gv2EPPtDep->SetMarkerColor(2);
    gv2EPPtDep->SetMarkerSize(0.8);
    for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        double middlepoint = (ptBins[iPtBin]+ptBins[iPtBin+1])/2.0;
        double binWidth = middlepoint-ptBins[iPtBin];
        gv2RealEPPtDep->SetPoint(iPtBin,middlepoint,meanvRealEPTruePtBins[1][iPtBin]); //only v2, which is i=1
        gv2RealEPPtDep->SetPointError(iPtBin,binWidth,meanErrorvRealEPTruePtBins[1][iPtBin]);
        gv2EPPtDep->SetPoint(iPtBin,middlepoint,vEPPtBins[1][iPtBin]); //only v2, which is i=1
        gv2EPPtDep->SetPointError(iPtBin,binWidth,vEPErrorPtBins[1][iPtBin]);
    }

    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(0);

    TH2F *hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -0.2, 6.0, 1, 0.0, 0.5);
    myhset( *hfr, "p_{T}", "v_{2}",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    gv2RealEPPtDep->Draw("same,p");
    gv2EPPtDep->Draw("same,p");
    fPtDep->Draw("same");
    
    
    // ================ FIGURE ==================
    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(0);
    //mpad->SetGridx(0); mpad->SetGridy(0);
    
    int iMiddleBin = int( hEta[iType]->GetNbinsX()/2 );
    double etaAtZero = hEta[iType]->GetBinContent(iMiddleBin);
    hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -1.0, 1.0, 10, etaAtZero-50., etaAtZero+50.);
    myhset( *hfr, "#eta", "1/N_{ev} dN/d#eta",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    hEta[iType]->Draw("same,p");
    
    TLegend *leg = new TLegend(0.20, 0.15, 0.80, 0.45,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->AddEntry(hEta[iType],"toy MC","p");
    leg->Draw();
    
    
    // ================ FIGURE ==================
    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(1);
    //mpad->SetGridx(0); mpad->SetGridy(0);
    
    hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, 0.0, 6.0, 10, 5.e-7, 5.e4);
    myhset( *hfr, "p_{T}  [GeV/c]", "1/N_{ev} dN/dp_{T}",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    hPt[iType]->Draw("same,p");
    fInputPt->Draw("same,l");
    
    leg = new TLegend(0.20, 0.15, 0.80, 0.45,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->AddEntry(hPt[iType],"toy MC","p");
    leg->AddEntry(hPt[iType],"Input","l");
    leg->Draw();
    
    
    // ================ FIGURE ==================
    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(0);
    //mpad->SetGridx(0); mpad->SetGridy(0);
    
    hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, -3.2, 3.2, 10, 0.0,1.1*fPhiInput->Eval(0.0));
    myhset( *hfr, "#phi  [rad/#pi]", "1/N_{ev} dN/d#phi",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    hPhi[iType]->Draw("same,p");
    fitPhi->Draw("same");
    fPhiInput->Draw("same");
    fAcceptanceFunc->Draw("same");
    
    leg = new TLegend(0.20, 0.15, 0.80, 0.45,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->AddEntry(hPhi[iType],"toy MC","p");
    leg->AddEntry(fPhiInput,Form("input v_{n}'s: %3.2f, %3.2f, %3.2f, %3.2f, %3.2f",inputFlow[0],inputFlow[1],inputFlow[2],inputFlow[3],inputFlow[4]),"l");
    leg->AddEntry(fAcceptanceFunc,"Acceptance function","l");
    leg->AddEntry(fitPhi,"fit, resulting parameters:","l");
    leg->AddEntry(fitPhi,Form("v_{1} = %4.3f +- %4.3f",fitPhi->GetParameter(1), TMath::Max(0.001,fitPhi->GetParError(1))),"");
    leg->AddEntry(fitPhi,Form("v_{2} = %4.3f +- %4.3f",fitPhi->GetParameter(2), TMath::Max(0.001,fitPhi->GetParError(2))),"");
    leg->AddEntry(fitPhi,Form("v_{3} = %4.3f +- %4.3f",fitPhi->GetParameter(3), TMath::Max(0.001,fitPhi->GetParError(3))),"");
    leg->AddEntry(fitPhi,Form("v_{4} = %4.3f +- %4.3f",fitPhi->GetParameter(4), TMath::Max(0.001,fitPhi->GetParError(4))),"");
    leg->AddEntry(fitPhi,Form("v_{5} = %4.3f +- %4.3f",fitPhi->GetParameter(5), TMath::Max(0.001,fitPhi->GetParError(5))),"");
    leg->Draw();
    
    
    // ================ FIGURE ==================
    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(0);
    //mpad->SetGridx(0); mpad->SetGridy(0);
    
    hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, 0.5, 5.5, 1, 0.0, 0.25);
    myhset( *hfr, "n", "v_{n}",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    hFlowIn->Draw("same");
    gv2with1cumul->Draw("same,p");
    gv2with4cumul->Draw("same,p");
    gv2with6cumul->Draw("same,p");
    gv2withEP->Draw("same,p");
    gv2withSP->Draw("same,p");
    gv2Obs->Draw("same,p");
    gv2withRealEP->Draw("same,p");
    gv2withRealEPTrue->Draw("same,p");
    
    leg = new TLegend(0.20 ,0.70,0.80,0.90,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->AddEntry(hFlowIn,"input flow","l");
    leg->AddEntry(gv2with1cumul,"mean of Q","p");
    leg->AddEntry(gv2with4cumul,"4-ple cumulants","p");
    leg->AddEntry(gv2with6cumul,"6-ple cumulants","p");
    leg->AddEntry(gv2withEP,"event plane method","p");
    leg->AddEntry(gv2withSP,"scalar product method","p");
    leg->AddEntry(gv2Obs,"real event plane method, v_{obs}","p");
    leg->AddEntry(gv2withRealEP,"real event plane method","p");
    leg->AddEntry(gv2withRealEPTrue,"real event plane method w/ true R","p");
    leg->Draw();
    
    
    // ================ FIGURE ==================
    
    TLine *lineZero = new TLine(0.5, 0.0, 5.5, 0.0);
    lineZero->SetLineStyle(2);
    lineZero->SetLineWidth(1);
    
    mc(iFig++);
    gStyle->SetOptStat(0); gStyle->SetOptFit(0);  gStyle->SetOptTitle(0);
    gPad->SetLogx(0); gPad->SetLogy(0);
    //mpad->SetGridx(0); mpad->SetGridy(0);
    
    hfr = new TH2F(Form("hfr%02d",iFig)," ", 1, 0.5, 5.5, 1, -0.05, 0.1);
    myhset( *hfr, "n", "measured - input",0.9,1.4, 0.06,0.05, 0.01,0.001, 0.04,0.05, 510,510);
    hfr->Draw();
    
    gv2ratioWith1cumul->Draw("same,p");
    gv2ratioWith4cumul->Draw("same,p");
    gv2ratioWith6cumul->Draw("same,p");
    gv2ratioWithEP->Draw("same,p");
    gv2ratioWithSP->Draw("same,p");
    gv2ratioObs->Draw("same,p");
    gv2ratioWithRealEP->Draw("same,p");
    gv2ratioWithRealEPTrue->Draw("same,p");
    lineZero->Draw("same");
    
    leg = new TLegend(0.20 ,0.70,0.80,0.90,"","brNDC");
    leg->SetTextSize(0.037);leg->SetBorderSize(0);leg->SetFillStyle(0);
    leg->AddEntry(gv2ratioWith1cumul,"mean of Q","p");
    leg->AddEntry(gv2ratioWith4cumul,"4-ple cumulants","p");
    leg->AddEntry(gv2ratioWith6cumul,"6-ple cumulants","p");
    leg->AddEntry(gv2ratioWithEP,"event plane method","p");
    leg->AddEntry(gv2ratioWithSP,"scalar product method","p");
    leg->AddEntry(gv2ratioObs,"real event plane method, v_{obs}","p");
    leg->AddEntry(gv2ratioWithRealEP,"real event plane method","p");
    leg->AddEntry(gv2ratioWithRealEPTrue,"real event plane method w/ true R","p");
    leg->Draw();
    

    // ================ UNTUNED FIGURES FOR CHECKS ONLY ==================
    if(bCheckAllHistos){
        mc(iFig++);
        hMultiplicity[iType]->Draw();
        mc(iFig++);
        hSqrtSumWeights[iType]->Draw();
        mc(iFig++);
        hSqrtSumWeightsA[iType]->Draw();
        mc(iFig++);
        hSqrtSumWeightsB[iType]->Draw();
        for(int i = 0; i < nCoef; i++){
            mc(iFig++);
            hQx[iType][i]->Draw();
            mc(iFig++);
            hQy[iType][i]->Draw();
            mc(iFig++);
            hQ[iType][i]->Draw();
            mc(iFig++);
            hQ2[iType][i]->Draw();
            mc(iFig++);
            hQ4[iType][i]->Draw();
            mc(iFig++);
            hQ6[iType][i]->Draw();
            mc(iFig++);
            hSPnominator[iType][i]->Draw();
            mc(iFig++);
            hSPdenominator[iType][i]->Draw();
            mc(iFig++);
            hEPnominator[iType][i]->Draw();
            mc(iFig++);
            hEPdenominator[iType][i]->Draw();
        }
    }
}

double CalculateV2Using4Cumul( double Q4, double Q2, double weight ){
    double arg = TMath::Abs( 2.*Q2*Q2 - Q4 );
    return TMath::Power( arg , 1./4 ) / weight;
}


double CalculateV2Using6Cumul( double Q6, double Q4, double Q2, double weight ){
    double arg = TMath::Abs( Q6 - 9.*Q4*Q2 + 12.*Q2*Q2*Q2 );
    return TMath::Power( arg / 4. , 1./6 ) / weight;
}


double CalculateV2UsingEPorSP( double up, double down, double weight ){
    return up / TMath::Sqrt( down ) / weight;
}

double CalculateV2ErrorUsing4Cumul( double Q4, double Q2, double eQ4, double eQ2, double weight, double weightError ){
    double vn = CalculateV2Using4Cumul(Q4, Q2, weight);
    return vn * ( eQ4/Q4 + 2.*eQ2/Q2 + weightError/weight );
}

double CalculateV2ErrorUsing6Cumul( double Q6, double Q4, double Q2, double eQ6, double eQ4, double eQ2, double weight, double weightError ){
    double vn = CalculateV2Using6Cumul(Q6, Q4, Q2, weight);
    return vn * ( eQ6/Q6 + eQ4/Q4 + 4.*eQ2/Q2 + weightError/weight );
}


double CalculateV2ErrorUsingEPorSP( double up, double down, double eUp, double eDown, double weight, double weightError ){
    double vn = CalculateV2UsingEPorSP(up,down,weight);
    return vn * ( eUp/up + 0.5*eDown/down + weightError/weight );
}

double fFitMeanQ( double *x, double *p ){
    double Q = x[0];
    double norm = p[0]; // should be order 1
    double meanQ = p[1];
    double width = p[2];
    double term1 = norm*2.*Q/width/width;
    double term2 = TMath::Exp( -1. * ( Q*Q + meanQ*meanQ ) / width/width );
    double term3 = TMath::BesselI0( 2.*Q*meanQ/width/width );
    return term1*term2*term3;
}

double CalculateResolutionKOne( double Xi ){
    double term1 = sqrt(pi)/2.0 * Xi;
    double term2 = TMath::Exp( -Xi*Xi/2.0 );
    double term3 = TMath::BesselI0( Xi*Xi/2.0 ) + TMath::BesselI1( Xi*Xi/2.0 );
    return term1*term2*term3;
}

void checkUnderOverFlow( TH1 *h ){
        if(h->GetBinContent(0)>0) cout << h->GetName() << " underflow bin not empty: " << h->GetBinContent(0) << endl;
        if(h->GetBinContent(h->GetXaxis()->GetNbins()+1)>0) cout << h->GetName() << " overflow bin not empty: " << h->GetBinContent(h->GetXaxis()->GetNbins()+1) << endl;
}

double AcceptanceFunc(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double effA = p[1];
    double effB = p[2];
    double detMaxA = p[3];
    double detMinA = p[4];
    double detMaxB = p[5];
    double detMinB = p[6];
    double value = 0.0;
    if(isInAcc(phi,detMaxA,detMinA)) {
        value = scaling*effA/(2.0*TMath::Pi());
    } else if(isInAcc(phi,detMaxB,detMinB)) {
        value = scaling*effB/(2.0*TMath::Pi());
    }
    return value;
}

// Test if in det acceptance:
bool isInAcc(double phi, double detMax, double detMin) {
    return phi <= detMax && phi > detMin;
}


double v2PtDependence(double *x, double *p){
    double pT = x[0];
    double alpha = p[0];
    double pTmax = p[1];
    double v2Max = p[2];
    return v2Max*v2PtDependenceFun(pT,pTmax,alpha);
}


// This is just a lookalike model, not based on any hard evidence or model.
double v2PtDependenceFun(double pT, double pTmax, double alpha){
    return TMath::Power(pT/pTmax,alpha)*TMath::Exp(-alpha*(pT/pTmax - 1.0));
}

