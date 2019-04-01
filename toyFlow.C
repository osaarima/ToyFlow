/*
 *  toyMCforHT.C
 *
 *  Created by Sami Rasanen on 2/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>

// OWN
#include "JToyMCTrack.h"  // TLorentzVector + isHT, isHard, isSoft, isCharged, isIsolated, isLeading + ID number
#include "JHistos.h"      // encapsulates all histograms

// ROOT 
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TClonesArray.h"

// MISC
#include "TStopwatch.h"

using namespace std;

struct Qvalues{
    double Qx;
    double Qy;
    double Q;
    double Q2;
    double Q4;
    double Q6;
    double psi;
    double vObs;
};

struct EffValues{
    double cBar;
    double sBar;
    double lambda2nSMinus;
    double lambda2nSPlus;
    double a2nPlus;
    double a2nMinus;
};

void GetEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, TClonesArray **listAllPtBins, int nPtBins, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, int nMult, JHistos *histos, bool bv2PtDep);
void GetDetectorParticles(TClonesArray *listAll, TClonesArray *listAllDet, TClonesArray *listSubADet, TClonesArray *listSubBDet, double detAMax, double detAMin, double detBMax, double detBMin, double detAEff, double detBEff, TRandom3 *rand, JHistos *histos);
void AnalyzeEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, TClonesArray **listAllPtBins, double *truePsi, bool bUseWeightning, bool bCorr, JHistos *histos, EffValues **effValues, int const nHisto);
double CalculateCumulants(TClonesArray *list, int n, bool bUseWeightning, Qvalues *qval, EffValues *effValues, JHistos *histos, int const nHisto, bool bCorr);
double CalculateWeight(JToyMCTrack *track, bool bUseWeightning);
double CalculateDotProduct(double ax, double ay, double bx, double by);
double SingleParticlePt(double *x, double *para);
double v2PtDependence(double *x, double *p);
double v2PtDependenceFun(double pT, double pTmax, double alpha);
double SingleParticlePhi(double *x, double *para);
double AcceptanceFunc(double *x, double *p);
double AcceptanceFuncSin(double *x, double *p);
double AcceptanceFuncCos(double *x, double *p);
bool isInAcc(double phi, double detMax, double detMin);
double DeltaPhi(double phi1, double phi2);
double calculatePsi(double Qy, double Qx, double n);
void efficiencyCalc(EffValues *effValues, int iHarm, double detAMax, double detAMin, double detBMax, double detBMin, double detAEff, double detBEff);
int getPtBin(double pT);



// ==============
// Main program
// ==============
int main(int argc, char **argv) {
	
    const int nEvents = argc>1 ? atoi(argv[1]) : 100;
    const double dNdeta = argc>2 ? atoi(argv[2]) : 1000;
    const double usePtDep = argc>3 ? atoi(argv[3]) : 0;
    const double useWeighting = argc>4 ? atoi(argv[4]) : 0;
    TString sFileText = argc>5 ? argv[5] : "Test";
    
    const double etaRange = 0.8;
    
    const int nMult = int( 2. * etaRange * dNdeta );
    cout << "nEvents = " << nEvents << ", dNdeta = " << dNdeta << " and total multiplicity = " << nMult << endl;
    
    const int nHarmonics = 5;
    const int nPtBins = 10;
    const double vn[nHarmonics] = {0.0, 0.15, 0.08, 0.03, 0.01}; // Peripheral
    //const double vn[nHarmonics] = {0.0, 0.01, 0.00, 0.00, 0.00}; // Very central
    
    // Effieciency
    //double const detAMax = 2.0*TMath::Pi()/3.0, detAMin = TMath::Pi()/3.0;
    //double const detBMax = 0.0,                 detBMin = -TMath::Pi()/2.0;
    //double const detAEff = 0.4,                 detBEff = 0.85;
    double const detAMax = TMath::Pi(),         detAMin = 2.0;
    double const detBMax = 2.0,                 detBMin = -TMath::Pi();
    double const detAEff = 0.6,                 detBEff = 1.0;

    
    bool bRandomPsi = true;
    bool bUseWeightning = useWeighting;
    bool bv2PtDep = usePtDep;
    cout << "bRandomPsi=" << bRandomPsi << ", bUseWeightning=" << bUseWeightning << ", bv2PtDep=" << bv2PtDep << endl;
    double alpha = 2.0;
    double defpTMax = 2.0;
    double const pi = TMath::Pi();
    
	TStopwatch timer;
	timer.Start();

    TString sUseWeightning;
    if(bUseWeightning) sUseWeightning="weight";
    else sUseWeightning="noWeight";

    TString sRandomPsi;
    if(bRandomPsi) sRandomPsi="randomPsi";
    else sRandomPsi="noRandomPsi";
	
    TString sv2PtDep;
    if(bv2PtDep) sv2PtDep="ptDep";
    else sv2PtDep="noptDep";
	
    TString outFileName = Form("toyFlow_%s_%s_%s_dNdeta-%.0f_nEvents-%d-%s.root",sUseWeightning.Data(),sRandomPsi.Data(),sv2PtDep.Data(),dNdeta,nEvents,sFileText.Data());
	TFile *fOut = TFile::Open( outFileName, "RECREATE" );
	
	TRandom3 *randomGenerator = new TRandom3();
	JHistos *histos = new JHistos();
    
    double Psi[nHarmonics] = {0};
    
    double Tdec = 0.12;
    double vr = 0.6;
    double Teff = Tdec * TMath::Sqrt( (1.+vr)/(1.-vr) );  //NOTE ======================================================================================
    cout << "Initial pT distribution : Tdec = " << Tdec << ", vr = " << vr << ", Teff = " << Teff << " and slope = 1/Teff = " << 1./Teff << endl;
    TF1 *fPtDistribution = new TF1("fPtDistribution", SingleParticlePt, 0.0, 10.0, 1);
    fPtDistribution->SetParameter(0,1./Teff);
    
    TF1 *fPhiDistribution = new TF1("fPhiDistribution", SingleParticlePhi, -pi, pi, 13);
    double params[13] = {vn[0], vn[1], vn[2], vn[3], vn[4], Psi[0], Psi[1], Psi[2], Psi[3], Psi[4], defpTMax,alpha,defpTMax};
    fPhiDistribution->SetParameters(params);
    
	TClonesArray *allHadrons = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventB = new TClonesArray("JToyMCTrack", nMult+1);
	TClonesArray *allHadronsDet = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventDetA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventDetB = new TClonesArray("JToyMCTrack", nMult+1);

	TClonesArray *allHadronsPtBins[nPtBins];
    for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
        allHadronsPtBins[iPtBin] = new TClonesArray("JToyMCTrack", nMult+1);
    }
	
    // save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",21, 0.5, 21.5);
    hInputNumbers->Fill(1, double(nEvents));
    hInputNumbers->Fill(2, double(dNdeta));
    hInputNumbers->Fill(3, etaRange);
    hInputNumbers->Fill(4, double(nMult));
    hInputNumbers->Fill(5, vn[0]);
    hInputNumbers->Fill(6, vn[1]);
    hInputNumbers->Fill(7, vn[2]);
    hInputNumbers->Fill(8, vn[3]);
    hInputNumbers->Fill(9, vn[4]);
    hInputNumbers->Fill(10, Tdec);
    hInputNumbers->Fill(11, vr);
    hInputNumbers->Fill(12, Teff);
    hInputNumbers->Fill(13, 1./Teff);
    hInputNumbers->Fill(14, detAMax);
    hInputNumbers->Fill(15, detAMin);
    hInputNumbers->Fill(16, detBMax);
    hInputNumbers->Fill(17, detBMin);
    hInputNumbers->Fill(18, detAEff);
    hInputNumbers->Fill(19, detBEff);
    hInputNumbers->Fill(20, alpha);
    hInputNumbers->Fill(21, defpTMax);

    // Calculate efficiency related values
    EffValues *effValues[nHarmonics+1];
    //cout << "effValues: " << effValues << endl;
    for (int iHarm=1; iHarm < nHarmonics+1; iHarm++) {
        effValues[iHarm] = new EffValues;
        efficiencyCalc(effValues[iHarm], iHarm, detAMax, detAMin, detBMax, detBMin, detAEff, detBEff);
        //cout << "effValues[" << iHarm << "]: " << effValues[iHarm] << endl;
    }
    
    int nOut = nEvents/20;
    for(int counter = 0; counter < nEvents; counter++){
        if(counter % nOut == 0) cout << 100*counter/nEvents << " % finished " << endl;
        /* Use Clear() or Clear("C") instead of Delete(). This will improve program execution time.
         *
         * TClonesArray object classes containing pointers allocate memory. To avoid causing memory leaks,
         * special Clear("C") must be used for clearing TClonesArray. When option "C" is specified, ROOT
         * automatically executes the Clear() method (by default it is empty contained in TObject). This
         * method must be overridden in the relevant TClonesArray object class, implementing the reset
         * procedure for pointer objects. */
        allHadrons->Clear("C");
        subEventA->Clear("C");
        subEventB->Clear("C");
        allHadronsDet->Clear("C");
        subEventDetA->Clear("C");
        subEventDetB->Clear("C");

        for(int iPtBin=0;iPtBin<nPtBins;iPtBin++) {
            allHadronsPtBins[iPtBin]->Clear("C");
        }
        
        if(bRandomPsi){
            for(int j = 0; j < 5; j++){
                Psi[j] = randomGenerator->Uniform(-pi,pi); // Should be pi/n_harmonic ?
                params[5+j] = Psi[j];
            }
            fPhiDistribution->SetParameters(params);
        }
        
        GetEvent(allHadrons, subEventA, subEventB, allHadronsPtBins, nPtBins, randomGenerator, fPtDistribution, fPhiDistribution, nMult, histos, bv2PtDep);
        GetDetectorParticles(allHadrons, allHadronsDet, subEventDetA, subEventDetB, 
                             detAMax, detAMin, detBMax, detBMin, detAEff, detBEff,
                             randomGenerator, histos);
        AnalyzeEvent(allHadrons, subEventA, subEventB, allHadronsPtBins, Psi, bUseWeightning, false, histos, effValues, 0); //True MC
        AnalyzeEvent(allHadronsDet, subEventDetA, subEventDetB, allHadronsPtBins, Psi, bUseWeightning, false, histos, effValues, 1); //Detector, no corr
        AnalyzeEvent(allHadronsDet, subEventDetA, subEventDetB, allHadronsPtBins, Psi, bUseWeightning, true, histos, effValues, 2); //Detector, corr
    }
                                    
	fOut->Write();
	timer.Print(); 
    //delete effValues;
	return 0;
}
                                    
                                    
// ============== END MAIN PROGRAM

void GetEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, TClonesArray **listAllPtBins, int nPtBins, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, int nMult, JHistos *histos, bool bv2PtDep){
    
    histos->hMultiplicity[0]->Fill(1.*nMult);
    
    JToyMCTrack track;
    TLorentzVector lVec;
    int iTrack[nPtBins] = {0};

    int nSub = int( nMult/2.0 );
    
    for(int i = 0; i < nMult; i++){
        double pT = fPt->GetRandom();
        if(bv2PtDep) fPhi->SetParameter(10,pT);
        double phi = fPhi->GetRandom();
        double eta = rand->Uniform(-0.8,0.8);
        
        histos->hPt[0]->Fill(pT);
        histos->hPhi[0]->Fill(phi);
        histos->hEta[0]->Fill(eta);
        histos->hPhiEta[0]->Fill(phi,eta);
        
        double px = pT*TMath::Cos(phi);
        double py = pT*TMath::Sin(phi);
        double pz = pT*TMath::SinH(eta);
        double  E = TMath::Sqrt(pT*pT+pz*pz); //massless
        lVec.SetPxPyPzE(px,py,pz,E);
        track.SetTrack( lVec, 0, 0, 0, 1, 0, 0, i); // Lorentz vector + isHT, isHard, isSoft, isCharged, isIsolated, isLeading, nID
        new((*listAll)[i]) JToyMCTrack( track );
        if( i < nSub ) {
            new((*listSubA)[i]) JToyMCTrack( track );
        } else {
            new((*listSubB)[i-nSub]) JToyMCTrack( track );
        }

        int nPtBin = getPtBin(pT);
        new((*listAllPtBins[nPtBin])[iTrack[nPtBin]++]) JToyMCTrack( track );
    }
    return;
}

void GetDetectorParticles(TClonesArray *listAll, 
                          TClonesArray *listAllDet, TClonesArray *listSubADet, TClonesArray *listSubBDet,
                          double detAMax, double detAMin, double detBMax, double detBMin, double detAEff, double detBEff,
                          TRandom3 *rand, JHistos *histos){

    int numberOfTracks = listAll->GetEntriesFast();

    double eta,phi,pT,px,py,pz,E;
    int detMultiplicity = 0;
    int iTrack=0;

    JToyMCTrack newTrack;
    TLorentzVector lVec;

    for(int i = 0; i < numberOfTracks; i++){
        
        JToyMCTrack *track = (JToyMCTrack*)listAll->At(i);
        
        pT = track->GetPt();
        phi = track->GetPhi();
        eta = track->GetEta();
        
        //newTrack.SetE(track.E());
        //newTrack = track;
        px = track->GetPx();
        py = track->GetPy();
        pz = track->GetPz();
        E = track->GetE(); //massless
        lVec.SetPxPyPzE(px,py,pz,E);
        newTrack.SetTrack( lVec, 0, 0, 0, 1, 0, 0, i); // Lorentz vector + isHT, isHard, isSoft, isCharged, isIsolated, isLeading, nID

        if((isInAcc(phi,detAMax,detAMin) && rand->Uniform(0.0,1.0) < detAEff) ||
           (isInAcc(phi,detBMax,detBMin) && rand->Uniform(0.0,1.0) < detBEff) ) {
            new((*listAllDet)[iTrack++]) JToyMCTrack( newTrack );
            detMultiplicity++;
            histos->hPt[1]->Fill(pT);
            histos->hPt[2]->Fill(pT);
            histos->hPhi[1]->Fill(phi);
            histos->hPhi[2]->Fill(phi);
            histos->hEta[1]->Fill(eta);
            histos->hEta[2]->Fill(eta);
            histos->hPhiEta[1]->Fill(phi,eta);
            histos->hPhiEta[2]->Fill(phi,eta);

            //newTrack.SetVect(track.Vect());
        }
    }
    histos->hMultiplicity[1]->Fill(1.*detMultiplicity);
    histos->hMultiplicity[2]->Fill(1.*detMultiplicity);

    numberOfTracks = listAllDet->GetEntriesFast();
    int nSub = int( numberOfTracks/2.0 );


    for(int i = 0; i < numberOfTracks; i++){

        JToyMCTrack *track = (JToyMCTrack*)listAllDet->At(i);

        px = track->GetPx();
        py = track->GetPy();
        pz = track->GetPz();
        E = track->GetE(); //massless
        lVec.SetPxPyPzE(px,py,pz,E);
        newTrack.SetTrack( lVec, 0, 0, 0, 1, 0, 0, i); // Lorentz vector + isHT, isHard, isSoft, isCharged, isIsolated, isLeading, nID
        if( i < nSub ) {
            new((*listSubADet)[i]) JToyMCTrack( newTrack );
        } else {
            new((*listSubBDet)[i-nSub]) JToyMCTrack( newTrack );
        }
    }


}



void AnalyzeEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, TClonesArray **listAllPtBins, double *truePsi, bool bUseWeightning, bool bCorr, JHistos *histos, EffValues **effValues, int const nHisto){
    double Qx, QxA, QxB, Qy, QyA, QyB, Q, QA, QB, Q2, Q2A, Q2B, Q4, Q4A, Q4B, Q6, Q6A, Q6B, psi, psiA, psiB, vObs, vObsA, vObsB, vnSPnominator,vnSPdenominator,vnEPnominator,vnEPdenominator;
    double sqrtSumWeights = 0.0, sqrtSumWeightsA = 0.0, sqrtSumWeightsB = 0.0;

    Qvalues *qval = new Qvalues;

    for(int i = 0; i < 5; i++){


        int n = i+1; // which harmonic

        sqrtSumWeights  = CalculateCumulants(listAll, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        Qx = qval->Qx; Qy = qval->Qy; Q = qval->Q; Q2 = qval->Q2; Q4 = qval->Q4; Q6 = qval->Q6; psi = qval->psi; vObs = qval->vObs;

        sqrtSumWeightsA  = CalculateCumulants(listSubA, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        QxA = qval->Qx; QyA = qval->Qy; QA = qval->Q; Q2A = qval->Q2; Q4A = qval->Q4; Q6A = qval->Q6; psiA = qval->psi; vObsA = qval->vObs;

        sqrtSumWeightsB  = CalculateCumulants(listSubB, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        QxB = qval->Qx; QyB = qval->Qy; QB = qval->Q; Q2B = qval->Q2; Q4B = qval->Q4; Q6B = qval->Q6; psiB = qval->psi; vObsB = qval->vObs;

        double rSub = TMath::Cos(n*(psiA-psiB));

        histos->hQx[nHisto][i]->Fill(Qx);
        histos->hQy[nHisto][i]->Fill(Qy);
        histos->hQ[nHisto][i]->Fill(Q);
        histos->hQ2[nHisto][i]->Fill(Q2);
        histos->hQ4[nHisto][i]->Fill(Q4);
        histos->hQ6[nHisto][i]->Fill(Q6);
        histos->hrSub[nHisto][i]->Fill(rSub);
        histos->hPsiAB[nHisto][i][0]->Fill(n*psiA,n*psiB);
        histos->hPsiDiff[nHisto][i][0]->Fill(psiA-psiB);
        histos->hPsiDiffN[nHisto][i][0]->Fill(n*(psiA-psiB));
        if(rSub<0.0) {
            histos->hPsiAB[nHisto][i][1]->Fill(n*psiA,n*psiB);
            histos->hPsiDiff[nHisto][i][1]->Fill(psiA-psiB);
            histos->hPsiDiffN[nHisto][i][1]->Fill(n*(psiA-psiB));
        }
        histos->hTrueReso[nHisto][i]->Fill(TMath::Cos(n*(psi-truePsi[i]))); //Note: truePsi goes from 0 to 4
        histos->hvObs[nHisto][i]->Fill(vObs);


        vnSPnominator = CalculateDotProduct(Qx, Qy, QxA, QyA);
        vnSPdenominator = CalculateDotProduct(QxA, QyA, QxB, QyB);

        vnEPnominator = vnSPnominator / QA;
        vnEPdenominator = vnSPdenominator / QA / QB;

        histos->hSPnominator[nHisto][i]->Fill(vnSPnominator);
        histos->hSPdenominator[nHisto][i]->Fill(vnSPdenominator);
        histos->hEPnominator[nHisto][i]->Fill(vnEPnominator);
        histos->hEPdenominator[nHisto][i]->Fill(vnEPdenominator);

        if(nHisto==0) { //Look only trueMC for now.
            for(int iPtBin=0;iPtBin<10;iPtBin++) {
                if(listAllPtBins[iPtBin]->GetEntriesFast()<3) continue; // If less than 2 particles, vObs will be nan.
                sqrtSumWeights  = CalculateCumulants(listAllPtBins[iPtBin], n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
                psi = qval->psi; vObs = qval->vObs;
                //cout << "n=" << n <<", iPtBin=" << iPtBin << ", psi=" << psi << ", vObs=" << vObs << endl;
                histos->hTrueResoPtBins[nHisto][i][iPtBin]->Fill(TMath::Cos(n*(psi-truePsi[i]))); //Note: truePsi goes from 0 to 4
                histos->hvObsPtBins[nHisto][i][iPtBin]->Fill(vObs);
            }
        }
    }

    histos->hSqrtSumWeights[nHisto]->Fill( sqrtSumWeights );
    histos->hSqrtSumWeightsA[nHisto]->Fill( sqrtSumWeightsA );
    histos->hSqrtSumWeightsB[nHisto]->Fill( sqrtSumWeightsB );

    delete qval;

    return;
}

double CalculateCumulants(TClonesArray *list, int n, bool bUseWeightning, Qvalues *qval, EffValues *effValues, JHistos *histos, int const nHisto, bool bCorr){

    int numberOfTracks = list->GetEntriesFast();
    //cout << "numberOfTracks=" << numberOfTracks << endl;

    double normalization = 0.0;
    double Qx = 0.0, Qy = 0.0;
    double xx = 0.0, yy = 0.0;
    // For corrections:
    double QxPrime, QyPrime, QxPPrime, QyPPrime;
    double weight, phi;

    //Note in AliRoot for eff corrections:
    // if (harmonic == 2) {

    for(int i = 0; i < numberOfTracks; i++){

        JToyMCTrack *track = (JToyMCTrack*)list->At(i);

        weight = CalculateWeight(track, bUseWeightning);
        normalization += weight*weight;

        phi = track->GetPhi();

        xx = TMath::Cos(n*phi);
        yy = TMath::Sin(n*phi);

        if(bCorr) {
            //From DOI: 10.1103/PhysRevC.77.034904
            //Recentering
            QxPrime = xx - effValues->cBar;
            QyPrime = yy - effValues->sBar;

            //Twisting
            QxPPrime = (QxPrime - effValues->lambda2nSMinus * QyPrime)/(1 - effValues->lambda2nSMinus * effValues->lambda2nSPlus);
            QyPPrime = (QyPrime - effValues->lambda2nSPlus * QxPrime)/(1 - effValues->lambda2nSMinus * effValues->lambda2nSPlus);

            //Rescaling
            xx = QxPPrime / effValues->a2nPlus;
            yy = QyPPrime / effValues->a2nMinus;
        }


        Qx += weight * xx;
        Qy += weight * yy;
    }

    double helperNormalization = normalization;
    normalization = TMath::Sqrt( normalization );

    double helperQx = Qx; double helperQy = Qy;
    Qx /= normalization; Qy /= normalization;

    double Q2 = Qx*Qx + Qy*Qy;
    double Q = TMath::Sqrt( Q2 );
    double Q4 = Q2*Q2;
    double Q6 = Q4*Q2;

    double psi = calculatePsi(Qy,Qx,n);
    double QyMod, QxMod;
    double vObs = 0.0;
    double helperPsi = 0.0;

    for(int i = 0; i < numberOfTracks; i++){
        JToyMCTrack *track = (JToyMCTrack*)list->At(i);
        weight = CalculateWeight(track, bUseWeightning);
        phi = track->GetPhi();
        xx = weight*TMath::Cos(n*phi);
        yy = weight*TMath::Sin(n*phi);
        QxMod = (helperQx-xx)/sqrt(helperNormalization-weight*weight);
        QyMod = (helperQy-yy)/sqrt(helperNormalization-weight*weight);
        helperPsi = calculatePsi(QyMod, QxMod ,n); // Remove autocorrelations.
        vObs += TMath::Cos(n*(phi-helperPsi));
    }

    (*qval).Qx = Qx;
    (*qval).Qy = Qy;
    (*qval).Q  = Q;
    (*qval).Q2 = Q2;
    (*qval).Q4 = Q4;
    (*qval).Q6 = Q6;
    (*qval).psi = psi;
    (*qval).vObs = vObs/((double)numberOfTracks);

    return normalization;
}

double CalculateWeight(JToyMCTrack *track, bool bUseWeightning){
    double weight;
    if(bUseWeightning) weight = track->GetPt();
    else weight = 1.;
    return weight;
}

double CalculateDotProduct(double ax, double ay, double bx, double by){
    return ax*bx + ay*by;
}


double SingleParticlePt(double *x, double *p){
    double pT = x[0];
    double slope = p[0];
    return TMath::Exp(-slope*pT);
}

// This is just a lookalike model, not based on any hard evidence or model.
double v2PtDependenceFun(double pT, double pTmax, double alpha){
    return TMath::Power(pT/pTmax,alpha)*TMath::Exp(-alpha*(pT/pTmax - 1.0));
}

double SingleParticlePhi(double *x, double *p){
    double phi   = x[0];
    double v1    = p[0];
    double v2max = p[1];
    double v3    = p[2];
    double v4    = p[3];
    double v5    = p[4];
    double Psi1  = p[5];
    double Psi2  = p[6];
    double Psi3  = p[7];
    double Psi4  = p[8];
    double Psi5  = p[9];
    double pT    = p[10];
    double alpha = p[11];
    double pTmax = p[12];

    return 1.0 + 2.*v1*TMath::Cos(phi-Psi1) + 2.*v2max*v2PtDependenceFun(pT,pTmax,alpha)*TMath::Cos(2.*(phi-Psi2))
        + 2.*v3*TMath::Cos(3.*(phi-Psi3)) + 2.*v4*TMath::Cos(4.*(phi-Psi4))
        + 2.*v5*TMath::Cos(5.*(phi-Psi5));
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

double AcceptanceFuncCos(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double effA = p[1];
    double effB = p[2];
    double detMaxA = p[3];
    double detMinA = p[4];
    double detMaxB = p[5];
    double detMinB = p[6];
    double multiply = p[7];
    double value = 0.0;
    if(isInAcc(phi,detMaxA,detMinA)) {
        value = scaling*effA/(2.0*TMath::Pi()) * TMath::Cos(multiply*phi);
    } else if(isInAcc(phi,detMaxB,detMinB)) {
        value = scaling*effB/(2.0*TMath::Pi()) * TMath::Cos(multiply*phi);
    }
    return value;
}

double AcceptanceFuncSin(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double effA = p[1];
    double effB = p[2];
    double detMaxA = p[3];
    double detMinA = p[4];
    double detMaxB = p[5];
    double detMinB = p[6];
    double multiply = p[7];
    double value = 0.0;
    if(isInAcc(phi,detMaxA,detMinA)) {
        value = scaling*effA/(2.0*TMath::Pi()) * TMath::Sin(multiply*phi);
    } else if(isInAcc(phi,detMaxB,detMinB)) {
        value = scaling*effB/(2.0*TMath::Pi()) * TMath::Sin(multiply*phi);
    }
    return value;
}


// Test if in det acceptance:
bool isInAcc(double phi, double detMax, double detMin) {
    return phi <= detMax && phi > detMin;
}

// Delta_phi on correlation function
double DeltaPhi(double phi1, double phi2) {
    double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
    return res>-TMath::Pi()/3.0 ? res : 2.0*TMath::Pi()+res ; 
}

double calculatePsi(double Qy, double Qx, double n) {
    return TMath::ATan2(Qy,Qx)/n;
}

void efficiencyCalc(EffValues *effValues, int iHarm, double detAMax, double detAMin, double detBMax, double detBMin, double detAEff, double detBEff) {

    double const pi = TMath::Pi();
    // u_n = x_n + iy_n ≡ cos nφ + isin nφ
    //
    // a_2n^+-= 1 ± \bar{c_2n} = 1 ± \bar{cos 2nφ}
    //
    // 19-02-- paper ref 31:
    // Let us denote by A(φ) the probability that a particle
    // with azimuthal angle φ be detected: A(φ) represents the
    // acceptance-efficiency profile of the detector. We choose
    // the normalization Integral A(φ) dφ/2π = 1.
    //
    // \bar{f} ≡  Integral A(φ) f(φ) dφ/2π
    //
    // 
    //
    //
    // /home/alidock/alice/AliRoot/STEER/STEERBase/AliEventplane.cxx
    // 
    //
    //
    TF1 *fAcceptanceFunc;
    TF1 *fAcceptanceFuncSin;
    TF1 *fAcceptanceFuncCos;
    double integral=0.0, cBar=0.0, sBar=0.0, c2nBar=0.0, s2nBar=0.0, a2nPlus=0.0, a2nMinus=0.0, lambda2nCPlus=0.0, lambda2nCMinus=0.0, lambda2nSPlus=0.0, lambda2nSMinus=0.0;
    fAcceptanceFunc    = new TF1("fAcceptanceFunc", AcceptanceFunc, -pi, pi, 10);
    fAcceptanceFuncSin = new TF1("fAcceptanceFuncSin", AcceptanceFuncSin, -pi, pi, 10);
    fAcceptanceFuncCos = new TF1("fAcceptanceFuncCos", AcceptanceFuncCos, -pi, pi, 10);
    fAcceptanceFunc->SetParameter(0,1.0);
    fAcceptanceFuncSin->SetParameter(0,1.0);
    fAcceptanceFuncCos->SetParameter(0,1.0);
    fAcceptanceFunc->SetParameter(1,detAEff);
    fAcceptanceFuncSin->SetParameter(1,detAEff);
    fAcceptanceFuncCos->SetParameter(1,detAEff);
    fAcceptanceFunc->SetParameter(2,detBEff);
    fAcceptanceFuncSin->SetParameter(2,detBEff);
    fAcceptanceFuncCos->SetParameter(2,detBEff);
    fAcceptanceFunc->SetParameter(3,detAMax);
    fAcceptanceFuncSin->SetParameter(3,detAMax);
    fAcceptanceFuncCos->SetParameter(3,detAMax);
    fAcceptanceFunc->SetParameter(4,detAMin);
    fAcceptanceFuncSin->SetParameter(4,detAMin);
    fAcceptanceFuncCos->SetParameter(4,detAMin);
    fAcceptanceFunc->SetParameter(5,detBMax);
    fAcceptanceFuncSin->SetParameter(5,detBMax);
    fAcceptanceFuncCos->SetParameter(5,detBMax);
    fAcceptanceFunc->SetParameter(6,detBMin);
    fAcceptanceFuncSin->SetParameter(6,detBMin);
    fAcceptanceFuncCos->SetParameter(6,detBMin);


    integral = fAcceptanceFunc->Integral(-pi,pi);
    //cout << "Acc func integral: " << integral << endl;
    fAcceptanceFunc->SetParameter(0,1.0/integral);
    fAcceptanceFuncSin->SetParameter(0,1.0/integral);
    fAcceptanceFuncCos->SetParameter(0,1.0/integral);
    //integral = fAcceptanceFunc->Integral(-pi,pi);
    //cout << "Acc func integral after scaling: " << integral << endl;

    fAcceptanceFuncSin->SetParameter(7,(double)iHarm);
    fAcceptanceFuncCos->SetParameter(7,(double)iHarm);
    sBar = fAcceptanceFuncSin->Integral(-pi,pi);
    cBar = fAcceptanceFuncCos->Integral(-pi,pi);
    //cout << "sBar=" << sBar << ", cBar=" << cBar << endl;

    fAcceptanceFuncSin->SetParameter(7,2.0*(double)iHarm);
    fAcceptanceFuncCos->SetParameter(7,2.0*(double)iHarm);
    s2nBar = fAcceptanceFuncSin->Integral(-pi,pi);
    c2nBar = fAcceptanceFuncCos->Integral(-pi,pi);

    a2nPlus  = 1.0 + c2nBar;
    a2nMinus = 1.0 - c2nBar;

    lambda2nCPlus  = c2nBar/a2nPlus;
    lambda2nCMinus = c2nBar/a2nMinus;
    lambda2nSPlus  = s2nBar/a2nPlus;
    lambda2nSMinus = s2nBar/a2nMinus;
    delete fAcceptanceFunc;
    fAcceptanceFunc = NULL;
    delete fAcceptanceFuncSin;
    fAcceptanceFuncSin = NULL;
    delete fAcceptanceFuncCos;
    fAcceptanceFuncCos = NULL;

    (*effValues).cBar = cBar;
    (*effValues).sBar = sBar;
    (*effValues).lambda2nSMinus = lambda2nSMinus;
    (*effValues).lambda2nSPlus = lambda2nSPlus;
    (*effValues).a2nPlus = a2nPlus;
    (*effValues).a2nMinus = a2nMinus;

    cout << "==================== For harmonic " << iHarm << " ====================" << endl;
    cout << "cBar           = " << cBar <<           ", sBar          = " << sBar << endl;
    cout << "lambda2nSMinus = " << lambda2nSMinus << ", lambda2nSPlus = " << lambda2nSPlus << endl;
    cout << "a2nPlus        = " << a2nPlus <<        ", a2nMinus      = " << a2nMinus << endl;
    
}

int getPtBin(double pT) {
    double const ptBins[11] = {0.0,0.2,0.6,1.0,1.5,1.8,2.2,2.8,3.5,4.5,6.0};
    int nBin;
    for(int iPtBin=0;iPtBin<10;iPtBin++) {
        if(pT>ptBins[iPtBin] && pT<ptBins[iPtBin+1]) nBin = iPtBin;
    }
    return nBin;
}

