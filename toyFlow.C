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
};

struct EffValues{
    double cBar;
    double sBar;
    double lambda2nSMinus;
    double lambda2nSPlus;
    double a2nPlus;
    double a2nMinus;
};

void GetEvent(TClonesArray *allHadrons, TClonesArray *subEventA, TClonesArray *subEventB, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, int nMulti, JHistos *histos);
void GetDetectorParticles(TClonesArray *listAll, TClonesArray *listAllDet, TClonesArray *listSubADet, TClonesArray *listSubBDet, TRandom3 *rand, JHistos *histos);
void AnalyzeEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, bool bUseWeightning, bool bCorr, JHistos *histos, EffValues **effValues, int const nHisto);
double CalculateCumulants(TClonesArray *list, int n, bool bUseWeightning, Qvalues *qval, EffValues *effValues, JHistos *histos, int const nHisto, bool bCorr);
double CalculateWeight(JToyMCTrack *track, bool bUseWeightning);
double CalculateDotProduct(double ax, double ay, double bx, double by);
double SingleParticlePt(double *x, double *para);
double SingleParticlePhi(double *x, double *para);
double AcceptanceFunc(double *x, double *p);
double AcceptanceFuncSin(double *x, double *p);
double AcceptanceFuncCos(double *x, double *p);
double DeltaPhi(double phi1, double phi2);
void efficiencyCalc(EffValues *effValues, int iHarm);

double const detAMax = 2.0*TMath::Pi()/3.0, detAMin = TMath::Pi()/3.0;
double const detBMax = 0.0,                 detBMin = -TMath::Pi()/2.0;
double const detAEff = 0.4,                 detBEff = 0.85;



// ==============
// Main program
// ==============
int main(int argc, char **argv) {
	
    const int nEvents = argc>1 ? atoi(argv[1]) : 100;
    const double dNdeta = argc>2 ? atoi(argv[2]) : 1000;
    TString sFileText = argc>3 ? argv[3] : "Test";
    
    const double etaRange = 0.8;
    
    const int nMult = int( 2. * etaRange * dNdeta );
    cout << "nEvents = " << nEvents << ", dNdeta = " << dNdeta << " and total multiplicity = " << nMult << endl;
    
    const int nHarmonics = 5;
    const double vn[nHarmonics] = {0.0, 0.15, 0.08, 0.03, 0.01}; // Peripheral
    //const double vn[nHarmonics] = {0.0, 0.01, 0.00, 0.00, 0.00}; // Very central
    
    bool bRandomPsi = true;
    bool bUseWeightning = false;
    double const pi = TMath::Pi();
    
	TStopwatch timer;
	timer.Start();

    TString sUseWeightning;
    if(bUseWeightning) sUseWeightning="weight";
    else sUseWeightning="noWeight";

    TString sRandomPsi;
    if(bRandomPsi) sRandomPsi="randomPsi";
    else sRandomPsi="noRandomPsi";
	
    TString outFileName = Form("toyFlow_%s_%s_dNdeta-%.0f_nEvents-%d-%s.root",sUseWeightning.Data(),sRandomPsi.Data(),dNdeta,nEvents,sFileText.Data());
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
    
    TF1 *fPhiDistribution = new TF1("fPhiDistribution", SingleParticlePhi, -pi, pi, 10);
    fPhiDistribution->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4], Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
    
	TClonesArray *allHadrons = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventB = new TClonesArray("JToyMCTrack", nMult+1);
	TClonesArray *allHadronsDet = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventDetA = new TClonesArray("JToyMCTrack", nMult+1);
    TClonesArray *subEventDetB = new TClonesArray("JToyMCTrack", nMult+1);
	
    // save input numbers
    TH1D *hInputNumbers = new TH1D("hInputNumbers","hInputNumbers",13, 0.5, 13.5);
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

    // Calculate efficiency related values
    EffValues *effValues[nHarmonics+1];
    for (int iHarm=1; iHarm < nHarmonics; iHarm++) {
        effValues[iHarm] = new EffValues;
        efficiencyCalc(effValues[iHarm], iHarm);
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
        
        if(bRandomPsi){
            for(int j = 0; j < 5; j++){
                Psi[j] = randomGenerator->Uniform(-pi,pi);
            }
            fPhiDistribution->SetParameters(vn[0], vn[1], vn[2], vn[3], vn[4], Psi[0], Psi[1], Psi[2], Psi[3], Psi[4]);
        }
        
        GetEvent(allHadrons, subEventA, subEventB, randomGenerator, fPtDistribution, fPhiDistribution, nMult, histos);
        GetDetectorParticles(allHadrons, allHadronsDet, subEventDetA, subEventDetB, randomGenerator, histos);
        AnalyzeEvent(allHadrons, subEventA, subEventB, bUseWeightning, false, histos, effValues, 0); //True MC
        AnalyzeEvent(allHadronsDet, subEventDetA, subEventDetB, bUseWeightning, false, histos, effValues, 1); //Detector, no corr
        AnalyzeEvent(allHadronsDet, subEventDetA, subEventDetB, bUseWeightning, true, histos, effValues, 2); //Detector, corr
    }
                                    
	fOut->Write();
	timer.Print(); 
    //delete effValues;
	return 0;
}
                                    
                                    
// ============== END MAIN PROGRAM

void GetEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, TRandom3 *rand, TF1 *fPt, TF1 *fPhi, int nMult, JHistos *histos){
    
    histos->hMultiplicity[0]->Fill(1.*nMult);
    
    JToyMCTrack track;
    TLorentzVector lVec;

    int nSub = int( nMult/2.0 );
    
    for(int i = 0; i < nMult; i++){
        double pT = fPt->GetRandom();
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
    }
    return;
}

void GetDetectorParticles(TClonesArray *listAll, 
                          TClonesArray *listAllDet, TClonesArray *listSubADet, TClonesArray *listSubBDet,
                          TRandom3 *rand, JHistos *histos){

    int numberOfTracks = listAll->GetEntriesFast();

    double eta,phi,pT,px,py,pz,E;
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

        if((phi < detAMax && phi > detAMin && rand->Uniform(0.0,1.0) > detAEff) ||
           (phi < detBMax && phi > detBMin && rand->Uniform(0.0,1.0) > detBEff) ) {
            new((*listAllDet)[iTrack++]) JToyMCTrack( newTrack );
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



void AnalyzeEvent(TClonesArray *listAll, TClonesArray *listSubA, TClonesArray *listSubB, bool bUseWeightning, bool bCorr, JHistos *histos, EffValues **effValues, int const nHisto){
    double Qx, QxA, QxB, Qy, QyA, QyB, Q, QA, QB, Q2, Q2A, Q2B, Q4, Q4A, Q4B, Q6, Q6A, Q6B, psi, psiA, psiB, vnSPnominator,vnSPdenominator,vnEPnominator,vnEPdenominator;
    double sqrtSumWeights = 0.0, sqrtSumWeightsA = 0.0, sqrtSumWeightsB = 0.0;

    Qvalues *qval = new Qvalues;

    for(int i = 0; i < 5; i++){


        int n = i+1; // which harmonic

        sqrtSumWeights  = CalculateCumulants(listAll, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        Qx = qval->Qx; Qy = qval->Qy; Q = qval->Q; Q2 = qval->Q2; Q4 = qval->Q4; Q6 = qval->Q6; psi = qval->psi;

        sqrtSumWeightsA  = CalculateCumulants(listSubA, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        QxA = qval->Qx; QyA = qval->Qy; QA = qval->Q; Q2A = qval->Q2; Q4A = qval->Q4; Q6A = qval->Q6; psiA = qval->psi;

        sqrtSumWeightsB  = CalculateCumulants(listSubB, n, bUseWeightning, qval, effValues[n], histos, nHisto, bCorr);
        QxB = qval->Qx; QyB = qval->Qy; QB = qval->Q; Q2B = qval->Q2; Q4B = qval->Q4; Q6B = qval->Q6; psiB = qval->psi;

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


        vnSPnominator = CalculateDotProduct(Qx, Qy, QxA, QyA);
        vnSPdenominator = CalculateDotProduct(QxA, QyA, QxB, QyB);

        vnEPnominator = vnSPnominator / QA;
        vnEPdenominator = vnSPdenominator / QA / QB;

        histos->hSPnominator[nHisto][i]->Fill(vnSPnominator);
        histos->hSPdenominator[nHisto][i]->Fill(vnSPdenominator);
        histos->hEPnominator[nHisto][i]->Fill(vnEPnominator);
        histos->hEPdenominator[nHisto][i]->Fill(vnEPdenominator);
    }

    histos->hSqrtSumWeights[nHisto]->Fill( sqrtSumWeights );
    histos->hSqrtSumWeightsA[nHisto]->Fill( sqrtSumWeightsA );
    histos->hSqrtSumWeightsB[nHisto]->Fill( sqrtSumWeightsB );

    delete qval;

    return;
}

double CalculateCumulants(TClonesArray *list, int n, bool bUseWeightning, Qvalues *qval, EffValues *effValues, JHistos *histos, int const nHisto, bool bCorr){

    int numberOfTracks = list->GetEntriesFast();

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
        normalization += weight * weight;

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

            // How about autocorrelation? Flow lecture p 15.
        }


        Qx += weight * xx;
        Qy += weight * yy;
    }

    normalization = TMath::Sqrt( normalization );

    Qx /= normalization; Qy /= normalization;

    double Q2 = Qx*Qx + Qy*Qy;
    double Q = TMath::Sqrt( Q2 );
    double Q4 = Q2*Q2;
    double Q6 = Q4*Q2;
    double psi = TMath::ATan2(Qy,Qx)/((double)n);
    double vObs;

    for(int i = 0; i < numberOfTracks; i++){

        JToyMCTrack *track = (JToyMCTrack*)list->At(i);
        phi = track->GetPhi();
        vObs = TMath::Cos(n*(phi-psi));
        histos->hvObs[nHisto][n-1]->Fill(vObs);
    }

    (*qval).Qx = Qx;
    (*qval).Qy = Qy;
    (*qval).Q  = Q;
    (*qval).Q2 = Q2;
    (*qval).Q4 = Q4;
    (*qval).Q6 = Q6;
    (*qval).psi = psi;

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

double SingleParticlePhi(double *x, double *p){
    double phi = x[0];
    double v1  = p[0];
    double v2  = p[1];
    double v3  = p[2];
    double v4  = p[3];
    double v5  = p[4];
    double Psi1 = p[5];
    double Psi2 = p[6];
    double Psi3 = p[7];
    double Psi4 = p[8];
    double Psi5 = p[9];

    return 1.0 + 2.*v1*TMath::Cos(phi-Psi1) + 2.*v2*TMath::Cos(2.*(phi-Psi2))
        + 2.*v3*TMath::Cos(3.*(phi-Psi3)) + 2.*v4*TMath::Cos(4.*(phi-Psi4))
        + 2.*v5*TMath::Cos(5.*(phi-Psi5));
}

double AcceptanceFunc(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double accA = detAEff;
    double accB = detBEff;
    double value = 0.0;
    if(phi < detAMax && phi > detAMin) {
        value = scaling*accA/(2.0*TMath::Pi());
    } else if(phi < detBMax && phi > detBMin) {
        value = scaling*accB/(2.0*TMath::Pi());
    }
    return value;
}

double AcceptanceFuncCos(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double multiply = p[1];
    double accA = detAEff;
    double accB = detBEff;
    double value = 0.0;
    if(phi < detAMax && phi > detAMin) {
        value = scaling*accA/(2.0*TMath::Pi()) * TMath::Cos(multiply*phi);
    } else if(phi < detBMax && phi > detBMin) {
        value = scaling*accB/(2.0*TMath::Pi()) * TMath::Cos(multiply*phi);
    }
    return value;
}

double AcceptanceFuncSin(double *x, double *p){
    double phi = x[0];
    double scaling = p[0];
    double multiply = p[1];
    double accA = detAEff;
    double accB = detBEff;
    double value = 0.0;
    if(phi < detAMax && phi > detAMin) {
        value = scaling*accA/(2.0*TMath::Pi()) * TMath::Sin(multiply*phi);
    } else if(phi < detBMax && phi > detBMin) {
        value = scaling*accB/(2.0*TMath::Pi()) * TMath::Sin(multiply*phi);
    }
    return value;
}



// Delta_phi on correlation function
double DeltaPhi(double phi1, double phi2) {
    double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
    return res>-TMath::Pi()/3.0 ? res : 2.0*TMath::Pi()+res ; 
}

void efficiencyCalc(EffValues *effValues, int iHarm) {

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

    integral = fAcceptanceFunc->Integral(-pi,pi);
    //cout << "Acc func integral: " << integral << endl;
    fAcceptanceFunc->SetParameter(0,1.0/integral);
    fAcceptanceFuncSin->SetParameter(0,1.0/integral);
    fAcceptanceFuncCos->SetParameter(0,1.0/integral);
    //integral = fAcceptanceFunc->Integral(-pi,pi);
    //cout << "Acc func integral after scaling: " << integral << endl;

    fAcceptanceFuncSin->SetParameter(1,(double)iHarm);
    fAcceptanceFuncCos->SetParameter(1,(double)iHarm);
    sBar = fAcceptanceFuncSin->Integral(-pi,pi);
    cBar = fAcceptanceFuncCos->Integral(-pi,pi);
    //cout << "sBar=" << sBar << ", cBar=" << cBar << endl;

    fAcceptanceFuncSin->SetParameter(1,2.0*(double)iHarm);
    fAcceptanceFuncCos->SetParameter(1,2.0*(double)iHarm);
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
    
}

