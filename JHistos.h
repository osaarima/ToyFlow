/*
 *  JHistos.h
 *  
 */

#include "TH1D.h"
#include "TH2D.h"

class JHistos {

public:
	JHistos();
	virtual ~JHistos(){;}
	
	// Iclusive
	TH1D *hPt[3];
	TH1D *hPhi[3];
	TH1D *hEta[3];
    TH2D *hPhiEta[3];
    TH1D *hMultiplicity[3];
    TH1D *hSqrtSumWeights[3];
    TH1D *hSqrtSumWeightsA[3];
    TH1D *hSqrtSumWeightsB[3];
    // integrated flow
    TH1D *hQx[3][5];
    TH1D *hQy[3][5];
    TH1D *hQ[3][5];
    TH1D *hQ2[3][5];
    TH1D *hQ4[3][5];
    TH1D *hQ6[3][5];
    TH1D *hvObs[3][5];
    TH1D *hrSub[3][5];
    TH2D *hPsiAB[3][5][2];
    TH1D *hPsiDiff[3][5][2];
    TH1D *hPsiDiffN[3][5][2];
    TH1D *hTrueReso[3][5];
    TH1D *hSPnominator[3][5];
    TH1D *hSPdenominator[3][5];
    TH1D *hEPnominator[3][5];
    TH1D *hEPdenominator[3][5];
};
