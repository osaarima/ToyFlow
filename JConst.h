#ifndef JCONST_H
#define JCONST_H

#define PI TMath::Pi()

const double etaRange = 0.8;

const int nHarmonics = 5;

// for pT-dependence
const double alpha = 2.0;
const double defpTMax = 2.0;
const double Tdec = 0.12;
const double vr = 0.6;

#define SECTORS_N 8
#define RINGS_N 5

//To be implemented:
//#define PTBINS_N 9
//static double pTBins[PTBINS_N+1] = {0.0, 0.2, 0.6, 1.2, 2.0, 3.0, 4.2, 5.6, 7.2, 9.0};

#endif
