/*
 *  JToyMCTrack.cxx
 *  
 *
 *  Created by Sami Rasanen on 2/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "JToyMCTrack.h"
#include "TMath.h"

JToyMCTrack::JToyMCTrack(){
	lVec.SetXYZM(0.0, 0.0, 0.0, 0.0);
	isHT = 0;
	isHard = 0;
	isSoft = 0;
	isCharged = 0;
	isIsolated = 0;
	isLeading = 0;
	nID = 0;
	coneActivity = 999999.;
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in){
	lVec.SetXYZM(lVec_in.Px(),lVec_in.Py(),lVec_in.Pz(),lVec_in.Mag());
	isHT = 0;
	isHard = 0;
	isSoft = 0;
	isCharged = 0;
	isIsolated = 0;
	isLeading = 0;
	nID = 0;
	coneActivity = 999999.;
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in){
	lVec.SetXYZM(lVec_in.Px(),lVec_in.Py(),lVec_in.Pz(),lVec_in.Mag());
	isHT = isHT_in;
	isHard = isHard_in;
	isSoft = isSoft_in;
	isCharged = 0;
	isIsolated = 0;
	isLeading = 0;
	nID = 0;
	coneActivity = 999999.;
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in, int isCharged_in){
	lVec.SetXYZM(lVec_in.Px(),lVec_in.Py(),lVec_in.Pz(),lVec_in.Mag());
	isHT = isHT_in;
	isHard = isHard_in;
	isSoft = isSoft_in;
	isCharged = isCharged_in;
	isIsolated = 0;
	isLeading = 0;
	nID = 0;
	coneActivity = 999999.;
}

JToyMCTrack::JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in, int isCharged_in, int isIsolated_in, int IsLeading_in){
	lVec.SetXYZM(lVec_in.Px(),lVec_in.Py(),lVec_in.Pz(),lVec_in.Mag());
	isHT = isHT_in;
	isHard = isHard_in;
	isSoft = isSoft_in;
	isCharged = isCharged_in;
	isIsolated = isIsolated_in;
	isLeading = IsLeading_in;
	nID = 0;
	coneActivity = 999999.;
}

void JToyMCTrack::SetMass(double m){
	double px = lVec.Px();
	double py = lVec.Py();
	double pz = lVec.Pz();
	double E_in = TMath::Sqrt( m*m + px*px + py*py + pz*pz );	
	lVec.SetE(E_in);
}

void JToyMCTrack::SetTrack(TLorentzVector lVec_in, int isHT_in = 0, int isHard_in = 0, int isSoft_in = 0, 
						   int isCharged_in = 0, int isIsolated_in = 0, int IsLeading_in = 0, int nID_in = 0){
	lVec = lVec_in;
	isHT = isHT_in;
	isHard = isHard_in;
	isSoft = isSoft_in;
	isCharged = isCharged_in;
	isIsolated = isIsolated_in;
	isLeading = IsLeading_in;
	nID = nID_in;	
}

ClassImp(JToyMCTrack)

