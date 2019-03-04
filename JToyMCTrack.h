/*
 *  JToyMCTrack.h
 *  
 *
 *  Created by Sami Rasanen on 2/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TLorentzVector.h"

class JToyMCTrack : public TObject {
	
	public:
		JToyMCTrack();
		JToyMCTrack(TLorentzVector lVec_in);
		JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in);
		JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in, int isCharged_in);
		JToyMCTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in, int isCharged_in, int isIsolated_in, int IsLeading_in);
		
		virtual ~JToyMCTrack(){;}
		
		double GetPx(){return lVec.Px();}
		double GetPy(){return lVec.Py();}
		double GetPz(){return lVec.Pz();}
		double GetPt(){return lVec.Perp();}
		double GetE(){return lVec.E();}
		double GetPhi(){return lVec.Phi();}
		double GetEta(){return lVec.Eta();}
		double GetMass(){return lVec.Mag();}
	
		TLorentzVector GetLVector(){return lVec;}
		
		int GetIsIsolated(){return isIsolated;}
		int GetIsLeading(){return isLeading;}
		int GetIsHT(){return isHT;}
		int GetIsHard(){return isHard;}
		int GetIsSoft(){return isSoft;}
		int GetIsCharged(){return isCharged;}
		double GetConeActivity(){return coneActivity;}
		
		int GetIDnumber(){return nID;}
		
		void SetPx(double px_in){lVec.SetPx(px_in);}
		void SetPy(double py_in){lVec.SetPy(py_in);}
		void SetPz(double pz_in){lVec.SetPz(pz_in);}
		void SetPt(double pt_in){lVec.SetPerp(pt_in);}
		void SetPhi(double phi_in){lVec.SetPhi(phi_in);}
		void SetMass(double m_in);
		
		void SetLVector(TLorentzVector lVec_in){lVec = lVec_in;}
		
		void SetTrack(TLorentzVector lVec_in, int isHT_in, int isHard_in, int isSoft_in, 
					  int isCharged_in, int isIsolated_in, int IsLeading_in, int nID);
		
		void SetIsIsolated(int isIsolated_in){isIsolated = isIsolated_in;}
		void SetIsLeading(int isLeading_in){isLeading = isLeading_in;}
		void SetIsHT(int isHT_in){isHT = isHT_in;}
		void SetIsHard(int isHard_in){isHard = isHard_in;}
		void SetIsSoft(int isSoft_in){isSoft = isSoft_in;}
		void SetIsCharged(int isCharged_in){isCharged = isCharged_in;}
		void SetConeActivity(double coneActivity_in)
			{coneActivity_in>0?coneActivity=coneActivity_in:coneActivity=999999.;}
		
		void SetIDnumber(int nID_in){nID = nID_in;}
		
	private:
		TLorentzVector lVec;
		int isHT;
		int isHard;
		int isSoft;
		int isIsolated;
		int isLeading;
		int isCharged;
		int nID;
		double coneActivity;

	protected:
		ClassDef(JToyMCTrack,1)

};

