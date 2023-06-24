#define electron_cxx
#include "electron.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>
#include <string>

using namespace std;

int main(int argc, const char*  argv[])
{
  string outfile = argv[2];
  Long64_t maxEvents = atof(argv[3]);
  TString p = (argv[1]);  
  electron t(p);
  t.Loop(outfile, maxEvents);
  return 0;
}

void electron::Loop(string outputfile, int maxEvents)
{
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntries();
  if(maxEvents >= 0) nentries = maxEvents;
  Long64_t nbytes = 0, nb = 0;
  int nTotal = 0;
  
  char* ouputname = const_cast<char*>(outputfile.c_str());
  TFile *output = new TFile(outputfile.c_str(),"RECREATE"); //Converting char into a string type

  //Declaration of pT histograms

  TH1F* Zee200_220 = new TH1F("Z_ee200_220","Z->ee Invariant_Mass200_220;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zee220_270 = new TH1F("Z_ee220_270","Z->ee Invariant_Mass220_270;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zee270_320 = new TH1F("Z_ee270_320","Z->ee Invariant_Mass270_320;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zee320     = new TH1F("Z_ee320","Z->ee Invariant_Mass320;Z_{mass}(GeV);No. of Events"    ,50,60,120);
  
  TH1F* Zeg200_220 = new TH1F("Z_eg200_220","Z->eg Invariant_Mass200_220;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zeg220_270 = new TH1F("Z_eg220_270","Z->eg Invariant_Mass220_270;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zeg270_320 = new TH1F("Z_eg270_320","Z->eg Invariant_Mass270_320;Z_{mass}(GeV);No. of Events",50,60,120);
  TH1F* Zeg320     = new TH1F("Z_eg320","Z->eg Invariant_Mass320;Z_{mass}(GeV);No. of Events"    ,50,60,120);
  
  //Declaration of eta histograms
  TH1F* Zee_eta0_2     = new TH1F("Z_ee eta 0-0.2","Z->ee Invariant_Mass_eta_0_0.2;Z_{mass}(GeV);No. of Events"      ,50,60,120);
  TH1F* Zee_eta2_4     = new TH1F("Z_ee eta 0.2-0.4","Z->ee Invariant_Mass_eta_0.2_0.4;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zee_eta4_8     = new TH1F("Z_ee eta 0.4-0.8","Z->ee Invariant_Mass_eta_0.4_0.8;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zee_eta8_12    = new TH1F("Z_ee eta 0.8-1.2","Z->ee Invariant_Mass_eta_0.8_1.2;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zee_eta12_14   = new TH1F("Z_ee eta 1.2-1.44","Z->ee Invariant_Mass_eta_1.2_1.44;Z_{mass}(GeV);No. of Events",50,60,120);
  
  TH1F* Zeg_eta0_2     = new TH1F("Z_eg eta 0-0.2","Z->eg Invariant_Mass_eta_0_0.2;Z_{mass}(GeV);No. of Events"      ,50,60,120);
  TH1F* Zeg_eta2_4     = new TH1F("Z_eg eta 0.2-0.4","Z->eg Invariant_Mass_eta_0.2_0.4;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zeg_eta4_8     = new TH1F("Z_eg eta 0.4-0.8","Z->eg Invariant_Mass_eta_0.4_0.8;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zeg_eta8_12    = new TH1F("Z_eg eta 0.8-1.2","Z->eg Invariant_Mass_eta_0.8_1.2;Z_{mass}(GeV);No. of Events"  ,50,60,120);
  TH1F* Zeg_eta12_14   = new TH1F("Z_eg eta 1.2-1.44","Z->eg Invariant_Mass_eta_1.2_1.44;Z_{mass}(GeV);No. of Events",50,60,120);
  
  int n_ele = 0, n_pho = 0, n_HLTpho = 0;
  
  for (Long64_t jentry=0; jentry < nentries; jentry++) 
    {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      double ele_M; // Invariant Mass
      
      //HLT condition
      if((HLTPho>>11) & 1 == 1)
	{
	  n_HLTpho++;

	  for (int i=0; i<nEle; i++)
	    {
	      TLorentzVector etag;      

	      //electron cuts
	      bool ele_cuts = (   elePt->at(i) > 30 && 
				  abs(eleEta->at(i)) < 2.5 &&
				  eleIDbit->at(i)>>2 & 1==1   );
	      if (ele_cuts)
		{
		  etag.SetPtEtaPhiE(elePt->at(i),eleEta->at(i),elePhi->at(i),eleE->at(i));
		} 
	      
	      for(int j=0; j<nPho; j++)
		{
		  n_ele++;
		  //	      cout<<"Tag mass: "<<etag.M()<<endl;
		  
		  //photon cuts
		  Float_t uncorrectedPhoEt = phoEt->at(j);
		  bool MediumphotonID = (
					 ( phoHoverE->at(j)                <  0.0109   ) &&
					 ( phoSigmaIEtaIEtaFull5x5->at(j)  <  0.026    ) &&
					 ( fabs(phoSCEta->at(j)) < 2.5)                  && 
					 ( TMath::Max( ( phoPFChWorstIso->at(j)  - rho*EAchargedworst(phoSCEta->at(j)) ), 0.0) < 4.98 )  &&
					 ( TMath::Max( ( phoPFNeuIso->at(j) - rho*EAneutral(phoSCEta->at(j)) ), 0.0) < (5.24 + (0.0072 * uncorrectedPhoEt) + (0.000017 * pow(uncorrectedPhoEt, 2.0))) )  &&
					 ( TMath::Max( ( phoPFPhoIso->at(j) - rho*EAphoton(phoSCEta->at(j))  ), 0.0) < (1.20 + (0.0019 * uncorrectedPhoEt)) )
					 );
		  
		  if(MediumphotonID)
		    {
		      n_pho++;
		      TLorentzVector probe;
		      probe.SetPtEtaPhiE(phoEt->at(j),phoEta->at(j),phoPhi->at(j),phoE->at(j));
		      //Invariant Mass
		      ele_M = (etag + probe).M();	  
		      
		      if(ele_M>60 && ele_M<120)
			{
			  if(phohasPixelSeed->at(j) == 1)
			    //photon pT for ee
			    {
			      if(phoEt->at(j) >= 200 && phoEt->at(j) < 220) {Zee200_220 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 220 && phoEt->at(j) < 270) {Zee220_270 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 270 && phoEt->at(j) < 320) {Zee270_320 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 320 )                      {Zee320     -> Fill(ele_M);}
			      
			      //eta 
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.0 && fabs(phoEta->at(j)) < 0.2)  {Zee_eta0_2->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.2 && fabs(phoEta->at(j)) < 0.4)  {Zee_eta2_4->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.4 && fabs(phoEta->at(j)) < 0.8)  {Zee_eta4_8->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.8 && fabs(phoEta->at(j)) < 1.2)  {Zee_eta8_12->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 1.2 && fabs(phoEta->at(j)) < 1.442){Zee_eta12_14->Fill(ele_M);}
			      
			    }
			  else if(phohasPixelSeed->at(j) == 0)
			    {
			      //photon pT for eg
			      if(phoEt->at(j) >= 200 && phoEt->at(j) < 220) {Zeg200_220 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 220 && phoEt->at(j) < 270) {Zeg220_270 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 270 && phoEt->at(j) < 320) {Zeg270_320 -> Fill(ele_M);}
			      if(phoEt->at(j) >= 320 )                      {Zeg320     -> Fill(ele_M);}
			      
			      //eta 			      
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.0 && fabs(phoEta->at(j)) < 0.2)  {Zeg_eta0_2->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.2 && fabs(phoEta->at(j)) < 0.4)  {Zeg_eta2_4->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.4 && fabs(phoEta->at(j)) < 0.8)  {Zeg_eta4_8->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 0.8 && fabs(phoEta->at(j)) < 1.2)  {Zeg_eta8_12->Fill(ele_M);}
			      if((phoEt->at(j)) > 225 && fabs(phoEta->at(j)) >= 1.2 && fabs(phoEta->at(j)) < 1.442){Zeg_eta12_14->Fill(ele_M);}
			      
			    }
			}
		    }
		  
		}
	    }
	}
    }
  cout<<"Total events: "<<nentries<<endl;
  cout<<"Number of ee pair: "<<n_ele<<endl;
  cout<<"Number of eg pair: "<<n_pho<<endl;
      
  TCanvas *c = new TCanvas();
  
  //************************** Set Stats() initialize ****************************
  //pT Plots
  Zee200_220 ->SetStats(0);      Zee220_270 ->SetStats(0);
  Zee270_320 ->SetStats(0);      Zee320   ->  SetStats(0);
  
  //eta Plots
  Zee_eta0_2 -> SetStats(0);     Zee_eta2_4 -> SetStats(0);
  Zee_eta4_8 -> SetStats(0);     Zee_eta8_12 -> SetStats(0);
  Zee_eta12_14 -> SetStats(0);
  
  //pT Plots 
  Zeg200_220 ->SetStats(0);     Zeg220_270 ->SetStats(0);
  Zeg270_320 ->SetStats(0);     Zeg320 ->SetStats(0);
  
  //eta Plots
  Zeg_eta0_2 -> SetStats(0);     Zeg_eta2_4 -> SetStats(0);
  Zeg_eta4_8 -> SetStats(0);     Zeg_eta8_12 -> SetStats(0);
  Zeg_eta12_14 -> SetStats(0);
  //****************************** Set Stats() ends ******************************
 
  Zee200_220->Write(); 
  Zee220_270->Write();
  Zee270_320->Write();
  Zee320->Write();
  
  Zeg200_220->Write();
  Zeg220_270->Write();
  Zeg270_320->Write();
  Zeg320->Write();
  
  Zee_eta0_2->Write();
  Zee_eta2_4->Write();
  Zee_eta4_8->Write();
  Zee_eta8_12->Write();
  Zee_eta12_14->Write();
  
  Zeg_eta0_2->Write();
  Zeg_eta2_4->Write();
  Zeg_eta4_8->Write();
  Zeg_eta8_12->Write();
  Zeg_eta12_14->Write();
  
  output -> Write();
}

Double_t electron::EAchargedworst(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.06118;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.05545;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.05000;
  return EffectiveArea;
}

Double_t electron::EAneutral(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.03279;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.05682;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.07969;
  return EffectiveArea;
}

Double_t electron::EAphoton(Double_t eta){
  Float_t EffectiveArea = 0.0;
  if(fabs(eta) >= 0.0   && fabs(eta) < 0.5   ) EffectiveArea = 0.1078;
  if(fabs(eta) >= 0.5   && fabs(eta) < 1.0   ) EffectiveArea = 0.1068;
  if(fabs(eta) >= 1.0   && fabs(eta) < 1.479 ) EffectiveArea = 0.08661;
  return EffectiveArea;
}
