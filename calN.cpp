#include <vector>
#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "functions.h"

using namespace std;


double getL(TLorentzVector *m) {
  double iHcal = 179.; //cm
  double fHcal = 295.; //cm
  return abs((fHcal-iHcal)*1./m->CosTheta());
}

int main() {

  TH1D *photonE = new TH1D("photonE", "Photon E [GeV]",1000, 0, 1000);
  TH1D *chi1E = new TH1D("chi1E", "#chi_{1} E [GeV]",1000, 0, 1000);
  TH1D *nChi1 = new TH1D("nChi1", "n #chi_{1}/event (all photon -> dark photon)",100,0,100);

  double mAp = 0.3; //GeV
  double mChi1 = 0.01; //GeV
  double mChi2 = 0.25; //GeV

  double hcalDensity = 8.83; //g/cm3 Cu:Zn = 7:3, cu (29) 63.546 g/mol, zn (30) 65.38 g/mol
  double nAvogadro = 6.0221409*pow(10,23);
  double lNp = hcalDensity*0.7/63.546*29.*nAvogadro + hcalDensity*0.3/65.38*30.*nAvogadro;

  double eps_ = 1.;
  double g12_ = 1.;
  double mT_ = 0.93827208816; //Gev

  TFile *tf = TFile::Open("QCDPhoton.root","READ");
  TTreeReader reader("QCDPhoton", tf);
  TTreeReaderValue<vector<vector<TLorentzVector>>> photon(reader, "photon");

  double N = 0.; 
  while (reader.Next()) {
    int cChi1 = 0;
    for (const auto cJet:*photon) {
      for (auto p:cJet) {
        // dark photon
        if ((p.E()*p.E() - mAp*mAp) < 0.) continue;
        photonE->Fill(p.E());
        double pp = sqrt(p.E()*p.E() - mAp*mAp)/p.P();
        auto d4v = TLorentzVector();
        d4v.SetPxPyPzE(p.Px()*pp, p.Py()*pp, p.Pz()*pp, p.E());
        //chi1
        auto tmp = TGenPhaseSpace();
        double mass1[2] = {mChi1, mChi1};
        tmp.SetDecay(d4v, 2, mass1);
        auto w = tmp.Generate();
        for (int x = 0; x < 2; x++) {
          auto tmpChi1 = tmp.GetDecay(x);
          if (abs(tmpChi1->Theta()) > TMath::Pi()/2.) continue;
          if (tmpChi1->E() <  ((mChi2*mChi2 - mChi1*mChi1)/(2.*mT_) + mChi2)) continue;
          chi1E->Fill(tmpChi1->E());
          cChi1 += 1;
          N += w*getL(tmpChi1)*lNp*get_xSec_proton(tmpChi1->E(), mChi1, mT_, mChi2, mAp, eps_, g12_);
        }
      }
    }
    nChi1->Fill(cChi1); 
  }
  cout << N << endl;
  TCanvas *c1 = new TCanvas("","",800,600);
  c1->SetLogy(1);
  photonE->Draw("Hist");
  c1->SaveAs("photonE.png");
  chi1E->Draw("Hist");
  c1->SaveAs("chi1E.png");
  nChi1->Draw("Hist");
  c1->SaveAs("nChi1.png");
}
