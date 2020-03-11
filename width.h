#include <iostream>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include "stdio.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TF3.h"
#include "TF1.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/Functor.h"
#include "Math/ParamFunctor.h"
#include "Math/WrappedParamFunction.h"
#include "Math/GaussIntegrator.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/WrappedTF1.h"

using namespace std;

double width_X_total (double eps, double mA) 
{
	double alpha = 1/137.;
	double me = 0.000511;
	double mmu = 0.1056584;
	double gEE = eps*eps/3*alpha*mA*(1+2*me*me/mA/mA)*sqrt(1-4*me*me/mA/mA);
	double gMuMu = 0.;
	if ( mA > 2*mmu ) gMuMu = eps*eps/3*alpha*mA*(1+2*mmu*mmu/mA/mA)*sqrt(1-4*mmu*mmu/mA/mA);
	return gEE + gMuMu;
}

double width_chi_total (double gOff, double eps, double mA, double m2, double m1) 
{
	double alpha = 1/137.;
	double me = 0.000511;
	double aVal = (m1+m2)*(m1+m2) + 4*me*me;
	double bVal = 2*m1*m2*(m1*m1+m2*m2+m1*m2)+2*me*me*(m1*m1+m2*m2+4*m1*m2+3*me*me);
	double cVal = m1*m1+m2*m2+2*me*me-mA*mA;
	double dVal = m1*m1+me*me;
	double eVal = m2*m2-me*me;
	double fVal = m1*m1-me*me;
	double gVal = m2*m2+me*me;

	TF1 f("ThreeBodyWidth", "([0]*[2]-[1]-[2]*[2]+2*[2]*x-2*x*x)/(x-[2]+[3]+0.5*(([4]-x)*([5]+x)-sqrt(x*x-2*x*[6]+[4]*[4])*sqrt(x*x-2*x*[3]+[5]*[5]))/x)-([0]*[2]-[1]-[2]*[2]+2*[2]*x-2*x*x)/(x-[2]+[3]+0.5*(([4]-x)*([5]+x)+sqrt(x*x-2*x*[6]+[4]*[4])*sqrt(x*x-2*x*[3]+[5]*[5]))/x)+([0]+2*x-2*[2])*log( (-[2]+x+ [3]+0.5*(([4]-x)*([5]+x)+sqrt(x*x-2*x*[6]+[4]*[4])*sqrt(x*x-2*x*[3]+[5]*[5]))/x )/(-[2]+x+ [3]+0.5*(([4]-x)*([5]+x)-sqrt(x*x-2*x*[6]+[4]*[4])*sqrt(x*x-2*x*[3]+[5]*[5]))/x ) ) - sqrt(x*x-2*x*[6]+[4]*[4])*sqrt(x*x-2*x*[3]+[5]*[5])/x", 0, 1000);
	f.SetParameter(0,aVal);
	f.SetParameter(1,bVal);
	f.SetParameter(2,cVal);
	f.SetParameter(3,dVal);
	f.SetParameter(4,eVal);
	f.SetParameter(5,fVal);
	f.SetParameter(6,gVal);

	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);
	
	double coeff = gOff*gOff * eps*eps * alpha / 16 / TMath::Pi()/TMath::Pi() / m2/m2/m2;
	double xmin = (m1+me)*(m1+me);
	double xmax = (m2-me)*(m2-me);
	return coeff * ig.Integral(xmin, xmax);
}
