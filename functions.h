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

double qSQ (double mT, double eT) {
	return 2*mT*(eT - mT);
}

double gE (double mT, double eT) {
	return pow(1+qSQ(mT, eT)/0.71,-2);
}

double gM (double mT, double eT) {
	double mUP = 2.79;
	return mUP*gE(mT, eT);
}

double f1Func (double mT, double eT) {
	return ( gE(mT, eT) + qSQ(mT, eT)/4/mT/mT )/(1 + qSQ(mT, eT)/4/mT/mT );
}

double kf2Func (double mT, double eT) {
	return ( gM(mT, eT) - gE(mT, eT) )/(1 + qSQ(mT, eT)/4/mT/mT );
}

double primaryMECalc (int mTSwitch, double eH, double eX, double eT, double mT, double mA, double mX, double mH, double eps, double g12) {
	double alpha = 1/137.;
	double f1 = 1.;
	double kf2 = 0.;
	double m0 = mT*(eX*eX + eH*eH) - 0.5*(mH - mX)*(mH - mX)*(eH - eX + mT) + mT*mT*(eH - eX) + mX*mX*eH - mH*mH*eX;
	double m1 = 0.;
	if (mTSwitch == 2) {	
		m1 = mT*( pow( eH + eX - 0.5*(mH*mH - mX*mX)/mT ,2) + (eX - eH + 2*mT)*(eH - eX - 0.5*(mH - mX)*(mH - mX)/mT) );
		f1 = f1Func(mT, eT);
		kf2 = kf2Func(mT, eT);
	}
	double mSQ = 8*4*TMath::Pi()*alpha*eps*eps*g12*g12*mT/pow(2*mT*(eH - eX) - mA*mA,2)*( m0*pow(f1+kf2,2) + m1*( -(f1+kf2)*kf2 + pow(kf2,2)/4/mT*(eX - eH + 2*mT) ) );
	return mSQ;
}

double secondaryMECalc(double s1, double s2, double m1, double m2, double mX, double me, double eps, double g12) {
	double alpha = 1/137.;
	double aSq = ( (s1+s2)*(pow(m1+m2,2)+4*me*me)-s1*s1-s2*s2-2*m1*m2*(m1*m1+m2*m2+m1*m2)-2*me*me*(m1*m1+m2*m2+4*m1*m2+3*me*me) )/pow(m1*m1+m2*m2+2*me*me-s1-s2-mX*mX,2)/pow(m2,2);
	return eps*eps*g12*g12*alpha*aSq;
}

double get_max_m2 (double m0, double mT, double m1) {
	return sqrt( mT*mT + 2*m0*mT + m1*m1 ) - mT;
}

double lam (double x, double y, double z) {
	return x*x + y*y + z*z -2*(x*y + y*z + z*x);
}

double get_s (double e1, double m1, double mT) {
	return mT*mT + m1*m1 + 2*e1*mT;
}

double get_PS (double e1, double m1, double mT) {
	return mT/8/TMath::Pi()/lam(get_s(e1,m1,mT),mT*mT,m1*m1);	
}

double get_ETRest(double e1, double m1, double mT, double m2) {
	return (get_s(e1,m1,mT)+mT*mT-m2*m2)/2/sqrt(get_s(e1,m1,mT));
}

double get_PTRest(double e1, double m1, double mT, double m2) {
	return sqrt( pow(get_ETRest(e1,m1,mT,m2),2)-mT*mT );
}

double get_EPsiRest(double e1, double m1, double mT, double m2) {
	return (get_s(e1,m1,mT)-mT*mT+m2*m2)/2/sqrt(get_s(e1,m1,mT));
}

double get_PPsiRest(double e1, double m1, double mT, double m2) {
	return sqrt( pow(get_EPsiRest(e1,m1,mT,m2),2)-m2*m2 );
}

double get_cosh(double e1, double m1, double mT) {
	return (e1+mT)/sqrt(get_s(e1,m1,mT));
}

double get_sinh(double e1, double m1, double mT) {
	return sqrt(get_cosh(e1,m1,mT)*get_cosh(e1,m1,mT) - 1);
}

double get_ETlabMax (double e1, double m1, double mT, double m2) {
	return get_ETRest(e1,m1,mT,m2)*get_cosh(e1,m1,mT) + get_PTRest(e1,m1,mT,m2)*get_sinh(e1,m1,mT);
}

double get_ETlabMin (double e1, double m1, double mT, double m2) {
	return get_ETRest(e1,m1,mT,m2)*get_cosh(e1,m1,mT) - get_PTRest(e1,m1,mT,m2)*get_sinh(e1,m1,mT);
}

double get_EPsilabMax (double e1, double m1, double mT, double m2) {
	return get_EPsiRest(e1,m1,mT,m2)*get_cosh(e1,m1,mT) + get_PPsiRest(e1,m1,mT,m2)*get_sinh(e1,m1,mT);
}

double get_ESecRest (double m1, double m2, double mX) {
	return (m2*m2 - m1*m1 + mX*mX)/2/m2;
}

double get_PSecRest (double m1, double m2, double mX) {
	return sqrt( pow( get_ESecRest(m1,m2,mX), 2) - mX*mX );
}

double get_xSec_electron(double e1, double m1, double mT, double m2, double mX, double eps, double g12) 
{
	double val_xSec = 0.;
	double xMax = get_ETlabMax(e1,m1,mT,m2);
	double xMin = get_ETlabMin(e1,m1,mT,m2);

	TF1 f("xSec_electron","8*[2]/pow(2*[2]*([2]-x)-[4]*[4],2)*([2]*([0]*[0]+pow([0]+[2]-x,2))-0.5*pow([3]-[1],2)*(2*[2]-x)+[2]*[2]*([2]-x)+[1]*[1]*([0]+[2]-x)-[3]*[3]*[0] )",xMin,xMax);
	f.SetParameter(0,e1);
	f.SetParameter(1,m1);
	f.SetParameter(2,mT);
	f.SetParameter(3,m2);
	f.SetParameter(4,mX);

	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);

	val_xSec = 4*TMath::Pi()/137 * eps * eps * g12 * g12 * get_PS (e1,m1,mT) * ig.Integral(xMin, xMax);
	return val_xSec*3.9057*pow(10,-28);
}

double get_ET_electron (double e1, double m1, double mT, double m2, double mX) 
{	
	double xMax = get_ETlabMax(e1,m1,mT,m2);
	double xMin = get_ETlabMin(e1,m1,mT,m2);

	TF1 *f = new TF1("ET_electron","8*[2]/pow(2*[2]*([2]-x)-[4]*[4],2)*([2]*([0]*[0]+pow([0]+[2]-x,2))-0.5*pow([3]-[1],2)*(2*[2]-x)+[2]*[2]*([2]-x)+[1]*[1]*([0]+[2]-x)-[3]*[3]*[0] )",xMin,xMax);
	f->SetParameter(0,e1);
	f->SetParameter(1,m1);
	f->SetParameter(2,mT);
	f->SetParameter(3,m2);
	f->SetParameter(4,mX);

	double r_num = f->GetRandom(xMin, xMax);

	return r_num;
}

double get_xSec_proton(double e1, double m1, double mT, double m2, double mX, double eps, double g12) 
{
	double val_xSec = 0.;
	double xMax = get_ETlabMax(e1,m1,mT,m2);
	double xMin = get_ETlabMin(e1,m1,mT,m2);

	TF1 f("xSec_proton","8*[2]/pow(2*[2]*([2]-x)-[4]*[4],2)*( ([2]*([0]*[0]+pow([0]+[2]-x,2))-0.5*pow([3]-[1],2)*(2*[2]-x)+[2]*[2]*([2]-x)+[1]*[1]*([0]+[2]-x)-[3]*[3]*[0])*pow( (1/pow(1+2*[2]*(x-[2])/0.71,2)+2*(x-[2])/4/[2])/(1+2*(x-[2])/4/[2]) + (1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]) ,2)+[2]*( pow(2*[0]+[2]-x-([3]*[3]-[1]*[1])/2/[2],2) + ([2]+x)*([2]-x-pow([3]-[1],2)/2/[2] ) ) * ( -((1/pow(1+2*[2]*(x-[2])/0.71,2)+2*(x-[2])/4/[2])/(1+2*(x-[2])/4/[2])+(1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]))*(1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2])+pow((1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]),2)/4/[2]*(x+[2]) ) )",xMin,xMax);
	f.SetParameter(0,e1);
	f.SetParameter(1,m1);
	f.SetParameter(2,mT);
	f.SetParameter(3,m2);
	f.SetParameter(4,mX);

	ROOT::Math::WrappedTF1 wf1(f);
	ROOT::Math::GaussIntegrator ig;
	ig.SetFunction(wf1);
	ig.SetRelTolerance(0.001);

	val_xSec = 4*TMath::Pi()/137 * eps * eps * g12 * g12 * get_PS (e1,m1,mT) * ig.Integral(xMin, xMax);
	return val_xSec*3.9057*pow(10,-28);
}

double get_ET_proton (double e1, double m1, double mT, double m2, double mX) 
{
	double xMax = get_ETlabMax(e1,m1,mT,m2);
	double xMin = get_ETlabMin(e1,m1,mT,m2);

//	cout << xMin << " " << xMax << endl;

	TF1 *f = new TF1("xET_proton","8*[2]/pow(2*[2]*([2]-x)-[4]*[4],2)*( ([2]*([0]*[0]+pow([0]+[2]-x,2))-0.5*pow([3]-[1],2)*(2*[2]-x)+[2]*[2]*([2]-x)+[1]*[1]*([0]+[2]-x)-[3]*[3]*[0])*pow( (1/pow(1+2*[2]*(x-[2])/0.71,2)+2*(x-[2])/4/[2])/(1+2*(x-[2])/4/[2]) + (1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]) ,2)+[2]*( pow(2*[0]+[2]-x-([3]*[3]-[1]*[1])/2/[2],2) + ([2]+x)*([2]-x-pow([3]-[1],2)/2/[2] ) ) * ( -((1/pow(1+2*[2]*(x-[2])/0.71,2)+2*(x-[2])/4/[2])/(1+2*(x-[2])/4/[2])+(1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]))*(1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2])+pow((1.79/pow(1+2*[2]*(x-[2])/0.71,2))/(1+2*(x-[2])/4/[2]),2)/4/[2]*(x+[2]) ) )",xMin,xMax);
	f->SetParameter(0,e1);
	f->SetParameter(1,m1);
	f->SetParameter(2,mT);
	f->SetParameter(3,m2);
	f->SetParameter(4,mX);

	double r_num = f->GetRandom(xMin, xMax);

	return r_num;
}

double get_flux (double e0) 
{
	return 6.4 * pow(10,-4) * pow(0.5/e0,2);
}

double get_Neve (double xs, double acc, double n_yr, double e0, int n_module)
{
	double atomU_Ar = 39.948;
	double n_proton_Ar = 18;
	double mass_DUNE = n_module * 10. * pow(10,9);
	double avoga_n = 6.02212 * pow(10,23);
	double n_electron = mass_DUNE * n_proton_Ar / atomU_Ar * avoga_n; // the number of protons = the number of electrons
	double yr_to_sec = 365 * 24 * 3600;
	
	return xs * acc * n_electron * yr_to_sec * n_yr * get_flux(e0);
}

double get_gammaH (double e1, double m1, double mT, double m2) {
	return get_EPsilabMax(e1,m1,mT,m2) / m2;
}

double get_gamma (double e1, double m1, double mT, double m2, double mX) {
	return ( get_ESecRest(m1,m2,mX) * get_EPsilabMax(e1,m1,mT,m2) / m2 + get_PSecRest(m1,m2,mX) * sqrt( pow(get_EPsilabMax(e1,m1,mT,m2)/m2,2) -1 ) )/mX;
}