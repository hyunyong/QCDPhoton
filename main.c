#include "selection.h"
#include "functions.h"

using namespace std;

int main(int argc, char** argv) 
{
	gRandom = new TRandom3(0); // Setting the random seed for TGenPhaseSpace

	TRandom3 *my_phi = new TRandom3();
	my_phi->SetSeed(0); // random seed

	// The following two random number generators will be used for defining the direction of imcoming chi1.
	TRandom3 *my_cos_p = new TRandom3();
	my_cos_p->SetSeed(0); // random seed
	TRandom3 *my_phi_p = new TRandom3();
	my_phi_p->SetSeed(0); // random seed

	// Reading information from the param_card. All mass quantities in the param_card should be in GeV.
	int d_var[6];
	float f_var[1];
	float dim_var[3];
	float cuts_e_var[7];
	float cuts_p_var[7];
	float trk_eff_var[2];
	float m_couple_var[5];
	FILE *infofile = fopen(argv[1],"r");
	fscanf(infofile,"%d %d %d %d %d %d", &d_var[0], &d_var[1], &d_var[2], &d_var[3], &d_var[4], &d_var[5]);
	fscanf(infofile,"%d %f", &d_var[6], &f_var[0]);
	fscanf(infofile,"%f %f %f", &dim_var[0], &dim_var[1], &dim_var[2]);
	fscanf(infofile,"%f %f %f %f %f %f %f", &cuts_e_var[0], &cuts_e_var[1], &cuts_e_var[2], &cuts_e_var[3], &cuts_e_var[4], &cuts_e_var[5], &cuts_e_var[6]);
	fscanf(infofile,"%f %f %f %f %f %f %f", &cuts_p_var[0], &cuts_p_var[1], &cuts_p_var[2], &cuts_p_var[3], &cuts_p_var[4], &cuts_p_var[5], &cuts_p_var[6]);
	fscanf(infofile,"%f %f", &trk_eff_var[0], &trk_eff_var[1]);
	fscanf(infofile,"%f %f %f %f %f", &m_couple_var[0], &m_couple_var[1], &m_couple_var[2], &m_couple_var[3], &m_couple_var[4]);
	fclose(infofile);

	// Defining the output file
	ofstream dataout;
	dataout.open(argv[2],ios::out);

	int det_type = d_var[0];
	int type_target = d_var[1];
	int n_eve = d_var[2];
	int inum_in = d_var[3];
	int inum = d_var[4];
	int jnum = d_var[5];

	int n_module = d_var[6];
	double n_yr = f_var[0];

	double det_dim[3];
	for (int i=0; i<3; i++) det_dim[i] = dim_var[i];
	double cuts_e[9];
	for (int i=0; i<7; i++) cuts_e[i] = cuts_e_var[i];
	for (int i=0; i<2; i++) cuts_e[i+7] = trk_eff_var[i];
	double cuts_p[9];
	for (int i=0; i<7; i++) cuts_p[i] = cuts_p_var[i];
	for (int i=0; i<2; i++) cuts_p[i+7] = trk_eff_var[i];

	double m0 = m_couple_var[0];
	double m1 = m_couple_var[1];
	double m2 = m_couple_var[2];
	double g11 = m_couple_var[3]; // coupling for elastic scattering
	double g12 = m_couple_var[4]; // coupling for inelastic scattering

	double mp = 0.938;
	double me = 0.000511;
	double mT;
	if (type_target == 1) mT = me;
	else mT = mp;

	double rm1 = m2/m1; // mass ratio of heavier state to ligher DM

//	int scan_counter = 1;
	double maxm2 = get_max_m2 (m0, mT, m1);

	if (m2 > maxm2) {
		cout << "m2 is not kinematically accessible. Scanning is not performed." << endl;
	}
	else 
	{
		cout << "Scan started." << endl;
		dataout << m0 << " " << m1 << " " << m2 << " " << mT << " " << g12 << endl;

		double xSec; // variable for total cross section
		double mX, eps;

		for (int i = inum_in; i<inum; i++) 
		{
//		for (int i = 100; i<101; i++) {
			mX = pow(10,-2.899 + i * 0.1); // mass of dark photon in GeV

			// Below line is for the visible decay case
//			if (mX > 2*m1) mX = 2*m1*0.999;
			
//			for (int j = 0; j<jnum; j++)  // Recommeded jnum = 40 so that you scan 10^-7 to 10^-3
			for (int j = jnum; j<jnum+20; j++)  // Recommeded jnum = 40 so that you scan 10^-7 to 10^-3
			{
//			for (int j = 200; j<201; j++) {

//				if (scan_counter%50 == 0) cout << scan_counter << "th scanning..." << endl;

				eps = pow(10,-7 + j * 0.1); // kinetic mixing parameter
				if ( (m2-m1) <= me*2 && rm1 > 1.0001) {
//					cout << m2 << " " << m1 << " " << mX << " " << m2-m1-2*me << " case1" << endl;
					continue;
				}
				if ( (m2-m1) > mX && mX <= me*2 && rm1 > 1.0001) {
//					cout << m2 << " " << m1 << " " << mX << " " << mX - m2 + m1 << " " << mX - 2*me << " case2" << endl;
					continue;
				}

				double w_tot, w_passed;
	
				// Event generation module start from here.
				w_tot = 0.; // weight total
				w_passed = 0.; // weight total for passed events
				bool passed_flag;

				for (int eve = 0; eve < n_eve; eve++) {

					if ((eve+1)%1000 == 0) cout << eve+1 << "th event generating" << endl;
	
					double w_prod = 1, // weight product
					passed_flag = false;

					// Primary scattering part
					double e1 = m0;
					double p1 = sqrt(e1*e1 - m1*m1);
					double pT_x, pT_y, pT_z, pT_e, pT_mag, pT_cos, pT_phi;

					if (type_target==1) pT_e = get_ET_electron(e1, m1, mT, m2, mX);
					else pT_e = get_ET_proton(e1, m1, mT, m2, mX);

					pT_mag = sqrt(pT_e*pT_e - mT*mT);
					pT_cos = (m2*m2-m1*m1-2*mT*mT-2*e1*mT+2*pT_e*mT+2*e1*pT_e)/2/p1/pT_mag;
					pT_phi = 2*TMath::Pi()*my_phi->Rndm();
					pT_x = pT_mag*sqrt(1-pT_cos*pT_cos)*cos(pT_phi);
					pT_y = pT_mag*sqrt(1-pT_cos*pT_cos)*sin(pT_phi);
					pT_z = pT_mag*pT_cos;

					// Rotating ourgoing particles
					double cos_rot_p = -1 + 2.* my_cos_p->Rndm();
					double sin_rot_p = sqrt(1 - pow(cos_rot_p,2));
					double phi_rot_p = 2 * TMath::Pi() * my_phi_p->Rndm();

					double p_rec[6];

					p_rec[0] = ( pT_x * cos_rot_p + pT_z * sin_rot_p ) * cos(phi_rot_p) - pT_y * sin(phi_rot_p);
					p_rec[1] = ( pT_x * cos_rot_p + pT_z * sin_rot_p ) * sin(phi_rot_p) + pT_y * cos(phi_rot_p);
					p_rec[2] = - pT_x * sin_rot_p + pT_z * cos_rot_p;
					p_rec[3] = pT_e;

					double pinv_rot_x = p1 * sin_rot_p * cos(phi_rot_p) - p_rec[0];
					double pinv_rot_y = p1 * sin_rot_p * sin(phi_rot_p) - p_rec[1];
					double pinv_rot_z = p1 * cos_rot_p - p_rec[2];

//					TLorentzVector visT(pT_x, pT_y, pT_z, pT_e);
//					TLorentzVector intX2(-pT_x,-pT_y,p1-pT_z,e1+mT-pT_e);
					TLorentzVector visT(p_rec[0], p_rec[1], p_rec[2], p_rec[3]);
					TLorentzVector intX2(pinv_rot_x, pinv_rot_y, pinv_rot_z, e1+mT-pT_e);
					
					p_rec[4] = visT.M();
					p_rec[5] = visT.P();

//					TLorentzVector *vis1;
//					TLorentzVector *vis2;

					bool flag_scenario;

					if (rm1 < 1.0001) {
						// cut checks
						if (type_target==1) passed_flag = test_eelas_event(p_rec, det_type, det_dim, cuts_e);
						else passed_flag = test_pelas_event(p_rec, cuts_p); // Probably, proton does not fly much.
					}
					else if (m2 >= m1 + mX) // two-body decay module
					{
						flag_scenario = false;

						double p1_sec[6];
						double p2_sec[6];
						double p_LLP[6];

						// Decay of heavier state
						TGenPhaseSpace sd;
						double ss[2] = {m1,mX};
						sd.SetDecay(intX2,2,ss,"");
						double wsd = sd.Generate();
						TLorentzVector *intDM = sd.GetDecay(0);
						TLorentzVector *intA = sd.GetDecay(1);
						w_prod = w_prod * wsd; 

						// Decay of on-shell dark photon
						TGenPhaseSpace td;
						double ts[2] = {me, me};
						td.SetDecay(*intA,2,ts,"");
						double wtd = td.Generate();
						TLorentzVector *v1 = td.GetDecay(0);
						TLorentzVector *v2 = td.GetDecay(1);
//						vis1 = v1;
//						vis2 = v2;
						w_prod = w_prod * wtd;

						p1_sec[0] = v1->Px();
						p1_sec[1] = v1->Py();
						p1_sec[2] = v1->Pz();
						p1_sec[3] = v1->E();
						p1_sec[4] = v1->M();
						p1_sec[5] = v1->P();

						p2_sec[0] = v2->Px();
						p2_sec[1] = v2->Py();
						p2_sec[2] = v2->Pz();
						p2_sec[3] = v2->E();
						p2_sec[4] = v2->M();
						p2_sec[5] = v2->P();

						p_LLP[0] = intA->Px();
						p_LLP[1] = intA->Py();
						p_LLP[2] = intA->Pz();
						p_LLP[3] = intA->E();
						p_LLP[4] = intA->M();
						p_LLP[5] = intA->P();

//						cout << "----- " << p_rec[3] << " " << p1_sec[3] << " " << p2_sec[3] << " " << intDM->E() << " -----" << endl;
						// cut checks
						if (type_target==1) passed_flag = test_erecoil_event(p_rec, p1_sec, p2_sec, p_LLP, det_type, det_dim, cuts_e, flag_scenario, mX, m1, m2, eps, g12);
						else passed_flag = test_precoil_event(p_rec, p1_sec, p2_sec, p_LLP, det_type, det_dim, cuts_p, flag_scenario, mX, m1, m2, eps, g12);

					}
					else // Decay of heavier state 
					{
						flag_scenario = true;

						double p1_sec[6];
						double p2_sec[6];
						double p_LLP[6];

						TGenPhaseSpace sd;
						double ss[3] = {m1, me, me};
						sd.SetDecay(intX2,3,ss,"");
						double wsd = sd.Generate();
						TLorentzVector *intDM = sd.GetDecay(0);
						TLorentzVector *v1 = sd.GetDecay(1);
						TLorentzVector *v2 = sd.GetDecay(2);
//						vis1 = v1;
//						vis2 = v2;
						double s1 = pow((*intDM + *v1).M(),2);
						double s2 = pow((*intDM + *v2).M(),2);
//						w_prod = w_prod * wsd;
						w_prod = w_prod * wsd * secondaryMECalc(s1, s2, m1, m2, mX, me, 1., 1.); // Couplings are set to be unity because we are interested in relative weights.

						p1_sec[0] = v1->Px();
						p1_sec[1] = v1->Py();
						p1_sec[2] = v1->Pz();
						p1_sec[3] = v1->E();
						p1_sec[4] = v1->M();
						p1_sec[5] = v1->P();

						p2_sec[0] = v2->Px();
						p2_sec[1] = v2->Py();
						p2_sec[2] = v2->Pz();
						p2_sec[3] = v2->E();
						p2_sec[4] = v2->M();
						p2_sec[5] = v2->P();

						p_LLP[0] = intX2.Px();
						p_LLP[1] = intX2.Py();
						p_LLP[2] = intX2.Pz();
						p_LLP[3] = intX2.E();
						p_LLP[4] = intX2.M();
						p_LLP[5] = intX2.P();

//						cout << "----- " << p_rec[3] << " " << p1_sec[3] << " " << p2_sec[3] << " " << intDM->E() << " -----" << endl;
						// cut checks
						if (type_target==1) passed_flag = test_erecoil_event(p_rec, p1_sec, p2_sec, p_LLP, det_type, det_dim, cuts_e, flag_scenario, mX, m1, m2, eps, g12);
						else passed_flag = test_precoil_event(p_rec, p1_sec, p2_sec, p_LLP, det_type, det_dim, cuts_p, flag_scenario, mX, m1, m2, eps, g12);
					}

					if (passed_flag) w_passed += w_prod;
					w_tot += w_prod;

				} // loop for eve

				double acceptance = w_passed / w_tot;
				if (rm1 < 1.0001) {
					// Module for calculating cross section, calculated later
					if (type_target == 1) xSec = get_xSec_electron(m0, m1, mT, m2, mX, eps, g11);
					else xSec = get_xSec_proton(m0, m1, mT, m2, mX, eps, g11);
					dataout << mX << " " << eps << " " << xSec << " " << acceptance << endl;
					cout << mX << " " << eps << " " << xSec << " " << acceptance << " " << get_Neve(xSec,acceptance,n_yr,m0,n_module) << endl;
				}
				else {
					// Module for calculating cross section, calculated later
					if (type_target == 1) xSec = get_xSec_electron(m0, m1, mT, m2, mX, eps, g12);
					else xSec = get_xSec_proton(m0, m1, mT, m2, mX, eps, g12);

					//////////////////////////////////////////////////////////
					/// Below is for test only                             ///
					//////////////////////////////////////////////////////////
/*
					// Module for calculating decay width
					double widthLonglived, maxGamma;
					if (m2 < (m1 + mX)) {
						widthLonglived = width_chi_total(1.,eps,mX,m2,m1);
						maxGamma = get_gammaH(m0, m1, mT, m2);
					}
					else {
						widthLonglived = width_X_total(eps, mX);
						maxGamma = get_gamma(m0, m1, mT, m2, mX);
					}
			
					// Module for calculating decay length
					double maxLength = maxGamma * 2.99792458*pow(10,10)*6.57895*pow(10,-25)/widthLonglived;
					double accept_lmax = get_accept_lmax(maxLength);

					cout << accept_lmax << " ";
*/
					//////////////////////////////////////////////////////////
					/// End of the test                                    ///
					//////////////////////////////////////////////////////////

					dataout << mX << " " << eps << " " << xSec << " " << acceptance << endl;
					cout << mX << " " << eps << " " << xSec << " " << acceptance << " " << get_Neve(xSec,acceptance,n_yr,m0,n_module) << endl;

				}
				
//				scan_counter++;

///////////////////////////////////////////////////////
// The default for get_Neve is 40 ton.yr.            //
// For other statistics, get_Meve should be changed. //
///////////////////////////////////////////////////////

				if (get_Neve(xSec,acceptance,n_yr,m0,n_module) > 2.3) 
				{
					jnum = j-1;
					break;
				}

			} // loop for variable j

		cout << "mX = " << mX << " scanning done!" << endl;

		} // loop for variable i

		cout << "Scan done!" << endl;
	}
	
	dataout.close();

	return 0;

}
