#include "width.h"

using namespace std;

double *get_primary_point (int i, double *dim) { // The function for getting a primary scattering point
	
	TRandom3 *my_random = new TRandom3();
	my_random->SetSeed(0);

	static double cor_p[3];

//	cout << "pointer value check: " << dim[0] << endl;

	if (i == 1) {
		double x_gen = dim[0] * my_random->Rndm();
		double y_gen = dim[1] * my_random->Rndm();
		double z_gen = dim[2] * my_random->Rndm();
		cor_p[0] = x_gen;
		cor_p[1] = y_gen;
		cor_p[2] = z_gen;
	}
	else if (i == 2) {
		double r_gen = dim[1] * my_random->Rndm();
		double a_gen = 2. * TMath::Pi() * my_random->Rndm();
		double z_gen = 2. * dim[0] * my_random->Rndm() - dim[0];
		cor_p[0] = r_gen * cos(a_gen);
		cor_p[1] = r_gen * sin(a_gen);
		cor_p[2] = z_gen;
	}
	else {
		double r_gen = dim[0] * my_random->Rndm();
		double c_gen = 2. * my_random->Rndm() - 1.;
		double a_gen = 2. * TMath::Pi() * my_random->Rndm();
		cor_p[0] = r_gen * sqrt(1 - c_gen * c_gen) * cos(a_gen);
		cor_p[1] = r_gen * sqrt(1 - c_gen * c_gen) * sin(a_gen);
		cor_p[2] = r_gen * c_gen;
	}

	return cor_p;
}

bool test_fid (double *cor_p, double *dir_l, int i, double *dim, double len) {
	bool test_result = false;

	double x_s = cor_p[0] + dir_l[0] * len;
	double y_s = cor_p[1] + dir_l[1] * len;
	double z_s = cor_p[2] + dir_l[2] * len;

	if (i == 1) {
		if ( x_s < dim[0] && x_s > 0. && y_s < dim[1] && y_s > 0. && z_s < dim[2] && z_s > 0. ) test_result = true;
	}
	else if (i == 2) {
		double r_s = sqrt( x_s*x_s + y_s*y_s );
		if ( abs(z_s)<dim[0] && r_s < dim[1] ) test_result = true;
	}
	else {
		double r_s = sqrt( x_s*x_s + y_s*y_s + z_s*z_s );
		if ( r_s < dim[0] ) test_result = true;
	}

	return test_result;
}

double *get_end (double *cor_p, double *dir_l, double len) 
{
	
	static double coord_s[3];

	coord_s[0] = cor_p[0] + dir_l[0] * len;
	coord_s[1] = cor_p[1] + dir_l[1] * len;
	coord_s[2] = cor_p[2] + dir_l[2] * len;

	return coord_s;
}

double get_ang (double *v1, double *v2)
{
	return acos ( (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]) / v1[5] / v2[5] ); 
}

double get_eerr_electron (double e)
{
	// Input value e should be in GeV.
//	return 0.01;
//	return 2. + 15/sqrt(e);
	double val;
	if (e < 0.4) val = 20.;
	if (e < 1.0) val = 10.;
	else val = 2. + 8/sqrt(e);
	return val;
//	val = 0.0000000001;
//	return val;
	// The return value in %.
}

double get_eerr_proton (double e)
{
	// Input value p should be in GeV.
	double val;

	if (e < 0.4) val = 10.;
//	else val = 5. + 30/sqrt(e);
	else val = 5. + 5/sqrt(e);
	
	return val;
	// The return value in %.
}

double get_travel_len (double e)
{
	// The following empirical model is from a fit. 
	return 13.0436 * log(1 + 32.8561 * e);
//	return 0.0001;
	// A length is returned in cm.
}

double get_travel_len_p (double e)
{
	// The following empirical model is from a fit. 
	return -5225.27*pow(e,11)+30369.7*pow(e,10)-77405.7*pow(e,9)+113788.*pow(e,8)-106827.*pow(e,7)+67149.*pow(e,6)-29026.9*pow(e,5)+8992.89*pow(e,4)-2366.26*pow(e,3)+842.307*pow(e,2)+5.73883*e;
//	return 0.0001;
	// A length is returned in cm.
}

bool test_eelas_event (double *visT, int det_type, double *det_dim, double *cuts)
{
	TRandom3 *my_err = new TRandom3();
	my_err->SetSeed(0); // random seed
	TRandom3 *my_eff = new TRandom3();
	my_eff->SetSeed(0); // random seed

	bool test = true;

	double e_rec = cuts[0]; // energy threshold for recoil proton
	double eff_rec = cuts[2]; // efficiency for recoil electron

	double par_e; // particle energy
	double err_e; // energy error
	double par_eff; // efficiency test number
	err_e = 0.01 * visT[3] * get_eerr_electron(visT[3]); // 0.01 converts % to fraction.
	par_e = my_err->Gaus(visT[3],err_e); // smeared energy
	if (eff_rec < 1) par_eff = my_eff->Rndm();
	else par_eff = 0.;

	if ( par_e < e_rec || par_eff > eff_rec ) test = false;
	else
	{
		// primary scattering coordinate
		double *coord_p;
		coord_p = get_primary_point(det_type,det_dim);

		double travel_len_p = get_travel_len( visT[3] );
		double dir_p[3] = {visT[0]/visT[5], visT[1]/visT[5], visT[2]/visT[5]};
		bool p_end_in_fid = test_fid(coord_p, dir_p, det_type, det_dim, travel_len_p);

//		p_end_in_fid = true; // This should be commened out if this condition is not tested. 
		if (!p_end_in_fid) test = false; // fully-contained?
	}
	
	return test;
}

bool test_pelas_event (double *visT, double *cuts)
{
	TRandom3 *my_err = new TRandom3();
	my_err->SetSeed(0); // random seed
	TRandom3 *my_eff = new TRandom3();
	my_eff->SetSeed(0); // random seed

	bool test = true;

	double e_rec = cuts[0]; // energy threshold for recoil proton
	double e_rec_max = 2.; // max energy for recoil proton
	double eff_rec = cuts[2]; // efficiency for recoil proton

	double par_e; // particle energy
	double err_e; // energy error
	double par_eff; // efficiency test number
	err_e = 0.01 * visT[3] * get_eerr_proton(visT[3]); // 0.01 converts % to fraction.
	par_e = my_err->Gaus(visT[3],err_e); // smeared energy
	if (eff_rec < 1) par_eff = my_eff->Rndm();
	else par_eff = 0.;

	if ( par_e < e_rec || visT[3] > e_rec_max || par_eff > eff_rec ) test = false;
	
	return test;
}

bool test_erecoil_event (double *visT, double *vis1, double *vis2, double *LLP, int det_type, double *det_dim, double *cuts, bool flag_scenario, double mX, double m1, double m2, double eps, double g12)
{
	double vc = 2.99792458*pow(10,10); // speed of light in cm/s

	TRandom3 *my_err = new TRandom3();
	my_err->SetSeed(0); // random seed
	TRandom3 *my_eff = new TRandom3();
	my_eff->SetSeed(0); // random seed
	TRandom3 *my_time = new TRandom3();
	my_time->SetSeed(0); // random seed
	TRandom3 *my_trk = new TRandom3();
	my_trk->SetSeed(0); // random seed

	bool test = true;

	double e_rec = cuts[0]; // energy threshold for recoil electron
	double e_sec = cuts[1]; // energy threshold for secondary electron
	double eff_rec = cuts[2]; // efficiency for recoil electron
	double eff_sec = cuts[3]; // efficiency for secondary electron
	double ang_rec = cuts[4]; // angular resolution for recoil electron
	double ang_sec = cuts[5]; // angular resolution for secondary electron
	double vert_res = cuts[6]; // vertex resolution

	double twotrk_miss = cuts[7]; // photon-like track miss identification probability
	double threetrk_miss = cuts[8]; // 3e-like track miss identification probability
	
	double par_e; // particle energy
	double err_e; // energy error
	double par_eff; // efficiency test number
	err_e = 0.01 * visT[3] * get_eerr_electron(visT[3]); // 0.01 converts % to fraction.
	par_e = my_err->Gaus(visT[3],err_e); // smeared energy
//	par_e = visT[3];
	if (eff_rec < 1) par_eff = my_eff->Rndm();
	else par_eff = 0.; 

	if ( par_e < e_rec || par_eff > eff_rec ) // To test whether or not the smeared recoil electron energy is greater than the threshold and also the electron is accepted. 
	{
//		cout << visT[3] << " " << par_e << " " << par_eff << " ";
//		cout << "1st requirement failed." << endl;
		test = false;
	}
	else
	{
		err_e = 0.01 * vis1[3] * get_eerr_electron(vis1[3]); // 0.01 converts % to fraction.
		par_e = my_err->Gaus(vis1[3],err_e); // smeared energy
//		par_e = vis1[3];
		if (eff_sec < 1) par_eff = my_eff->Rndm();
		else par_eff = 0.;

		if ( par_e < e_sec || par_eff > eff_sec ) // To test the smeared energy and the acceptance of the 1st decay product.
		{
//			cout << vis1[3] << " " << par_e << " " << par_eff << " ";
//			cout << "2nd requirement failed." << endl;
			test = false;
		}
		else
		{
			err_e = 0.01 * vis2[3] * get_eerr_electron(vis2[3]); // 0.01 converts % to fraction.
			par_e = my_err->Gaus(vis2[3],err_e); // smeared energy
//			par_e = vis2[3];
			if (eff_sec < 1) par_eff = my_eff->Rndm();
			else par_eff = 0.;

			if ( par_e < e_sec || par_eff > eff_sec ) // To test the smeared energy and the acceptance of the 2nd decay product.
			{
//				cout << vis2[3] << " " << par_e << " " << par_eff << " ";
//				cout << "3rd requirement failed." << endl;
				test = false;
			}
			else
			{

				// Module for calculating decay width
				double width_LLP; 
				if (flag_scenario) width_LLP = width_chi_total(g12,eps,mX,m2,m1); // flag_scenario = true -> chi2 is long-lived.
				else width_LLP = width_X_total(eps, mX);
				
				double gamma_LLP = LLP[3]/LLP[4];
				double tau_LLP_lab = gamma_LLP * 6.57895*pow(10,-25) / width_LLP;
				double t_LLP = my_time->Exp(tau_LLP_lab);
				double len_disp = vc * t_LLP * sqrt(pow(gamma_LLP,2)-1) / gamma_LLP; // decay length of long-lived particle in cm 
			
				// primary scattering coordinate
				double *coord_p;
				coord_p = get_primary_point(det_type,det_dim);

				double dir_long[3] = {LLP[0]/LLP[5], LLP[1]/LLP[5], LLP[2]/LLP[5]};
				bool sec_in_fid = true;
				if ( len_disp > vert_res ) sec_in_fid = test_fid(coord_p, dir_long, det_type, det_dim, len_disp);

				if (!sec_in_fid) // To test whether or not the secondary vertex is in the detector fiducial volume.
				{
//					cout << flag_scenario << ", " << gamma_LLP << ", " << width_LLP << ", " << tau_LLP_lab << ", " << t_LLP << ", " << len_disp << ", " << coord_p[0] << " " << coord_p[1] << " " << coord_p[2] << " ";
//					cout << "4th requirement failed." << endl;
					test = false; // whether or not the displaced vertex is within the fiducial volume
				}
				else 
				{

					double travel_len_p = get_travel_len( visT[3] );
					double dir_p[3] = {visT[0]/visT[5], visT[1]/visT[5], visT[2]/visT[5]};
					bool p_end_in_fid = test_fid(coord_p, dir_p, det_type, det_dim, travel_len_p);

//					p_end_in_fid = true; // This should be commened out if this condition is not tested. 
					if (!p_end_in_fid) 
					{
//						cout << "5th requirement failed." << endl;
						test = false; // fully-contained?
					}
					else
					{
						double *coord_sec;
						coord_sec = get_end(coord_p, dir_long, len_disp);

						double travel_len_vis1 = get_travel_len(vis1[3]);
						double dir_vis1[3] = {vis1[0]/vis1[5], vis1[1]/vis1[5], vis1[2]/vis1[5]};
						bool vis1_end_in_fid = test_fid(coord_sec, dir_vis1, det_type, det_dim, travel_len_vis1);

//						vis1_end_in_fid = true; // This should be commened out if this condition is not tested. 
						if (!vis1_end_in_fid) 
						{
//							cout << "6th requirement failed." << endl;
							test = false; // fully-contained?
						}
						else
						{
							double travel_len_vis2 = get_travel_len(vis2[3]);
							double dir_vis2[3] = {vis2[0]/vis2[5], vis2[1]/vis2[5], vis2[2]/vis2[5]};
							bool vis2_end_in_fid = test_fid(coord_sec, dir_vis2, det_type, det_dim, travel_len_vis2);

//							vis2_end_in_fid = true; // This should be commened out if this condition is not tested. 
							if (!vis2_end_in_fid) 
							{
//								cout << "7th requirement failed." << endl;
								test = false; // fully-contained?
							}
							else
							{
								test = true;

/*
								//////////////////////////////////////////////////////////
								// Simplified version of angular resolution criteria    //
								//////////////////////////////////////////////////////////
								double ang_ee = get_ang(vis1, vis2);
								double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) );
								if (ang_ee < ang_sec || ang_pe < ang_rec ) test = false;
*/

								//////////////////////////////////////////////////////////
								// Angular resolution criteria should be written below. //
								//////////////////////////////////////////////////////////
								double ang_ee = get_ang(vis1, vis2);

								if ( len_disp > vert_res ) // If the displaced vertex is really tagged as displaced!
								{ 
									if ( len_disp > 1.0 * travel_len_p ) // This implies that if the length of displaced vertex is greater than 0.9 of recoil electron, the secondary vertex is identified.
									{
										double twotrk_iden = my_trk->Rndm();
//										if (twotrk_iden < twotrk_miss )
										if (twotrk_iden < 0. )
										{
//											cout << "8th requirement failed." << endl;
											test = false;
										}
									}
									else
									{
										double ang_pLLP = get_ang(visT, LLP); // angle between recoil electron and LLP
										if ( ang_pLLP > 5 * ang_rec ) // This implies that if the secondary vertex is >5 deg away from the recoil electron, the secondary vertex is identified.
										{
											double twotrk_iden = my_trk->Rndm();
											if (twotrk_iden < twotrk_miss ) 
											{
//												cout << "9th requirement failed." << endl;
												test = false;
											}
										}
										else if ( len_disp > 0.2 * travel_len_p ) 
										{
//											cout << "10th requirement failed." << endl;
											test = false; // If the secondary vertex is <5 deg from the recoil electron and >20% of recoil electron track, the secondary vertex is buried, and three tracks are not resolved.
										}
										else
										{
											double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
											if ( ang_pe < ang_rec && ang_ee < ang_sec ) // three particles merged
											{
												double threetrk_iden = my_trk->Rndm();
												if (threetrk_iden < threetrk_miss ) 
												{
//													cout << "11th requirement failed." << endl;
													test = false;
												}
											}
											else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
											{
												double twotrk_iden = my_trk->Rndm();
												if (twotrk_iden < twotrk_miss ) 
												{
//													cout << "12th requirement failed." << endl;
													test = false;
												}
											}
											else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil electron
											{
												double twotrk_iden = my_trk->Rndm();
												if (twotrk_iden < twotrk_miss ) 
												{
//													cout << "13th requirement failed." << endl;
													test = false;
												}
											}
										}
									}
								}
								else // If three partciles are merged! 
								{
									double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
									if ( ang_pe < ang_rec && ang_ee < ang_sec ) // three particles merged
									{
										double threetrk_iden = my_trk->Rndm();
										if (threetrk_iden < threetrk_miss ) 
										{
//											cout << "14th requirement failed." << endl;
											test = false;
										}
									}
									else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "15th requirement failed." << endl;
											test = false;
										}
									}
									else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil electron
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "16th requirement failed." << endl;
											test = false;
										}
									}

								}
								//////////////////////////////////////////////////////////
								// Up to here.                                          //
								//////////////////////////////////////////////////////////

							}

						}
					}

				}

			}
		}
	}

	return test;
}

/*
bool test_erecoil_event (double *visT, double *vis1, double *vis2, double *LLP, int det_type, double *det_dim, double *cuts, bool flag_scenario, double mX, double m1, double m2, double eps, double g12)
{
	double vc = 2.99792458*pow(10,10); // speed of light in cm/s

	TRandom3 *my_err = new TRandom3();
	my_err->SetSeed(0); // random seed
	TRandom3 *my_eff = new TRandom3();
	my_eff->SetSeed(0); // random seed
	TRandom3 *my_time = new TRandom3();
	my_time->SetSeed(0); // random seed
	TRandom3 *my_trk = new TRandom3();
	my_trk->SetSeed(0); // random seed

	bool test = true;

	double e_rec = cuts[0]; // energy threshold for recoil electron
	double e_sec = cuts[1]; // energy threshold for secondary electron
	double eff_rec = cuts[2]; // efficiency for recoil electron
	double eff_sec = cuts[3]; // efficiency for secondary electron
	double ang_rec = cuts[4]; // angular resolution for recoil electron
	double ang_sec = cuts[5]; // angular resolution for secondary electron
	double vert_res = cuts[6]; // vertex resolution

	double twotrk_miss = cuts[7]; // photon-like track miss identification probability
	double threetrk_miss = cuts[8]; // 3e-like track miss identification probability
	
	double par_e; // particle energy
	double err_e; // energy error
	double par_eff; // efficiency test number
	err_e = 0.01 * visT[3] * get_eerr_electron(visT[3]); // 0.01 converts % to fraction.
	par_e = my_err->Gaus(visT[3],err_e); // smeared energy
	if (eff_rec < 1) par_eff = my_eff->Rndm();
	else par_eff = 0.; 

	if ( par_e < e_rec || par_eff > eff_rec ) // To test whether or not the smeared recoil electron energy is greater than the threshold and also the electron is accepted. 
	{
//		cout << visT[3] << " " << par_e << " " << par_eff << " ";
//		cout << "1st requirement failed." << endl;
		test = false;
	}
	else
	{
		err_e = 0.01 * vis1[3] * get_eerr_electron(vis1[3]); // 0.01 converts % to fraction.
		par_e = my_err->Gaus(vis1[3],err_e); // smeared energy
		if (eff_sec < 1) par_eff = my_eff->Rndm();
		else par_eff = 0.;

		if ( par_e < e_sec || par_eff > eff_sec ) // To test the smeared energy and the acceptance of the 1st decay product.
		{
//			cout << vis1[3] << " " << par_e << " " << par_eff << " ";
//			cout << "2nd requirement failed." << endl;
			test = false;
		}
		else
		{
			err_e = 0.01 * vis2[3] * get_eerr_electron(vis2[3]); // 0.01 converts % to fraction.
			par_e = my_err->Gaus(vis2[3],err_e); // smeared energy
			if (eff_sec < 1) par_eff = my_eff->Rndm();
			else par_eff = 0.;

			if ( par_e < e_sec || par_eff > eff_sec ) // To test the smeared energy and the acceptance of the 2nd decay product.
			{
//				cout << vis2[3] << " " << par_e << " " << par_eff << " ";
//				cout << "3rd requirement failed." << endl;
				test = false;
			}
			else
			{
				// Module for calculating decay width
				double width_LLP; 
				if (flag_scenario) width_LLP = width_chi_total(g12,eps,mX,m2,m1); // flag_scenario = true -> chi2 is long-lived.
				else width_LLP = width_X_total(eps, mX);
				
				double gamma_LLP = LLP[3]/LLP[4];
				double tau_LLP_lab = gamma_LLP * 6.57895*pow(10,-25) / width_LLP;
				double t_LLP = my_time->Exp(tau_LLP_lab);
				double len_disp = vc * t_LLP * sqrt(pow(gamma_LLP,2)-1) / gamma_LLP; // decay length of long-lived particle in cm 
			
				// primary scattering coordinate
				double *coord_p;
				coord_p = get_primary_point(det_type,det_dim);

				double dir_long[3] = {LLP[0]/LLP[5], LLP[1]/LLP[5], LLP[2]/LLP[5]};
				bool sec_in_fid = true;
				if ( len_disp > vert_res ) sec_in_fid = test_fid(coord_p, dir_long, det_type, det_dim, len_disp);

				if (!sec_in_fid) // To test whether or not the secondary vertex is in the detector fiducial volume.
				{
//					cout << flag_scenario << ", " << gamma_LLP << ", " << width_LLP << ", " << tau_LLP_lab << ", " << t_LLP << ", " << len_disp << ", " << coord_p[0] << " " << coord_p[1] << " " << coord_p[2] << " ";
//					cout << "4th requirement failed." << endl;
					test = false; // whether or not the displaced vertex is within the fiducial volume
				}
				else 
				{
					double travel_len_p = get_travel_len( visT[3] );
					double dir_p[3] = {visT[0]/visT[5], visT[1]/visT[5], visT[2]/visT[5]};
					bool p_end_in_fid = test_fid(coord_p, dir_p, det_type, det_dim, travel_len_p);

//					p_end_in_fid = true; // This should be commened out if this condition is not tested. 
					if (!p_end_in_fid) 
					{
//						cout << "5th requirement failed." << endl;
						test = false; // fully-contained?
					}
					else
					{
						double *coord_sec;
						coord_sec = get_end(coord_p, dir_long, len_disp);

						double travel_len_vis1 = get_travel_len(vis1[3]);
						double dir_vis1[3] = {vis1[0]/vis1[5], vis1[1]/vis1[5], vis1[2]/vis1[5]};
						bool vis1_end_in_fid = test_fid(coord_sec, dir_vis1, det_type, det_dim, travel_len_vis1);

//						vis1_end_in_fid = true; // This should be commened out if this condition is not tested. 
						if (!vis1_end_in_fid) 
						{
//							cout << "6th requirement failed." << endl;
							test = false; // fully-contained?
						}
						else
						{
							double travel_len_vis2 = get_travel_len(vis2[3]);
							double dir_vis2[3] = {vis2[0]/vis2[5], vis2[1]/vis2[5], vis2[2]/vis2[5]};
							bool vis2_end_in_fid = test_fid(coord_sec, dir_vis2, det_type, det_dim, travel_len_vis2);

//							vis2_end_in_fid = true; // This should be commened out if this condition is not tested. 
							if (!vis2_end_in_fid) 
							{
//								cout << "7th requirement failed." << endl;
								test = false; // fully-contained?
							}
							else
							{
								test = true;

								//////////////////////////////////////////////////////////
								// Angular resolution criteria should be written below. //
								//////////////////////////////////////////////////////////
								double ang_ee = get_ang(vis1, vis2);
								if ( len_disp > vert_res ) 
								{ 
									if ( len_disp > 0.9 * travel_len_p ) // This implies that if the length of displaced vertex is greater than 0.9 of recoil electron, the secondary vertex is identified.
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "8th requirement failed." << endl;
											test = false;
										}
									}
									else
									{
										double ang_pLLP = get_ang(visT, LLP); // angle between recoil electron and LLP
										if ( ang_pLLP > 5 * ang_rec ) // This implies that if the secondary vertex is >5 deg away from the recoil electron, the secondary vertex is identified.
										{
											double twotrk_iden = my_trk->Rndm();
											if (twotrk_iden < twotrk_miss ) 
											{
//												cout << "9th requirement failed." << endl;
												test = false;
											}
										}
										else if ( len_disp > 0.2 * travel_len_p ) 
										{
//											cout << "10th requirement failed." << endl;
											test = false; // If the secondary vertex is <5 deg from the recoil electron and >20% of recoil electron track, the secondary vertex is buried, and three tracks are not resolved.
										}
										else
										{
											double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
											if ( ang_pe < ang_rec && ang_ee < ang_sec ) // three particles merged
											{
												double threetrk_iden = my_trk->Rndm();
												if (threetrk_iden < threetrk_miss ) 
												{
//													cout << "11th requirement failed." << endl;
													test = false;
												}
											}
											else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
											{
												double twotrk_iden = my_trk->Rndm();
												if (twotrk_iden < twotrk_miss ) 
												{
//													cout << "12th requirement failed." << endl;
													test = false;
												}
											}
											else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil electron
											{
												double twotrk_iden = my_trk->Rndm();
												if (twotrk_iden < twotrk_miss ) 
												{
//													cout << "13th requirement failed." << endl;
													test = false;
												}
											}
										}
									}
								}
								else 
								{
									double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
									if ( ang_pe < ang_rec && ang_ee < ang_sec ) // three particles merged
									{
										double threetrk_iden = my_trk->Rndm();
										if (threetrk_iden < threetrk_miss ) 
										{
//											cout << "14th requirement failed." << endl;
											test = false;
										}
									}
									else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "15th requirement failed." << endl;
											test = false;
										}
									}
									else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil electron
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "16th requirement failed." << endl;
											test = false;
										}
									}

								}
								//////////////////////////////////////////////////////////
								// Up to here.                                          //
								//////////////////////////////////////////////////////////

							}

						}
					}

				}

			}
		}
	}

	return test;
}
*/

bool test_precoil_event (double *visT, double *vis1, double *vis2, double *LLP, int det_type, double *det_dim, double *cuts, bool flag_scenario, double mX, double m1, double m2, double eps, double g12)
{
	double vc = 2.99792458*pow(10,10); // speed of light in cm/s

	TRandom3 *my_err = new TRandom3();
	my_err->SetSeed(0); // random seed
	TRandom3 *my_eff = new TRandom3();
	my_eff->SetSeed(0); // random seed
	TRandom3 *my_time = new TRandom3();
	my_time->SetSeed(0); // random seed
	TRandom3 *my_trk = new TRandom3();
	my_trk->SetSeed(0); // random seed

	bool test = true;

	double e_rec = cuts[0]; // energy threshold for recoil proton
	double e_rec_max = 2.; // max energy for recoil proton
	double e_sec = cuts[1]; // energy threshold for secondary electron
	double eff_rec = cuts[2]; // efficiency for recoil proton
	double eff_sec = cuts[3]; // efficiency for secondary electron
	double ang_rec = cuts[4]; // angular resolution for recoil proton
	double ang_sec = cuts[5]; // angular resolution for secondary electron
	double vert_res = cuts[6]; // vertex resolution

	double twotrk_miss = cuts[7]; // photon-like track miss id probability
	
	double e_temp = pow(visT[5],2)/2/0.938;
	double par_e; // proton kinetic energy
	double err_e; // kinetic energy error
	double par_eff; // efficiency test number
	err_e = 0.01 * e_temp * get_eerr_proton(e_temp); // 0.01 converts % to fraction.
	par_e = my_err->Gaus(e_temp,err_e); // smeared energy
	if (eff_rec < 1) par_eff = my_eff->Rndm();
	else par_eff = 0.; 

	if ( par_e < e_rec || visT[3] > e_rec_max || par_eff > eff_rec ) test = false;
	else
	{
		err_e = 0.01 * vis1[3] * get_eerr_electron(vis1[3]); // 0.01 converts % to fraction.
		par_e = my_err->Gaus(vis1[3],err_e); // smeared energy
		if (eff_sec < 1) par_eff = my_eff->Rndm();
		else par_eff = 0.;

		if ( par_e < e_sec || par_eff > eff_sec ) test = false;
		else
		{
			err_e = 0.01 * vis2[3] * get_eerr_electron(vis2[3]); // 0.01 converts % to fraction.
			par_e = my_err->Gaus(vis2[3],err_e); // smeared energy
			if (eff_sec < 1) par_eff = my_eff->Rndm();
			else par_eff = 0.;

			if ( par_e < e_sec || par_eff > eff_sec ) test = false;
			else
			{
				// Module for calculating decay width
				double width_LLP; 
				if (flag_scenario) width_LLP = width_chi_total(g12,eps,mX,m2,m1); // flag_scenario = true -> chi2 is long-lived.
				else width_LLP = width_X_total(eps, mX);
				
				double gamma_LLP = LLP[3]/LLP[4];
				double tau_LLP_lab = gamma_LLP * 6.57895*pow(10,-25) / width_LLP;
				double t_LLP = my_time->Exp(tau_LLP_lab);
				double len_disp = vc * t_LLP * sqrt(pow(gamma_LLP,2)-1) / gamma_LLP; // decay length of long-lived particle in cm without taking epsilon into account
			
				// primary scattering coordinate
				double *coord_p;
				coord_p = get_primary_point(det_type,det_dim);

				double dir_long[3] = {LLP[0]/LLP[5], LLP[1]/LLP[5], LLP[2]/LLP[5]};
				bool sec_in_fid = true;
				if ( len_disp > vert_res ) sec_in_fid = test_fid(coord_p, dir_long, det_type, det_dim, len_disp);

				if (!sec_in_fid) test = false; // whether or not the displaced vertex is within the fiducial volume
				else 
				{
					double travel_len_p = get_travel_len_p( e_temp );
					double dir_p[3] = {visT[0]/visT[5], visT[1]/visT[5], visT[2]/visT[5]};
					bool p_end_in_fid = test_fid(coord_p, dir_p, det_type, det_dim, travel_len_p);

//					p_end_in_fid = true; // This should be commened out if this condition is not tested. 
					if (!p_end_in_fid) test = false; // fully-contained?
					else
					{
						double *coord_sec;
						coord_sec = get_end(coord_p, dir_long, len_disp);

						double travel_len_vis1 = get_travel_len(vis1[3]);
						double dir_vis1[3] = {vis1[0]/vis1[5], vis1[1]/vis1[5], vis1[2]/vis1[5]};
						bool vis1_end_in_fid = test_fid(coord_sec, dir_vis1, det_type, det_dim, travel_len_vis1);

//						vis1_end_in_fid = true; // This should be commened out if this condition is not tested. 
						if (!vis1_end_in_fid) test = false; // fully-contained?
						else
						{
							double travel_len_vis2 = get_travel_len(vis2[3]);
							double dir_vis2[3] = {vis2[0]/vis2[5], vis2[1]/vis2[5], vis2[2]/vis2[5]};
							bool vis2_end_in_fid = test_fid(coord_sec, dir_vis2, det_type, det_dim, travel_len_vis2);

//							vis2_end_in_fid = true; // This should be commened out if this condition is not tested. 
							if (!vis2_end_in_fid) test = false; // fully-contained?
							else
							{
								//////////////////////////////////////////////////////////
								// Angular resolution criteria should be written below. //
								//////////////////////////////////////////////////////////
								double ang_ee = get_ang(vis1, vis2);

								if ( len_disp > vert_res ) // If the displaced vertex is really tagged as displaced!
								{ 
									if ( len_disp > 1.0 * travel_len_p ) // This implies that if the length of displaced vertex is greater than 0.9 of recoil electron, the secondary vertex is identified.
									{
										double twotrk_iden = my_trk->Rndm();
//										if (twotrk_iden < twotrk_miss )
										if (twotrk_iden < 0. )
										{
//											cout << "8th requirement failed." << endl;
											test = false;
										}
									}
									else
									{
										double ang_pLLP = get_ang(visT, LLP); // angle between recoil electron and LLP
										if ( ang_pLLP > 5 * ang_rec ) // This implies that if the secondary vertex is >5 deg away from the recoil electron, the secondary vertex is identified.
										{
											double twotrk_iden = my_trk->Rndm();
											if (twotrk_iden < twotrk_miss ) 
											{
//												cout << "9th requirement failed." << endl;
												test = false;
											}
										}
										else if ( len_disp > 0.2 * travel_len_p ) 
										{
//											cout << "10th requirement failed." << endl;
											test = false; // If the secondary vertex is <5 deg from the recoil electron and >20% of recoil electron track, the secondary vertex is buried, and three tracks are not resolved.
										}
										else
										{
											double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
											if ( ang_pe < ang_rec && ang_ee < ang_sec ) // proton and two electrons merged
											{
												test = false;
											}
											else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
											{
												double twotrk_iden = my_trk->Rndm();
												if (twotrk_iden < twotrk_miss ) 
												{
//													cout << "12th requirement failed." << endl;
													test = false;
												}
											}
											else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil proton
											{
												test = false;
											}
										}
									}
								}
								else // If three partciles are merged! 
								{
									double ang_pe = min( get_ang(visT, vis2), get_ang(vis1, visT) ); // smaller of angles between recoil electron and secondary electrons
									if ( ang_pe < ang_rec && ang_ee < ang_sec ) // proton and two electrons merged
									{
										test = false;
									}
									else if ( ang_pe > ang_rec && ang_ee < ang_sec ) // two secondary particles merged
									{
										double twotrk_iden = my_trk->Rndm();
										if (twotrk_iden < twotrk_miss ) 
										{
//											cout << "15th requirement failed." << endl;
											test = false;
										}
									}
									else if ( ang_pe < ang_rec && ang_ee > ang_sec ) // one secondary particles merged with recoil proton
									{
										test = false;
									}

								}
								//////////////////////////////////////////////////////////
								// Up to here.                                          //
								//////////////////////////////////////////////////////////
							}

						}
					}

				}

			}
		}
	}

	return test;
}

////////////////////////////////////////////
/// Below is for test only               ///
////////////////////////////////////////////
double get_accept_lmax (double lmax) {
	double acc;
	if (lmax > 2500) acc = 0.00138474*(lmax - 2500) + 4.32541;
	else acc = 0.0000000386333 * lmax*lmax + ((4.32541 - 1)/2500 - 0.0000000386333 * 2500)*lmax + 1; 
	return 1/acc;
}
