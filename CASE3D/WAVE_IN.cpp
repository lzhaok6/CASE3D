//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include "data.h"

int loc;
double** WAVE_IN(int NNODE, double** GCOORD, double* T, int TIME, double** PIN, double DT, double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO, std::vector<int> shadow_pts) {
	//In this code, we've only coded the spherical wave case since you don't have plane wave in real application
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern STRU_WET_SURF ss[ssnumber];

	int i, j, k, z;
	double R = 0.0; 
	double RO = 0.0; 
	int elenode2D_gs;
	if (mappingalgo == 4 || mappingalgo == 5) {
		if (element_type == 0) {
			elenode2D_gs = (hprefg_flu + 1)*(hprefg_flu + 1);
		}
		else if (element_type == 1) {
			elenode2D_gs = 3;
		}
	}
	else {
		if (element_type == 0) {
			elenode2D_gs = NINT*NINT;
		}
		else if (element_type == 1) {
			elenode2D_gs = 3;
		}
	}

	if (WAVE == 1) {
		if (step == 0) {
			if (TIME == 1) {
				for (j = 0; j < NNODE; j++) {
					//if ((T[TIME - 1] + (abs(GCOORD[j][1]) - xo) / C) >= 0) {
					//if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) >= 0) {
					if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) > -1e-6) {
						//PIN[j][0] = PPEAK*exp(-(T[TIME - 1] + ((abs(GCOORD[j][1]) - xo) / C)) / TAU);
						PIN[j][0] = PPEAK*exp(-(T[TIME - 1] + (((-xo) - GCOORD[j][1]) / C)) / TAU);

					}
					else {
						PIN[j][0] = 0.0;
					}
				}
			}
			else { //TIME>0
				for (j = 0; j < NNODE; j++) {
					//if ((T[TIME - 1] + (abs(GCOORD[j][1]) - xo) / C) >= 0) {
					//if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) >= 0) {
					if ((T[TIME - 1] + ((-xo) - GCOORD[j][1]) / C) > -1e-6) {
						//PIN[j][1] = PPEAK*exp(-(T[TIME - 1] + ((abs(GCOORD[j][1]) - xo) / C)) / TAU);
						PIN[j][1] = PPEAK*exp(-(T[TIME - 1] + (((-xo) - GCOORD[j][1]) / C)) / TAU);
						//std::cout << " " << std::endl;
					}
					else {
						PIN[j][1] = 0.0;
					}
				}
			}
		}
		else { //step wave with no exponential decay
			if (TIME == 1) {
				for (j = 0; j < NNODE; j++) {
					if ((T[TIME - 1] + (-xo + GCOORD[j][0]) / C) > -1e-6) {
						PIN[j][0] = PPEAK;
					}
					else {
						PIN[j][0] = 0.0;
					}
				}
			}
			else { //TIME>0
				for (j = 0; j < NNODE; j++) {
					if ((T[TIME - 1] + (-xo + GCOORD[j][0]) / C) > -1e-6) {
						PIN[j][1] = PPEAK;
					}
					else {
						PIN[j][1] = 0.0;
					}
				}
			}
		}
	}
	else if (WAVE == 2) { //Spherical wave
		if (TIME == 1) {
			RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
			//the distance between reference point and explosion point
			//However, in spherical wave case, SONODE should be the point that spherical wave reach first when contacting with the free surface.  
			for (j = 0; j < NNODE; j++) {
				R = sqrt(pow((GCOORD[j][2] - ZC), 2) + pow((GCOORD[j][0] - XC), 2) + pow((GCOORD[j][1] - YC), 2));
				//if ((T[TIME - 1]) - (R - RO) / C >= 0) { //error prone (please refer to the comment in PWA)
				if (T[TIME - 1] - (R - RO) / C > -1e-6) {
					PIN[j][0] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
				}
				else {
					PIN[j][0] = 0.0;
				}
			}
			//std::cout << " " << std::endl;
		}
		else { //TIME>1
			RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
			for (j = 0; j < NNODE; j++) {
				R = sqrt(pow((GCOORD[j][2] - ZC), 2) + pow((GCOORD[j][0] - XC), 2) + pow((GCOORD[j][1] - YC), 2));
				//if ((T[TIME - 1]) - (R - RO) / C >= 0) {
				if (T[TIME - 1] - (R - RO) / C > -1e-6) {
					PIN[j][1] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
				}
				else {
					PIN[j][1] = 0.0;
				}
			}
		}

		if (tfm == 0) { //Scattered field model. Redefine the PIN in the shadow region
			for (i = 0; i < shadow_pts.size(); i++) {
				PIN[shadow_pts[i] - 1][0] = 0.0;
				PIN[shadow_pts[i] - 1][1] = 0.0;
			}
		}


		//if 
		if (tfm == 0 && (mappingalgo == 4 || mappingalgo == 5)) {
			if (incidentdisp_on_fluid == 1) {
				if (TIME == 1) {
					RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
					for (z = 0; z < owsfnumber; z++) {
						for (j = 0; j < ol[z].FSNEL; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
							for (k = 0; k < elenode2D_gs; k++) {
								R = sqrt(pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][2] - ZC), 2) + pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][0] - XC), 2) + pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][1] - YC), 2));
								if (T[TIME - 1] - (R - RO) / C > -1e-6) {
									ol[z].PIN_gs[j*elenode2D_gs + k][0] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
								}
								else {
									ol[z].PIN_gs[j*elenode2D_gs + k][0] = 0.0;
								}
							}
						}
					}
				}
				else { //TIME>1
					RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
					for (z = 0; z < owsfnumber; z++) {
						for (j = 0; j < ol[z].FSNEL; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
							for (k = 0; k < elenode2D_gs; k++) {
								R = sqrt(pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][2] - ZC), 2) + pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][0] - XC), 2) + pow((ol[z].GCOORD_flu_gs[j*elenode2D_gs + k][1] - YC), 2));
								if (T[TIME - 1] - (R - RO) / C > -1e-6) {
									ol[z].PIN_gs[j*elenode2D_gs + k][1] = PPEAK*(RO / R)*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
								}
								else {
									ol[z].PIN_gs[j*elenode2D_gs + k][1] = 0.0;
								}
							}
						}
					}
				}
			}
			else {
				double A = 0.0;
				double B = 0.0;
				double K = 0.0;
				if (TNT == 1) { //TNT (1.60 g/cc) Ref: Geers and Hunter (2002)
					A = 0.13;
					B = 0.23;
					K = 52.4e6;
				}
				else { //HBX-1 (1.72 g/cc) Ref: Ref: Geers and Hunter (2002)
					A = 0.15;
					B = 0.29;
					K = 56.7e6;
				}
				if (TIME == 1) {
					RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
					for (z = 0; z < ssnumber; z++) {
						for (j = 0; j < ss[z].Node_stru; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
							R = sqrt(pow((wsflist[z]->nodecoord[3 * j + 2] - ZC), 2) + pow((wsflist[z]->nodecoord[3 * j + 0] - XC), 2) + pow((wsflist[z]->nodecoord[3 * j + 1] - YC), 2));
							if (Colewave == 1) {
								PPEAK = K*pow(R / pow(W*0.453592, 1.0 / 3.0), -(1 + A)); //pa
								TAU = pow(W*0.453592, 1.0 / 3.0)*0.084*pow(pow(W*0.453592, 1.0 / 3.0) / R, -B) / 1000; //sec
							}
							if (T[TIME - 1] - (R - RO) / C > -1e-6) {
								ss[z].PIN[j][0] = PPEAK*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
							}
							else {
								ss[z].PIN[j][0] = 0.0;
							}
						}
					}
				}
				else { //TIME>1
					RO = sqrt(pow((ZC - ZO), 2) + pow((YC - YO), 2) + pow((XC - XO), 2));
					for (z = 0; z < ssnumber; z++) {
						for (j = 0; j < ss[z].Node_stru; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
							R = sqrt(pow((wsflist[z]->nodecoord[3 * j + 2] - ZC), 2) + pow((wsflist[z]->nodecoord[3 * j + 0] - XC), 2) + pow((wsflist[z]->nodecoord[3 * j + 1] - YC), 2));
							if (Colewave == 1) {
								PPEAK = K*pow(R / pow(W*0.453592, 1.0 / 3.0), -(1 + A)); //pa
								TAU = pow(W*0.453592, 1.0 / 3.0)*0.084*pow(pow(W*0.453592, 1.0 / 3.0) / R, -B) / 1000; //sec
							}
							if (T[TIME - 1] - (R - RO) / C > -1e-6) {
								ss[z].PIN[j][1] = PPEAK*exp(-((T[TIME - 1]) - (R - RO) / C) / TAU);
							}
							else {
								ss[z].PIN[j][1] = 0.0;
							}
						}
					}
				}
				std::cout << " " << std::endl; 
			}
		}
	}
	return (PIN);
}