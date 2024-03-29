#include "header.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "data.h"
#include "adapter.h"
#include "mpcci.h"
#include <vector>
#include <string>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <algorithm> 


//NRB determines the NRB local node numbering and the associated NRB arrays
void TIME_INT(int NNODE, double** GCOORD, int***LNA_3D, int**IEN, int NEL, int TIME, double *T, double DT, int NDT,
	double*** HMASTER, double* Q, double KAPPA, double PPEAK, double TAU, double XC, double YC,
	double ZC, double XO, double YO, double ZO, double ***SHOD, /*double** gamman, double** gamma_tn*/ double gamman[], double gamma_tn[],
	double****Gn, double****SHG, double****gamma_t, double ****gamma, double*****G) {
	int h, i, j, k, q, z, ii, jj, kk, m;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber];

	int* counter1;
	int* counter2;
	int* counter3;
	int* counter4;
	int* counter5;
	int* counter6;
	int* LNAct1; int* LNAct2; int* LNAct3;
	double* SHOD1;
	int*IENct1; int*IENct2; int*IENct3;

	if (tensorfactorization == 1) {
		int ctt1 = 0;
		counter1 = new int[NINT*NINT*NINT];
		counter2 = new int[NINT*NINT*NINT];
		counter3 = new int[NINT*NINT*NINT];
		counter4 = new int[NINT*NINT*NINT];
		counter5 = new int[NINT*NINT*NINT];
		counter6 = new int[NINT*NINT*NINT];
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					counter1[ctt1] = i*NINT*NINT + j*NINT + k;
					counter2[ctt1] = j*NINT*NINT + i*NINT + k;
					counter3[ctt1] = j*NINT*NINT + k*NINT + i;
					counter4[ctt1] = k*NINT*NINT + j*NINT + i;
					counter5[ctt1] = k*NINT*NINT + i*NINT + j;
					counter6[ctt1] = i*NINT*NINT + k*NINT + j;
					ctt1 += 1;
				}
			}
		}
		int ctt2 = 0;
		LNAct1 = new int[NINT*NINT*NINT];
		LNAct2 = new int[NINT*NINT*NINT];
		LNAct3 = new int[NINT*NINT*NINT];
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					LNAct1[ctt2] = LNA_3D[k][i][j] - 1;
					LNAct2[ctt2] = LNA_3D[i][k][j] - 1;
					LNAct3[ctt2] = LNA_3D[i][j][k] - 1;
					ctt2 += 1;
				}
			}
		}
		int ctt3 = 0;
		int ctt4 = 0;
		int ctt5 = 0;
		SHOD1 = new double[NINT*NINT];
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				SHOD1[i*NINT + j] = SHOD[1][i][j];
			}
		}
		ctt1 = 0;
		IENct1 = new int[NINT*NINT*NINT*NEL];
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				IENct1[ctt1] = IEN[LNAct1[j]][i] - 1;
				ctt1 += 1;
			}
		}
		ctt2 = 0;
		IENct2 = new int[NINT*NINT*NINT*NEL];
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				IENct2[ctt2] = IEN[LNAct2[j]][i] - 1;
				ctt2 += 1;
			}
		}
		ctt3 = 0;
		IENct3 = new int[NINT*NINT*NINT*NEL];
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				IENct3[ctt3] = IEN[LNAct3[j]][i] - 1;
				ctt3 += 1;
			}
		}
	}

	NRBstruct nrb;
	nrb = NRB(NNODE, GCOORD, LNA_3D);
	//clean nr[z].Jacob_2D
	for (z = 0; z < nrbsurfnumber; z++) {
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			delete[] nr[z].Jacob_2D[i];
		}
		delete[] nr[z].Jacob_2D;
	}

	//time history record
	std::clock_t start;
	double duration;
	//double duration_int;
	start = std::clock();
	printf("3DFrigate - computation\n");
	fflush(stdout);
	//Obtain system time information:
	auto ti = std::time(nullptr);
	auto tm = *std::localtime(&ti);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
	auto timestr = oss.str();

	int elenode3D = 0;
	int elenode2D = 0;
	if (element_type == 0) { //hex element
		elenode3D = NINT*NINT*NINT;
		elenode2D = NINT*NINT;
	}
	if (element_type == 1) { //tet element
		if (N == 1) {
			elenode3D = 4;
			elenode2D = 3;
		}
		else {
			std::cout << "High-order tet element is not supported yet" << std::endl;
			system("PAUSE ");
		}
	}

	//dynamic variables initiation
	double *DSDOT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF CONDENSATION
	double **ds;   //SOLUTION ARRAY FOR CONDENSATION
	double **FEEDOT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF DISP. POTENTIAL
	double **FEEDOT_inc;
	double **FEE;  //SOLUTION ARRAY FOR DISP. POTENTIAL
	double **P; //SOLUTION ARRAY FOR DYNAMIC/SCATTERED PRESSURE
	double **PT; //SOLUTION ARRAY FOR TOTAL PRESSURE
	double *PH;
	double** PIN;
	double *FFORCE; //ARRAY OF INTERNAL + EXTERNAL FORCE ON FLUID
	double *BNRB; //SOLUTION ARRAY FOR FORCE OF NRBC ON FLUID
	double *HF; //GLOBAL REACTANCE MATIRX
	double **HFTEMP;
	HFTEMP = new double*[NEL];
	for (i = 0; i < NEL; i++) {
		HFTEMP[i] = new double[elenode3D];
	}

	double **HFTEMPn;
	if (tensorfactorization == 1) {
		HFTEMPn = new double*[NEL];
		for (i = 0; i < NEL; i++) {
			HFTEMPn[i] = new double[elenode3D];
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < elenode3D; j++) {
				HFTEMPn[i][j] = 0.0;
			}
		}
	}

	//double FEETEMP[elenode3D]; //LOCAL DISP. POTENTIAL
	double* FEETEMP;
	FEETEMP = new double[elenode3D];

	//double *DPS_ukn;
	double *XNRBORG; //SOLUTION ARRAY FOR NRBC DISPLACEMENT AT T=0 (ORG MEANS ORIGIN)
	//double *XCOR; //SOLUTION ARRAY FOR NRB
	double *XCOR_kn;
	//double *XCOR_ukn;
	double** PSI_inc;
	//double *DPS;
	//==============================================================================//
	DSDOT = new double[NNODE];
	ds = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		ds[i] = new double[3];
	}
	FEEDOT = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEEDOT[i] = new double[2];
	}
	FEEDOT_inc = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEEDOT_inc[i] = new double[2];
	}
	FEE = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEE[i] = new double[2];
	}
	P = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		P[i] = new double[2];
	}
	PT = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		PT[i] = new double[2];
	}
	PH = new double[NNODE];
	PIN = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		PIN[i] = new double[2];
	}
	FFORCE = new double[NNODE];
	BNRB = new double[NNODE];
	HF = new double[NNODE];
	//DPS_ukn = new double[NNODE];

	//XEST_kn = new double[NNODE];
	for (z = 0; z < nrbsurfnumber; z++) {
		nr[z].XEST_kn = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XEST_kn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XEST = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XEST[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XEST_ukn = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XEST_ukn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].angle_disp1 = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].angle_disp1[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].angle_disp2 = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].angle_disp2[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].disp_mag = new double**[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].disp_mag[i] = new double*[nr[z].NEL_nrb];
			for (j = 0; j < nr[z].NEL_nrb; j++) {
				nr[z].disp_mag[i][j] = new double[2];
			}
		}
		nr[z].XNRBORG2 = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XNRBORG2[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XNRB_kn = new double**[elenode2D];
		nr[z].XNRB_ukn = new double**[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XNRB_kn[i] = new double*[nr[z].NEL_nrb];
			nr[z].XNRB_ukn[i] = new double*[nr[z].NEL_nrb];
			for (j = 0; j < nr[z].NEL_nrb; j++) {
				nr[z].XNRB_kn[i][j] = new double[2];
				nr[z].XNRB_ukn[i][j] = new double[2];
			}
		}
		nr[z].XCOR_ukn = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XCOR_ukn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XCOR = new double*[elenode2D];
		for (i = 0; i < elenode2D; i++) {
			nr[z].XCOR[i] = new double[nr[z].NEL_nrb];
		}
	}

	XNRBORG = new double[NNODE];
	//XCOR = new double[NNODE];
	XCOR_kn = new double[NNODE];
	//XCOR_ukn = new double[NNODE];

	//=======================================================================//
	for (i = 0; i < NNODE; i++) {
		DSDOT[i] = 0.0;
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 3; j++) {
			ds[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			FEEDOT[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			FEEDOT_inc[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			FEE[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			P[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			PT[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		PH[i] = 0.0;
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {   //X DIRECTION AND Y DIRECTION 
			PIN[i][j] = 0.0;
		}
	}
	for (i = 0; i < NNODE; i++) {
		FFORCE[i] = 0.0;
	}

	for (i = 0; i < NNODE; i++) {
		BNRB[i] = 0.0;
	}

	for (i = 0; i < NNODE; i++) {
		HF[i] = 0.0;
	}

	for (i = 0; i < elenode3D; i++) {
		FEETEMP[i] = 0.0;
	}

	for (i = 0; i < NNODE; i++) {
		//DPS_ukn[i] = 0.0;
		XNRBORG[i] = 0.0;
		XCOR_kn[i] = 0.0;
	}

	for (z = 0; z < nrbsurfnumber; z++) {
		for (i = 0; i < elenode2D; i++) {
			for (j = 0; j < nr[z].NEL_nrb; j++) {
				nr[z].XEST[i][j] = 0.0;
				nr[z].XEST_kn[i][j] = 0.0;
				nr[z].XEST_ukn[i][j] = 0.0;
				nr[z].angle_disp1[i][j] = 1.0;
				nr[z].angle_disp2[i][j] = 1.0;
				nr[z].XNRBORG2[i][j] = 0.0;
				for (k = 0; k < 2; k++) {
					nr[z].disp_mag[i][j][k] = 0.0;
				}
			}
		}
	}

	for (z = 0; z < nrbsurfnumber; z++) {
		for (i = 0; i < elenode2D; i++) {
			for (j = 0; j < nr[z].NEL_nrb; j++) {
				for (k = 0; k < 2; k++) {
					nr[z].XNRB_kn[i][j][k] = 0.0;
					nr[z].XNRB_ukn[i][j][k] = 0.0;
				}
				nr[z].XCOR_ukn[i][j] = 0.0;
				nr[z].XCOR[i][j] = 0.0;
			}
		}
	}

	if (tfm == 0) {
		//declarition
		PSI_inc = new double*[NNODE];
		for (i = 0; i < NNODE; i++) {
			PSI_inc[i] = new double[2];
		}
		for (z = 0; z < owsfnumber; z++) { //this memory allocation scheme could have been improved
			ol[z].PSI = new double[ol[z].GIDNct];
			ol[z].DI = new double[ol[z].GIDNct];
			ol[z].DISPI = new double*[ol[z].GIDNct];
			for (j = 0; j < ol[z].GIDNct; j++) {
				ol[z].DISPI[j] = new double[2];
			}
		}
		if (debug == 0) {
			for (z = 0; z < nrbsurfnumber; z++) {
				nr[z].XNRB = new double**[elenode2D];
				for (i = 0; i < elenode2D; i++) {
					nr[z].XNRB[i] = new double*[nr[z].NEL_nrb];
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						nr[z].XNRB[i][j] = new double[2];
					}
				}
				for (i = 0; i < elenode2D; i++) {
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < 2; k++) {
							nr[z].XNRB[i][j][k] = 0.0;
						}
					}
				}
			}
		}
		else {
			for (z = 0; z < nrbsurfnumber; z++) {
				nr[z].XNRB = new double**[elenode2D];
				for (i = 0; i < elenode2D; i++) {
					nr[z].XNRB[i] = new double*[nr[z].NEL_nrb];
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						nr[z].XNRB[i][j] = new double[2];
					}
				}
				for (i = 0; i < elenode2D; i++) {
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < 2; k++) {
							nr[z].XNRB[i][j][k] = 0.0;
						}
					}
				}
			}
		}

		//DPS = new double[NNODE];
		//initialization
		for (i = 0; i < NNODE; i++) {
			for (j = 0; j < 2; j++) {
				PSI_inc[i][j] = 0.0;
			}
		}
		for (z = 0; z < owsfnumber; z++) {
			for (j = 0; j < ol[z].GIDNct; j++) {
				ol[z].PSI[j] = 0.0;
				ol[z].DI[j] = 0.0;
				for (k = 0; k < 2; k++) {
					ol[z].DISPI[j][k] = 0.0;
				}
			}
		}
	}

	double *WBS; //wet surface structure force derived from displacement sent back from Nastran 
	WBS = new double[NNODE];
	for (z = 0; z < owsfnumber; z++) { //this memory allocation scheme could have been improved
		if (mappingalgo == 4 || mappingalgo == 5) {
			ol[z].DISP_gs = new double*[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
			ol[z].DISP_norm = new double*[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
			for (j = 0; j < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); j++) {
				ol[z].DISP_norm[j] = new double[2];
				ol[z].DISP_gs[j] = new double[3];
			}
		}
		else if (mappingalgo == 2) {
			ol[z].DISP = new double*[ol[z].FSNEL*NINT*NINT];
			for (j = 0; j < ol[z].FSNEL*NINT*NINT; j++) {
				ol[z].DISP[j] = new double[3]; //defined in 3 directions 
			}
			ol[z].DISP_norm = new double*[ol[z].FSNEL*elenode2D];
			for (j = 0; j < ol[z].FSNEL*elenode2D; j++) {
				ol[z].DISP_norm[j] = new double[2];
			}
		}
	}

	for (z = 0; z < owsfnumber; z++) {
		for (j = 0; j < NNODE; j++) {
			WBS[j] = 0.0;
		}
		if (mappingalgo == 4 || mappingalgo == 5) {
			for (j = 0; j < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); j++) {
				for (k = 0; k < 2; k++) {
					ol[z].DISP_norm[j][k] = 0.0;
				}
				for (k = 0; k < 3; k++) {
					ol[z].DISP_gs[j][k] = 0.0;
				}
			}
		}
		else if (mappingalgo == 2) {
			for (j = 0; j < ol[z].FSNEL*NINT*NINT; j++) {
				for (k = 0; k < 3; k++) {
					ol[z].DISP[j][k] = 0.0;
				}
			}
			for (j = 0; j < ol[z].FSNEL*elenode2D; j++) {
				for (k = 0; k < 2; k++) {
					ol[z].DISP_norm[j][k] = 0.0;
				}
			}
		}
	}
	double* WP;
	WP = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		WP[i] = 0.0;
	}
	double* WPIN;
	WPIN = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		WPIN[i] = 0.0;
	}

	//Initialize the incident pressure and total pressure along with hydrostatic pressure for the first time step
	TIME = 1;
	for (i = 0; i < NNODE; i++) {
		PH[i] = abs(GCOORD[i][1]) * RHO * grav + PATM; //absolute pressure
	}
	PIN = WAVE_IN(NNODE, GCOORD, T, TIME, PIN, DT, PPEAK, TAU, XC, YC, ZC, XO, YO, ZO);

	if (tfm == 1) {
		for (j = 0; j < NNODE; j++) {
			P[j][0] += PIN[j][0]; //initial condition 
		}
		for (i = 0; i < NNODE; i++) {
			PT[i][0] = P[i][0] + PH[i];   //TOTAL PRESSURE (correct one)
		}
	}
	else {
		for (i = 0; i < NNODE; i++) {
			PT[i][0] = PIN[i][0] + PH[i];   //TOTAL PRESSURE (correct one)
		}
	}

	//P and PIN does not directly participate in the calculation (the necessary initial conditions are integrated analytically)
	//the loop below is only effective when the wave front touches the structural bottom at the first time step (6.17)

	//Initial conditions special for TFM
	//The below initial condition calculation is prepared for exponential decay wave form rather than the Abaqus smoothed wave form
	//initialize FEEDOT (initialized to -(1/2)dt, correspond to staggered integration (a.k.a leap frog integration))
	//currently, the initial condition of Abaqus wave profile is not coded

	if (tfm == 1) {
		double rd = 0.0; double ro = 0.0;
		if (WAVE == 1) {
			for (i = 0; i < NNODE; i++) {
				if (-(2 * xo - 2 * abs(GCOORD[i][1]) + C*DT) / (2 * C) >= 0.0) {
					FEEDOT[i][0] = -PPEAK*TAU*1.0*(exp((DT / 2 + (xo - abs(GCOORD[i][1])) / C) / TAU) - 1);
				}
				else {
					FEEDOT[i][0] = -PPEAK*TAU*0.0*(exp((DT / 2 + (xo - abs(GCOORD[i][1])) / C) / TAU) - 1);
				}
			}
			//prototype: -PPEAK*TAU*heaviside(-(2*xo - 2*x + C*dt)/(2*C))*(exp((dt/2 + (xo - x)/C)/TAU) - 1)

			//initialize FEE (initialized to 0)
			for (i = 0; i < NNODE; i++) {
				if (-(xo - abs(GCOORD[i][1])) / C >= 0.0) {
					FEE[i][0] = -(PPEAK*TAU*exp(-abs(GCOORD[i][1]) / (C*TAU))*1.0*(xo*exp(abs(GCOORD[i][1]) / (C*TAU)) - abs(GCOORD[i][1])*exp(abs(GCOORD[i][1]) / (C*TAU)) - C*TAU*exp(xo / (C*TAU)) + C*TAU*exp(abs(GCOORD[i][1]) / (C*TAU)))) / C;
				}
				else {
					FEE[i][0] = -(PPEAK*TAU*exp(-abs(GCOORD[i][1]) / (C*TAU))*0.0*(xo*exp(abs(GCOORD[i][1]) / (C*TAU)) - abs(GCOORD[i][1])*exp(abs(GCOORD[i][1]) / (C*TAU)) - C*TAU*exp(xo / (C*TAU)) + C*TAU*exp(abs(GCOORD[i][1]) / (C*TAU)))) / C;
				}
			}
			//prototype: -(PPEAK*TAU*exp(-x/(C*TAU))*heaviside(-(XO - x)/C)*(XO*exp(x/(C*TAU)) - x*exp(x/(C*TAU)) - C*TAU*exp(XO/(C*TAU)) + C*TAU*exp(x/(C*TAU))))/C

			//initialize XNBORG (initialize to 0)
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NRBNODE; j++) {
					if (-(2 * xo - 2 * abs(GCOORD[nr[z].NRBA[j] - 1][1]) + C*DT) / (2 * C) >= 0.0) {
						XNRBORG[nr[z].NRBA[j] - 1] = -(PPEAK*TAU*1.0*(exp((DT / 2 + (xo - abs(GCOORD[nr[z].NRBA[j] - 1][1])) / C) / TAU) - 1)) / (C*RHO);
						//XNRBORG[j] = -(PPEAK*TAU*1.0*(exp((DT / 2 + (xo - abs(GCOORD[nr[z].NRBA[j] - 1][1])) / C) / TAU) - 1)) / (C*RHO);
					}
					else {
						XNRBORG[nr[z].NRBA[j] - 1] = -(PPEAK*TAU*0.0*(exp((DT / 2 + (xo - abs(GCOORD[nr[z].NRBA[j] - 1][1])) / C) / TAU) - 1)) / (C*RHO);
						//XNRBORG[j] = -(PPEAK*TAU*0.0*(exp((DT / 2 + (xo - abs(GCOORD[nr[z].NRBA[j] - 1][1])) / C) / TAU) - 1)) / (C*RHO);
					}

				}
			}
			//prototype: -(PPEAK*TAU*heaviside(-(2*XO - 2*x + C*dt)/(2*C))*(exp((dt/2 + (XO - x)/C)/TAU) - 1))/(C*RHO)
		}
		else { //spherical wave
			for (i = 0; i < NNODE; i++) {
				rd = pow((pow((GCOORD[i][0] - XC), 2) + pow((GCOORD[i][1] - YC), 2) + pow((GCOORD[i][2] - ZC), 2)), 0.5);
				ro = pow((pow((XC - XO), 2) + pow((YC - YO), 2) + pow((ZC - ZO), 2)), 0.5);
				if (-(2 * rd - 2 * ro + C*DT) / (2 * C) >= 0) {
					FEEDOT[i][0] = -(PPEAK*TAU*ro*1.0*(exp((DT / 2 + (rd - ro) / C) / TAU) - 1)) / rd;
					FEEDOT_inc[i][0] = FEEDOT[i][0];
				}
				else {
					FEEDOT[i][0] = -(PPEAK*TAU*ro*0.0*(exp((DT / 2 + (rd - ro) / C) / TAU) - 1)) / rd;
					FEEDOT_inc[i][0] = FEEDOT[i][0];
				}
				//-(PPEAK*TAU*ro*heaviside(-(2 * rd - 2 * ro + C*dt) / (2 * C))*(exp((dt / 2 + (rd - ro) / C) / TAU) - 1)) / rd
				if (-(rd - ro) / C >= 0) {
					FEE[i][0] = -(PPEAK*TAU*ro*1.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*rd);
				}
				else {
					FEE[i][0] = -(PPEAK*TAU*ro*0.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*rd);
				}
				//-(PPEAK*TAU*ro*heaviside(-(rd - ro) / C)*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*rd)
			}
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NRBNODE; j++) {
					rd = pow((pow((GCOORD[nr[z].NRBA[j] - 1][0] - XC), 2) + pow((GCOORD[nr[z].NRBA[j] - 1][1] - YC), 2) + pow((GCOORD[nr[z].NRBA[j] - 1][2] - ZC), 2)), 0.5);
					ro = pow((pow((XC - XO), 2) + pow((YC - YO), 2) + pow((ZC - ZO), 2)), 0.5);
					//initialize to 0 
					if (-(rd - ro) / C >= 0) {
						//XNRBORG[j] = -(PPEAK*TAU*ro*1.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*1.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
						XNRBORG[nr[z].NRBA[j] - 1] = -(PPEAK*TAU*ro*1.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*1.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
					}
					else {
						//XNRBORG[j] = -(PPEAK*TAU*ro*0.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*0.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
						XNRBORG[nr[z].NRBA[j] - 1] = -(PPEAK*TAU*ro*0.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*0.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
					}
					// - (PPEAK*TAU*ro*heaviside(-(rd - ro)/C)*(rd - ro + C*TAU - C*TAU*exp(rd/(C*TAU) - ro/(C*TAU))))/(C*RHO*rd^2) - (PPEAK*TAU*ro*heaviside(-(rd - ro)/C)*(exp((rd - ro)/(C*TAU)) - 1))/(C*RHO*rd)
				}
			}
		}
	}

	//=========================Get free surface points=========================//
	int fspt_num = 0;
	double searchrange = 1e-2;
	for (j = 0; j < NNODE; j++) {
		if (fsdebug == 0) {
			if (abs(GCOORD[j][1]) <= searchrange) {
				fspt_num += 1;
			}
		}
		else {
			/*
			if (abs(GCOORD[j][1] + SY) <= 1e-6) {
				fspt_num += 1;
			}
			*/
		}
	}
	int *fspt;
	fspt = new int[fspt_num];
	for (j = 0; j < fspt_num; j++) {
		fspt[j] = 0;
	}
	int count = 0;
	for (j = 0; j < NNODE; j++) {
		if (fsdebug == 0) {
			if (abs(GCOORD[j][1]) <= searchrange) {
				fspt[count] = j + 1;
				count += 1;
			}
		}
		else {
			/*
			if (abs(GCOORD[j][1] + SY) <= 1e-6) {
				fspt[count] = j + 1;
				count += 1;
			}
			*/
		}
	}
	if (count != fspt_num || count == 0) {
		std::cout << "the free surface point recognization is wrong or there's no free surface. Continue?" << std::endl;
		system("PAUSE ");
	}

	//Initiate the first coupling
	initcoupling();

	//Update the current MpCCI time
	current_time = T[0];

	interface_mappingstruct in;

	in = interface_mapping(1, GCOORD, WP, IEN, LNA_3D);
	dotransfer();
	in = interface_mapping(0, GCOORD, WP, IEN, LNA_3D);

	//Generate the information file
	std::string name3 = "parameters_" + timestr + ".txt";
	std::ofstream myfile2;
	myfile2.open(name3);
	//General information:
	myfile2 << "Code name: 3Dbarge_TFM_Abaqus_sym_mesh" << std::endl;
	myfile2 << "Simulation date and time: " << timestr << std::endl;
	myfile2 << "Mesh information:" << std::endl;
	myfile2 << "N: " << N << " mesh size: " << " wetted surface number: " << owsfnumber << " Nodenumber: " << NNODE << " Element number: " << NEL << " h/p refinement rate: " << refine << std::endl;
	myfile2 << "Explosive information: " << std::endl;
	myfile2 << "stdoff point (spherical): " << XO << " " << YO << " " << ZO << " explosion center: " << XC << " " << YC << " " << ZC << " peak pressure (Mpa): " << PPEAK << " decay rate: " << TAU << " WAVE (1: plane, 2: spherical): " << WAVE << std::endl;
	myfile2 << "Time integration information" << std::endl;
	myfile2 << "CFL: " << CFLFRAC << " dt: " << DT << " dt scale factor: " << dtscale << " total time: " << TTERM << " damping: " << BETA << " explicit central difference (Leap frog)" << std::endl;
	myfile2 << "Fluid properties: " << std::endl;
	myfile2 << "CAV: " << CAV << " C: " << C << " RHO: " << RHO << " atmospheric pressure: " << PATM << " saturated pressure: " << PSAT << std::endl;
	myfile2 << "FSI coupling" << std::endl;
	myfile2 << "mapping algorithm: " << mappingalgo << std::endl;
	myfile2 << "debug mode: " << debug << std::endl;

	double angle = 0.0; double R = 0.0;

	double* WPTEMP;
	double* BFTEMP;
	double* DISPTEMP;
	double* WBSTEMP;
	double* NRBDISPTEMP;
	double* BNRBTEMP;
	WPTEMP = new double[elenode2D];
	BFTEMP = new double[elenode2D];
	DISPTEMP = new double[elenode2D];
	WBSTEMP = new double[elenode2D];
	NRBDISPTEMP = new double[elenode2D];
	BNRBTEMP = new double[elenode2D];
	for (i = 0; i < elenode2D; i++) {
		WPTEMP[i] = 0.0;
		BFTEMP[i] = 0.0;
		NRBDISPTEMP[i] = 0.0;
		BNRBTEMP[i] = 0.0;
	}
	for (i = 0; i < elenode2D; i++) {
		DISPTEMP[i] = 0.0;
		WBSTEMP[i] = 0.0;
	}

	double BNRBt = 0.0;
	double dst = 0.0;
	//=================================Start the time integration routine=================================//
	int oc = 0;
	double lastDT = 0.0;

	int NDT_out = 0;
	std::vector <int> T_out;
	std::vector <int> sampline;
	std::vector <double> sampline_gc; //coordinate of the searched nodes
	std::vector <int> sampline_found;
	std::ofstream *outline;
	int* sampline2;
	std::string tecplotfile = "tecplot.dat";
	std::ofstream tecplotfilehd;
	tecplotfilehd.open(tecplotfile);
	//std::string energyfile = "history.txt";
	//std::ofstream energyfilehd;
	//energyfilehd.open(energyfile);
	//Get the sample points on a line to observe the wave propagation pressure distribution
	count = 0;
	for (i = 0; i < NNODE; i++) {
		if (abs(GCOORD[i][0] - 0.0) < 1e-6 && abs(GCOORD[i][2] - 0.0) < 1e-6) {
			sampline.push_back(i + 1);
			sampline_gc.push_back(GCOORD[i][1]); 
			count += 1;
		}
	}
	std::sort(sampline_gc.rbegin(), sampline_gc.rend()); //arrange the coordinate in decending order

	for (i = 0; i < sampline_gc.size(); i++) {
		for (j = 0; j < sampline.size(); j++) {
			if (GCOORD[sampline[j] - 1][1] == sampline_gc[i]) {
				sampline_found.push_back(sampline[j]);
			}
		}
	}
	for (i = 0; i < NDT; i++) {
		if (T[i] >= output_int*NDT_out) {
			T_out.push_back(i);
			NDT_out += 1;
		}
	}
	std::string name2 = "PT_line_result_";
		
	outline = new std::ofstream[NDT_out];
	std::string filename2;
	if (output == 1) {
		for (i = 0; i < NDT_out; i++) {
			filename2 = name2 + std::to_string(T_out[i] * DT * 1000) + "ms " + timestr + ".txt";
			outline[i].open(filename2);
		}
	}

	double hd = 0;
	//for (i = 0; i < NDT - 1; i++) { //error prone: i other than time should not present in this loop

	if (improvednrb == 1) {
		//Initialize the pressure gradient vector
		for (z = 0; z < nrbsurfnumber; z++) {
			nr[z].P_dev = new double**[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].P_dev[i] = new double*[elenode2D];
				for (j = 0; j < elenode2D; j++) {
					nr[z].P_dev[i][j] = new double[3];
				}
			}
		}
		for (z = 0; z < nrbsurfnumber; z++) {
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				for (j = 0; j < elenode2D; j++) {
					for (k = 0; k < 3; k++) {
						nr[z].P_dev[i][j][k] = 0.0;
					}
				}
			}
		}
	}

	if (Bleich == 1) {
		int wavdirc[3] = { 0,1,0 };
	}

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

	for (i = 0; i < NDT; i++) {
		TIME = i + 2;
		std::cout << "time is:   " << TIME << std::endl;
		/*current_time is send to MpCCI server in function dotransfer*/
		current_time = T[TIME - 1]; //start loop from "first time step"
		//current_time = T[TIME - 2];
		/*
		//The last time step size would be adjusted (stretch or shrink) to meet the target termination time
		if (i == NDT - 2) { //the last time step
			DT = TTERM - DT*(NDT - 2);
			time_step_size = DT;
			current_time = TTERM;
		}
		*/

		if (i == NDT - 1) { //the last time step
			time_step_size = TTERM - DT*(NDT - 1);
			current_time = TTERM;
		}

		PIN = WAVE_IN(NNODE, GCOORD, T, TIME, PIN, DT, PPEAK, TAU, XC, YC, ZC, XO, YO, ZO); //USED TO UPDATE PIN IN THIS SUBROUTINE

		for (k = 0; k < NNODE; k++) {
			WBS[k] = 0.0;
		}

		if (tfm == 1) {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					if (nodeforcemap2 == 1) {
						if (debug_hydro == 0) {
							WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM; //correct version with structural gravity
						}
						else {
							//WP[ol[z].GIDN[j] - 1] = PH[ol[z].GIDN[j] - 1] - PATM;
							WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1];
						}
						//ol[z].WP[ol[z].GIDN[j] - 1] = PH[ol[z].GIDN[j] - 1] - PATM; //pure hydrostatic pressure
					//ol[z].WP[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][0]; //incident pressure
					}
					else { //absolute pressure
						WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM;
						//ol[z].WP[ol[z].GIDN[j] - 1] = PH[ol[z].GIDN[j] - 1] - PATM; //pure hydrostatic pressure
						//ol[z].WP[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][0]; //incident pressure
					}
					//ol[z].WP[ol[z].GIDN[j] - 1]
					//= PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]; //correct version if structural gravity is not specified in Abaqus
					if (Bleich == 1) {
						if (nodeforcemap2 == 1) {
							WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
						}
						else {
							WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
						}
					}
				}
			}
		}
		else {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM; // //correct version with structural gravity
					if (Bleich == 1) {
						WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
					}
					WPIN[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][1] + PIN[ol[z].GIDN[j] - 1][0];
				}
			}
		}

		//start = std::clock();
		//=======================define double* nodeforce in fluid code==========================//
		//mapping the fluid force ABF from user defined mesh to MpCCI defined mesh on coupling surface using interpolation
		in = interface_mapping(1, GCOORD, WP, IEN, LNA_3D);
		//after this subroutine, the nodeforce should already be mapped onto coupling surface (data.h)
		//int fluid2structure, int**IEN_3D, int***LNA_3D, int**LNA_2D, int NNODE, double *Z
		dotransfer();
		//=============map nodal displacement from coupled surface to fluid mesh===================//
		in = interface_mapping(0, GCOORD, WP, IEN, LNA_3D);
		if (tfm == 0) { //Scattered field model
			//double angle = 0.0; //cos value
			double r = 0.0;
			double angle = 0.0;
			//===============Incident fluid displacement on FSI boundary (used to derive structure force)================//
			//Change the definition of DISPI to node based rather than element based. 
			for (z = 0; z < owsfnumber; z++) {
				angle = 0.0;
				for (j = 0; j < ol[z].FSNEL; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
					for (k = 0; k < elenode2D; k++) {
						if (WAVE == 1) { //plane wave
							angle = ol[z].norm[j][0] * wavdirc[0] + ol[z].norm[j][1] * wavdirc[1] + ol[z].norm[j][2] * wavdirc[2];
							ol[z].PSI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1]
								= (DT / 2.0)*(WPIN[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1]); //INCIDENT DISP. PREDICTOR by time integration	
							ol[z].DI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1]
								= ol[z].PSI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1] / (RHO*C);     //trapezoidal integration of the double integrator
							ol[z].DISPI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][1]
								= ol[z].DISPI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][0] + angle*ol[z].DI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1];
						}
						else { //Spherical wave 
							r = sqrt(pow((GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][0] - XC), 2) + pow((GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][1] - YC), 2) + pow((GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][2] - ZC), 2));
							ol[z].PSI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1]
								= (DT / 2.0)*(WPIN[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1]); //INCIDENT DISP. PREDICTOR by time integration	
							PSI_inc[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][1] = ol[z].PSI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1];
							angle = (ol[z].norm[j][0] * (GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][0] - XC) / r + ol[z].norm[j][1] * (GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][1] - YC) / r + ol[z].norm[j][2] * (GCOORD[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][2] - ZC) / r);
							ol[z].DI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1]
								= ol[z].PSI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1] / (RHO*C) + (PSI_inc[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][0] + PSI_inc[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1][1]) * (DT / 2) / (RHO * r);  //trapezoidal integration of the double integrator
							ol[z].DISPI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][1]
								= ol[z].DISPI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][0] + angle*ol[z].DI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1]; //UPDATED INCIDENT STRUCUTRE DISPLACEMENT 
						}
					}
				}
			}
			//SET VALUES AT T(i+1) TO VALUES AT T(i)
			/*
			for (j = 0; j < NNODE; j++) {
				for (z = 0; z < owsfnumber; z++) {
					ol[z].DISPI[j][0] = ol[z].DISPI[j][1];  //SOLUTION ARRAY FOR INCIDENT STRUCUTRE DISPLACEMENT
				}
				PSI_inc[j][0] = PSI_inc[j][1];
			}
			*/
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					ol[z].DISPI[j][0] = ol[z].DISPI[j][1];
					for (k = 0; k < elenode2D; k++) {
						PSI_inc[ol[z].IEN_gb[k][j] - 1][0] = PSI_inc[ol[z].IEN_gb[k][j] - 1][1];
					}
				}
			}
		}

		if (tfm == 1) { //Total field model
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].FSNEL; j++) { //for every fluid element linked to structure
					for (h = 0; h < elenode2D; h++) {
						WBSTEMP[h] = 0.0;
						for (k = 0; k < elenode2D_gs; k++) { //For mappingalgo 5, k stands for the quadrature nodes; For mappingalgo2, k stands for the interpolation nodes
							if (mappingalgo == 5 || mappingalgo == 4) {
								ol[z].DISP_norm[j*elenode2D_gs + k][1] = ol[z].norm[j][0] * ol[z].DISP_gs[j*elenode2D_gs + k][0] + ol[z].norm[j][1] * ol[z].DISP_gs[j*elenode2D_gs + k][1] + ol[z].norm[j][2] * ol[z].DISP_gs[j*elenode2D_gs + k][2];
							}
							else {
								ol[z].DISP_norm[j*elenode2D_gs + k][1] = ol[z].norm[j][0] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][0] + ol[z].norm[j][1] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][1] + ol[z].norm[j][2] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][2];
							}
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * (-1) * RHO * (ol[z].DISP_norm[j*elenode2D_gs + k][1] + (ol[z].DISP_norm[j*elenode2D_gs + k][1] - ol[z].DISP_norm[j*elenode2D_gs + k][0]));
						}
					}
					for (k = 0; k < elenode2D; k++) {
						WBS[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1] += WBSTEMP[k];
					}
				}
			}
		}
		else {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].FSNEL; j++) { //for every fluid element linked to structure
					for (h = 0; h < elenode2D; h++) {
						WBSTEMP[h] = 0.0;
						for (k = 0; k < elenode2D_gs; k++) {
							if (mappingalgo == 5 || mappingalgo == 4) {
								ol[z].DISP_norm[j*elenode2D_gs + k][1] = ol[z].norm[j][0] * ol[z].DISP_gs[j*elenode2D_gs + k][0] + ol[z].norm[j][1] * ol[z].DISP_gs[j*elenode2D_gs + k][1] + ol[z].norm[j][2] * ol[z].DISP_gs[j*elenode2D_gs + k][2];
							}
							else {
								//ol[z].DISP_norm[j*elenode2D_gs + k][1] = ol[z].norm[j][0] * ol[z].DISP[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][0] + ol[z].norm[j][1] * ol[z].DISP[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][1] + ol[z].norm[j][2] * ol[z].DISP[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][2];
								ol[z].DISP_norm[j*elenode2D_gs + k][1] = ol[z].norm[j][0] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][0] + ol[z].norm[j][1] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][1] + ol[z].norm[j][2] * ol[z].DISP[j*elenode2D_gs + ol[z].FP_2D[k] - 1][2];
							}
							std::cout << "We need to derive the DISPI on gauss point!!! (not done yet)" << std::endl;
							system("PAUSE ");
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * (-1) * RHO * (ol[z].DISP_norm[j*elenode2D_gs + k][1] + (ol[z].DISP_norm[j*elenode2D_gs + k][1] - ol[z].DISP_norm[j*elenode2D_gs + k][0]) - ol[z].DISPI[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1][1]);
						}
					}
					for (k = 0; k < elenode2D; k++) {
						//ol[z].WBS[ol[z].IEN_lc[ol[z].FP_2D[k] - 1][j] - 1] += WBSTEMP[k];
						WBS[ol[z].IEN_gb[ol[z].FP_2D[k] - 1][j] - 1] += WBSTEMP[k];
					}
				}
			}
		}

		for (j = 0; j < NNODE; j++) {
			BNRB[j] = 0.0;
		}

		//NRB PREDICTOR
		//NRBA was originally designed to store all the NRB nodes on all NRB surfaces in FSP code
		//We want to make NRBA local within each NRB surface in the current code 
		//IEN_gb is a 2D connectivity matrix
		//Derive the angle between pressure gradient (the same direction with displacement/velocity gradient) and the normal direction of the NRB surface elements
		if (improvednrb == 1) {
			for (z = 0; z < nrbsurfnumber; z++) {
				for (ii = 0; ii < nr[z].NEL_nrb; ii++) {
					for (jj = 0; jj < elenode2D; jj++) {
						for (kk = 0; kk < 3; kk++) {
							nr[z].P_dev[ii][jj][kk] = 0.0;
						}
					}
				}
			}
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						//loop through each point in element and accumulate the value
						for (ii = 0; ii < NINT; ii++) {
							for (jj = 0; jj < NINT; jj++) {
								for (kk = 0; kk < NINT; kk++) {
									for (h = 0; h < 3; h++) { //the three directions (x y z )
										if (debug == 0) {
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[j][k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
										}
										else {
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * PH[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
											//although we've devided the RHO, the value is still quite large.
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * P[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * (FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] - FEE_inc[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0]) / RHO;
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[j][k] - 1] * FEEDOT[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		//NRB predictor
		if (tfm == 1) { //TFM
			if (WAVE == 1) { //Plane wave
				for (z = 0; z < nrbsurfnumber; z++) {
					angle = -1.0;
					for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
						for (k = 0; k < elenode2D; k++) {
							if (improvednrb == 1) {
								//modify the value base on the angle between normal displacement direction and normal surface direction
								if (pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5) != 0) {
									if (debug == 0) {
										nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									}
									else {
										//nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
										nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (nr[z].norm[j][0]) + (nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (nr[z].norm[j][1]) + (nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									}
									//std::cout << " " << std::endl;
								}
								else {
									nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = 1.0;
								}
							}
							if (debug == 0) {
								nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*(0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1])
									- (0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
							}
							else {
								//nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*(0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) 
									//- (angle) * (0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
								nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*(0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
							}
							if (debug == 0) {
								nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0];
							}
							else {
								//nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] 
									//- nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
								//nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*(P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] - PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0]);
								//std::cout << " " << std::endl;
								nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0]; //plane wave approximation (PWA)
							}
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
						for (k = 0; k < elenode2D; k++) {
							nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j];
							nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j];
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < elenode2D; k++) {
							//if (debug == 0) {
							NRBDISPTEMP[k] = -RHO* (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							//}
							//else {
								//NRBDISPTEMP[k] = -RHO* (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1]);
							//}
							//std::cout << " " << std::endl;
						}
						for (h = 0; h < elenode2D; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < elenode2D; k++) {
								BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * NRBDISPTEMP[k];
							}
						}
						for (k = 0; k < elenode2D; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
			else { //Spherical wave
				for (z = 0; z < nrbsurfnumber; z++) { //from bottom to front surface (k=0 means all node case) 
					for (j = 0; j < nr[z].NRBNODE; j++) {
						FEEDOT_inc[nr[z].NRBA[j] - 1][1] = FEEDOT_inc[nr[z].NRBA[j] - 1][0] + (0.5 * DT)*(PIN[nr[z].NRBA[j] - 1][0] + PIN[nr[z].NRBA[j] - 1][1]);
						//prototype: FEEDOT_inc[nr[z].NRBA[j][0] - 1][1] = (PPEAK*TAU - PPEAK*TAU*exp(pow((pow(XC, 2) - 2 * XC*d.GCOORD[nr[z].NRBA[j][0] - 1][0] + pow(YC, 2) - 2 * YC*d.GCOORD[nr[z].NRBA[j][0] - 1][1] + pow(ZC, 2) - 2 * ZC*d.GCOORD[nr[z].NRBA[j][0] - 1][2] + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][0], 2) + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][1], 2) + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][2], 2)), 0.5) / (C*TAU))*exp(-(DT*i) / TAU)*exp(-pow((pow(XC, 2) - 2 * XC*XO + pow(XO, 2) + pow(YC, 2) - 2 * YC*YO + pow(YO, 2) + pow(ZC, 2) - 2 * ZC*ZO + pow(ZO, 2)), 0.5) / (C*TAU))*exp(-DT / TAU))*(sign((C*DT - pow((pow(XC, 2) - 2 * XC*d.GCOORD[nr[z].NRBA[j][0] - 1][0] + pow(YC, 2) - 2 * YC*d.GCOORD[nr[z].NRBA[j][0] - 1][1] + pow(ZC, 2) - 2 * ZC*d.GCOORD[nr[z].NRBA[j][0] - 1][2] + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][0], 2) + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][1], 2) + pow(d.GCOORD[nr[z].NRBA[j][0] - 1][2], 2)), 0.5) + pow((pow(XC, 2) - 2 * XC*XO + pow(XO, 2) + pow(YC, 2) - 2 * YC*YO + pow(YO, 2) + pow(ZC, 2) - 2 * ZC*ZO + pow(ZO, 2)), 0.5) + C*DT*i) / C) / 2.0 + 0.5);
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
						for (k = 0; k < elenode2D; k++) {
							if (improvednrb == 1) {
								//modify the value base on the angle between normal displacement direction and normal surface direction
								if (pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5) != 0) {
									if (debug == 0) {
										nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									}
									else {
										nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									}
									//std::cout << " " << std::endl;
								}
								else {
									nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = 1.0;
								}
							}
							R = pow(pow(GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] - XC, 2) + pow(GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1] - YC, 2) + pow(GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][2] - ZC, 2), 0.5);
							angle = (nr[z].norm[j][0] * (GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] - XC) + nr[z].norm[j][1] * (GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1] - YC) + nr[z].norm[j][2] * (GCOORD[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][2] - ZC)) / R;
							if (debug == 0) {
								nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*((0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) +
									(0.5 * DT / (RHO*R))*(FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1])) -
									(0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
							}
							else {
								//nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*((0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) +
									//(0.5 * DT / (RHO*R))*(FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]))
									//- (angle) * (0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
								nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j] = angle*((0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) +
									(0.5 * DT / (RHO*R))*(FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + FEEDOT_inc[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]));
							}
							if (debug == 0) {
								nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0]; //plane wave approximation (PWA)
							}
							else {
								//nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0]
									//- nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (0.5 * DT / (RHO*C))*(PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]);
								nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0];
								//nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								//nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0]; //plane wave approximation (PWA)
							}
							nr[z].XNRBORG2[nr[z].DP_2D[k] - 1][j] = angle*XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1];
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
						for (k = 0; k < elenode2D; k++) {
							nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j];
							nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j];
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) { //IEN_gb needs to be changed (NINT*N)
						for (h = 0; h < elenode2D; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < elenode2D; k++) {
								if (debug_algo5 == 0) {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (-RHO) * (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRBORG2[nr[z].DP_2D[k] - 1][j]);
								}
								else {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (0);
								}
							}
						}
						for (k = 0; k < elenode2D; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
		}
		else { //SFM
			for (z = 0; z < nrbsurfnumber; z++) { //from bottom to front surface (k=0 means all node case)
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						if (improvednrb == 1) {
							//modify the value base on the angle between normal displacement direction and normal surface direction
							if (pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5) != 0) {
								if (debug == 0) {
									nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
								else {
									nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									//nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
							}
							else {
								nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] = 1.0;
							}
						}
						if (debug == 0) {
							nr[z].XEST[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0];
						}
						else {
							//nr[z].XEST[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*abs(P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0]);
							nr[z].XEST[nr[z].DP_2D[k] - 1][j] = nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * (DT / (RHO*C))*P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0];
						}
					}
				}
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						if (debug == 0) {
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST[nr[z].DP_2D[k] - 1][j]; //0 IS BEFORE MODIFICATION, 1 IS AFTER MODIFICATION
						}
						else {
							//nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST[nr[z].DP_2D[k] - 1][j]; //0 IS BEFORE MODIFICATION, 1 IS AFTER MODIFICATION
							nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] + DT * nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0];
							//nr[z].XNRB[nr[z].DP_2D[k] - 1][j][2] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] + (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0]); //0 IS BEFORE MODIFICATION, 1 IS AFTER MODIFICATION
						}
					}
				}
			}
			if (WAVE == 1) { //Plane wave
				for (z = 0; z < nrbsurfnumber; z++) {
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < elenode2D; k++) {
							if (debug == 0) {
								NRBDISPTEMP[k] = -RHO* (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							}
							else {
								NRBDISPTEMP[k] = -RHO* (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							}
						}
						for (h = 0; h < elenode2D; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < elenode2D; k++) {
								BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * NRBDISPTEMP[k];
							}
						}
						for (k = 0; k < elenode2D; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
			else { //Spherical wave
				for (z = 0; z < nrbsurfnumber; z++) { //from bottom to front surface (k=0 means all node case)
					//XNRB is the predicted displacement normal to the structure (since P is normal to the surface)
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (h = 0; h < elenode2D; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < elenode2D; k++) {
								if (debug_algo5 == 0) {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (-RHO) * (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
								}
								else {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (1);
								}
							}
						}
						for (k = 0; k < elenode2D; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
		}

		//time integration
		//start = std::clock();
		#pragma omp parallel for num_threads(6)
		for (int j = 0; j < NNODE; j++) {
			DSDOT[j] = (ds[j][1] - ds[j][0]) / DT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF CONDENSATION
			FEEDOT[j][1] = FEEDOT[j][0] + DT*(P[j][0] + (BETA*DT*(pow(C, 2))*DSDOT[j]));
			FEE[j][1] = FEE[j][0] + DT*FEEDOT[j][1];
			FEEDOT[j][0] = FEEDOT[j][1];  //USE TWO VALUE TO STORE THE N TIME STEP AND N+1 TIME STEP
			HF[j] = 0.0;
		} //1 MEANS N+1 TIME STEP; 0 MEANS N TIME STEP
		//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
		//std::cout << "total CPU time (ms): " << duration << std::endl;
		//std::cout << " " << std::endl;

		//start = std::clock();
		//ctt3 = 0;
		if (tensorfactorization == 0) {
			//start = std::clock();
			#pragma omp parallel for num_threads(6)
			for (int j = 0; j < NEL; j++) { //the loop takes 1+NEL+1+NEL=2*NEL+2 operations
				//-------------------------matrix multiplication-------------------------//
				for (int z = 0; z < elenode3D; z++) { //takes NEL*(2*NINT^3+2) operations
					HFTEMP[j][z] = 0.0; //takes NEL*NINT^3 operations
					for (int k = 0; k < elenode3D; k++) { //takes NEL*NINT^3*(2*NINT^3+2) operaitons
						HFTEMP[j][z] = HFTEMP[j][z] + HMASTER[j][z][k] * FEE[IEN[k][j] - 1][1];
					} //The premise is that HMASTER is the same for every element (only true for equally shaped element!!)
				}
			}
			for (int j = 0; j < NEL; j++) {
				for (int k = 0; k < elenode3D; k++) { //takes NEL*(2*NINT^3+2) operations
					HF[IEN[k][j] - 1] = HF[IEN[k][j] - 1] + HFTEMP[j][k];   //global level reactance matrix //takes 2*NINT^3*NEL operations
				} //HF is like global stiffness matrix (include variable in it)
			}
			//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
			//std::cout << "total CPU time (ms): " << duration << std::endl;
			//std::cout << " " << std::endl;
		}
		else {
			//#pragma omp parallel for num_threads(6) /*private(gamman, gamma_tn)*/
			for (int j = 0; j < NEL; j++) { //takes 2*NEL+2 operations 
				double gamman[3 * NINT*NINT*NINT];
				double gamma_tn[3 * NINT*NINT*NINT];
				for (int ii = 0; ii < NINT; ii++) { //p NEL*(2*NINT+2) operations
					for (int h = 0; h < NINT; h++) { //j,i,i //takes NEL*NINT*(2*NINT+2) operations
						for (int z = 0; z < NINT; z++) { //k,k,j //takes NEL*NINT^2*(2*NINT+2) operations
							gamman[3 * counter1[ii*NINT*NINT + h*NINT + z] + 0] = 0.0;
							gamman[3 * counter2[ii*NINT*NINT + h*NINT + z] + 1] = 0.0;
							gamman[3 * counter3[ii*NINT*NINT + h*NINT + z] + 2] = 0.0;
							for (int q = 0; q < NINT; q++) { //q (recurrent addition in this dimension) //takes NEL*NINT^3*(2*NINT+2) operations	
								gamman[3 * counter1[ii*NINT*NINT + h*NINT + z] + 0] += SHOD1[q*NINT + ii] * FEE[IENct1[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								gamman[3 * counter2[ii*NINT*NINT + h*NINT + z] + 1] += SHOD1[q*NINT + ii] * FEE[IENct2[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								gamman[3 * counter3[ii*NINT*NINT + h*NINT + z] + 2] += SHOD1[q*NINT + ii] * FEE[IENct3[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
							}
							//ctt1 += 1;
						}
					}
				}
				//takes 15*NINT^3*NEL FLOP
				for (int ii = 0; ii < NINT*NINT*NINT; ii++) { //takes NEL*(2 * NINT + 2) operations
					gamma_tn[3 * ii + 0] = Gn[j][0][0][ii] * gamman[3 * ii + 0] + Gn[j][0][1][ii] * gamman[3 * ii + 1] + Gn[j][0][2][ii] * gamman[3 * ii + 2];
					//ipk
					gamma_tn[3 * ii + 1] = Gn[j][1][1][ii] * gamman[3 * ii + 1] + Gn[j][0][1][ii] * gamman[3 * ii + 0] + Gn[j][1][2][ii] * gamman[3 * ii + 2];
					//ijp
					gamma_tn[3 * ii + 2] = Gn[j][2][2][ii] * gamman[3 * ii + 2] + Gn[j][0][2][ii] * gamman[3 * ii + 0] + Gn[j][1][2][ii] * gamman[3 * ii + 1];
					//oc += 15;
				}

				//takes 6*NINT^4*NEL FLOP
				//ctt2 = 0;
				for (int h = 0; h < NINT; h++) {  //i //takes NEL*NINT*(2*NINT+2) operations
					for (int k = 0; k < NINT; k++) { //j  //takes NEL*NINT^2*(2*NINT+2) operations
						for (int z = 0; z < NINT; z++) { //k //takes NEL*NINT^3*(2*NINT+2) operations
							HFTEMPn[j][LNAct3[h*NINT*NINT + k*NINT + z]] = 0.0;
							for (int ii = 0; ii < NINT; ii++) { //p (recurrent addition) //p NEL*(2*NINT+2) operations
								HFTEMPn[j][LNAct3[h*NINT*NINT + k*NINT + z]] += SHOD1[h*NINT + ii] * gamma_tn[3 * (ii*NINT*NINT + k*NINT + z) + 0] + SHOD1[k*NINT + ii] * gamma_tn[3 * (h*NINT*NINT + ii*NINT + z) + 1] + SHOD1[z*NINT + ii] * gamma_tn[3 * (h*NINT*NINT + k*NINT + ii) + 2];
							}
							//ctt2 += 1;
						}
					}
				}
			}
			for (j = 0; j < NEL; j++) {
				//takes 1*NINT^3*NEL FLOP
				for (h = 0; h < NINT*NINT*NINT; h++) { //takes NEL*(2 * NINT + 2) operations
					HF[IENct3[j*NINT*NINT*NINT + h]] += HFTEMPn[j][LNAct3[h]]; //takes NEL*2*NINT^3 operations
					//ctt3 += 1;
				}
			}
			/*
			for (ii = 0; ii < NINT; ii++) { //p NEL*(2*NINT+2) operations
					for (h = 0; h < NINT; h++) { //j,i,i //takes NEL*NINT*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k,k,j //takes NEL*NINT^2*(2*NINT+2) operations
							gamma[ii][h][z][0] = 0.0;
							gamma[h][ii][z][1] = 0.0;
							gamma[h][z][ii][2] = 0.0;
							for (q = 0; q < NINT; q++) { //q (recurrent addition in this dimension) //takes NEL*NINT^3*(2*NINT+2) operation
								gamma[ii][h][z][0] += SHOD[1][q][ii] * FEE[IEN[LNA_3D[q][h][z] - 1][j] - 1][1]; //takes 9*NEL*NINT^4 operations
																												//ipk
								gamma[h][ii][z][1] += SHOD[1][q][ii] * FEE[IEN[LNA_3D[h][q][z] - 1][j] - 1][1];
								//ijp
								gamma[h][z][ii][2] += SHOD[1][q][ii] * FEE[IEN[LNA_3D[h][z][q] - 1][j] - 1][1];
								//oc += 6;
							}
						}
					}
				}
				for (ii = 0; ii < NINT; ii++) { //takes NEL*(2 * NINT + 2) operations
					for (h = 0; h < NINT; h++) { //takes NEL*NINT*(2*NINT + 2) operations
						for (z = 0; z < NINT; z++) { //takes NEL*NINT^2*(2*NINT+2) operations
							gamma_t[ii][h][z][0] = G[0][0][ii][h][z] * gamma[ii][h][z][0] + G[0][1][ii][h][z] * gamma[ii][h][z][1] + G[0][2][ii][h][z] * gamma[ii][h][z][2]; //takes 18*NINT^3*NEL operations
																																											 //ipk
							gamma_t[ii][h][z][1] = G[1][1][ii][h][z] * gamma[ii][h][z][1] + G[0][1][ii][h][z] * gamma[ii][h][z][0] + G[1][2][ii][h][z] * gamma[ii][h][z][2];
							//ijp
							gamma_t[ii][h][z][2] = G[2][2][ii][h][z] * gamma[ii][h][z][2] + G[0][2][ii][h][z] * gamma[ii][h][z][0] + G[1][2][ii][h][z] * gamma[ii][h][z][1];
							oc += 15;
						}
					}
				}
				for (h = 0; h < NINT; h++) {  //i //takes NEL*NINT*(2*NINT+2) operations
					for (k = 0; k < NINT; k++) { //j  //takes NEL*NINT^2*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k //takes NEL*NINT^3*(2*NINT+2) operations
							HFTEMP[j][LNA_3D[h][k][z] - 1] = 0.0;
							for (ii = 0; ii < NINT; ii++) { //p (recurrent addition) //p NEL*(2*NINT+2) operations
								HFTEMP[j][LNA_3D[h][k][z] - 1] += SHOD[1][h][ii] * gamma_t[ii][k][z][0] + SHOD[1][k][ii] * gamma_t[h][ii][z][1] + SHOD[1][z][ii] * gamma_t[h][k][ii][2]; //takes 7*NINT^4*NEL
							}
						}
					}
				}
				for (h = 0; h < NINT; h++) { //takes NEL*(2 * NINT + 2) operations
					for (q = 0; q < NINT; q++) { //takes NEL*NINT*(2*NINT + 2) operations
						for (z = 0; z < NINT; z++) { //takes NEL*NINT^2*(2*NINT+2) operations
							HF[IEN[LNA_3D[h][q][z] - 1][j] - 1] += HFTEMP[j][LNA_3D[h][q][z] - 1]; //takes NEL*2*NINT^3 operations
						}
					}
				}
			*/
		}
		//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
		//std::cout << "total CPU time (ms): " << duration << std::endl;
		//std::cout << " " << std::endl;


		for (j = 0; j < NNODE; j++) {
			FFORCE[j] = -HF[j];
		}

		for (j = 0; j < NNODE; j++) { //for every fluid element linked to structure
			FFORCE[j] += WBS[j];
		}

		for (j = 0; j < nrb.NNODE_nrb; j++) {
			FFORCE[nrb.NRBA_t[j] - 1] += BNRB[nrb.NRBA_t[j] - 1];
		}
		//The combination of FFORCE passes the test

		hd = 0.0;
		for (j = 0; j < NNODE; j++) {
			hd += FFORCE[j];
		}

		if (tfm == 1) {
			for (j = 0; j < fspt_num; j++) {
				Q[fspt[j] - 1] = 1.0;
				FFORCE[fspt[j] - 1] = 0.0 / pow(C, 2.0);
			}
		}
		else {
			for (j = 0; j < fspt_num; j++) {
				Q[fspt[j] - 1] = 1.0;
				FFORCE[fspt[j] - 1] = -PIN[fspt[j] - 1][1] / pow(C, 2.0); //PIN[][1] error prone
			}
		}

		for (j = 0; j < NNODE; j++) {
			ds[j][2] = FFORCE[j] / Q[j];
		}
		//The combination of ds[j][2] passes the test
		//The combination of nr[z].ADMASTERG passes the tests

		//PRESSURE CORRECTION ON NRB NODES
		for (j = 0; j < nrb.NNODE_nrb; j++) {
			KAPPA = (DT*C*nrb.ADMASTERG[nrb.NRBA_t[j] - 1]) / (2 * Q[nrb.NRBA_t[j] - 1]);
			//PRESSURE CORRECTION FACTOR
			ds[nrb.NRBA_t[j] - 1][2] = ds[nrb.NRBA_t[j] - 1][1] + ((ds[nrb.NRBA_t[j] - 1][2] - ds[nrb.NRBA_t[j] - 1][1]) / (1 + KAPPA));
		} //after ds is updated, P can be updated.

		//======================NODE-BY-NODE CAVITATION CHECK=========================//
		if (tfm == 1) {
			if (CAV == 1) {
				for (j = 0; j < NNODE; j++) {
					if ((pow(C, 2.0))*ds[j][2] > PSAT - PH[j]) {
						P[j][1] = pow(C, 2.0)*ds[j][2];
					}
					else {
						P[j][1] = PSAT - PH[j];
					}
				}
			}
			else if (CAV == 0) {
				for (j = 0; j < NNODE; j++) {
					P[j][1] = pow(C, 2.0)*ds[j][2];
				}
			}
		}
		else {
			if (CAV == 1) {
				for (j = 0; j < NNODE; j++) {
					if ((pow(C, 2.0))*ds[j][2] > PSAT - (PH[j] + PIN[j][1])) {
						P[j][1] = pow(C, 2.0)*ds[j][2];
					}
					else {
						P[j][1] = PSAT - (PH[j] + PIN[j][1]);
					}
				}
			}
			else if (CAV == 0) {
				for (j = 0; j < NNODE; j++) {
					P[j][1] = pow(C, 2.0)*ds[j][2];
				}
			}
		}

		if (improvednrb == 1) {
			for (z = 0; z < nrbsurfnumber; z++) {
				for (ii = 0; ii < nr[z].NEL_nrb; ii++) {
					for (jj = 0; jj < elenode2D; jj++) {
						for (kk = 0; kk < 3; kk++) {
							nr[z].P_dev[ii][jj][kk] = 0.0;
						}
					}
				}
			}
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						//loop through each point in element and accumulate the value
						for (ii = 0; ii < NINT; ii++) {
							for (jj = 0; jj < NINT; jj++) {
								for (kk = 0; kk < NINT; kk++) {
									for (h = 0; h < 3; h++) { //the three directions (x y z )
										if (debug == 0) {
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[j][k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
										}
										else {
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * P[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * (FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] - FEE_inc[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1]) / RHO;
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[j][k] - 1] * FEEDOT[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		if (tfm == 1) {
			for (z = 0; z < nrbsurfnumber; z++) {
				/*
				//XNRB corrector
				for (j = 0; j < nr[z].NRBNODE; j++) {
					DPS_ukn[nr[z].NRBA[j] - 1] = 0.5*DT*(P[nr[z].NRBA[j] - 1][1] + P[nr[z].NRBA[j] - 1][0]);
					XCOR_ukn[nr[z].NRBA[j] - 1] = DPS_ukn[nr[z].NRBA[j] - 1] / (RHO*C);
				}
				*/
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						if (improvednrb == 1) {
							//modify the value base on the angle between normal displacement direction and normal surface direction
							if (pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5) != 0) {
								if (debug == 0) {
									nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
								else {
									nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
								//std::cout << " " << std::endl;
							}
							else {
								nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = 1.0;
							}
						}
						if (debug == 0) {
							nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C);
						}
						else {
							//nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C)
								//- 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * PIN[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C);
							//nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C);
							nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
							//nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j] = 0.5 * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1]);
							nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1]);
						}
					}
				}
				for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
					for (k = 0; k < elenode2D; k++) {
						//nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + XCOR_ukn[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1];
						nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j];
					}
				}
				for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
					for (k = 0; k < elenode2D; k++) {
						nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][0] = nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1];
						nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1];
					}
				}
			}
		}
		else { //SFM
			for (z = 0; z < nrbsurfnumber; z++) {
				/*
				for (j = 0; j < nr[z].NRBNODE; j++) {
					DPS[nr[z].NRBA[j] - 1] = 0.5*DT*(P[nr[z].NRBA[j] - 1][1] + P[nr[z].NRBA[j] - 1][0]); // PRESSURE CORRECTOR FOR NRB
				}
				*/
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						if (improvednrb == 1) {
							//modify the value base on the angle between normal displacement direction and normal surface direction
							if (pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5) != 0) {
								if (debug == 0) {
									nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
								else {
									nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
									//nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = ((-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0]) * (-nr[z].norm[j][0]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1]) * (-nr[z].norm[j][1]) + (-nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2]) * (-nr[z].norm[j][2])) / pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
								}
							}
							else {
								nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] = 1.0;
							}
						}
						if (debug == 0) {
							nr[z].XCOR[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C);
						}
						else {
							//nr[z].XCOR[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * abs(P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0]) + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * abs(P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1])) / (RHO*C);
							//nr[z].XCOR[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * P[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1][1]) / (RHO*C);
							nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
							nr[z].XCOR[nr[z].DP_2D[k] - 1][j] = 0.5 * DT * (nr[z].angle_disp1[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][0] + nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1]);
						}
					}
				}
				/*
				for (j = 0; j < nr[z].NRBNODE; j++) {
					XCOR[nr[z].NRBA[j] - 1] = DPS[nr[z].NRBA[j] - 1] / (RHO*C);  //DISPLACEMENT CORRECTOR
				}
				*/
				for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
					for (k = 0; k < elenode2D; k++) {
						if (debug == 0) {
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] + nr[z].XCOR[nr[z].DP_2D[k] - 1][j];
						}
						else {
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] + nr[z].XCOR[nr[z].DP_2D[k] - 1][j];
							//nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1] = pow(pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][0], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][1], 2) + pow(nr[z].P_dev[j][nr[z].DP_2D[k] - 1][2], 2), 0.5);
							//nr[z].XNRB[nr[z].DP_2D[k] - 1][j][2] = nr[z].angle_disp2[nr[z].DP_2D[k] - 1][j] * nr[z].disp_mag[nr[z].DP_2D[k] - 1][j][1];
						}
					}
				}
			}
		}

		//UPDATED TOTAL PRESSURE 
		if (tfm == 1) {
			for (j = 0; j < NNODE; j++) {
				PT[j][1] = PH[j] + P[j][1];
			}
		}
		else {
			for (j = 0; j < NNODE; j++) {
				PT[j][1] = PH[j] + PIN[j][1] + P[j][1];
			}
		}
		
		if (output == 1) {
			for (k = 0; k < NDT_out; k++) {
				if (i == T_out[k]) {
					for (j = 0; j < sampline_found.size(); j++) {
						outline[k] << GCOORD[sampline_found[j] - 1][1] << " " << PT[sampline_found[j] - 1][0] << " " << PT[sampline_found[j] - 1][1] << std::endl;
					}
				}
			}
			if (i == 0) {
				//output the tecplot data file
				for (j = 0; j < NNODE; j++) {
					tecplotfilehd << GCOORD[j][0] << " " << GCOORD[j][1] << " " << GCOORD[j][2] << " " << PIN[j][0] << " " << PIN[j][1] << " " << PT[j][1] << " " << ds[j][2] << " " << BNRB[j] << " " << FEE[j][1] << " " << HF[j] << " " << FFORCE[j] << std::endl;
				}
				tecplotfilehd << std::endl;
				for (j = 0; j < NEL; j++) {
					for (k = 0; k < elenode3D; k++) {
						tecplotfilehd << IEN[k][j] << " ";
					}
					tecplotfilehd << std::endl;
				}
			}
		}
	
		for (j = 0; j < NNODE; j++) {
			ds[j][0] = ds[j][1];
			ds[j][1] = ds[j][2];
			FEE[j][0] = FEE[j][1];
			PIN[j][0] = PIN[j][1];
			P[j][0] = P[j][1];
			PT[j][0] = PT[j][1];
		} //PIN+P is the dynamic pressure

		if (tfm == 0) {
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < elenode2D; k++) {
						if (debug == 0) {
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1];
						}
						else {
							nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1];
							//nr[z].XNRB[nr[z].DP_2D[k] - 1][j][0] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1]; //U(n-1)=U(n)
							//nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB[nr[z].DP_2D[k] - 1][j][2]; //U(n)=U(n+1)
						}
					}
				}
			}
		}

		for (z = 0; z < owsfnumber; z++) {
			if (mappingalgo == 4 || mappingalgo == 5) {
				for (j = 0; j < ol[z].FSNEL; j++) {
					for (k = 0; k < elenode2D_gs; k++) {
						ol[z].DISP_norm[j*elenode2D_gs + k][0] = ol[z].DISP_norm[j*elenode2D_gs + k][1];
					}
				}
			}
			else {
				for (j = 0; j < ol[z].FSNEL; j++) {
					for (k = 0; k < elenode2D; k++) {
						ol[z].DISP_norm[j*elenode2D + k][0] = ol[z].DISP_norm[j*elenode2D + k][1];
					}
				}
			}
		}
		//Output the pressure history under a specified point
		//extern double BF_val[4];
		//energyfilehd << current_time << " " << nr[1].angle_disp1[nr[1].DP_2D[2] - 1][291] << " " << nr[1].angle_disp2[nr[1].DP_2D[2] - 1][291] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][0] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][1] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][2] << std::endl;
		//energyfilehd << current_time << " " << nr[0].angle_disp1[nr[0].DP_2D[2] - 1][94] << " " << nr[0].angle_disp2[nr[0].DP_2D[2] - 1][94] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][0] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][1] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][2] << " " << P[nr[0].IEN_gb[nr[0].DP_2D[2] - 1][94] - 1][1] << std::endl;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	myfile2 << "total CPU time: " << duration << std::endl;
	return;
}