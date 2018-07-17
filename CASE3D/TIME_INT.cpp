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


//NRB determines the NRB local node numbering and the associated NRB arrays
void TIME_INT(int NNODE, double** GCOORD, int***LNA_3D, int**IEN, int NEL, int TIME, double *T, double DT, int NDT, 
	double*** HMASTER, double* Q, double KAPPA, double PPEAK, double TAU, double XC, double YC, 
	double ZC, double XO, double YO, double ZO, double ***SHOD, double** gamman, double** gamma_tn, double****Gn, double****SHG,
	double****gamma_t, double ****gamma, double*****G) {

	int h, i, j, k, q, z, ii, jj, kk, m;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber]; 

	int ctt1 = 0;
	int* counter1 = new int[NINT*NINT*NINT];
	int* counter2 = new int[NINT*NINT*NINT];
	int* counter3 = new int[NINT*NINT*NINT];
	int* counter4 = new int[NINT*NINT*NINT];
	int* counter5 = new int[NINT*NINT*NINT];
	int* counter6 = new int[NINT*NINT*NINT];
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
	int* LNAct1 = new int[NINT*NINT*NINT];
	int* LNAct2 = new int[NINT*NINT*NINT];
	int* LNAct3 = new int[NINT*NINT*NINT];
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

	double* SHOD1 = new double[NINT*NINT];
	for (i = 0; i < NINT; i++) {
		for (j = 0; j < NINT; j++) {
			SHOD1[i*NINT + j] = SHOD[1][i][j];
		}
	}

	ctt1 = 0;
	int*IENct1 = new int[NINT*NINT*NINT*NEL];
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			IENct1[ctt1] = IEN[LNAct1[j]][i] - 1;
			ctt1 += 1;
		}
	}

	ctt2 = 0;
	int*IENct2 = new int[NINT*NINT*NINT*NEL];
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			IENct2[ctt2] = IEN[LNAct2[j]][i] - 1;
			ctt2 += 1;
		}
	}

	ctt3 = 0;
	int*IENct3 = new int[NINT*NINT*NINT*NEL];
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			IENct3[ctt3] = IEN[LNAct3[j]][i] - 1;
			ctt3 += 1;
		}
	}

	NRBstruct nrb; 
	nrb = NRB(NNODE, GCOORD, LNA_3D);
	//NRB(NNODE, GCOORD, LNA_3D);

	//time history record
	std::clock_t start;
	double duration;
	double duration_int;
	start = std::clock();
	printf("3DFrigate - computation\n");
	fflush(stdout);
	//Obtain system time information:
	auto ti = std::time(nullptr);
	auto tm = *std::localtime(&ti);
	std::ostringstream oss;
	oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
	auto timestr = oss.str();

	//dynamic variables initiation
	double *DSDOT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF CONDENSATION
	double **ds;   //SOLUTION ARRAY FOR CONDENSATION
	double **FEEDOT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF DISP. POTENTIAL
	double **FEEDOT_inc;
	double **FEEDOT_INC;
	double **FEE;  //SOLUTION ARRAY FOR DISP. POTENTIAL
	double **FEE_inc;
	double **P; //SOLUTION ARRAY FOR DYNAMIC/SCATTERED PRESSURE
	double **PT; //SOLUTION ARRAY FOR TOTAL PRESSURE
	double *PH;
	double** PIN;
	double *BS; //STRUCTURE FORCE TO FLUID 
	double *FFORCE; //ARRAY OF INTERNAL + EXTERNAL FORCE ON FLUID
	double *BNRB; //SOLUTION ARRAY FOR FORCE OF NRBC ON FLUID
	double *HF; //GLOBAL REACTANCE MATIRX
	double *HFn;
	//double HFTEMP[(N + 1)*(N + 1)*(N + 1)];  //REACTANCE MATIRX PRODUCT VARAIBLE
	double **HFTEMP;
	HFTEMP = new double*[NEL];
	for (i = 0; i < NEL; i++) {
		HFTEMP[i] = new double[NINT*NINT*NINT];
	}
	double HFTEMPn[(N + 1)*(N + 1)*(N + 1)];
	double FEETEMP[(N + 1)*(N + 1)*(N + 1)]; //LOCAL DISP. POTENTIAL
	double *DPS_ukn;
	double *XNRBORG; //SOLUTION ARRAY FOR NRBC DISPLACEMENT AT T=0 (ORG MEANS ORIGIN)
	//double *XCOR; //SOLUTION ARRAY FOR NRB
	double *XCOR_kn;
	//double *XCOR_ukn;
	double** PSI_inc;
	double *DPS;
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
	FEEDOT_INC = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEEDOT_INC[i] = new double[2];
	}
	FEEDOT_inc = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEEDOT_inc[i] = new double[2];
	}
	FEE = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEE[i] = new double[2];
	}
	FEE_inc = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		FEE_inc[i] = new double[2];
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
	BS = new double[NNODE];
	FFORCE = new double[NNODE];
	BNRB = new double[NNODE];
	//BNRB = new double[u.NRBNODE];
	HF = new double[NNODE];
	HFn = new double[NNODE];
	//DPS_ukn = new double[u.NRBNODE];
	DPS_ukn = new double[NNODE];

	//XEST_kn = new double[NNODE];
	for (z = 0; z < nrbsurfnumber; z++) {
		nr[z].XEST_kn = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XEST_kn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XEST = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XEST[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XEST_ukn = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XEST_ukn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].angle_disp1 = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].angle_disp1[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].angle_disp2 = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].angle_disp2[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].disp_mag = new double**[NINT*NINT];
		for (i = 0; i < NINT*NINT;i++) {
			nr[z].disp_mag[i] = new double*[nr[z].NEL_nrb];
			for (j = 0; j < nr[z].NEL_nrb;j++) {
				nr[z].disp_mag[i][j] = new double[2];
			}
		}
		nr[z].XNRBORG2 = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XNRBORG2[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XNRB_kn = new double**[NINT*NINT];
		nr[z].XNRB_ukn = new double**[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XNRB_kn[i] = new double*[nr[z].NEL_nrb];
			nr[z].XNRB_ukn[i] = new double*[nr[z].NEL_nrb];
			for (j = 0; j < nr[z].NEL_nrb; j++) {
				nr[z].XNRB_kn[i][j] = new double[2];
				nr[z].XNRB_ukn[i][j] = new double[2];
			}
		}
		nr[z].XCOR_ukn = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].XCOR_ukn[i] = new double[nr[z].NEL_nrb];
		}
		nr[z].XCOR = new double*[NINT*NINT];
		for (i = 0; i < NINT*NINT; i++) {
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
			FEEDOT_INC[i][j] = 0.0;
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
			FEE_inc[i][j] = 0.0;
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
		BS[i] = 0.0;
	}
	for (i = 0; i < NNODE; i++) {
		FFORCE[i] = 0.0;
	}

	for (i = 0; i < NNODE; i++) {
		BNRB[i] = 0.0;
	}
	
	for (i = 0; i < NNODE; i++) {
		HF[i] = 0.0;
		HFn[i] = 0.0;
	}
	for (i = 0; i < (N + 1)*(N + 1)*(N + 1); i++) {
		//HFTEMP[i] = 0.0;
		HFTEMPn[i] = 0.0;
		FEETEMP[i] = 0.0;
	}
	
	for (i = 0; i < NNODE; i++) {
		DPS_ukn[i] = 0.0;
		XNRBORG[i] = 0.0;
		//XCOR[i] = 0.0;
		XCOR_kn[i] = 0.0;
		//XCOR_ukn[i] = 0.0;
	}

	for (z = 0; z < nrbsurfnumber; z++) {
		for (i = 0; i < NINT*NINT; i++) {
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
		for (i = 0; i < NINT*NINT; i++) {
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
			ol[z].PSI = new double[NNODE];
			ol[z].DI = new double[NNODE];
			ol[z].DISPI = new double*[NNODE];
			ol[z].WPIN = new double[NNODE];
			for (j = 0; j < NNODE; j++) {
				ol[z].DISPI[j] = new double[2];
			}
		}
		if (debug == 0) {
			for (z = 0; z < nrbsurfnumber; z++) {
				nr[z].XNRB = new double**[NINT*NINT];
				for (i = 0; i < NINT*NINT; i++) {
					nr[z].XNRB[i] = new double*[nr[z].NEL_nrb];
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						nr[z].XNRB[i][j] = new double[2];
					}
				}
				for (i = 0; i < NINT*NINT; i++) {
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
				nr[z].XNRB = new double**[NINT*NINT];
				for (i = 0; i < NINT*NINT; i++) {
					nr[z].XNRB[i] = new double*[nr[z].NEL_nrb];
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						nr[z].XNRB[i][j] = new double[2];
					}
				}
				for (i = 0; i < NINT*NINT; i++) {
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < 2; k++) {
							nr[z].XNRB[i][j][k] = 0.0;
						}
					}
				}
			}
		}

		DPS = new double[NNODE];
		//initialization
		for (i = 0; i < NNODE; i++) {
			for (j = 0; j < 2; j++) {
				PSI_inc[i][j] = 0.0;
			}
		}
		for (z = 0; z < owsfnumber; z++) {
			for (j = 0; j < NNODE; j++) {
				ol[z].PSI[j] = 0.0;
				ol[z].DI[j] = 0.0;
				ol[z].WPIN[j] = 0.0;
				for (k = 0; k < 2; k++) {
					ol[z].DISPI[j][k] = 0.0;
				}
			}
		}

		for (i = 0; i < NNODE; i++) {
			DPS[i] = 0.0;
		}
	}

	for (z = 0; z < owsfnumber; z++) { //this memory allocation scheme could have been improved
		ol[z].WBS = new double[NNODE];
		ol[z].DISP = new double*[NNODE];
		ol[z].DISP_norm = new double*[NNODE];
		ol[z].WP = new double[NNODE];
		for (j = 0; j < NNODE; j++) {
			ol[z].DISP[j] = new double[3]; //defined in 3 directions 
			ol[z].DISP_norm[j] = new double[2]; //norm displacement for two concecutive time steps
		}
	}
	for (z = 0; z < owsfnumber; z++) {
		for (j = 0; j < NNODE; j++) {
			ol[z].WBS[j] = 0.0;
			ol[z].WP[j] = 0.0;
			for (k = 0; k < 3; k++) {
				ol[z].DISP[j][k] = 0.0;
			}
			for (k = 0; k < 2; k++) {
				ol[z].DISP_norm[j][k] = 0.0;
			}
		}
	}

	//Initialize the incident pressure and total pressure along with hydrostatic pressure for the first time step
	TIME = 1;
	for (i = 0; i < NNODE; i++) {
		PH[i] = abs(GCOORD[i][1]) * RHO * grav + PATM; //absolute pressure
	}
	PIN = WAVE_IN(NNODE, GCOORD, T, TIME, PIN, DT, PPEAK, TAU, XC, YC, ZC, XO, YO, ZO);
	
	double hd2 = 0.0; 
	for (i = 0; i < NNODE;i++) {
		hd2 += PH[i];
	}

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
	for (i = 0; i < NNODE; i++) {
		hd2 += PT[i][0];
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
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < nr[z].NRBNODE; j++) {
					for (k = 0; k < nrbsurfnumber; k++) {
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
	}

	for (i = 0; i < NNODE; i++) {
		FEE_inc[i][0] = FEE[i][0];
	}
	for (i = 0; i < NNODE; i++) {
		FEEDOT_INC[i][0] = FEEDOT[i][0];
	}
	//=========================Get free surface points=========================//
	int fspt_num = 0;
	for (j = 0; j < NNODE; j++) {
		if (fsdebug == 0) {
			if (abs(GCOORD[j][1]) <= 1e-6) {
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
			if (abs(GCOORD[j][1]) <= 1e-6) {
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
		std::cout << "the free surface point recognization is wrong or there's no free surface" << std::endl;
	}

	//Initiate the first coupling
	initcoupling();

	//Update the current MpCCI time
	current_time = T[0];

	interface_mappingstruct in;
	in = interface_mapping(1, GCOORD);
	dotransfer();
	in = interface_mapping(0, GCOORD);

	//Generate the information file
	std::string name3 = "parameters_" + timestr + ".txt";
	std::ofstream myfile;
	myfile.open(name3);
	//General information:
	myfile << "Code name: 3Dbarge_TFM_Abaqus_sym_mesh" << std::endl;
	myfile << "Simulation date and time: " << timestr << std::endl;
	myfile << "Mesh information:" << std::endl;
	myfile << "N: " << N << " mesh size: "  << " wetted surface number: " << owsfnumber << " Nodenumber: " << NNODE << " Element number: " << NEL << " h/p refinement rate: " << refine << std::endl;
	myfile << "Explosive information: " << std::endl;
	myfile << "stdoff point (spherical): " << XO << " " << YO << " " << ZO << " explosion center: " << XC << " " << YC << " " << ZC << " peak pressure (Mpa): " << PPEAK << " decay rate: " << TAU << " WAVE (1: plane, 2: spherical): " << WAVE << std::endl;
	myfile << "Time integration information" << std::endl;
	myfile << "CFL: " << CFLFRAC << " dt: " << DT << " dt scale factor: " << dtscale << " total time: " << TTERM << " damping: " << BETA << " explicit central difference (Leap frog)" << std::endl;
	myfile << "Fluid properties: " << std::endl;
	myfile << "CAV: " << CAV << " C: " << C << " RHO: " << RHO << " atmospheric pressure: " << PATM << " saturated pressure: " << PSAT << std::endl;
	myfile << "FSI coupling" << std::endl;
	myfile << "mapping algorithm: " << mappingalgo << std::endl;
	myfile << "debug mode: " << debug << std::endl;

	double angle = 0.0; double R = 0.0;

	double* WPTEMP;
	double* BFTEMP;
	double* DISPTEMP;
	double* WBSTEMP; 
	double* NRBDISPTEMP;
	double* BNRBTEMP;
	WPTEMP = new double[NINT*NINT];
	BFTEMP = new double[NINT*NINT];
	DISPTEMP = new double[NINT*NINT];
	WBSTEMP = new double[NINT*NINT];
	NRBDISPTEMP = new double[NINT*NINT];
	BNRBTEMP = new double[NINT*NINT];
	for (i = 0; i < NINT*NINT; i++) {
		WPTEMP[i] = 0.0;
		BFTEMP[i] = 0.0;
		NRBDISPTEMP[i] = 0.0;
		BNRBTEMP[i] = 0.0;
	}
	for (i = 0; i < NINT*NINT; i++) {
		DISPTEMP[i] = 0.0;
		WBSTEMP[i] = 0.0;
	}
	
	double BNRBt = 0.0; 
	double dst = 0.0;
	//=================================Start the time integration routine=================================//
	int oc = 0;
	double lastDT = 0.0;

	std::string energyfile = "history.txt";
	std::ofstream energyfilehd;
	energyfilehd.open(energyfile);

	//Get the sample points on a line to observe the wave propagation pressure distribution
	std::vector<int> sampline;
	count = 0;
	for (i = 0; i < NNODE; i++) {
		if (abs(GCOORD[i][0] - 0.0) < 1e-6 && abs(GCOORD[i][2] - 0.0) < 1e-6) {
			sampline.push_back(i + 1);
			count += 1;
		}
	}
	//sort the arrary sampline and store it in sampline2
	int* hold;
	hold = new int[sampline.size()];
	for (i = 0; i < sampline.size(); i++) {
		hold[i] = round(abs(GCOORD[sampline[i] - 1][1]) / YHE) - round(SY / YHE);
	}
	int* sampline2;
	sampline2 = new int[sampline.size()];
	for (i = 0; i < sampline.size(); i++) {
		sampline2[hold[i]] = sampline[i];
	}
	/*
	int* sampline2;
	int** samplinec2;
	sampline2 = new int[count];
	samplinec2 = new int*[4];
	for (i = 0; i < 4; i++) {
		samplinec2[i] = new int[count];
	}
	for (i = 0; i < count; i++) {
		sampline2[i] = 0;
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < count; j++) {
			samplinec2[i][j] = 0;
		}
	}
	int YNEL = round((SY + DY) / YHE);
	double** AYIN;
	AYIN = new double*[YNEL];
	for (i = 0; i < YNEL; i++) {
		AYIN[i] = new double[N + 1];
	}
	LOBATTOstruct b;
	if (N > 1) {
		b = LOBATTO(N);
	}
	double*YE;
	YE = new double[YNEL + 1];
	for (i = 0; i < YNEL + 1; i++) {
		YE[i] = i*YHE;
	}
	double MY; double dY;
	for (i = 0; i < YNEL; i++) {
		MY = 0.5*(YE[i + 1] + YE[i]);
		dY = 0.5*(YE[i + 1] - YE[i]);
		AYIN[i][0] = YE[i];
		if (N > 1) {
			for (j = 0; j < N - 1; j++) {
				AYIN[i][j + 1] = MY + b.Z[j] * dY;
			}
		}
		AYIN[i][N] = YE[i + 1];
	}
	double* ypt;
	ypt = new double[count];
	int cnt = 0;
	for (i = 0; i < round(DY / YHE); i++) {
		for (j = 0; j < NINT; j++) {
			ypt[cnt] = AYIN[SYNEL + i][j];
			cnt += 1;
		}
		cnt -= 1;
	}
	if (cnt + 1 != count) {
		std::cout << "ypt is wrong" << std::endl;
		//system("PAUSE ");
	}
	int count2 = 0;
	while (count2 != count) {
		for (i = 0; i < count; i++) {
			if (abs(GCOORD[sampline[i] - 1][1] + ypt[count2]) < 1e-5) {
				sampline2[count2] = sampline[i];
				count2 += 1;
			}
		}
	}
	*/

	std::vector <int> T_out;
	int NDT_out = 0;
	for (i = 0; i < NDT; i++) {
		if (T[i] > output_int*NDT_out) {
			T_out.push_back(i);
			NDT_out += 1;
		}
	}
	std::string name2 = "PT_line_result_";
	std::ofstream *outline;
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

	//Initialize the pressure gradient vector
	for (z = 0; z < nrbsurfnumber; z++) {
		nr[z].P_dev = new double**[nr[z].NEL_nrb];
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			nr[z].P_dev[i] = new double*[NINT*NINT];
			for (j = 0; j < NINT*NINT; j++) {
				nr[z].P_dev[i][j] = new double[3];
			}
		}
	}
	for (z = 0; z < nrbsurfnumber; z++) {
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			for (j = 0; j < NINT*NINT; j++) {
				for (k = 0; k < 3; k++) {
					nr[z].P_dev[i][j][k] = 0.0;
				}
			}
		}
	}

	if (Bleich == 1) {
		int wavdirc[3] = { 0,1,0 };
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

		for (z = 0; z < owsfnumber; z++) {
			//for (k = 0; k < NNODE; k++) {
			for (k = 0; k < ol[z].GIDNct; k++) {
				//ol[j].WP[k] = 0.0;
				ol[z].WBS[ol[z].GIDN[k] - 1] = 0.0;
			}
		}

		if (tfm == 1) {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					if (nodeforcemap2 == 1) {
						ol[z].WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM; //correct version with structural gravity
						//ol[z].WP[ol[z].GIDN[j] - 1] = PH[ol[z].GIDN[j] - 1] - PATM; //pure hydrostatic pressure
						//ol[z].WP[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][0]; //incident pressure
					}
					else { //absolute pressure
						ol[z].WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM;
						//ol[z].WP[ol[z].GIDN[j] - 1] = PH[ol[z].GIDN[j] - 1] - PATM; //pure hydrostatic pressure
						//ol[z].WP[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][0]; //incident pressure
					}
					//ol[z].WP[ol[z].GIDN[j] - 1]
					//= PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]; //correct version if structural gravity is not specified in Abaqus
					if (Bleich == 1) {
						if (nodeforcemap2 == 1) {
							ol[z].WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
						}
						else {
							ol[z].WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
						}
					}
				}
			}
		}
		else {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					ol[z].WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM; // //correct version with structural gravity
					if (Bleich == 1) {
						ol[z].WP[ol[z].GIDN[j] - 1] = (PT[ol[z].GIDN[j] - 1][0] - PH[ol[z].GIDN[j] - 1]) / SX / SZ;
					}
					ol[z].WPIN[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][1] + PIN[ol[z].GIDN[j] - 1][0];
				}
			}
		}

		//start = std::clock();
		//=======================define double* nodeforce in fluid code==========================//
		//mapping the fluid force ABF from user defined mesh to MpCCI defined mesh on coupling surface using interpolation
		in = interface_mapping(1, GCOORD);
		//after this subroutine, the nodeforce should already be mapped onto coupling surface (data.h)
		//int fluid2structure, int**IEN_3D, int***LNA_3D, int**LNA_2D, int NNODE, double *Z

		dotransfer();

		//=============map nodal displacement from coupled surface to fluid mesh===================//
		in = interface_mapping(0, GCOORD);

		if (tfm == 0) { //Scattered field model
			//double angle = 0.0; //cos value
			double r = 0.0;
			double angle = 0.0;
			//===============Incident fluid displacement on FSI boundary (used to derive structure force)================//
			//Change the definition of DISPI to node based rather than element based. 
			for (z = 0; z < owsfnumber; z++) {
				angle = 0.0;
				for (j = 0; j < ol[z].FSNEL; j++) {  //FOR STRUCTURE ELEMENT ON THE FSI BOUNDARY
					for (k = 0; k < NINT*NINT; k++) {
						if (WAVE == 1) {
							angle = ol[z].norm[j][0] * wavdirc[0] + ol[z].norm[j][1] * wavdirc[1] + ol[z].norm[j][2] * wavdirc[2];
							ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
								= (DT / 2.0)*(ol[z].WPIN[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]); //INCIDENT DISP. PREDICTOR by time integration	
							ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
								= ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] / (RHO*C);     //trapezoidal integration of the double integrator
							ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]
								= ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + angle*ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1];
						}
						else { //Spherical wave 
							r = sqrt(pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] - XC), 2) + pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - YC), 2) + pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2] - ZC), 2));
							ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
								= (DT / 2.0)*(ol[z].WPIN[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]); //INCIDENT DISP. PREDICTOR by time integration	
							PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] = ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1];
							angle = (ol[z].norm[j][0] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] - XC) / r + ol[z].norm[j][1] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - YC) / r + ol[z].norm[j][2] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2] - ZC) / r);
							ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
								= ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] / (RHO*C) + (PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]) * (DT / 2) / (RHO * r);  //trapezoidal integration of the double integrator
							ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]
								= ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + angle*ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]; //UPDATED INCIDENT STRUCUTRE DISPLACEMENT 
						}
					}
				}
			}
			//SET VALUES AT T(i+1) TO VALUES AT T(i)
			for (j = 0; j < NNODE; j++) {
				for (z = 0; z < owsfnumber; z++) {
					ol[z].DISPI[j][0] = ol[z].DISPI[j][1];  //SOLUTION ARRAY FOR INCIDENT STRUCUTRE DISPLACEMENT
				} //DISPI[1] = DISPI[0] + DI
				PSI_inc[j][0] = PSI_inc[j][1];
			}
		}

		if (tfm == 1) { //Total field model
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].FSNEL; j++) { //for every fluid element linked to structure
					for (h = 0; h < NINT*NINT; h++) {
						WBSTEMP[h] = 0.0;
						for (k = 0; k < NINT*NINT; k++) {
							ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] = ol[z].norm[j][0] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + ol[z].norm[j][1] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + ol[z].norm[j][2] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2];
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * (-1) * RHO * (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0]));
						}
					}
					for (k = 0; k < NINT*NINT; k++) {
						ol[z].WBS[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] += WBSTEMP[k];
						//std::cout << "" << std::endl;
					}
				}
			}
		}
		else {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].FSNEL; j++) { //for every fluid element linked to structure
					for (h = 0; h < NINT*NINT; h++) {
						WBSTEMP[h] = 0.0;
						for (k = 0; k < NINT*NINT; k++) {
							ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] = ol[z].norm[j][0] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + ol[z].norm[j][1] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + ol[z].norm[j][2] * ol[z].DISP[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2];
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * (-1) * RHO * (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0]) - ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]); //Error prone: Get the normal vector for DISPI??? 
						}
					}
					for (k = 0; k < NINT*NINT; k++) {
						ol[z].WBS[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] += WBSTEMP[k];
					}
				}
			}
		}

		hd = 0.0;
		for (z = 0; z < owsfnumber; z++) {
			for (j = 0; j < NNODE; j++) {
				hd += ol[z].WBS[j];
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
					for (jj = 0; jj < NINT*NINT; jj++) {
						for (kk = 0; kk < 3; kk++) {
							nr[z].P_dev[ii][jj][kk] = 0.0;
						}
					}
				}
			}

			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < NINT*NINT; k++) {
						//loop through each point in element and accumulate the value
						for (ii = 0; ii < NINT; ii++) {
							for (jj = 0; jj < NINT; jj++) {
								for (kk = 0; kk < NINT; kk++) {
									for (h = 0; h < 3; h++) { //the three directions (x y z )
										if (debug == 0) {
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
										}
										else {
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * PH[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
											//although we've devided the RHO, the value is still quite large.
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * P[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * (FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] - FEE_inc[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0]) / RHO;
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEEDOT[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][0] / RHO;
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
						for (k = 0; k < NINT*NINT; k++) {
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
						for (k = 0; k < NINT*NINT; k++) {
							nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j];
							nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j];
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (k = 0; k < NINT*NINT; k++) {
							//if (debug == 0) {
							NRBDISPTEMP[k] = -RHO* (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							//}
							//else {
								//NRBDISPTEMP[k] = -RHO* (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1]);
							//}
							//std::cout << " " << std::endl;
						}
						for (h = 0; h < NINT*NINT; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < NINT*NINT; k++) {
								BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * NRBDISPTEMP[k];
							}
						}
						for (k = 0; k < NINT*NINT; k++) {
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
						for (k = 0; k < NINT*NINT; k++) {
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
						for (k = 0; k < NINT*NINT; k++) {
							nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_kn[nr[z].DP_2D[k] - 1][j];
							nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XEST_ukn[nr[z].DP_2D[k] - 1][j];
						}
					}
					for (j = 0; j < nr[z].NEL_nrb; j++) { //IEN_gb needs to be changed (NINT*N)
						for (h = 0; h < NINT*NINT; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < NINT*NINT; k++) {
								BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (-RHO) * (nr[z].XNRB_kn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] + nr[z].XNRBORG2[nr[z].DP_2D[k] - 1][j]);
							}
						}
						for (k = 0; k < NINT*NINT; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
		}
		else { //SFM
			for (z = 0; z < nrbsurfnumber; z++) { //from bottom to front surface (k=0 means all node case)
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < NINT*NINT; k++) {
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
					for (k = 0; k < NINT*NINT; k++) {
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
						for (k = 0; k < NINT*NINT; k++) {
							if (debug == 0) {
								NRBDISPTEMP[k] = -RHO* (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							}
							else {
								NRBDISPTEMP[k] = -RHO* (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
							}
						}
						for (h = 0; h < NINT*NINT; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < NINT*NINT; k++) {
								BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * NRBDISPTEMP[k];
							}
						}
						for (k = 0; k < NINT*NINT; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
			else { //Spherical wave
				for (z = 0; z < nrbsurfnumber; z++) { //from bottom to front surface (k=0 means all node case)
					//XNRB is the predicted displacement normal to the structure (since P is normal to the surface)
					for (j = 0; j < nr[z].NEL_nrb; j++) {
						for (h = 0; h < NINT*NINT; h++) {
							BNRBTEMP[h] = 0.0;
							for (k = 0; k < NINT*NINT; k++) {
								if (debug == 0) {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (-RHO) * (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
								}
								else {
									BNRBTEMP[h] += nr[z].ADMASTER[j][h][k] * (-RHO) * (nr[z].XNRB[nr[z].DP_2D[k] - 1][j][1] - XNRBORG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1]);
								}
							}
						}
						for (k = 0; k < NINT*NINT; k++) {
							BNRB[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1] += BNRBTEMP[k];
						}
					}
				}
			}
		}
		/*
		hd = 0.0;
		for (j = 0; j < NNODE; j++) {
			hd += ol[0].WBS[j];
		}
		*/	
		//time integration
		//start = std::clock();
		#pragma omp parallel for num_threads(6)
		for (int j = 0; j < NNODE; j++) {
			DSDOT[j] = (ds[j][1] - ds[j][0]) / DT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF CONDENSATION
			FEEDOT[j][1] = FEEDOT[j][0] + DT*(P[j][0] + (BETA*DT*(pow(C, 2))*DSDOT[j]));
			FEE[j][1] = FEE[j][0] + DT*FEEDOT[j][1];
			FEEDOT[j][0] = FEEDOT[j][1];  //USE TWO VALUE TO STORE THE N TIME STEP AND N+1 TIME STEP
			if (debug == 1) {
				FEEDOT_INC[j][1] = FEEDOT_INC[j][0] + 0.5*DT*(PIN[j][0] + PIN[j][1]);
				FEE_inc[j][1] = FEE_inc[j][0] + DT*FEEDOT_INC[j][1];
				FEEDOT_INC[j][0] = FEEDOT_INC[j][1];  //USE TWO VALUE TO STORE THE N TIME STEP AND N+1 TIME STEP
			}
			HF[j] = 0.0;
			HFn[j] = 0.0;
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
				for (int z = 0; z < NINT*NINT*NINT; z++) { //takes NEL*(2*NINT^3+2) operations
					HFTEMP[j][z] = 0.0; //takes NEL*NINT^3 operations
					for (int k = 0; k < NINT*NINT*NINT; k++) { //takes NEL*NINT^3*(2*NINT^3+2) operaitons
						HFTEMP[j][z] = HFTEMP[j][z] + HMASTER[j][z][k] * FEE[IEN[k][j] - 1][1];
					} //The premise is that HMASTER is the same for every element (only true for equally shaped element!!)
				}
			}
			for (int j = 0; j < NEL; j++) {
				for (int k = 0; k < NINT*NINT*NINT; k++) { //takes NEL*(2*NINT^3+2) operations
					HF[IEN[k][j] - 1] = HF[IEN[k][j] - 1] + HFTEMP[j][k];   //global level reactance matrix //takes 2*NINT^3*NEL operations
				} //HF is like global stiffness matrix (include variable in it)
			}
			//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
			//std::cout << "total CPU time (ms): " << duration << std::endl;
			//std::cout << " " << std::endl;
		}
		else {
			//#pragma omp parallel for num_threads(6)
			for (j = 0; j < NEL; j++) { //takes 2*NEL+2 operations 
				//takes 6*NEL*NINT^4 FLOP
				//ctt1 = 0;
				for (ii = 0; ii < NINT; ii++) { //p NEL*(2*NINT+2) operations
					for (h = 0; h < NINT; h++) { //j,i,i //takes NEL*NINT*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k,k,j //takes NEL*NINT^2*(2*NINT+2) operations
							gamman[counter1[ii*NINT*NINT + h*NINT + z]][0] = 0.0;
							gamman[counter2[ii*NINT*NINT + h*NINT + z]][1] = 0.0;
							gamman[counter3[ii*NINT*NINT + h*NINT + z]][2] = 0.0;
							for (q = 0; q < NINT; q++) { //q (recurrent addition in this dimension) //takes NEL*NINT^3*(2*NINT+2) operations	
								gamman[counter1[ii*NINT*NINT + h*NINT + z]][0] += SHOD1[q*NINT + ii] * FEE[IENct1[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								gamman[counter2[ii*NINT*NINT + h*NINT + z]][1] += SHOD1[q*NINT + ii] * FEE[IENct2[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								gamman[counter3[ii*NINT*NINT + h*NINT + z]][2] += SHOD1[q*NINT + ii] * FEE[IENct3[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
							}
							//ctt1 += 1;
						}
					}
				}

				//takes 15*NINT^3*NEL FLOP
				for (ii = 0; ii < NINT*NINT*NINT; ii++) { //takes NEL*(2 * NINT + 2) operations
					gamma_tn[ii][0] = Gn[j][0][0][ii] * gamman[ii][0] + Gn[j][0][1][ii] * gamman[ii][1] + Gn[j][0][2][ii] * gamman[ii][2];
					//ipk
					gamma_tn[ii][1] = Gn[j][1][1][ii] * gamman[ii][1] + Gn[j][0][1][ii] * gamman[ii][0] + Gn[j][1][2][ii] * gamman[ii][2];
					//ijp
					gamma_tn[ii][2] = Gn[j][2][2][ii] * gamman[ii][2] + Gn[j][0][2][ii] * gamman[ii][0] + Gn[j][1][2][ii] * gamman[ii][1];
					//oc += 15;
				}

				//takes 6*NINT^4*NEL FLOP
				//ctt2 = 0;
				for (h = 0; h < NINT; h++) {  //i //takes NEL*NINT*(2*NINT+2) operations
					for (k = 0; k < NINT; k++) { //j  //takes NEL*NINT^2*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k //takes NEL*NINT^3*(2*NINT+2) operations
							HFTEMPn[LNAct3[h*NINT*NINT + k*NINT + z]] = 0.0;
							for (ii = 0; ii < NINT; ii++) { //p (recurrent addition) //p NEL*(2*NINT+2) operations
								HFTEMPn[LNAct3[h*NINT*NINT + k*NINT + z]] += SHOD1[h*NINT + ii] * gamma_tn[ii*NINT*NINT + k*NINT + z][0] + SHOD1[k*NINT + ii] * gamma_tn[h*NINT*NINT + ii*NINT + z][1] + SHOD1[z*NINT + ii] * gamma_tn[h*NINT*NINT + k*NINT + ii][2];
							}
							//ctt2 += 1;
						}
					}
				}
				//takes 1*NINT^3*NEL FLOP
				for (h = 0; h < NINT*NINT*NINT; h++) { //takes NEL*(2 * NINT + 2) operations
					HF[IENct3[j*NINT*NINT*NINT + h]] += HFTEMPn[LNAct3[h]]; //takes NEL*2*NINT^3 operations
					//ctt3 += 1;
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
		}

		for (j = 0; j < NNODE; j++) {
			FFORCE[j] = -HF[j];
		}

		for (z = 0; z < owsfnumber; z++) {
			for (j = 0; j < ol[z].GIDNct; j++) { //for every fluid element linked to structure
				FFORCE[ol[z].GIDN[j] - 1] += ol[z].WBS[ol[z].GIDN[j] - 1];
				//count += 1;
			}
		}

		for (j = 0; j < nrb.NNODE_nrb; j++) {
			FFORCE[nrb.NRBA_t[j] - 1] += BNRB[nrb.NRBA_t[j] - 1];
		}
		//The combination of FFORCE passes the test

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
					for (jj = 0; jj < NINT*NINT; jj++) {
						for (kk = 0; kk < 3; kk++) {
							nr[z].P_dev[ii][jj][kk] = 0.0;
						}
					}
				}
			}
			for (z = 0; z < nrbsurfnumber; z++) {
				for (j = 0; j < nr[z].NEL_nrb; j++) {
					for (k = 0; k < NINT*NINT; k++) {
						//loop through each point in element and accumulate the value
						for (ii = 0; ii < NINT; ii++) {
							for (jj = 0; jj < NINT; jj++) {
								for (kk = 0; kk < NINT; kk++) {
									for (h = 0; h < 3; h++) { //the three directions (x y z )
										if (debug == 0) {
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
										}
										else {
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * P[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1];
											//nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * (FEE[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] - FEE_inc[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1]) / RHO;
											nr[z].P_dev[j][nr[z].DP_2D[k] - 1][h] += SHG[nr[z].NRBELE_ARR[j] - 1][h][LNA_3D[ii][jj][kk] - 1][nr[z].DP[k] - 1] * FEEDOT[IEN[LNA_3D[ii][jj][kk] - 1][nr[z].NRBELE_ARR[j] - 1] - 1][1] / RHO;
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
					for (k = 0; k < NINT*NINT; k++) {
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
					for (k = 0; k < NINT*NINT; k++) {
						//nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + XCOR_ukn[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][j] - 1];
						nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][1] = nr[z].XNRB_ukn[nr[z].DP_2D[k] - 1][j][0] + nr[z].XCOR_ukn[nr[z].DP_2D[k] - 1][j];
					}
				}
				for (j = 0; j < nr[z].NEL_nrb; j++) { //Need to be changed
					for (k = 0; k < NINT*NINT; k++) {
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
					for (k = 0; k < NINT*NINT; k++) {
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
					for (k = 0; k < NINT*NINT; k++) {
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
					for (k = 0; k < NINT*NINT; k++) {
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
			for (j = 0; j < ol[z].GIDNct; j++) {
				ol[z].DISP_norm[ol[z].GIDN[j] - 1][0] = ol[z].DISP_norm[ol[z].GIDN[j] - 1][1];
			}
		}


		if (output == 1) {
			for (k = 0; k < NDT_out; k++) {
				if (i == T_out[k]) {
					for (j = 0; j < sampline.size(); j++) {
						outline[k] << GCOORD[sampline2[j] - 1][1] << " " << P[sampline2[j] - 1][0] << " " << P[sampline2[j] - 1][1] << " " << PT[sampline2[j] - 1][1] << std::endl;
					}
				}
			}
		}

		//Output the pressure history under a specified point
		//extern double BF_val[4];
		//energyfilehd << current_time << " " << nr[1].angle_disp1[nr[1].DP_2D[2] - 1][291] << " " << nr[1].angle_disp2[nr[1].DP_2D[2] - 1][291] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][0] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][1] << " " << nr[1].P_dev[291][nr[1].DP_2D[2] - 1][2] << std::endl;
		energyfilehd << current_time << " " << nr[0].angle_disp1[nr[0].DP_2D[2] - 1][94] << " " << nr[0].angle_disp2[nr[0].DP_2D[2] - 1][94] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][0] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][1] << " " << nr[0].P_dev[94][nr[0].DP_2D[2] - 1][2] << " " << P[nr[0].IEN_gb[nr[0].DP_2D[2] - 1][94] - 1][1] << std::endl;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	myfile << "total CPU time: " << duration << std::endl;
	return;
}