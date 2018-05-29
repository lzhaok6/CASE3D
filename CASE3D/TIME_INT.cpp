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
void TIME_INT(int NNODE, double** GCOORD, double* W, int**LNA_2D, int***LNA_3D, int**IEN, int NEL, double* S, double***SHL,
	int TIME, double *T, double DT, int NDT, double* Z, double** AYIN, double*** HMASTER, double* Q, double***phi_fem, double* timer, double* ampt, double KAPPA,
	double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO, double ***SHOD, double ****gamma, double****gamma_t, double*****G,
	double** gamman, double** gamma_tn, double****Gn, double*** phi_femg, double*** phi_sem2) {
	int h, i, j, k, q, z, ii, jj, m;
	//extern int owsfnumber;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber];
	//extern int output;
	/*
	int*** counter = new int**[NINT];
	for (i = 0; i < NINT;i++) {
	counter[i] = new int*[NINT];
	for (j = 0; j < NINT;j++) {
	counter[i][j] = new int[NINT];
	}
	}

	int ctt = 0;
	for (i = 0; i < NINT; i++) {
	for (j = 0; j < NINT; j++) {
	for (k = 0; k < NINT; k++) {
	counter[i][j][k] = ctt;
	ctt += 1;
	}
	}
	}
	*/

	/*
	int ctt = 0;
	int** counter = new int*[NINT*NINT*NINT];
	for (i = 0; i < NINT*NINT*NINT; i++) {
	counter[i] = new int[3];
	}
	for (i = 0; i < NINT;i++) {
	for (j = 0; j < NINT;j++) {
	for (k = 0; k < NINT;k++) {
	counter[ctt][0] = i*NINT*NINT + j*NINT + k;
	counter[ctt][1] = j*NINT*NINT + i*NINT + k;
	counter[ctt][2] = j*NINT*NINT + k*NINT + i;
	ctt += 1;
	}
	}
	}
	*/

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

	NRB(NNODE, GCOORD, W, LNA_3D, IEN, NEL, SHL, SHOD);

	//time history record
	std::clock_t start;
	double duration;
	double duration_int;
	start = std::clock();
	printf("3Dbarge - computation\n");
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
	double **FEE;  //SOLUTION ARRAY FOR DISP. POTENTIAL
	double **P; //SOLUTION ARRAY FOR DYNAMIC/SCATTERED PRESSURE
	double **PT; //SOLUTION ARRAY FOR TOTAL PRESSURE
	double *PH;
	double** PIN;
	double *BS; //STRUCTURE FORCE TO FLUID 
	double *FFORCE; //ARRAY OF INTERNAL + EXTERNAL FORCE ON FLUID
	double *BNRB; //SOLUTION ARRAY FOR FORCE OF NRBC ON FLUID
	double *HF; //GLOBAL REACTANCE MATIRX
	double *HFn;
	double HFTEMP[(N + 1)*(N + 1)*(N + 1)];  //REACTANCE MATIRX PRODUCT VARAIBLE
	double HFTEMPn[(N + 1)*(N + 1)*(N + 1)];
	double FEETEMP[(N + 1)*(N + 1)*(N + 1)]; //LOCAL DISP. POTENTIAL
	double *DPS_ukn;
	double *XEST; //SOLUTION ARRAY FOR NRB  (ESTIMATED DISPLACEMENT???)
	double *XEST_kn;
	double *XEST_ukn;
	double *XNRBORG; //SOLUTION ARRAY FOR NRBC DISPLACEMENT AT T=0 (ORG MEANS ORIGIN)
	double *XNRBORG2;
	double *XCOR; //SOLUTION ARRAY FOR NRB
	double *XCOR_kn;
	double *XCOR_ukn;
	double **XNRB_kn;
	double **XNRB_ukn;
	double** PSI_inc;
	double **XNRB; 
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
	BS = new double[NNODE];
	FFORCE = new double[NNODE];
	//BNRB = new double[u.NRBNODE];
	BNRB = new double[NNODE];
	HF = new double[NNODE];
	HFn = new double[NNODE];
	//DPS_ukn = new double[u.NRBNODE];
	DPS_ukn = new double[NNODE];
	//XEST = new double[u.NRBNODE];
	XEST = new double[NNODE];
	//XEST_kn = new double[u.NRBNODE];
	XEST_kn = new double[NNODE];
	//XEST_ukn = new double[u.NRBNODE];
	XEST_ukn = new double[NNODE];
	//XNRBORG = new double[u.NRBNODE];
	XNRBORG = new double[NNODE];
	//XNRBORG2 = new double[u.NRBNODE];
	XNRBORG2 = new double[NNODE];
	//XCOR = new double[u.NRBNODE];
	XCOR = new double[NNODE];
	//XCOR_kn = new double[u.NRBNODE];
	XCOR_kn = new double[NNODE];
	//XCOR_ukn = new double[u.NRBNODE];
	XCOR_ukn = new double[NNODE];
	//XNRB_kn = new double*[u.NRBNODE];
	XNRB_kn = new double*[NNODE];
	/*
	for (i = 0; i < u.NRBNODE; i++) {
		XNRB_kn[i] = new double[2];
	}
	XNRB_ukn = new double*[u.NRBNODE];
	for (i = 0; i < u.NRBNODE; i++) {
		XNRB_ukn[i] = new double[2];
	}
	*/
	
	for (i = 0; i < NNODE; i++) {
		XNRB_kn[i] = new double[2];
	}
	XNRB_ukn = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		XNRB_ukn[i] = new double[2];
	}
	
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
		BS[i] = 0.0;
	}
	for (i = 0; i < NNODE; i++) {
		FFORCE[i] = 0.0;
	}
	/*
	for (i = 0; i < u.NRBNODE; i++) {
		BNRB[i] = 0.0;
	}
	*/
	
	for (i = 0; i < NNODE; i++) {
		BNRB[i] = 0.0;
	}
	
	for (i = 0; i < NNODE; i++) {
		HF[i] = 0.0;
		HFn[i] = 0.0;
	}
	for (i = 0; i < (N + 1)*(N + 1)*(N + 1); i++) {
		HFTEMP[i] = 0.0;
		HFTEMPn[i] = 0.0;
		FEETEMP[i] = 0.0;
	}
	/*
	for (i = 0; i < u.NRBNODE; i++) {
		DPS_ukn[i] = 0.0;
		XEST[i] = 0.0;
		XEST_kn[i] = 0.0;
		XEST_ukn[i] = 0.0;
		XNRBORG[i] = 0.0;
		XNRBORG2[i] = 0.0;
		XCOR[i] = 0.0;
		XCOR_kn[i] = 0.0;
		XCOR_ukn[i] = 0.0;
	}
	for (i = 0; i < u.NRBNODE; i++) {
		for (j = 0; j < 2; j++) {
			XNRB_kn[i][j] = 0.0;
			XNRB_ukn[i][j] = 0.0;
		}
	}
	*/
	for (i = 0; i < NNODE; i++) {
		DPS_ukn[i] = 0.0;
		XEST[i] = 0.0;
		XEST_kn[i] = 0.0;
		XEST_ukn[i] = 0.0;
		XNRBORG[i] = 0.0;
		XNRBORG2[i] = 0.0;
		XCOR[i] = 0.0;
		XCOR_kn[i] = 0.0;
		XCOR_ukn[i] = 0.0;
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < 2; j++) {
			XNRB_kn[i][j] = 0.0;
			XNRB_ukn[i][j] = 0.0;
		}
	}
	
	if (tfm == 0) {
		//declarition
		PSI_inc = new double*[NNODE];
		for (i = 0; i < NNODE; i++) {
			PSI_inc[i] = new double[2];
		}
		for (i = 0; i < owsfnumber; i++) { //this memory allocation scheme could have been improved
			ol[i].PSI = new double[NNODE];
			ol[i].DI = new double[NNODE];
			ol[i].DISPI = new double*[NNODE];
			ol[i].WPIN = new double[NNODE];
			for (j = 0; j < NNODE; j++) {
				ol[i].DISPI[j] = new double[2];
			}
		}
		XNRB = new double*[NNODE];
		for (i = 0; i < NNODE; i++) {
			XNRB[i] = new double[2];
		}
		DPS = new double[NNODE];
		//initialization
		for (i = 0; i < NNODE; i++) {
			for (j = 0; j < 2; j++) {
				PSI_inc[i][j] = 0.0;
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < NNODE; j++) {
				ol[i].PSI[j] = 0.0;
				ol[i].DI[j] = 0.0;
				ol[i].WPIN[j] = 0.0;
				for (k = 0; k < 2; k++) {
					ol[i].DISPI[j][k] = 0.0;
				}
			}
		}
		for (i = 0; i < NNODE; i++) {
			for (j = 0; j < 2; j++) {
				XNRB[i][j] = 0.0;
			}
		}
		for (i = 0; i < NNODE; i++) {
			DPS[i] = 0.0;
		}
	}

	for (i = 0; i < owsfnumber; i++) { //this memory allocation scheme could have been improved
		ol[i].WBS = new double[NNODE];
		ol[i].DISP = new double*[NNODE];
		ol[i].DISP_norm = new double*[NNODE];
		ol[i].WP = new double[NNODE];
		for (j = 0; j < NNODE; j++) {
			ol[i].DISP[j] = new double[3]; //defined in 3 directions 
			ol[i].DISP_norm[j] = new double[2]; //norm displacement for two concecutive time steps
		}
	}
	for (i = 0; i < owsfnumber; i++) {
		for (j = 0; j < NNODE; j++) {
			ol[i].WBS[j] = 0.0;
			ol[i].WP[j] = 0.0;
			for (k = 0; k < 3; k++) {
				ol[i].DISP[j][k] = 0.0;
			}
			for (k = 0; k < 2; k++) {
				ol[i].DISP_norm[j][k] = 0.0;
			}
		}
	}

	//The variable used to keep track of the nodecoord from MpCCI during the last time step
	for (m = 0; m < owsfnumber; m++) {
		if (ol[m].FSNEL > 0) {
			ol[m].nodecoord_mpcci = new double[3 * ol[m].GIDNct_st];
		}
	}

	for (m = 0; m < owsfnumber; m++) {
		if (ol[m].FSNEL > 0) {
			for (i = 0; i < 3 * ol[m].GIDNct_st; i++) {
				ol[m].nodecoord_mpcci[i] = ol[m].location;
			}
		}
	}

	//=======================================================================================//
	//Dynamic variables for fluid nodal force
	if (mappingalgo == 1 || mappingalgo == 3) {
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BF1 = new double*[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				ol[i].BF1[j] = new double[ol[i].FSNEL]; //BF is not evenly distributed in this case
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				for (k = 0; k < ol[i].FSNEL; k++) {
					ol[i].BF1[j][k] = 0.0;
				}
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BF3 = new double*[(hpref + 1)*(hpref + 1)];
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				ol[i].BF3[j] = new double[ol[i].FSNEL / refine / refine];
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				for (k = 0; k < ol[i].FSNEL / refine / refine; k++) {
					ol[i].BF3[j][k] = 0.0;
				}
			}
		}
	}
	else if (mappingalgo == 2) {
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BF2 = new double[NNODE];
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < NNODE; j++) {
				ol[i].BF2[j] = 0.0;
			}
		}
	}
	else if (mappingalgo == 4) {
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BP = new double*[(hpref + 1)*(hpref + 1)];
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				ol[i].BP[j] = new double[ol[i].FSNEL / refine / refine];
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				for (k = 0; k < ol[i].FSNEL / refine / refine; k++) {
					ol[i].BP[j][k] = 0.0;
				}
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BF1 = new double*[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				ol[i].BF1[j] = new double[ol[i].FSNEL]; //BF is not evenly distributed in this case
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				for (k = 0; k < ol[i].FSNEL; k++) {
					ol[i].BF1[j][k] = 0.0;
				}
			}
		}
	}
	else if (mappingalgo == 5) {
		for (i = 0; i < owsfnumber; i++) {
			ol[i].BPG = new double*[(hprefg + 1)*(hprefg + 1)];
			for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {
				ol[i].BPG[j] = new double[ol[i].FSNEL / refine / refine];
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {
				for (k = 0; k < ol[i].FSNEL / refine / refine; k++) {
					ol[i].BPG[j][k] = 0.0;
				}
			}
		}

		for (i = 0; i < owsfnumber; i++) {
			ol[i].BF1 = new double*[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				ol[i].BF1[j] = new double[ol[i].FSNEL]; //BF is not evenly distributed in this case
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				for (k = 0; k < ol[i].FSNEL; k++) {
					ol[i].BF1[j][k] = 0.0;
				}
			}
		}
	}


	//Might not be necessary
	//delete[] W;
	//delete[] S;

	//Initialize the incident pressure and total pressure along with hydrostatic pressure for the first time step
	TIME = 1;
	for (i = 0; i < NNODE; i++) {
		if (Bleich == 1) {
			PH[i] = (abs(GCOORD[i][1]) + (DRAFT - abs(GCOORD[0][1]))) * RHO * grav + PATM;
		}
		else {
			PH[i] = abs(GCOORD[i][1]) * RHO * grav + PATM; //absolute pressure
		}
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
		//spherical wave
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
		for (j = 0; j < nr[0].NRBNODE; j++) {
			for (k = 0; k < nrbsurfnumber; k++) {
				rd = pow((pow((GCOORD[nr[k].NRBA[j] - 1][0] - XC), 2) + pow((GCOORD[nr[k].NRBA[j] - 1][1] - YC), 2) + pow((GCOORD[nr[k].NRBA[j] - 1][2] - ZC), 2)), 0.5);
				ro = pow((pow((XC - XO), 2) + pow((YC - YO), 2) + pow((ZC - ZO), 2)), 0.5);
				//initialize to 0 
				if (-(rd - ro) / C >= 0) {
					//XNRBORG[j] = -(PPEAK*TAU*ro*1.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*1.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
					XNRBORG[nr[k].NRBA[j] - 1] = -(PPEAK*TAU*ro*1.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*1.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
				}
				else {
					//XNRBORG[j] = -(PPEAK*TAU*ro*0.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*0.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
					XNRBORG[nr[k].NRBA[j] - 1] = -(PPEAK*TAU*ro*0.0*(rd - ro + C*TAU - C*TAU*exp(rd / (C*TAU) - ro / (C*TAU)))) / (C*RHO*pow(rd, 2)) - (PPEAK*TAU*ro*0.0*(exp((rd - ro) / (C*TAU)) - 1)) / (C*RHO*rd);
				}
				// - (PPEAK*TAU*ro*heaviside(-(rd - ro)/C)*(rd - ro + C*TAU - C*TAU*exp(rd/(C*TAU) - ro/(C*TAU))))/(C*RHO*rd^2) - (PPEAK*TAU*ro*heaviside(-(rd - ro)/C)*(exp((rd - ro)/(C*TAU)) - 1))/(C*RHO*rd)
			}
		}
	}

	//=========================Get free surface points=========================//
	int fspt_num = 0;
	for (j = 0; j < NNODE; j++) {
		if (abs(GCOORD[j][1]) <= 1e-6) {
			fspt_num += 1;
		}
	}
	int *fspt;
	fspt = new int[fspt_num];
	for (j = 0; j < fspt_num; j++) {
		fspt[j] = 0;
	}
	int count = 0;
	for (j = 0; j < NNODE; j++) {
		if (abs(GCOORD[j][1]) <= 1e-6) {
			fspt[count] = j + 1;
			count += 1;
		}
	}
	if (count != fspt_num || count == 0) {
		std::cout << "the free surface point recognization is wrong or there's no free surface" << std::endl;
	}

	//Initiate the first coupling
	initcoupling();

	//Update the current MpCCI time
	current_time = T[0];

	//Obtain the local node numbering for the base fluid mesh
	TD_LOCAL_NODEstruct ctbase;
	ctbase = TD_LOCAL_NODE(hpref);
	TD_LOCAL_NODEstruct ctbase_algo5;
	ctbase_algo5 = TD_LOCAL_NODE(hprefg);

	//Perform the initial interface mapping
	LOBATTOstruct b;
	b = LOBATTO(hprefg);
	GLLQUADstruct f;
	f = GLLQUAD(b.Z, b.WL, hprefg, 0); //obtain Gauss-Legendre nodes (Not Gauss-Lobatto nodes)

	interface_mappingstruct in;
	in = interface_mapping(1, IEN, LNA_3D, LNA_2D, ctbase.LNA, ctbase_algo5.LNA, Z, TIME, GCOORD, phi_fem, f.W, phi_femg, phi_sem2);
	dotransfer();
	in = interface_mapping(0, IEN, LNA_3D, LNA_2D, ctbase.LNA, ctbase_algo5.LNA, Z, TIME, GCOORD, phi_fem, f.W, phi_femg, phi_sem2);

	//Get the sample points on a line to observe the wave propagation pressure distribution
	std::vector<int> sampline;
	count = 0;
	if (Bleich == 1) {
		for (i = 0; i < NNODE; i++) {
			if (abs(GCOORD[i][0] - 0.0) < 1e-6 && abs(GCOORD[i][2] - SZ / 2) < 1e-6) {
				sampline.push_back(i + 1);
				count += 1;
			}
		}
	}
	else {
		for (i = 0; i < NNODE; i++) {
			if (abs(GCOORD[i][0] - 0.0) < 1e-6 && abs(GCOORD[i][2] - 0.0) < 1e-6) {
				sampline.push_back(i + 1);
				count += 1;
			}
		}
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
		system("PAUSE ");
	}
	//get the points under the structural bottom corners
	int** samplinec;
	samplinec = new int*[4];
	for (i = 0; i < 4; i++) {
		samplinec[i] = new int[count];
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < count; j++) {
			samplinec[i][j] = 0;
		}
	}
	int ct1 = 0; int ct2 = 0; int ct3 = 0; int ct4 = 0;
	for (i = 0; i < NNODE; i++) {
		if (GCOORD[i][1] <= -SY && abs(GCOORD[i][0] + SX) < 1e-6 && abs(GCOORD[i][2] - SZ / 2) < 1e-6) {
			samplinec[0][ct1] = i + 1;
			ct1 += 1;
		}
		if (GCOORD[i][1] <= -SY && abs(GCOORD[i][0] + SX) < 1e-6 && abs(GCOORD[i][2] + SZ / 2) < 1e-6) {
			samplinec[1][ct2] = i + 1;
			ct2 += 1;
		}
		if (GCOORD[i][1] <= -SY && abs(GCOORD[i][0]) < 1e-6 && abs(GCOORD[i][2] - SZ / 2) < 1e-6) {
			samplinec[2][ct3] = i + 1;
			ct3 += 1;
		}
		if (GCOORD[i][1] <= -SY && abs(GCOORD[i][0]) < 1e-6 && abs(GCOORD[i][2] + SZ / 2) < 1e-6) {
			samplinec[3][ct4] = i + 1;
			ct4 += 1;
		}
	}
	//sort the arrary
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
	int count2 = 0; int count3 = 0; int count4 = 0; int count5 = 0; int count6 = 0;
	while (count2 != count || count3 != count || count4 != count || count5 != count || count6 != count) {
		for (i = 0; i < count; i++) {
			if (abs(GCOORD[sampline[i] - 1][1] + ypt[count2]) < 1e-5) {
				sampline2[count2] = sampline[i];
				count2 += 1;
			}
			if (abs(GCOORD[samplinec[0][i] - 1][1] + ypt[count3]) < 1e-5) {
				samplinec2[0][count3] = samplinec[0][i];
				count3 += 1;
			}
			if (abs(GCOORD[samplinec[1][i] - 1][1] + ypt[count4]) < 1e-5) {
				samplinec2[1][count4] = samplinec[1][i];
				count4 += 1;
			}
			if (abs(GCOORD[samplinec[2][i] - 1][1] + ypt[count5]) < 1e-5) {
				samplinec2[2][count5] = samplinec[2][i];
				count5 += 1;
			}
			if (abs(GCOORD[samplinec[3][i] - 1][1] + ypt[count6]) < 1e-5) {
				samplinec2[3][count6] = samplinec[3][i];
				count6 += 1;
			}
		}
	}

	
		//Generate the output files
		std::vector <int> T_out;
		int NDT_out = 0;
		for (i = 0; i < NDT; i++) {
			if (T[i] > output_int*NDT_out) {
				T_out.push_back(i);
				NDT_out += 1;
			}
		}
		/*
		std::string name = "PT_result_";
		//std::ofstream outfstr[Count];
		std::ofstream *outfstr;
		outfstr = new std::ofstream[NDT_out];
		std::string filename;
		if (output == 1) {
			for (i = 0; i < NDT_out; i++) {
				filename = name + std::to_string(T_out[i] * DT * 1000) + "ms " + timestr + ".txt";
				outfstr[i].open(filename);
			}
		}
		*/

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
		

	//Generate the information file
	std::string name3 = "parameters_" + timestr + ".txt";
	std::ofstream myfile;
	myfile.open(name3);
	//General information:
	myfile << "Code name: 3Dbarge_TFM_Abaqus_sym_mesh" << std::endl;
	myfile << "Simulation date and time: " << timestr << std::endl;
	myfile << "Mesh information:" << std::endl;
	myfile << "N: " << N << " mesh size: " << " XHE: " << XHE << " YHE: " << YHE << " ZHE: " << ZHE << " AX: " << AX << " BZ: " << BZ << " DY: " << DY << " SX: " << SX << " SY: " << SY << " SZ: " << SZ << " SXNEL: " << SXNEL << " SYNEL: " << SYNEL << " SZNEL: " << SZNEL << " wetted surface number: " << owsfnumber << " Nodenumber: " << NNODE << " Element number: " << NEL << " h/p refinement rate: " << refine << std::endl;
	myfile << "Explosive information: " << std::endl;
	myfile << "stdoff point (spherical): " << XO << " " << YO << " " << ZO << " explosion center: " << XC << " " << YC << " " << ZC << " stdoff (plane wave): " << xo << " peak pressure (Mpa): " << PPEAK << " decay rate: " << TAU << " WAVE (1: plane, 2: spherical): " << WAVE << " Abaquswaveform (1: smoothed, 0: sharp exponential decay): " << Abaquswaveform << std::endl;
	myfile << "Time integration information" << std::endl;
	myfile << "CFL: " << CFLFRAC << " dt: " << DT << " dt scale factor: " << dtscale << " total time: " << TTERM << " damping: " << BETA << " explicit central difference (Leap frog)" << std::endl;
	myfile << "Fluid properties: " << std::endl;
	myfile << "CAV: " << CAV << " C: " << C << " RHO: " << RHO << " atmospheric pressure: " << PATM << " saturated pressure: " << PSAT << std::endl;
	myfile << "FSI coupling" << std::endl;
	myfile << "mapping algorithm: " << mappingalgo << std::endl;
	myfile << "debug mode: " << debug << std::endl;

	double angle = 0.0; double R = 0.0;
	double norm[4][3];
	norm[0][0] = 0.0; norm[0][1] = -1.0; norm[0][2] = 0.0;  //bottom
	norm[1][0] = 0.0; norm[1][1] = 0.0; norm[1][2] = -1.0;  //back
	norm[2][0] = -1.0; norm[2][1] = 0.0; norm[2][2] = 0.0;  //left
	norm[3][0] = 0.0; norm[3][1] = 0.0; norm[3][2] = 1.0;   //front

	/*
	//generate a debug file
	std::string name4 = "pressure history at point 272572.txt";
	std::ofstream myfiledebug;
	myfiledebug.open(name4);
	*/

	//int*Wetnodes = new int[];
	int nx = 0; int ny = 0; double lcx = 0.0; double lcy = 0.0; double nomx = 0.0; double denomx = 0.0; double nomy = 0.0; double denomy = 0.0; 
	double phig = 0.0; 
	//double basep[2] = { -1,1 };
	double* basep; //points on base fluid mesh 
	basep = new double[NINT];
	if (FEM == 0) {
		for (i = 0; i < NINT; i++) {
			basep[i] = S[i];
		}
	}
	else {
		basep[0] = -1;
		basep[1] = 1; 
	}
	
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

	int** wetnorm;
	if (tfm == 0) {
		wetnorm = new int*[4];
		for (i = 0; i < 4; i++) {
			wetnorm[i] = new int[3];
		}
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 3; j++) {
				wetnorm[i][j] = 0;
			}
		}
		wetnorm[0][0] = 0; wetnorm[0][1] = 0; wetnorm[0][2] = 1;
		wetnorm[1][0] = 1; wetnorm[1][1] = 0; wetnorm[1][2] = 0;
		wetnorm[2][0] = 0; wetnorm[2][1] = 1; wetnorm[2][2] = 0;
		wetnorm[3][0] = 0; wetnorm[3][1] = 0; wetnorm[3][2] = -1;
	}

	std::string energyfile = "energy_history.txt";
	std::ofstream energyfilehd;
	energyfilehd.open(energyfile);

	//for (i = 0; i < NDT - 1; i++) { //error prone: i other than time should not present in this loop
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
		//int NNODE, double** GCOORD, double* T, int TIME, double** PIN, double DT, double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO
		for (j = 0; j < owsfnumber; j++) {
			//for (k = 0; k < NNODE; k++) {
			for (k = 0; k < ol[j].GIDNct; k++) {
				//ol[j].WP[k] = 0.0;
				ol[j].WBS[ol[j].GIDN[k] - 1] = 0.0;
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
				}
			}
		}
		else {
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < ol[z].GIDNct; j++) {
					ol[z].WP[ol[z].GIDN[j] - 1] = PT[ol[z].GIDN[j] - 1][0] - PATM; // //correct version with structural gravity 
					ol[z].WPIN[ol[z].GIDN[j] - 1] = PIN[ol[z].GIDN[j] - 1][1] + PIN[ol[z].GIDN[j] - 1][0];
				}
			}
		}

		//start = std::clock();
		//=======================define double* nodeforce in fluid code==========================//
		//mapping the fluid force ABF from user defined mesh to MpCCI defined mesh on coupling surface using interpolation
		in = interface_mapping(1, IEN, LNA_3D, LNA_2D, ctbase.LNA, ctbase_algo5.LNA, Z, TIME, GCOORD, phi_fem, f.W, phi_femg, phi_sem2);
		//after this subroutine, the nodeforce should already be mapped onto coupling surface (data.h)
		//int fluid2structure, int**IEN_3D, int***LNA_3D, int**LNA_2D, int NNODE, double *Z
	
		dotransfer();

		//=============map nodal displacement from coupled surface to fluid mesh===================//
		in = interface_mapping(0, IEN, LNA_3D, LNA_2D, ctbase.LNA, ctbase_algo5.LNA, Z, TIME, GCOORD, phi_fem, f.W, phi_femg, phi_sem2);
		
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
						r = sqrt(pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] - XC), 2) + pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - YC), 2) + pow((GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2] - ZC), 2));
						ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
							= (DT / 2.0)*(ol[z].WPIN[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]); //INCIDENT DISP. PREDICTOR by time integration	
						PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] = ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1];
						angle = (wetnorm[z][0] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] - XC) / r + wetnorm[z][1] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - YC) / r + wetnorm[z][2] * (GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][2] - ZC) / r);
						ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]
							= ol[z].PSI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] / (RHO*C) + (PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + PSI_inc[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]) * (DT / 2) / (RHO * r);  //trapezoidal integration of the double integrator
						ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]
							= ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0] + angle*ol[z].DI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1]; //UPDATED INCIDENT STRUCUTRE DISPLACEMENT 
						if (abs(GCOORD[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][ol[z].dir] - ol[z].location) > 1e-5) {
							std::cout << "initial displacement DISPI is given to wrong points on fluid surface: " << z << std::endl;
						}
					}
				}
			}
			//SET VALUES AT T(i+1) TO VALUES AT T(i)
			for (j = 0; j < NNODE; j++) {
				for (k = 0; k < owsfnumber; k++) {
					ol[k].DISPI[j][0] = ol[k].DISPI[j][1];  //SOLUTION ARRAY FOR INCIDENT STRUCUTRE DISPLACEMENT
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
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * RHO * (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0]));
						}
					}
					for (k = 0; k < NINT*NINT; k++) {
						ol[z].WBS[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] += WBSTEMP[k];
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
							WBSTEMP[h] += ol[z].FPMASTER[j][h][k] * RHO * (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] + (ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1] - ol[z].DISP_norm[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][0]) - ol[z].DISPI[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1][1]); //Error prone: Get the normal vector for DISPI??? 
						}
					}
					for (k = 0; k < NINT*NINT; k++) {
						ol[z].WBS[IEN[ol[z].FP[k] - 1][ol[z].GIDF[j] - 1] - 1] += WBSTEMP[k];
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
		if (tfm == 1) {
			for (k = 0; k < nrbsurfnumber; k++) { //from bottom to front surface (k=0 means all node case) 
				for (j = 0; j < nr[0].NRBNODE; j++) {
					FEEDOT_inc[nr[k].NRBA[j] - 1][1] = FEEDOT_inc[nr[k].NRBA[j] - 1][0] + (0.5 * DT)*(PIN[nr[k].NRBA[j] - 1][0] + PIN[nr[k].NRBA[j] - 1][1]);
					//prototype: FEEDOT_inc[u.NRBA[j][0] - 1][1] = (PPEAK*TAU - PPEAK*TAU*exp(pow((pow(XC, 2) - 2 * XC*d.GCOORD[u.NRBA[j][0] - 1][0] + pow(YC, 2) - 2 * YC*d.GCOORD[u.NRBA[j][0] - 1][1] + pow(ZC, 2) - 2 * ZC*d.GCOORD[u.NRBA[j][0] - 1][2] + pow(d.GCOORD[u.NRBA[j][0] - 1][0], 2) + pow(d.GCOORD[u.NRBA[j][0] - 1][1], 2) + pow(d.GCOORD[u.NRBA[j][0] - 1][2], 2)), 0.5) / (C*TAU))*exp(-(DT*i) / TAU)*exp(-pow((pow(XC, 2) - 2 * XC*XO + pow(XO, 2) + pow(YC, 2) - 2 * YC*YO + pow(YO, 2) + pow(ZC, 2) - 2 * ZC*ZO + pow(ZO, 2)), 0.5) / (C*TAU))*exp(-DT / TAU))*(sign((C*DT - pow((pow(XC, 2) - 2 * XC*d.GCOORD[u.NRBA[j][0] - 1][0] + pow(YC, 2) - 2 * YC*d.GCOORD[u.NRBA[j][0] - 1][1] + pow(ZC, 2) - 2 * ZC*d.GCOORD[u.NRBA[j][0] - 1][2] + pow(d.GCOORD[u.NRBA[j][0] - 1][0], 2) + pow(d.GCOORD[u.NRBA[j][0] - 1][1], 2) + pow(d.GCOORD[u.NRBA[j][0] - 1][2], 2)), 0.5) + pow((pow(XC, 2) - 2 * XC*XO + pow(XO, 2) + pow(YC, 2) - 2 * YC*YO + pow(YO, 2) + pow(ZC, 2) - 2 * ZC*ZO + pow(ZO, 2)), 0.5) + C*DT*i) / C) / 2.0 + 0.5);
				}
				count = 0;
				for (j = 0; j < nr[k].NRBNODE; j++) {
					R = pow(pow(GCOORD[nr[k].NRBA[j] - 1][0] - XC, 2) + pow(GCOORD[nr[k].NRBA[j] - 1][1] - YC, 2) + pow(GCOORD[nr[k].NRBA[j] - 1][2] - ZC, 2), 0.5);
					angle = (norm[k][0] * (GCOORD[nr[k].NRBA[j] - 1][0] - XC) + norm[k][1] * (GCOORD[nr[k].NRBA[j] - 1][1] - YC) + norm[k][2] * (GCOORD[nr[k].NRBA[j] - 1][2] - ZC)) / R;
					XEST_kn[nr[k].NRBA[j] - 1] = angle*((0.5 * DT / (RHO*C))*(PIN[nr[k].NRBA[j] - 1][0] + PIN[nr[k].NRBA[j] - 1][1]) +
						(0.5 * DT / (RHO*R))*(FEEDOT_inc[nr[k].NRBA[j] - 1][0] + FEEDOT_inc[nr[k].NRBA[j] - 1][1])) -
						(0.5 * DT / (RHO*C))*(PIN[nr[k].NRBA[j] - 1][0] + PIN[nr[k].NRBA[j] - 1][1]);
					XEST_ukn[nr[k].NRBA[j] - 1] = (DT / (RHO*C))*P[nr[k].NRBA[j] - 1][0]; //plane wave approximation (PWA)
					XNRBORG2[nr[k].NRBA[j] - 1] = angle*XNRBORG[nr[k].NRBA[j] - 1];
					count += 1;
				}
				for (j = 0; j < nr[k].NRBNODE; j++) {
					XNRB_kn[nr[k].NRBA[j] - 1][1] = XNRB_kn[nr[k].NRBA[j] - 1][0] + XEST_kn[nr[k].NRBA[j] - 1];
					XNRB_ukn[nr[k].NRBA[j] - 1][1] = XNRB_ukn[nr[k].NRBA[j] - 1][0] + XEST_ukn[nr[k].NRBA[j] - 1];
				}
				for (j = 0; j < nr[k].NEL_nrb; j++) {
					for (h = 0; h < NINT*NINT; h++) {
						BNRBTEMP[h] = 0.0;
						for (z = 0; z < NINT*NINT; z++) {
							BNRBTEMP[h] += nr[k].ADMASTER[j][h][z] * (-RHO) * (XNRB_kn[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1][1] + XNRB_ukn[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1][1] + XNRBORG2[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1]);
						}
					}
					for (z = 0; z < NINT*NINT; z++) {
						BNRB[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1] += BNRBTEMP[z];
					}
				}
			}
		}
		else {
			for (k = 0; k < nrbsurfnumber; k++) { //from bottom to front surface (k=0 means all node case)
				//NRB PREDICTOR
				for (j = 0; j < nr[0].NRBNODE; j++) {
					XEST[nr[k].NRBA[j] - 1] = (DT / (RHO*C))*P[nr[k].NRBA[j] - 1][0];  //SOLUTION ARRAY FOR NRB (DX)
				}
				for (j = 0; j < nr[0].NRBNODE; j++) {
					XNRB[nr[k].NRBA[j] - 1][1] = XNRB[nr[k].NRBA[j] - 1][0] + XEST[nr[k].NRBA[j] - 1]; //0 IS BEFORE MODIFICATION, 1 IS AFTER MODIFICATION
				}
				//XNRB is the predicted displacement normal to the structure (since P is normal to the surface)
				for (j = 0; j < nr[k].NEL_nrb; j++) {
					for (h = 0; h < NINT*NINT; h++) {
						BNRBTEMP[h] = 0.0;
						for (z = 0; z < NINT*NINT; z++) {
							BNRBTEMP[h] += nr[k].ADMASTER[j][h][z] * (-RHO) * (XNRB[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1][1] - XNRBORG[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1]);
						}
					}
					for (z = 0; z < NINT*NINT; z++) {
						BNRB[nr[k].IEN_gb[nr[k].DP[z] - 1][j] - 1] += BNRBTEMP[z];
					}
				}
			}
		}

		//time integration
		for (j = 0; j < NNODE; j++) {
			DSDOT[j] = (ds[j][1] - ds[j][0]) / DT; //SOLUTION ARRAY FOR FIRST TIME DERIVATIVE OF CONDENSATION
			FEEDOT[j][1] = FEEDOT[j][0] + DT*(P[j][0] + (BETA*DT*(pow(C, 2))*DSDOT[j]));
			FEE[j][1] = FEE[j][0] + DT*FEEDOT[j][1];
			FEEDOT[j][0] = FEEDOT[j][1];  //USE TWO VALUE TO STORE THE N TIME STEP AND N+1 TIME STEP
			HF[j] = 0.0;
			HFn[j] = 0.0;
		} //1 MEANS N+1 TIME STEP; 0 MEANS N TIME STEP

		//start = std::clock();
		ctt3 = 0;
		if (tensorfactorization == 0) {
			for (j = 0; j < NEL; j++) { //the loop takes 1+NEL+1+NEL=2*NEL+2 operations
				//oc += 2;
				for (k = 0; k < NINT*NINT*NINT; k++) { //the loop takes NEL*(2*NINT^3+2) operations
					//oc += 2;
					FEETEMP[k] = FEE[IEN[k][j] - 1][1]; //LOCAL DISP. POTENTIAL //takes NEL*NINT^3 operations
					oc += 1;
				} //Take the FEE for every element
				  //-------------------------matrix multiplication-------------------------//
				for (z = 0; z < NINT*NINT*NINT; z++) { //takes NEL*(2*NINT^3+2) operations
					//oc += 2;
					HFTEMP[z] = 0.0; //takes NEL*NINT^3 operations
					//oc += 1;
					for (k = 0; k < NINT*NINT*NINT; k++) { //takes NEL*NINT^3*(2*NINT^3+2) operaitons
						//oc += 2;
						HFTEMP[z] = HFTEMP[z] + HMASTER[j][z][k] * FEETEMP[k]; //elemental level reactance matrix //takes 3*NEL*NINT^6 operations
						oc += 2;
					} //The premise is that HMASTER is the same for every element (only true for equally shaped element!!)
				}
				for (k = 0; k < NINT*NINT*NINT; k++) { //takes NEL*(2*NINT^3+2) operations
					//oc += 2;
					HF[IEN[k][j] - 1] = HF[IEN[k][j] - 1] + HFTEMP[k];   //global level reactance matrix //takes 2*NINT^3*NEL operations
					oc += 1;
				} //HF is like global stiffness matrix (include variable in it)
			}
		}
		else {
			for (j = 0; j < NEL; j++) { //takes 2*NEL+2 operations 
				//takes 6*NEL*NINT^4 FLOP
				ctt1 = 0;
				for (ii = 0; ii < NINT; ii++) { //p NEL*(2*NINT+2) operations
					ctt2 = 0;
					for (h = 0; h < NINT; h++) { //j,i,i //takes NEL*NINT*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k,k,j //takes NEL*NINT^2*(2*NINT+2) operations
							gamman[counter1[ctt1]][0] = 0.0;
							gamman[counter2[ctt1]][1] = 0.0;
							gamman[counter3[ctt1]][2] = 0.0;
							for (q = 0; q < NINT; q++) { //q (recurrent addition in this dimension) //takes NEL*NINT^3*(2*NINT+2) operations	
								gamman[counter1[ctt1]][0] += SHOD1[q*NINT + ii] * FEE[IENct1[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1]; 
								gamman[counter2[ctt1]][1] += SHOD1[q*NINT + ii] * FEE[IENct2[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								gamman[counter3[ctt1]][2] += SHOD1[q*NINT + ii] * FEE[IENct3[j*NINT*NINT*NINT + h*NINT*NINT + z*NINT + q]][1];
								ctt2 += 1;
							}
							ctt1 += 1;
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
				ctt2 = 0;
				for (h = 0; h < NINT; h++) {  //i //takes NEL*NINT*(2*NINT+2) operations
					for (k = 0; k < NINT; k++) { //j  //takes NEL*NINT^2*(2*NINT+2) operations
						for (z = 0; z < NINT; z++) { //k //takes NEL*NINT^3*(2*NINT+2) operations
							HFTEMPn[LNAct3[ctt2]] = 0.0;
							for (ii = 0; ii < NINT; ii++) { //p (recurrent addition) //p NEL*(2*NINT+2) operations
								HFTEMPn[LNAct3[ctt2]] += SHOD1[h*NINT + ii] * gamma_tn[ii*NINT*NINT + k*NINT + z][0] + SHOD1[k*NINT + ii] * gamma_tn[h*NINT*NINT + ii*NINT + z][1] + SHOD1[z*NINT + ii] * gamma_tn[h*NINT*NINT + k*NINT + ii][2]; 
							}
							ctt2 += 1;
						}
					}
				}

				//takes 1*NINT^3*NEL FLOP
				for (h = 0; h < NINT*NINT*NINT; h++) { //takes NEL*(2 * NINT + 2) operations
					HF[IENct3[ctt3]] += HFTEMPn[LNAct3[h]]; //takes NEL*2*NINT^3 operations
					ctt3 += 1;
				}
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

		for (k = 0; k < nrbsurfnumber; k++) {
			for (j = 0; j < nr[0].NRBNODE; j++) {
				FFORCE[nr[k].NRBA[j] - 1] += BNRB[nr[k].NRBA[j] - 1];
			}
		}

		for (j = 0; j < fspt_num; j++) {
			Q[fspt[j] - 1] = 1.0;
			FFORCE[fspt[j] - 1] = 0.0 / pow(C, 2.0);
		}

		for (j = 0; j < NNODE; j++) {
			ds[j][2] = FFORCE[j] / Q[j];
		}

		//PRESSURE CORRECTION ON NRB NODES
		for (k = 0; k < nrbsurfnumber; k++) {
			for (j = 0; j < nr[0].NRBNODE; j++) { //UPDATE CONDENSATION AT TIME T=0
				KAPPA = (DT*C*nr[k].ADMASTERG[nr[k].NRBA[j] - 1]) / (2 * Q[nr[k].NRBA[j] - 1]);
				//PRESSURE CORRECTION FACTOR
				ds[nr[k].NRBA[j] - 1][2] = ds[nr[k].NRBA[j] - 1][1] + ((ds[nr[k].NRBA[j] - 1][2] - ds[nr[k].NRBA[j] - 1][1]) / (1 + KAPPA));
			} //after ds is updated, P can be updated.
		}

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

		if (tfm == 1) {
			for (k = 0; k < nrbsurfnumber; k++) {
				//XNRB corrector
				for (j = 0; j < nr[0].NRBNODE; j++) {
					DPS_ukn[nr[k].NRBA[j] - 1] = 0.5*DT*(P[nr[k].NRBA[j] - 1][1] + P[nr[k].NRBA[j] - 1][0]);
					XCOR_ukn[nr[k].NRBA[j] - 1] = DPS_ukn[nr[k].NRBA[j] - 1] / (RHO*C);
					XNRB_ukn[nr[k].NRBA[j] - 1][1] = XNRB_ukn[nr[k].NRBA[j] - 1][0] + XCOR_ukn[nr[k].NRBA[j] - 1];
				}

				for (j = 0; j < nr[0].NRBNODE; j++) {
					XNRB_kn[nr[k].NRBA[j] - 1][0] = XNRB_kn[nr[k].NRBA[j] - 1][1];
					XNRB_ukn[nr[k].NRBA[j] - 1][0] = XNRB_ukn[nr[k].NRBA[j] - 1][1];
				}
			}
		}
		else {
			for (k = 0; k < nrbsurfnumber; k++) {
				for (j = 0; j < nr[0].NRBNODE; j++) {
					DPS[nr[k].NRBA[j] - 1] = 0.5*DT*(P[nr[k].NRBA[j] - 1][1] + P[nr[k].NRBA[j] - 1][0]); // PRESSURE CORRECTOR FOR NRB
				}

				for (j = 0; j < nr[0].NRBNODE; j++) {
					XCOR[nr[k].NRBA[j] - 1] = DPS[nr[k].NRBA[j] - 1] / (RHO*C);  //DISPLACEMENT CORRECTOR
				}

				for (j = 0; j < nr[0].NRBNODE; j++) {
					XNRB[nr[k].NRBA[j] - 1][1] = XNRB[nr[k].NRBA[j] - 1][0] + XCOR[nr[k].NRBA[j] - 1]; //CORRECTED DISPLACEMENT ON NRB
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
			for (k = 0; k < nrbsurfnumber; k++) {
				for (j = 0; j < nr[0].NRBNODE; j++) {
					XNRB[nr[k].NRBA[j] - 1][0] = XNRB[nr[k].NRBA[j] - 1][1];
				}
			}
		}
		
		for (z = 0; z < owsfnumber; z++) {
			for (j = 0; j < ol[z].GIDNct; j++) {
				ol[z].DISP_norm[ol[z].GIDN[j] - 1][0] = ol[z].DISP_norm[ol[z].GIDN[j] - 1][1];
			}
			for (j = 0; j < 3 * ol[z].GIDNct_st; j++) {
				ol[z].nodecoord_mpcci[j] = wsflist[z]->nodecoord[j]; //pass the value of node displacement in all 3 directions
			}
		}

		if (output == 1) {
			for (k = 0; k < NDT_out; k++) {
				if (i == T_out[k]) {
					for (j = 0; j < count2; j++) {
						outline[k] << GCOORD[sampline2[j] - 1][1] << " " << PT[sampline2[j] - 1][1] << " " << ds[sampline2[j] - 1][2] << " " << PT[samplinec2[0][j] - 1][1] << " " << ds[samplinec2[0][j] - 1][2] << " " << PT[samplinec2[1][j] - 1][1] << " " << ds[samplinec2[1][j] - 1][2] << " " << PT[samplinec2[2][j] - 1][1] << " " << ds[samplinec2[2][j] - 1][2] << " " << PT[samplinec2[3][j] - 1][1] << " " << ds[samplinec2[3][j] - 1][2] << std::endl;
					}
				}
			}
		}
		//Output the pressure history under a specified point
		extern double BF_val[4];
		energyfilehd << current_time << " " << in.energy_sent << " " << in.energy_rec << " " << BF_val[0] << " " << BF_val[1] << " " << BF_val[2] << " " << BF_val[3] << " " << ol[0].OBF_val << " " << ol[1].OBF_val << " " << ol[2].OBF_val << " " << ol[3].OBF_val << std::endl;
	}
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	myfile << "total CPU time: " << duration << std::endl;
	return;
}