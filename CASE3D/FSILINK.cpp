//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>

//construct the link between structure and fluid module A C D separately
void FSILINK(double* W, int*** LNA, int**IEN, double***SHL, double**GCOORD, int NNODE, double***SHOD) {
	int i, j, k, l, m, n;
	int IC; //counter 
	extern int owsfnumber;
	extern OWETSURF ol[4]; //defined in FSILINK 

	//weighted integrator
	for (i = 0; i < owsfnumber; i++) {
		ol[i].SF = new double*[NINT*NINT*NINT];
		for (j = 0; j < NINT*NINT*NINT; j++) {
			ol[i].SF[j] = new double[ol[i].FSNEL];
		}
	}
	//Initialize SF here!!!!!!!!!!!
	for (i = 0; i < owsfnumber; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			for (k = 0; k < ol[i].FSNEL; k++) {
				ol[i].SF[j][k] = 0.0; //not all element in SF is given value
			}
		}
	}

	//Define the normal direction of surfaces 
	//If the direction pointing inward structure is identical with positive axis, then it is positive. 
	ol[0].NORM = 1.0; //surface 1
	ol[1].NORM = 1.0; //surface 2
	ol[2].NORM = 1.0; //surface 3
	ol[3].NORM = -1.0; //surface 4

	//============================determine secondary variable weights SF on coupling boundary===================//
	//the definition of SF is changed in order to map from denser fluid mesh to looser coupling mesh in interface_mapping subroutine
	LOBATTOstruct bq;
	bq = LOBATTO(Nq); //one unit higher than the interpolation order
	GLLQUADstruct gq;
	gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
	LOCAL_SHAPEstruct ls;
	ls = LOCAL_SHAPE(LNA, Nq); //one order higher than the interpolation order
	double FUNC;
	ol[0].FP = new int[NINT*NINT];
	for (i = 0; i < NINT*NINT; i++) {
		ol[0].FP[i] = 0;
	}
	int count;
	count = 0;
	for (m = 0; m < NINT; m++) {
		for (l = 0; l < NINT; l++) {
			ol[0].FP[count] = LNA[m][l][N]; //collect the surface node local numbering in one array
			//ol[0].FP[count] = LNA[l][m][N];
			count += 1;
		}
	}
	if (count > NINT*NINT) {
		std::cout << "FP out of bound" << std::endl;
		system("PAUSE ");
	}

	for (n = 0; n < ol[0].FSNEL; n++) { //same SF for all the FSNEL element 
		for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
			for (l = 0; l < NINT*NINT; l++) {
				FUNC = 0.0;
				if (m == l) {  //collect the diagonal term
					for (j = 0; j < NINT; j++) { //j and k are integration points in two dimensions 
						for (k = 0; k < NINT; k++) {
							FUNC += W[j] * W[k] * (SHL[3][ol[0].FP[m] - 1][j*NINT*NINT + k*NINT + N] / SHOD[0][N][N]) * (SHL[3][ol[0].FP[l] - 1][j*NINT*NINT + k*NINT + N] / SHOD[0][N][N]) * (XHE / 2.0)*(YHE / 2.0);
							//FUNC += W[j] * W[k] * SHL[3][ol[0].FP[m] - 1][j][k][N] * SHL[3][ol[0].FP[l] - 1][j][k][N] * (1 / 2.0)*(1 / 2.0);
						}
					}
					ol[0].SF[ol[0].FP[m] - 1][n] = FUNC;
					if (ol[0].FP[m] > NINT*NINT*NINT) {
						std::cout << "SF out of bound" << std::endl;
						system("PAUSE ");
					}
				}
			}
			//std::cout << ol[0].SF[ol[0].FP[m] - 1][n] << std::endl;
		}
		//std::cout << " " << std::endl;
	}
	//std::cout << " " << std::endl;
	//SF should be the same for all elements on a specific coupling surface. 

	for (i = 0; i < owsfnumber; i++) {
		ol[i].FPMASTER = new double*[NINT*NINT];
		for (j = 0; j < NINT*NINT; j++) {
			ol[i].FPMASTER[j] = new double[NINT*NINT];
		}
		for (j = 0; j < NINT*NINT; j++) {
			for (k = 0; k < NINT*NINT; k++) {
				ol[i].FPMASTER[j][k] = 0.0;
			}
		}
	}

	//surface 1 and 4 (the opposite surfaces have the same ol[0].FPMASTER)
	for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
		for (l = 0; l < NINT*NINT; l++) {
			FUNC = 0.0;
			for (j = 0; j < NqINT; j++) { //j and k are for integration point (one unit higher than the interpolation point)
				for (k = 0; k < NqINT; k++) {
					FUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[0].FP[m] - 1][j*NqINT*NqINT + k*NqINT + Nq] / ls.SHOD[0][N][Nq]) * (ls.SHL[3][ol[0].FP[l] - 1][j*NqINT*NqINT + k*NqINT + Nq] / ls.SHOD[0][N][Nq]) * (XHE / 2.0)*(YHE / 2.0);
					std::cout << " " << std::endl;
				}
			}
			ol[0].FPMASTER[m][l] += FUNC;
			//ol[3].FPMASTER[m][l] += FUNC;
		}
	}
	std::cout << " " << std::endl;
	

	ol[1].FP = new int[NINT*NINT];
	for (i = 0; i < NINT*NINT; i++) {
		ol[1].FP[i] = 0;
	}
	count = 0;
	for (m = 0; m < NINT; m++) {
		for (l = 0; l < NINT; l++) {
			ol[1].FP[count] = LNA[N][m][l];
			//ol[1].FP[count] = LNA[N][l][m]; //FP=[2,3,6,7]
			count += 1;
		}
	}
	if (count > NINT*NINT) {
		std::cout << "FP out of bound" << std::endl;
		system("PAUSE ");
	}
	//ATTENTION: For surface 2 and surface 5, "l" should come before "m" (m represents x axis and l represents y axis)
	//since later on phi_fem_alias use this sequence. 

	for (n = 0; n < ol[1].FSNEL; n++) {
		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				FUNC = 0.0;
				if (m == l) {
					for (j = 0; j < NINT; j++) {
						for (k = 0; k < NINT; k++) {
							FUNC += W[j] * W[k] * (SHL[3][ol[1].FP[m] - 1][N*NINT*NINT + j*NINT + k] / SHOD[0][N][N]) * (SHL[3][ol[1].FP[l] - 1][N*NINT*NINT + j*NINT + k] / SHOD[0][N][N]) * (YHE / 2.0)*(ZHE / 2.0);
						}
					}
					ol[1].SF[ol[1].FP[m] - 1][n] = FUNC;
					if (ol[1].FP[m] > NINT*NINT*NINT) {
						std::cout << "SF out of bound" << std::endl;
						system("PAUSE ");
					}
				}
			}
		}
	}

	for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
		for (l = 0; l < NINT*NINT; l++) {
			FUNC = 0.0;
			for (j = 0; j < NqINT; j++) { //j and k are for integration point (one unit higher than the interpolation point)
				for (k = 0; k < NqINT; k++) {
					FUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[1].FP[m] - 1][Nq*NqINT*NqINT + j*NqINT + k] / ls.SHOD[0][N][Nq]) * (ls.SHL[3][ol[1].FP[l] - 1][Nq*NqINT*NqINT + j*NqINT + k] / ls.SHOD[0][N][Nq]) * (ZHE / 2.0)*(YHE / 2.0);
				}
			}
			ol[1].FPMASTER[m][l] += FUNC;
		}
	}

	ol[2].FP = new int[NINT*NINT];
	for (i = 0; i < NINT*NINT; i++) {
		ol[2].FP[i] = 0;
	}
	count = 0;
	for (m = 0; m < NINT; m++) {
		for (l = 0; l < NINT; l++) {
			ol[2].FP[count] = LNA[m][N][l];
			//ol[2].FP[count] = LNA[l][N][m];
			count += 1;
		}
	}
	if (count > NINT*NINT) {
		std::cout << "FP out of bound" << std::endl;
		system("PAUSE ");
	}

	for (n = 0; n < ol[2].FSNEL; n++) {
		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				FUNC = 0.0;
				if (m == l) {
					for (j = 0; j < NINT; j++) {
						for (k = 0; k < NINT; k++) {
							FUNC += W[j] * W[k] * (SHL[3][ol[2].FP[m] - 1][j*NINT*NINT + N*NINT + k] / SHOD[0][N][N]) * (SHL[3][ol[2].FP[l] - 1][j*NINT*NINT + N*NINT + k] / SHOD[0][N][N]) * (XHE / 2.0)*(ZHE / 2.0);
						}
					}
					ol[2].SF[ol[2].FP[m] - 1][n] = FUNC;
					if (ol[2].FP[m] > NINT*NINT*NINT) {
						std::cout << "SF out of bound" << std::endl;
						system("PAUSE ");
					}
				}
			}
		}
	}
	
	for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
		for (l = 0; l < NINT*NINT; l++) {
			FUNC = 0.0;
			for (j = 0; j < NqINT; j++) { //j and k are for integration point (one unit higher than the interpolation point)
				for (k = 0; k < NqINT; k++) {
					FUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[2].FP[m] - 1][j*NqINT*NqINT + Nq*NqINT + k] / ls.SHOD[0][N][Nq]) * (ls.SHL[3][ol[2].FP[l] - 1][j*NqINT*NqINT + Nq*NqINT + k] / ls.SHOD[0][N][Nq]) * (XHE / 2.0)*(ZHE / 2.0);
				}
			}
			ol[2].FPMASTER[m][l] += FUNC;
		}
	}
	std::cout << " " << std::endl;

	ol[3].FP = new int[NINT*NINT];
	for (i = 0; i < NINT*NINT; i++) {
		ol[3].FP[i] = 0;
	}
	count = 0;
	for (m = 0; m < NINT; m++) {
		for (l = 0; l < NINT; l++) {
			ol[3].FP[count] = LNA[m][l][0]; //collect the surface node local numbering
			count += 1;
		}
	}
	if (count > NINT*NINT) {
		std::cout << "FP out of bound" << std::endl;
		system("PAUSE ");
	}

	for (n = 0; n < ol[3].FSNEL; n++) {
		for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
			for (l = 0; l < NINT*NINT; l++) {
				FUNC = 0.0;
				if (m == l) {  //collect the diagonal term
					for (j = 0; j < NINT; j++) { //j and k are for integration point
						for (k = 0; k < NINT; k++) {
							FUNC += W[j] * W[k] * (SHL[3][ol[3].FP[m] - 1][j*NINT*NINT + k*NINT + 0] / SHOD[0][0][0]) * (SHL[3][ol[3].FP[l] - 1][j*NINT*NINT + k*NINT + 0] / SHOD[0][0][0]) * (XHE / 2.0)*(YHE / 2.0);
						}
					}
					ol[3].SF[ol[3].FP[m] - 1][n] = FUNC;
					if (ol[3].FP[m] > NINT*NINT*NINT) {
						std::cout << "SF out of bound" << std::endl;
						system("PAUSE ");
					}
				}
			}
		}
	}

	for (m = 0; m < NINT*NINT; m++) { //m and l are for interpolation point
		for (l = 0; l < NINT*NINT; l++) {
			FUNC = 0.0;
			for (j = 0; j < NqINT; j++) { //j and k are for integration point (one unit higher than the interpolation point)
				for (k = 0; k < NqINT; k++) {
					FUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[3].FP[m] - 1][j*NqINT*NqINT+k*NqINT+0] / ls.SHOD[0][0][0]) * (ls.SHL[3][ol[3].FP[l] - 1][j*NqINT*NqINT + k*NqINT + 0] / ls.SHOD[0][0][0]) * (XHE / 2.0)*(YHE / 2.0);
				}
			}
			ol[3].FPMASTER[m][l] += FUNC;
		}
	}

	TD_LOCAL_NODEstruct ct;
	ct = TD_LOCAL_NODE(hpref);

	//if (mappingalgo == 3 || (mappingalgo == 4 && refine > 1)) {
	//if (mappingalgo == 1 || mappingalgo==2 || mappingalgo == 3 || mappingalgo == 4) {
	//total wetted node container 	
		ol[0].GIDN = new int[(N*SXNEL + 1)*(N*SYNEL + 1)]; 
		ol[1].GIDN = new int[(N*SZNEL + 1)*(N*SYNEL + 1)];
		ol[2].GIDN = new int[(N*SXNEL + 1)*(N*SZNEL + 1)];
		ol[3].GIDN = new int[(N*SXNEL + 1)*(N*SYNEL + 1)];

		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (GCOORD[i][0] + SX > -1e-6 && GCOORD[i][1] + SY > -1e-6 && abs(GCOORD[i][2] + SZ / 2) < 1e-6) {
				ol[0].GIDN[count] = i + 1;
				count += 1;
			}
		}
		ol[0].GIDNct = count;
		if (AX > 0 && BZ > 0) {
			if (count != (N*SXNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "surface 1 point number is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (GCOORD[i][2] + SZ / 2 >= -1e-6 && GCOORD[i][2] - SZ / 2 <= 1e-6 && GCOORD[i][1] + SY >= -1e-6 && abs(GCOORD[i][0] + SX) < 1e-6) {
				ol[1].GIDN[count] = i + 1;
				count += 1;
			}
		}
		ol[1].GIDNct = count;
		if (AX > 0 && BZ > 0) {
			if (count != (N*SZNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "surface 2 point number is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (GCOORD[i][2] + SZ / 2 >= -1e-6 && GCOORD[i][2] - SZ / 2 <= 1e-6 && GCOORD[i][0] + SX >= -1e-6 && abs(GCOORD[i][1] + SY) < 1e-6) {
				ol[2].GIDN[count] = i + 1;
				count += 1;
			}
		}
		ol[2].GIDNct = count;
		if (count != (N*SXNEL + 1)*(N*SZNEL + 1)) {
			std::cout << "surface 3 point number is wrong" << std::endl;
			system("PAUSE ");
		}

		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (GCOORD[i][0] + SX >= -1e-6 && GCOORD[i][1] + SY >= -1e-6 && abs(GCOORD[i][2] - SZ / 2) < 1e-6) {
				ol[3].GIDN[count] = i + 1;
				count += 1;
			}
		}
		ol[3].GIDNct = count;

		if (AX > 0 && BZ > 0) {
			if (count != (N*SXNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "surface 4 point number is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		//TD_LOCAL_NODEstruct ct;
		//ct = TD_LOCAL_NODE(hpref);
		for (i = 0; i < owsfnumber; i++) {
			ol[i].IEN_base = new int*[(hpref + 1)*(hpref + 1)];
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				ol[i].IEN_base[j] = new int[ol[i].FSNEL / refine / refine];
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				for (k = 0; k < ol[i].FSNEL / refine / refine; k++) {
					ol[i].IEN_base[j][k] = 0;
				}
			}
		}
		if (ol[0].FSNEL > 0) {
			for (i = 0; i < SYNEL / refine; i++) { //row
				for (j = 0; j < SXNEL / refine; j++) { //column
					for (k = 0; k < hpref + 1; k++) {
						for (l = 0; l < hpref + 1; l++) {
							ol[0].IEN_base[ct.LNA[l][k] - 1][j + i*(SXNEL / refine)] = ol[0].GIDN[(i*hpref + (hpref - k))*(N*SXNEL + 1) + (j*hpref + l + 1) - 1];
						}
					}
				}
			}
			if (((SYNEL / refine - 1)*hpref + (hpref - 0))*(N*SXNEL + 1) + ((SXNEL / refine - 1)*hpref + hpref + 1) != (N*SXNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "ol[0].IEN_base is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		if (ol[1].FSNEL > 0) {
			for (i = 0; i < SYNEL / refine; i++) { //row
				for (j = 0; j < SZNEL / refine; j++) { //column
					for (k = 0; k < hpref + 1; k++) {
						for (l = 0; l < hpref + 1; l++) {
							ol[1].IEN_base[ct.LNA[l][k] - 1][j + i*(SZNEL / refine)] = ol[1].GIDN[(i*hpref + (hpref - k))*(N*SZNEL + 1) + (j*hpref + l + 1) - 1];
						}
					}
				}
			}
			if (((SYNEL / refine - 1)*hpref + (hpref - 0))*(N*SZNEL + 1) + ((SZNEL / refine - 1)*hpref + hpref + 1) != (N*SZNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "ol[1].IEN_base is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		if (ol[2].FSNEL > 0) {
			for (i = 0; i < SZNEL / refine; i++) { //row
				for (j = 0; j < SXNEL / refine; j++) { //column
					for (k = 0; k < hpref + 1; k++) {
						for (l = 0; l < hpref + 1; l++) {
							//row number: N*SZNEL - (i*hpref + (hpref - k))
							//column number: j*hpref + l + 1
							ol[2].IEN_base[ct.LNA[l][k] - 1][j + i*(SXNEL / refine)] = ol[2].GIDN[(N*SZNEL - (i*hpref + (hpref - k)))*(N*SXNEL + 1) + (j*hpref + l + 1) - 1];
						}
					}
				}
			}
			if ((N*SZNEL - (0 * hpref + (hpref - hpref)))*(N*SXNEL + 1) + ((SXNEL / refine - 1)*hpref + hpref + 1) != (N*SXNEL + 1)*(N*SZNEL + 1)) {
				std::cout << "ol[2].IEN_base is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		if (ol[3].FSNEL > 0) {
			for (i = 0; i < SYNEL / refine; i++) { //row
				for (j = 0; j < SXNEL / refine; j++) { //column
					for (k = 0; k < hpref + 1; k++) {
						for (l = 0; l < hpref + 1; l++) {
							ol[3].IEN_base[ct.LNA[l][k] - 1][j + i*(SXNEL / refine)] = ol[3].GIDN[(i*hpref + (hpref - k))*(N*SXNEL + 1) + (j*hpref + l + 1) - 1];
						}
					}
				}
			}
			if (((SYNEL / refine - 1)*hpref + (hpref - 0))*(N*SXNEL + 1) + ((SXNEL / refine - 1)*hpref + hpref + 1) != (N*SXNEL + 1)*(N*SYNEL + 1)) {
				std::cout << "ol[3].IEN_base is wrong" << std::endl;
				system("PAUSE ");
			}
			//ol[2].IEN_base[ct.LNA[l][k] - 1][j + i*SXNEL] = ol[2].GIDN[(N*SXNEL + 1)*(N*SZNEL + 1) - (i*(refine)+(refine - k))*(j*(refine)+l + 1) - 1]
		}

		//define connectivity matrix between coarse and fine mesh (used for h refinement)
		for (i = 0; i < owsfnumber; i++) {
			ol[i].eletran = new int*[refine*refine];
			for (j = 0; j < refine*refine; j++) {
				ol[i].eletran[j] = new int[ol[i].FSNEL / refine / refine];
			}
		}
		//The distribution of eletran is from left to right, from top to bottom
		//surface 1
		if (ol[0].FSNEL > 0) {
			count = 0;
			for (i = 0; i < SYNEL / refine; i++) {
				for (j = 0; j < SXNEL / refine; j++) {
					for (k = 0; k < refine; k++) {
						for (l = 0; l < refine; l++) {
							ol[0].eletran[l + k*refine][j + (SXNEL / refine)*i] = SXNEL*(refine*i + k) + refine*j + l + 1;
							count += 1;
						}
					}
				}
			}
			if (ol[0].eletran[refine*refine - 1][ol[0].FSNEL / refine / refine - 1] != ol[0].FSNEL) {
				std::cout << "eletran for surface 1 is wrong" << std::endl;
				system("PAUSE ");
			}
		}

		//surface 2
		if (ol[1].FSNEL > 0) {
			count = 0;
			for (i = 0; i < SYNEL / refine; i++) {
				for (j = 0; j < SZNEL / refine; j++) {
					for (k = 0; k < refine; k++) {
						for (l = 0; l < refine; l++) {
							ol[1].eletran[l + k*refine][j + (SZNEL / refine)*i] = SZNEL*(refine*i + k) + refine*j + l + 1;
							count += 1;
						}
					}
				}
			}
			if (ol[1].eletran[refine*refine - 1][ol[1].FSNEL / refine / refine - 1] != ol[1].FSNEL) {
				std::cout << "eletran for surface 2 is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		//surface 3
		if (ol[2].FSNEL > 0) {
			count = 0;
			for (i = 0; i < SZNEL / refine; i++) {
				for (j = 0; j < SXNEL / refine; j++) {
					for (k = 0; k < refine; k++) {
						for (l = 0; l < refine; l++) {
							ol[2].eletran[l + k*refine][j + (SXNEL / refine)*i] = SXNEL*(refine*i + k) + refine*j + l + 1;
							//ol[2].eletran[l + k*refine][j + (SXNEL / refine)*i] = SXNEL*(i + k) + refine*j + l + 1;
							count += 1;
						}
					}
				}
			}
			if (ol[2].eletran[refine*refine - 1][ol[2].FSNEL / refine / refine - 1] != ol[2].FSNEL) {
				std::cout << "eletran for surface 3 is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		//surface 4
		if (ol[3].FSNEL > 0) {
			count = 0;
			for (i = 0; i < SYNEL / refine; i++) {
				for (j = 0; j < SXNEL / refine; j++) {
					for (k = 0; k < refine; k++) {
						for (l = 0; l < refine; l++) {
							ol[3].eletran[l + k*refine][j + (SXNEL / refine)*i] = SXNEL*(refine*i + k) + refine*j + l + 1;
							count += 1;
						}
					}
				}
			}
			if (ol[3].eletran[refine*refine - 1][ol[3].FSNEL / refine / refine - 1] != ol[3].FSNEL) {
				std::cout << "eletran for surface 4 is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		//std::cout << " " << std::endl;
	//}

	//Formulate the model file which will be passed to MpCCI adapter
	TD_ELE_GENstruct a, b, c, d, e;
	std::ofstream myfile;
	myfile.open("model.txt");
	int BCXNEL = 0;
	int BCYNEL = 0;
	/*
	//surface 1 calculation
	if (mappingalgo == 1 || (mappingalgo == 4 && N > 1)) {
		BCXNEL = SXNEL;
		BCYNEL = SYNEL;
	}
	else if (mappingalgo == 2) {
		BCXNEL = SXNEL*N;
		BCYNEL = SYNEL*N;
	}
	else if (mappingalgo == 3 || (mappingalgo == 4 && refine > 1)) {
		BCXNEL = SXNEL / refine;
		BCYNEL = SYNEL / refine;
	}*/
	if (mappingalgo == 1 || mappingalgo == 3 || mappingalgo == 4 || mappingalgo == 5) {
		BCXNEL = SXNEL / refine;
		BCYNEL = SYNEL / refine;
	}
	else if (mappingalgo == 2) {
		BCXNEL = SXNEL*N;
		BCYNEL = SYNEL*N;
	}
	int BCNEL = BCXNEL*BCYNEL;
	int BCXNODE = BCXNEL + 1;
	int BCYNODE = BCYNEL + 1;
	int BCNODE = BCXNODE*BCYNODE;
	ol[0].GIDNct_st = BCNODE;
	a = TD_ELE_GEN(BCNEL, BCXNODE, BCNODE, BCXNEL, BCYNEL, SX, SY);
	//int NEL, int XNODE, int NNODE, int XNEL, int YNEL, double XL, double YL
	//build the relationship between the SEM mesh and FEM mesh (associate the FEM node numbering with SEM node numbering) 
	//get the node coordinate. 
	int ele1 = 0;
	int line = 0;
	if (mappingalgo == 2) {
		ol[0].IEN_gb = new int*[4];
		for (i = 0; i < 4; i++) {
			ol[0].IEN_gb[i] = new int[BCNEL];
		}
		for (i = 0; i < 4; i++) {
			for (j = 0; j < BCNEL; j++) {
				ol[0].IEN_gb[i][j] = 0;
			}
		}

		ol[0].XYHE_gb = new double*[2];
		for (i = 0; i < 2; i++) {
			ol[0].XYHE_gb[i] = new double[BCNEL];
		}
		for (i = 0; i < 2; i++) {
			for (j = 0; j < BCNEL; j++) {
				ol[0].XYHE_gb[i][j] = 0;
			}
		}

		for (i = 0; i < ol[0].FSNEL; i++) {
			for (j = 0; j < NINT - 1; j++) { //line
				for (k = 0; k < NINT - 1; k++) { //colomn
					if (k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line > BCNEL - 1) {
						std::cout << "out of bound!!!" << std::endl;
						system("PAUSE ");
					}
					ol[0].IEN_gb[0][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k][N - 1 - j][N] - 1][ol[0].GIDF[i] - 1];
					ol[0].IEN_gb[1][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k + 1][N - 1 - j][N] - 1][ol[0].GIDF[i] - 1];
					ol[0].IEN_gb[2][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k + 1][N - j][N] - 1][ol[0].GIDF[i] - 1];
					ol[0].IEN_gb[3][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k][N - j][N] - 1][ol[0].GIDF[i] - 1];
					ol[0].XYHE_gb[0][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = GCOORD[IEN[LNA[k + 1][N - 1 - j][N] - 1][ol[0].GIDF[i] - 1] - 1][0] - GCOORD[IEN[LNA[k][N - 1 - j][N] - 1][ol[0].GIDF[i] - 1] - 1][0];
					ol[0].XYHE_gb[1][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = GCOORD[IEN[LNA[k][N - j][N] - 1][ol[0].GIDF[i] - 1]][1] - GCOORD[IEN[LNA[k][N - 1 - j][N] - 1][ol[0].GIDF[i] - 1] - 1][1];
					ele1 += 1;
					if (ele1 % (BCXNEL*N) == 0) {
						line += 1;
					}
				}
			}
		} 
		if (ele1 != BCNEL && AX > 0 && BZ > 0) {
			std::cout << "the element count of surface1 is wrong" << std::endl;
			system("PAUSE ");
		}

		//change the definition of a.GCOORD to SEM coordinate
		for (i = 0; i < ele1; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 2; k++) {
					a.GCOORD[a.IEN[j][i] - 1][k] = GCOORD[ol[0].IEN_gb[j][i] - 1][k];
				}
			}
		}
		std::cout << " " << std::endl;
	}
	else {
		for (i = 0; i < BCNODE; i++) {
			a.GCOORD[i][0] += -SX;
		}
	}

	//if (BCNEL > 0) {
	if (BZ > 0) {
		//myfile << "EF wetsurface1 3 2" << std::endl;
		myfile << "EF wetsurface1 3 2" << std::endl;
		myfile << "NODES " << BCNODE << std::endl;
		for (i = 0; i < BCNODE; i++) {
			myfile << i << " " << a.GCOORD[i][0] << " " << a.GCOORD[i][1] << " " << -SZ / 2 << " " << std::endl;
		}
		//int BFSNEL = (SX / ol[0].XHE)*(SY / ol[0].YHE);
		myfile << "ELEMENTS " << BCNEL << std::endl;
		//output connectivity matrix
		for (i = 0; i < BCNEL; i++) {
			myfile << i;
			for (j = 0; j < 4; j++) {
				myfile << " " << a.IEN[j][i] - 1; //node numbering starts from 0 in model file
			}
			myfile << std::endl;
		}
	}
	//surface 2 calculation
	//model file output
	int ACXNEL;
	int ACYNEL;
	/*
	if (mappingalgo==1 || (mappingalgo == 4 && N > 1)) {
		ACXNEL = SZNEL;
		ACYNEL = SYNEL;
	}
	else if (mappingalgo==2) {
		ACXNEL = (SZNEL)*N;
		ACYNEL = SYNEL*N;
	}
	else if (mappingalgo == 3 || (mappingalgo == 4 && refine > 1)) {
		ACXNEL = SZNEL / refine;
		ACYNEL = SYNEL / refine;
	}
	*/
	if (mappingalgo == 1 || mappingalgo == 3 || mappingalgo == 4 || mappingalgo == 5) {
		ACXNEL = SZNEL / refine;
		ACYNEL = SYNEL / refine;
	}
	else if (mappingalgo == 2) {
		ACXNEL = SZNEL*N;
		ACYNEL = SYNEL*N;
	}
	int ACNEL = ACXNEL*ACYNEL;
	int ACXNODE = ACXNEL + 1;
	int ACYNODE = ACYNEL + 1;
	int ACNODE = ACXNODE*ACYNODE;
	ol[1].GIDNct_st = ACNODE;
	b = TD_ELE_GEN(ACNEL, ACXNODE, ACNODE, ACXNEL, ACYNEL, SZ, SY);
	ele1 = 0;
	line = 0;
	if (mappingalgo == 2) {
		ol[1].IEN_gb = new int*[4];
		for (i = 0; i < 4; i++) {
			ol[1].IEN_gb[i] = new int[ACNEL];
		}
		for (i = 0; i < 4; i++) {
			for (j = 0; j < ACNEL; j++) {
				ol[1].IEN_gb[i][j] = 0;
			}
		}
		ol[1].XYHE_gb = new double*[2];
		for (i = 0; i < 2; i++) {
			ol[1].XYHE_gb[i] = new double[ACNEL];
		}
		for (i = 0; i < 2; i++) {
			for (j = 0; j < ACNEL; j++) {
				ol[1].XYHE_gb[i][j] = 0;
			}
		}
		for (i = 0; i < ol[1].FSNEL; i++) {
			for (j = 0; j < NINT - 1; j++) { //line
				for (k = 0; k < NINT - 1; k++) { //colomn
					if (k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line > ACNEL - 1) {
						std::cout << "out of bound!!!" << std::endl;
						system("PAUSE ");
					}
					ol[1].IEN_gb[0][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = IEN[LNA[N][N - 1 - j][k] - 1][ol[1].GIDF[i] - 1];
					ol[1].IEN_gb[1][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = IEN[LNA[N][N - 1 - j][k + 1] - 1][ol[1].GIDF[i] - 1];
					ol[1].IEN_gb[2][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = IEN[LNA[N][N - j][k + 1] - 1][ol[1].GIDF[i] - 1];
					ol[1].IEN_gb[3][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = IEN[LNA[N][N - j][k] - 1][ol[1].GIDF[i] - 1];
					ol[1].XYHE_gb[0][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = GCOORD[IEN[LNA[N][N - 1 - j][k + 1] - 1][ol[1].GIDF[i] - 1] - 1][2] - GCOORD[IEN[LNA[N][N - 1 - j][k] - 1][ol[1].GIDF[i] - 1] - 1][2];
					ol[1].XYHE_gb[1][k + j*ACXNEL + (i - line*(SZNEL))*N + ACXNEL*N*line] = GCOORD[IEN[LNA[N][N - j][k] - 1][ol[1].GIDF[i] - 1] - 1][1] - GCOORD[IEN[LNA[N][N - 1 - j][k] - 1][ol[1].GIDF[i] - 1] - 1][1];
					ele1 += 1;
					if (ele1 % (ACXNEL*N) == 0) {
						line += 1;
					}
				}
			}
			//std::cout << " " << std::endl;
		}
		if (ele1 != ACNEL && AX > 0 && BZ > 0) {
			std::cout << "the element count of surface2 is wrong" << std::endl;
			system("PAUSE ");
		}

		//change the definition of a.GCOORD to SEM coordinate
		for (i = 0; i < ele1; i++) {
			for (j = 0; j < 4; j++) {
				b.GCOORD[b.IEN[j][i] - 1][0] = GCOORD[ol[1].IEN_gb[j][i] - 1][2];
				b.GCOORD[b.IEN[j][i] - 1][1] = GCOORD[ol[1].IEN_gb[j][i] - 1][1];
			}
		}
	}
	else {
		for (i = 0; i < ACNODE; i++) {
			b.GCOORD[i][0] += -SZ / 2;
		}
	}

	if (AX > 0) {
		//myfile << "EF wetsurface2 3 2" << std::endl;
		myfile << "EF wetsurface2 3 0" << std::endl;
		myfile << "NODES " << ACNODE << std::endl;
		for (i = 0; i < ACNODE; i++) {
			myfile << i << " " << -SX << " " << b.GCOORD[i][1] << " " << b.GCOORD[i][0] << " " << std::endl;
		}
		//int BFSNEL = (SX / ol[0].XHE)*(SY / ol[0].YHE);
		myfile << "ELEMENTS " << ACNEL << std::endl;
		//output connectivity matrix
		for (i = 0; i < ACNEL; i++) {
			myfile << i;
			for (j = 0; j < 4; j++) {
				myfile << " " << b.IEN[j][i] - 1; //node numbering starts from 0 in model file
			}
			myfile << std::endl;
		}
	}

	//surface 3 calculation
	int DCXNEL = 0;
	int DCYNEL = 0;
	/*
	if (mappingalgo==1 || (mappingalgo == 4 && N > 1)) {
		DCXNEL = SXNEL;
		DCYNEL = SZNEL;
	}
	else if (mappingalgo==2) {
		DCXNEL = SXNEL*N;
		DCYNEL = SZNEL*N;
	}
	else if (mappingalgo == 3 || (mappingalgo == 4 && refine > 1)) {
		DCXNEL = SXNEL / refine;
		DCYNEL = SZNEL / refine;
	}
	*/
	if (mappingalgo == 1 || mappingalgo == 3 || mappingalgo == 4 || mappingalgo == 5) {
		DCXNEL = SXNEL / refine;
		DCYNEL = SZNEL / refine;
	}
	else if (mappingalgo == 2) {
		DCXNEL = SXNEL*N;
		DCYNEL = SZNEL*N;
	}
	int DCNEL = DCXNEL*DCYNEL;
	int DCXNODE = DCXNEL + 1;
	int DCYNODE = DCYNEL + 1;
	int DCNODE = DCXNODE*DCYNODE;
	ol[2].GIDNct_st = DCNODE;
	ol[2].IEN_gb = new int*[4];
	for (i = 0; i < 4; i++) {
		ol[2].IEN_gb[i] = new int[DCNEL];
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < DCNEL; j++) {
			ol[2].IEN_gb[i][j] = 0;
		}
	}
	ol[2].XYHE_gb = new double*[2];
	for (i = 0; i < 2; i++) {
		ol[2].XYHE_gb[i] = new double[DCNEL];
	}
	for (i = 0; i < 2; i++) {
		for (j = 0; j < DCNEL; j++) {
			ol[2].XYHE_gb[i][j] = 0;
		}
	}
	c = TD_ELE_GEN(DCNEL, DCXNODE, DCNODE, DCXNEL, DCYNEL, SX, SZ);
	//TD_ELE_GEN(int NEL, int XNODE, int NNODE,int XNEL,int YNEL,double XL, double YL)
	//input value error prone! 
	ele1 = 0;
	line = 0;
	if (mappingalgo == 2) {
		for (i = 0; i < ol[2].FSNEL; i++) {
			for (j = 0; j < NINT - 1; j++) { //line
				for (k = 0; k < NINT - 1; k++) { //colomn
					if (k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line > DCNEL - 1) {
						std::cout << "out of bound!!!" << std::endl;
						system("PAUSE ");
					}
					ol[2].IEN_gb[0][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = IEN[LNA[k][N][N - 1 - j] - 1][ol[2].GIDF[i] - 1];
					ol[2].IEN_gb[1][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = IEN[LNA[k + 1][N][N - 1 - j] - 1][ol[2].GIDF[i] - 1];
					ol[2].IEN_gb[2][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = IEN[LNA[k + 1][N][N - j] - 1][ol[2].GIDF[i] - 1];
					ol[2].IEN_gb[3][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = IEN[LNA[k][N][N - j] - 1][ol[2].GIDF[i] - 1];
					ol[2].XYHE_gb[0][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = GCOORD[IEN[LNA[k + 1][N][N - 1 - j] - 1][ol[2].GIDF[i] - 1] - 1][0] - GCOORD[IEN[LNA[k][N][N - 1 - j] - 1][ol[2].GIDF[i] - 1] - 1][0];
					ol[2].XYHE_gb[1][k + j*DCXNEL + (i - line*(SXNEL))*N + DCXNEL*N*line] = GCOORD[IEN[LNA[k][N][N - j] - 1][ol[2].GIDF[i] - 1] - 1][2] - GCOORD[IEN[LNA[k][N][N - 1 - j] - 1][ol[2].GIDF[i] - 1] - 1][2];
					ele1 += 1;
					if (ele1 % (DCXNEL*N) == 0) {
						line += 1;
					}
				}
			}
			//std::cout << " " << std::endl;
		}
		if (ele1 != DCNEL) {
			std::cout << "the element count of surface3 is wrong" << std::endl;
			system("PAUSE ");
		}

		//change the definition of a.GCOORD to SEM coordinate
		for (i = 0; i < ele1; i++) {
			for (j = 0; j < 4; j++) {
				c.GCOORD[c.IEN[j][i] - 1][0] = GCOORD[ol[2].IEN_gb[j][i] - 1][0];
				c.GCOORD[c.IEN[j][i] - 1][1] = GCOORD[ol[2].IEN_gb[j][i] - 1][2];
			}
		}
	}
	else {
		for (i = 0; i < DCNODE; i++) {
			c.GCOORD[i][0] += -SX;
			c.GCOORD[i][1] += SZ / 2;
		}
	}
	if (SX > 0 && SZ>0) {
		//myfile << "EF wetsurface3 3 2" << std::endl;
		myfile << "EF wetsurface3 3 1" << std::endl;
		myfile << "NODES " << DCNODE << std::endl;
		for (i = 0; i < DCNODE; i++) {
			myfile << i << " " << c.GCOORD[i][0] << " " << -SY << " " << c.GCOORD[i][1] << " " << std::endl;
		}
		//int BFSNEL = (SX / ol[0].XHE)*(SY / ol[0].YHE);
		myfile << "ELEMENTS " << DCNEL << std::endl;
		//output connectivity matrix
		for (i = 0; i < DCNEL; i++) {
			myfile << i;
			for (j = 0; j < 4; j++) {
				myfile << " " << c.IEN[j][i] - 1; //node numbering starts from 0 in model file
			}
			myfile << std::endl;
		}
	}

	//surface 4 calculation (everything same with surface 1 except for coordinate)
	//a = TD_ELE_GEN(BCNEL, BCXNODE, BCNODE, BCXNEL, BCYNEL, BX, BY - DY);
	d = TD_ELE_GEN(BCNEL, BCXNODE, BCNODE, BCXNEL, BCYNEL, SX, SY);
	//int NEL, int XNODE, int NNODE, int XNEL, int YNEL, double XL, double YL
	//build the relationship between the SEM mesh and FEM mesh (associate the FEM node numbering with SEM node numbering) 
	//get the node coordinate. 
	ol[3].IEN_gb = new int*[4];
	for (i = 0; i < 4; i++) {
		ol[3].IEN_gb[i] = new int[BCNEL];
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < BCNEL; j++) {
			ol[3].IEN_gb[i][j] = 0;
		}
	}
	ol[3].XYHE_gb = new double*[2];
	for (i = 0; i < 2; i++) {
		ol[3].XYHE_gb[i] = new double[BCNEL];
	}
	for (i = 0; i < 2; i++) {
		for (j = 0; j < BCNEL; j++) {
			ol[3].XYHE_gb[i][j] = 0;
		}
	}
	ele1 = 0;
	line = 0;
	if (mappingalgo == 2) {
		for (i = 0; i < ol[3].FSNEL; i++) {
			for (j = 0; j < NINT - 1; j++) { //line
				for (k = 0; k < NINT - 1; k++) { //colomn
					if (k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line > BCNEL - 1) {
						std::cout << "out of bound!!!" << std::endl;
						system("PAUSE ");
					}
					ol[3].IEN_gb[0][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k][N - 1 - j][0] - 1][ol[3].GIDF[i] - 1];
					ol[3].IEN_gb[1][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k + 1][N - 1 - j][0] - 1][ol[3].GIDF[i] - 1];
					ol[3].IEN_gb[2][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k + 1][N - j][0] - 1][ol[3].GIDF[i] - 1];
					ol[3].IEN_gb[3][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = IEN[LNA[k][N - j][0] - 1][ol[3].GIDF[i] - 1];
					ol[3].XYHE_gb[0][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = GCOORD[IEN[LNA[k + 1][N - 1 - j][0] - 1][ol[3].GIDF[i] - 1] - 1][0] - GCOORD[IEN[LNA[k][N - 1 - j][0] - 1][ol[3].GIDF[i] - 1] - 1][0];
					ol[3].XYHE_gb[1][k + j*BCXNEL + (i - line*SXNEL)*N + BCXNEL*N*line] = GCOORD[IEN[LNA[k][N - j][0] - 1][ol[3].GIDF[i] - 1] - 1][1] - GCOORD[IEN[LNA[k][N - 1 - j][0] - 1][ol[3].GIDF[i] - 1] - 1][1];
					ele1 += 1;
					if (ele1 % (BCXNEL*N) == 0) {
						line += 1;
					}
				}
			}
		}
		if (ele1 != BCNEL && AX > 0 && BZ > 0) {
			std::cout << "the element count of surface4 is wrong" << std::endl;
			system("PAUSE ");
		}

		//change the definition of a.GCOORD to SEM coordinate
		for (i = 0; i < ele1; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 2; k++) {
					d.GCOORD[d.IEN[j][i] - 1][k] = GCOORD[ol[3].IEN_gb[j][i] - 1][k];
				}
			}
		}
		std::cout << " " << std::endl;
	}
	else {
		for (i = 0; i < BCNODE; i++) {
			d.GCOORD[i][0] += -SX;
		}
	}
	ol[3].GIDNct_st = BCNODE;
	if (BZ > 0) {
		//myfile << "EF wetsurface4 3 2" << std::endl;
		myfile << "EF wetsurface4 3 2" << std::endl;
		myfile << "NODES " << BCNODE << std::endl;
		for (i = 0; i < BCNODE; i++) {
			myfile << i << " " << d.GCOORD[i][0] << " " << d.GCOORD[i][1] << " " << SZ / 2 << " " << std::endl;
		}
		//int BFSNEL = (SX / ol[0].XHE)*(SY / ol[0].YHE);
		myfile << "ELEMENTS " << BCNEL << std::endl;
		//output connectivity matrix
		for (i = 0; i < BCNEL; i++) {
			myfile << i;
			for (j = 0; j < 4; j++) {
				myfile << " " << d.IEN[j][i] - 1; //node numbering starts from 0 in model file
			}
			myfile << std::endl;
		}
	}

	//Node weight on base fluid mesh for displacement mapping (algorithm 3 and 4) 
	if (mappingalgo == 3 || mappingalgo == 4) {
		for (i = 0; i < owsfnumber; i++) {
			ol[i].NW = new double*[(hpref + 1)*(hpref + 1)];
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				ol[i].NW[j] = new double[ol[i].FSNEL / refine / refine];
			}
		}
		for (i = 0; i < owsfnumber; i++) {
			for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
				for (k = 0; k < ol[i].FSNEL / refine / refine; k++) {
					ol[i].NW[j][k] = 0.0;
				}
			}
		}
		int* NWflag = new int[NNODE];
		for (m = 0; m < owsfnumber; m++) {
			for (k = 0; k < NNODE; k++) {
				NWflag[k] = 0;
			}
			for (i = 0; i < ol[m].FSNEL / refine / refine; i++) {
				for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
					NWflag[ol[m].IEN_base[j][i] - 1] += 1;
				}
			}
			for (i = 0; i < ol[m].FSNEL / refine / refine; i++) {
				for (j = 0; j < (hpref + 1)*(hpref + 1); j++) {
					ol[m].NW[j][i] = 1.0 / NWflag[ol[m].IEN_base[j][i] - 1];
				}
			}
		}
		//std::cout << " " << std::endl;

		//===========================end of the definition of SF==========================//
		delete[] NWflag;
	}

	for (i = 0; i < BCNODE; i++) {
		delete[] a.GCOORD[i];
	}
	delete[] a.GCOORD;
	for (i = 0; i < (NC + 1)*(NC + 1);i++) {
		delete[] a.IEN[i];
	}
	delete[] a.IEN;

	for (i = 0; i < BCNODE; i++) {
		delete[] d.GCOORD[i];
	}
	delete[] d.GCOORD;
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		delete[] d.IEN[i];
	}
	delete[] d.IEN;

	for (i = 0; i < ACNODE; i++) {
		delete[] b.GCOORD[i];
	}
	delete[] b.GCOORD;
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		delete[] b.IEN[i];
	}
	delete[] b.IEN;

	for (i = 0; i < DCNODE; i++) {
		delete[] c.GCOORD[i];
	}
	delete[] c.GCOORD;
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		delete[] c.IEN[i];
	}
	delete[] c.IEN;
	
	std::cout << " " << std::endl;
	return;
}
