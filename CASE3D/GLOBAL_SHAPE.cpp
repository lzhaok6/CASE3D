//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <ctime>

//GLOBAL_SHAPE determines the global shape function derivative at the Gauss quad. points.

struct GLOBAL_SHAPEstruct GLOBAL_SHAPE(int NEL, double***SHL, double**GCOORD, int**IEN, double* JACOB_tet, int*** LNA) {
	int i, j, k, l, m, n;
	GLOBAL_SHAPEstruct t;
	t.dummy = 1; 

	//Since storing SHG is too expansive, we decide to define SHG on the fly in the function that needs it. 
	/*
	std::clock_t start;
	start = std::clock();
	if (element_type == 0) {
		t.SHG = new double***[NEL];
		for (i = 0; i < NEL; i++) {
			t.SHG[i] = new double**[4];
			for (j = 0; j < 4; j++) {
				t.SHG[i][j] = new double*[NINT*NINT*NINT];
				for (k = 0; k < NINT*NINT*NINT; k++) {
					t.SHG[i][j][k] = new double[NINT*NINT*NINT];
				}
			}
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < NINT*NINT*NINT; k++) {
					for (l = 0; l < NINT*NINT*NINT; l++) {
						t.SHG[i][j][k][l] = 0.0;
					}
				}
			}
		}

		double* JACOB = new double[NINT*NINT*NINT]; 
		double** XS; //define XS to be a local variable which will be cleaned at the end of function
		XS = new double*[9];
		for (i = 0; i < 9; i++) {
			XS[i] = new double[NINT*NINT*NINT];
		}
		for (i = 0; i < 9; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				XS[i][j] = 0.0;
			}
		}
		LOBATTOstruct bq;
		GLLQUADstruct gq;
		LOCAL_GSHAPEstruct lg;
		bq = LOBATTO(N);
		gq = GLLQUAD(bq.Z, bq.WL, N, !FEM);
		lg = LOCAL_GSHAPE(gq.S, LNA, NINT);
		int*cn = new int[8];
		cn[0] = LNA[0][0][0]; cn[1] = LNA[N][0][0];
		cn[2] = LNA[N][N][0]; cn[3] = LNA[0][N][0];
		cn[4] = LNA[0][0][N]; cn[5] = LNA[N][0][N];
		cn[6] = LNA[N][N][N]; cn[7] = LNA[0][N][N];
		for (m = 0; m < NEL; m++) {
			//derive XS first
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = 0;
						for (l = 0; l < 8; l++) {
							XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
							XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta 
							XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta  
							XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi 
							XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
							XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta 
							XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi 
							XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta 
							XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
						}
						JACOB[i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k])
							- XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k])
							+ XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k]);
					}
				}
			}
			for (i = 0; i < NINT*NINT*NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						for (l = 0; l < NINT; l++) {
							t.SHG[m][0][i][j*NINT*NINT + k*NINT + l] = ((XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dx
							t.SHG[m][1][i][j*NINT*NINT + k*NINT + l] = ((-XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dy
							t.SHG[m][2][i][j*NINT*NINT + k*NINT + l] = ((XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] - XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dz
							t.SHG[m][3][i][j*NINT*NINT + k*NINT + l] = SHL[3][i][j*NINT*NINT + k*NINT + l]; //phi=phi 
						}
					}
				}
			}
		}
		for (i = 0; i < 9; i++) {
			delete[] XS[i];
		}
		delete[] XS; 
		delete[] JACOB; 
	}
	*/

	if (element_type == 1) {
		t.SHG_tet = new double**[NEL];
		for (i = 0; i < NEL; i++) {
			t.SHG_tet[i] = new double*[3]; // x y z
			for (j = 0; j < 3; j++) {
				t.SHG_tet[i][j] = new double[4]; //4 nodes
			}
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 4; k++) {
					t.SHG_tet[i][j][k] = 0.0;
				}
			}
		}

		for (i = 0; i < NEL; i++) {
			double x21 = GCOORD[IEN[1][i] - 1][0] - GCOORD[IEN[0][i] - 1][0];
			double y21 = GCOORD[IEN[1][i] - 1][1] - GCOORD[IEN[0][i] - 1][1];
			double z21 = GCOORD[IEN[1][i] - 1][2] - GCOORD[IEN[0][i] - 1][2];
			double x31 = GCOORD[IEN[2][i] - 1][0] - GCOORD[IEN[0][i] - 1][0];
			double y31 = GCOORD[IEN[2][i] - 1][1] - GCOORD[IEN[0][i] - 1][1];
			double z31 = GCOORD[IEN[2][i] - 1][2] - GCOORD[IEN[0][i] - 1][2];
			double x41 = GCOORD[IEN[3][i] - 1][0] - GCOORD[IEN[0][i] - 1][0];
			double y41 = GCOORD[IEN[3][i] - 1][1] - GCOORD[IEN[0][i] - 1][1];
			double z41 = GCOORD[IEN[3][i] - 1][2] - GCOORD[IEN[0][i] - 1][2];
			double x32 = GCOORD[IEN[2][i] - 1][0] - GCOORD[IEN[1][i] - 1][0];
			double y32 = GCOORD[IEN[2][i] - 1][1] - GCOORD[IEN[1][i] - 1][1];
			double z32 = GCOORD[IEN[2][i] - 1][2] - GCOORD[IEN[1][i] - 1][2];
			double x42 = GCOORD[IEN[3][i] - 1][0] - GCOORD[IEN[1][i] - 1][0];
			double y42 = GCOORD[IEN[3][i] - 1][1] - GCOORD[IEN[1][i] - 1][1];
			double z42 = GCOORD[IEN[3][i] - 1][2] - GCOORD[IEN[1][i] - 1][2];
			//First node
			t.SHG_tet[i][0][0] = (-y32*z42 + z32*y42) / JACOB_tet[i]; //x 
			t.SHG_tet[i][1][0] = (x32*z42 - z32*x42) / JACOB_tet[i]; //y 
			t.SHG_tet[i][2][0] = (-x32*y42 + y32*x42) / JACOB_tet[i]; //z
			//Second node
			t.SHG_tet[i][0][1] = (y31*z41 - z31*y41) / JACOB_tet[i];
			t.SHG_tet[i][1][1] = (-x31*z41 + z31*x41) / JACOB_tet[i];
			t.SHG_tet[i][2][1] = (x31*y41 - y31*x41) / JACOB_tet[i];
			//Third node
			t.SHG_tet[i][0][2] = (-y21*z41 + z21*y41) / JACOB_tet[i];
			t.SHG_tet[i][1][2] = (x21*z41 - z21*x41) / JACOB_tet[i];
			t.SHG_tet[i][2][2] = (-x21*y41 + y21*x41) / JACOB_tet[i];
			//Fourth node
			t.SHG_tet[i][0][3] = (y21*z31 - z21*y31) / JACOB_tet[i];
			t.SHG_tet[i][1][3] = (-x21*z31 + z21*x31) / JACOB_tet[i];
			t.SHG_tet[i][2][3] = (x21*y31 - y21*x31) / JACOB_tet[i];
		}

	}

	return t;

}
