//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>

struct JACOBIANstruct JACOBIAN(int NEL, double*****GSHL, double****GSHL_2D, double **GCOORD, int **IEN, int*** LNA) {
	int i, j, k, l, m, n;
	JACOBIANstruct t;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	//initialize t.XS
	t.XS = new double***[NEL];
	for (i = 0; i < NEL; i++) {
		t.XS[i] = new double**[3];
		for (j = 0; j < 3; j++) {
			t.XS[i][j] = new double*[3];
			for (k = 0; k < 3; k++) {
				t.XS[i][j][k] = new double[NINT*NINT*NINT];
			}
		}
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < NINT*NINT*NINT; l++) {
					t.XS[i][j][k][l] = 0.0;
				}
			}
		}
	}

	t.XS_2D = new double***[ol[0].NEL];
	for (i = 0; i < ol[0].NEL; i++) {
		t.XS_2D[i] = new double**[2];
		for (j = 0; j < 2; j++) {
			t.XS_2D[i][j] = new double*[2];
			for (k = 0; k < 2; k++) {
				t.XS_2D[i][j][k] = new double[NINT*NINT];
			}
		}
	}

	t.JACOB = new double*[NEL];
	for (i = 0; i < NEL; i++) {
		t.JACOB[i] = new double[NINT*NINT*NINT];
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			t.JACOB[i][j] = 0.0;
		}
	}

	//Jacobian matrix is evaluated on integration points
	//XS is the Jacobian matrix
	//corner nodes arrary

	int*cn = new int[8];
	cn[0] = LNA[0][0][0];
	cn[1] = LNA[N][0][0];
	cn[2] = LNA[N][N][0];
	cn[3] = LNA[0][N][0];
	cn[4] = LNA[0][0][N];
	cn[5] = LNA[N][0][N];
	cn[6] = LNA[N][N][N];
	cn[7] = LNA[0][N][N];

	for (m = 0; m < NEL; m++) {
		std::cout << m << std::endl;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < 8; l++) {
						//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
						t.XS[m][0][0][i*NINT*NINT + j*NINT + k] = t.XS[m][0][0][i*NINT*NINT + j*NINT + k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
						t.XS[m][1][0][i*NINT*NINT + j*NINT + k] = t.XS[m][1][0][i*NINT*NINT + j*NINT + k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi
						t.XS[m][2][0][i*NINT*NINT + j*NINT + k] = t.XS[m][2][0][i*NINT*NINT + j*NINT + k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi
						t.XS[m][0][1][i*NINT*NINT + j*NINT + k] = t.XS[m][0][1][i*NINT*NINT + j*NINT + k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta
						t.XS[m][1][1][i*NINT*NINT + j*NINT + k] = t.XS[m][1][1][i*NINT*NINT + j*NINT + k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
						t.XS[m][2][1][i*NINT*NINT + j*NINT + k] = t.XS[m][2][1][i*NINT*NINT + j*NINT + k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta
						t.XS[m][0][2][i*NINT*NINT + j*NINT + k] = t.XS[m][0][2][i*NINT*NINT + j*NINT + k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta
						t.XS[m][1][2][i*NINT*NINT + j*NINT + k] = t.XS[m][1][2][i*NINT*NINT + j*NINT + k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta
						t.XS[m][2][2][i*NINT*NINT + j*NINT + k] = t.XS[m][2][2][i*NINT*NINT + j*NINT + k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
					}
					//The determinant of jacobian matrix
					t.JACOB[m][i*NINT*NINT + j*NINT + k] = t.XS[m][0][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][1][1][i*NINT*NINT + j*NINT + k] * t.XS[m][2][2][i*NINT*NINT + j*NINT + k] - t.XS[m][2][1][i*NINT*NINT + j*NINT + k] * t.XS[m][1][2][i*NINT*NINT + j*NINT + k])
						- t.XS[m][1][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][0][1][i*NINT*NINT + j*NINT + k] * t.XS[m][2][2][i*NINT*NINT + j*NINT + k] - t.XS[m][2][1][i*NINT*NINT + j*NINT + k] * t.XS[m][0][2][i*NINT*NINT + j*NINT + k])
						+ t.XS[m][2][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][0][1][i*NINT*NINT + j*NINT + k] * t.XS[m][1][2][i*NINT*NINT + j*NINT + k] - t.XS[m][1][1][i*NINT*NINT + j*NINT + k] * t.XS[m][0][2][i*NINT*NINT + j*NINT + k]);
				}
			}
		}
	}
	//Initialize Jacob_2D!!!!!!!!!!!!!!!!!!!1
	ol[0].Jacob_2D = new double*[ol[0].NEL];
	for (i = 0; i < ol[0].NEL; i++) {
		ol[0].Jacob_2D[i] = new double[NINT*NINT];
	}
	//We can determine the 2D Jacobian determinant for boundary condition here. 
	for (m = 0; m < ol[0].NEL; m++) { //loop through each element on wetted surface
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (l = 0; l < 4; l++) {
					t.XS_2D[m][0][0][i*NINT + j] = t.XS_2D[m][0][0][i*NINT + j] + GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][0]; //dx/dxi
					t.XS_2D[m][1][0][i*NINT + j] = t.XS_2D[m][1][0][i*NINT + j] + GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][1]; //dy/dxi
					t.XS_2D[m][0][1][i*NINT + j] = t.XS_2D[m][0][1][i*NINT + j] + GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][0]; //dx/deta
					t.XS_2D[m][1][1][i*NINT + j] = t.XS_2D[m][1][1][i*NINT + j] + GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][1]; //dy/deta
				}
				ol[0].Jacob_2D[m][i*NINT + j] = t.XS_2D[m][0][0][i*NINT + j] * t.XS_2D[m][1][1][i*NINT + j] - t.XS_2D[m][1][0][i*NINT + j] * t.XS_2D[m][0][1][i*NINT + j];
			}
		}
	}

	std::cout << " " << std::endl;
	return t;
}