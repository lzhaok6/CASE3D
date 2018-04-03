//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

struct JACOBIANstruct JACOBIAN(int NEL, double*****GSHL, double **GCOORD, int **IEN, int*** LNA) {
	int i, j, k, l, m, n;
	JACOBIANstruct t; 
	
	t.XS = new double*****[NEL];   //dx/ds    3*3 matrix
	for (i = 0; i < NEL; i++) { //for x y z direction
		t.XS[i] = new double****[3];
		for (j = 0; j < 3; j++) { //for x y z direction
			t.XS[i][j] = new double***[3];
			for (k = 0; k < 3; k++) {
				t.XS[i][j][k] = new double**[NINT];
				for (l = 0; l < NINT; l++) {
					t.XS[i][j][k][l] = new double*[NINT];
					for (m = 0; m < NINT; m++) {
						t.XS[i][j][k][l][m] = new double[NINT];
					}
				}
			}
		}
	}
	//INITIALIZE t.XS
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < NINT; l++) {
					for (m = 0; m < NINT; m++) {
						for (n = 0; n < NINT; n++) {
							t.XS[i][j][k][l][m][n] = 0.0;
						}
					}
				}
			}
		}
	}

	t.JACOB = new double***[NEL];
	for (i = 0; i < NEL; i++) {
		t.JACOB[i] = new double**[NINT];
		for (j = 0; j < NINT; j++) {
			t.JACOB[i][j] = new double*[NINT];
			for (k = 0; k < NINT; k++) {
				t.JACOB[i][j][k] = new double[NINT];
			}
		}
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT; l++) {
					t.JACOB[i][j][k][l] = 0.0;
				}
			}
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

	//New Jacobian matrix definition
	for (m = 0; m < NEL; m++) {
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < 8; l++) {
						//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
						t.XS[m][0][0][i][j][k] = t.XS[m][0][0][i][j][k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
						t.XS[m][1][0][i][j][k] = t.XS[m][1][0][i][j][k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi
						t.XS[m][2][0][i][j][k] = t.XS[m][2][0][i][j][k] + GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi
						t.XS[m][0][1][i][j][k] = t.XS[m][0][1][i][j][k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta
						t.XS[m][1][1][i][j][k] = t.XS[m][1][1][i][j][k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
						t.XS[m][2][1][i][j][k] = t.XS[m][2][1][i][j][k] + GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta
						t.XS[m][0][2][i][j][k] = t.XS[m][0][2][i][j][k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta
						t.XS[m][1][2][i][j][k] = t.XS[m][1][2][i][j][k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta
						t.XS[m][2][2][i][j][k] = t.XS[m][2][2][i][j][k] + GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
					}
					//The determinant of Jacobian matrix
					t.JACOB[m][i][j][k] = t.XS[m][0][0][i][j][k] * (t.XS[m][1][1][i][j][k] * t.XS[m][2][2][i][j][k] - t.XS[m][2][1][i][j][k] * t.XS[m][1][2][i][j][k])
						- t.XS[m][1][0][i][j][k] * (t.XS[m][0][1][i][j][k] * t.XS[m][2][2][i][j][k] - t.XS[m][2][1][i][j][k] * t.XS[m][0][2][i][j][k])
						+ t.XS[m][2][0][i][j][k] * (t.XS[m][0][1][i][j][k] * t.XS[m][1][2][i][j][k] - t.XS[m][1][1][i][j][k] * t.XS[m][0][2][i][j][k]);
				}
			}
		}
	}

	std::cout << " " << std::endl;
	return t;
}