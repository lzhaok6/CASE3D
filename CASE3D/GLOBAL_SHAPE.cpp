//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

//GLOBAL_SHAPE determines the global shape function derivative at the Gauss quad. points.

struct GLOBAL_SHAPEstruct GLOBAL_SHAPE(int NEL, double******XS, double****JACOB, double*****SHL) {
	int i, j, k, l, m, n;
	GLOBAL_SHAPEstruct t;
	t.SHG = new double*****[NEL];
	for (m = 0; m < NEL; m++) {
		t.SHG[m] = new double****[4];
		for (i = 0; i < 4; i++) {
			t.SHG[m][i] = new double***[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				t.SHG[m][i][j] = new double**[NINT];
				for (k = 0; k < NINT; k++) {
					t.SHG[m][i][j][k] = new double*[NINT];
					for (l = 0; l < NINT; l++) {
						t.SHG[m][i][j][k][l] = new double[NINT];
					}
				}
			}
		}
	}

	for (i = 0; i < NEL; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < NINT*NINT*NINT; k++) {
				for (l = 0; l < NINT; l++) {
					for (m = 0; m < NINT; m++) {
						for (n = 0; n < NINT; n++) {
							t.SHG[i][j][k][l][m][n] = 0.0;
						}
					}
				}
			}
		}
	}

	//---------------BEGIN GLOBAL_SHAPE COMP.----------------------------//
	//NINT is integration point number 
	for (m = 0; m < NEL; m++) {
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < NINT; l++) {
						t.SHG[m][0][i][j][k][l] = ((XS[m][1][1][j][k][l] * XS[m][2][2][j][k][l] - XS[m][1][2][j][k][l] * XS[m][2][1][j][k][l])*SHL[0][i][j][k][l] + (-XS[m][1][0][j][k][l] * XS[m][2][2][j][k][l] + XS[m][2][0][j][k][l] * XS[m][1][2][j][k][l])*SHL[1][i][j][k][l] + (XS[m][1][0][j][k][l] * XS[m][2][1][j][k][l] - XS[m][2][0][j][k][l] * XS[m][1][1][j][k][l])*SHL[2][i][j][k][l]) / JACOB[m][j][k][l];
						t.SHG[m][1][i][j][k][l] = ((-XS[m][0][1][j][k][l] * XS[m][2][2][j][k][l] + XS[m][2][1][j][k][l] * XS[m][0][2][j][k][l])*SHL[0][i][j][k][l] + (XS[m][0][0][j][k][l] * XS[m][2][2][j][k][l] - XS[m][2][0][j][k][l] * XS[m][0][2][j][k][l])*SHL[1][i][j][k][l] + (-XS[m][0][0][j][k][l] * XS[m][2][1][j][k][l] + XS[m][2][0][j][k][l] * XS[m][0][1][j][k][l])*SHL[2][i][j][k][l]) / JACOB[m][j][k][l];
						t.SHG[m][2][i][j][k][l] = ((XS[m][0][1][j][k][l] * XS[m][1][2][j][k][l] - XS[m][1][1][j][k][l] * XS[m][0][2][j][k][l])*SHL[0][i][j][k][l] + (-XS[m][0][0][j][k][l] * XS[m][1][2][j][k][l] + XS[m][1][0][j][k][l] * XS[m][0][2][j][k][l])*SHL[1][i][j][k][l] + (XS[m][0][0][j][k][l] * XS[m][1][1][j][k][l] - XS[m][1][0][j][k][l] * XS[m][0][1][j][k][l])*SHL[2][i][j][k][l]) / JACOB[m][j][k][l];
						t.SHG[m][3][i][j][k][l] = SHL[3][i][j][k][l]; //phi=phi 
					}
				}
			}
		}
	}
	return t;
}
