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

struct GLOBAL_SHAPEstruct GLOBAL_SHAPE(int NEL, double***SHL, double****XS, double**JACOB) {
	int i, j, k, l, m, n, e;
	GLOBAL_SHAPEstruct t;

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

	std::clock_t start;
	start = std::clock();
	//---------------BEGIN GLOBAL_SHAPE COMP.----------------------------//
	for (m = 0; m < NEL; m++) {
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < NINT; l++) {
						t.SHG[m][0][i][j*NINT*NINT + k*NINT + l] = ((XS[m][1][1][j*NINT*NINT + k*NINT + l] * XS[m][2][2][j*NINT*NINT + k*NINT + l] - XS[m][1][2][j*NINT*NINT + k*NINT + l] * XS[m][2][1][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[m][1][0][j*NINT*NINT + k*NINT + l] * XS[m][2][2][j*NINT*NINT + k*NINT + l] + XS[m][2][0][j*NINT*NINT + k*NINT + l] * XS[m][1][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[m][1][0][j*NINT*NINT + k*NINT + l] * XS[m][2][1][j*NINT*NINT + k*NINT + l] - XS[m][2][0][j*NINT*NINT + k*NINT + l] * XS[m][1][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[m][j*NINT*NINT + k*NINT + l]; //d(phi)/dx
						t.SHG[m][1][i][j*NINT*NINT + k*NINT + l] = ((-XS[m][0][1][j*NINT*NINT + k*NINT + l] * XS[m][2][2][j*NINT*NINT + k*NINT + l] + XS[m][2][1][j*NINT*NINT + k*NINT + l] * XS[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (XS[m][0][0][j*NINT*NINT + k*NINT + l] * XS[m][2][2][j*NINT*NINT + k*NINT + l] - XS[m][2][0][j*NINT*NINT + k*NINT + l] * XS[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (-XS[m][0][0][j*NINT*NINT + k*NINT + l] * XS[m][2][1][j*NINT*NINT + k*NINT + l] + XS[m][2][0][j*NINT*NINT + k*NINT + l] * XS[m][0][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[m][j*NINT*NINT + k*NINT + l]; //d(phi)/dy
						t.SHG[m][2][i][j*NINT*NINT + k*NINT + l] = ((XS[m][0][1][j*NINT*NINT + k*NINT + l] * XS[m][1][2][j*NINT*NINT + k*NINT + l] - XS[m][1][1][j*NINT*NINT + k*NINT + l] * XS[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[m][0][0][j*NINT*NINT + k*NINT + l] * XS[m][1][2][j*NINT*NINT + k*NINT + l] + XS[m][1][0][j*NINT*NINT + k*NINT + l] * XS[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[m][0][0][j*NINT*NINT + k*NINT + l] * XS[m][1][1][j*NINT*NINT + k*NINT + l] - XS[m][1][0][j*NINT*NINT + k*NINT + l] * XS[m][0][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[m][j*NINT*NINT + k*NINT + l]; //d(phi)/dz
						t.SHG[m][3][i][j*NINT*NINT + k*NINT + l] = SHL[3][i][j*NINT*NINT + k*NINT + l]; //phi=phi 
					}
				}
			}
		}
	}

	double duration;
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
	std::cout << "total CPU time (ms): " << duration << std::endl;
	std::cout << " " << std::endl;
	return t;
}
