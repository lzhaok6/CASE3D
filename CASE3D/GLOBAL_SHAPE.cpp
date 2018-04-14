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

struct GLOBAL_SHAPEstruct GLOBAL_SHAPE(int NEL, double******XS, double****JACOB, double*****SHL, double****xs, double**jacob) {
	int i, j, k, l, m, n, e;
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

	t.shg = new double***[NEL];
	for (i = 0; i < NEL; i++) {
		t.shg[i] = new double**[4];
		for (j = 0; j < 4; j++) {
			t.shg[i][j] = new double*[NINT*NINT*NINT];
			for (k = 0; k < NINT*NINT*NINT; k++) {
				t.shg[i][j][k] = new double[NINT*NINT*NINT];
			}
		}
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < NINT*NINT*NINT; k++) {
				for (l = 0; l < NINT*NINT*NINT; l++) {
					t.shg[i][j][k][l] = 0.0;
				}
			}
		}
	}

	double **shl;
	shl = new double*[NINT*NINT*NINT];
	for (i = 0; i < NINT*NINT*NINT; i++) {
		shl[i] = new double[NINT*NINT*NINT];
	}
	for (i = 0; i < NINT*NINT*NINT; i++) {
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT; l++) {
					shl[i][j*NINT*NINT + k*NINT + l] = SHL[3][i][j][k][l];
				}
			}
		}
	}

	std::clock_t start;
	start = std::clock();
	//---------------BEGIN GLOBAL_SHAPE COMP.----------------------------//
	//NINT is integration point number 
	//store XS for j, k, l as one dimension ???
	/*
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
	*/
	
	for (m = 0; m < NEL; m++) {
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < NINT; l++) {
						t.shg[m][0][i][j*NINT*NINT + k*NINT + l] = ((xs[m][1][1][j*NINT*NINT + k*NINT + l] * xs[m][2][2][j*NINT*NINT + k*NINT + l] - xs[m][1][2][j*NINT*NINT + k*NINT + l] * xs[m][2][1][j*NINT*NINT + k*NINT + l])*SHL[0][i][j][k][l] + (-xs[m][1][0][j*NINT*NINT + k*NINT + l] * xs[m][2][2][j*NINT*NINT + k*NINT + l] + xs[m][2][0][j*NINT*NINT + k*NINT + l] * xs[m][1][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j][k][l] + (xs[m][1][0][j*NINT*NINT + k*NINT + l] * xs[m][2][1][j*NINT*NINT + k*NINT + l] - xs[m][2][0][j*NINT*NINT + k*NINT + l] * xs[m][1][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j][k][l]) / jacob[m][j*NINT*NINT + k*NINT + l];
						t.shg[m][1][i][j*NINT*NINT + k*NINT + l] = ((-xs[m][0][1][j*NINT*NINT + k*NINT + l] * xs[m][2][2][j*NINT*NINT + k*NINT + l] + xs[m][2][1][j*NINT*NINT + k*NINT + l] * xs[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j][k][l] + (xs[m][0][0][j*NINT*NINT + k*NINT + l] * xs[m][2][2][j*NINT*NINT + k*NINT + l] - xs[m][2][0][j*NINT*NINT + k*NINT + l] * xs[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j][k][l] + (-xs[m][0][0][j*NINT*NINT + k*NINT + l] * xs[m][2][1][j*NINT*NINT + k*NINT + l] + xs[m][2][0][j*NINT*NINT + k*NINT + l] * xs[m][0][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j][k][l]) / jacob[m][j*NINT*NINT + k*NINT + l];
						t.shg[m][2][i][j*NINT*NINT + k*NINT + l] = ((xs[m][0][1][j*NINT*NINT + k*NINT + l] * xs[m][1][2][j*NINT*NINT + k*NINT + l] - xs[m][1][1][j*NINT*NINT + k*NINT + l] * xs[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j][k][l] + (-xs[m][0][0][j*NINT*NINT + k*NINT + l] * xs[m][1][2][j*NINT*NINT + k*NINT + l] + xs[m][1][0][j*NINT*NINT + k*NINT + l] * xs[m][0][2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j][k][l] + (xs[m][0][0][j*NINT*NINT + k*NINT + l] * xs[m][1][1][j*NINT*NINT + k*NINT + l] - xs[m][1][0][j*NINT*NINT + k*NINT + l] * xs[m][0][1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j][k][l]) / jacob[m][j*NINT*NINT + k*NINT + l];
						t.shg[m][3][i][j*NINT*NINT + k*NINT + l] = SHL[3][i][j][k][l]; //phi=phi 
						//t.shg[m][3][i][j*NINT*NINT + k*NINT + l] = shl[i][j*NINT*NINT + k*NINT + l]; //phi=phi 
					}
				}
			}
		}
	}
	
	/*
	Eigen::MatrixXd JB(3, 3);
	Eigen::VectorXd shl(3);
	Eigen::VectorXd shg(3);
	for (e = 0; e < NEL; e++) {
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < NINT; j++) { //j,k,l are points 
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < NINT; l++) {
						//define the Eigen matrix to facilitate the calculation
						for (m = 0; m < 3; m++) {
							for (n = 0; n < 3; n++) {
								JB(m, n) = XS[e][m][n][j][k][l];
							}
							shl(m) = SHL[m][i][j][k][l];
						}
						shg = JB.inverse()*shl;
						t.SHG[e][0][i][j][k][l] = shg(0);
						t.SHG[e][1][i][j][k][l] = shg(1);
						t.SHG[e][2][i][j][k][l] = shg(2);
						t.SHG[e][3][i][j][k][l] = SHL[3][i][j][k][l]; //phi=phi 
					}
				}
			}
		}
		std::cout << e << std::endl;
	}
	*/
	double duration;
	duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
	std::cout << "total CPU time (ms): " << duration << std::endl;
	std::cout << " " << std::endl;
	return t;
}
