//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>

struct JACOBIANstruct JACOBIAN(int NEL, double **GCOORD, int **IEN, int*** LNA) {
	int i, j, k, l, m, n, z;
	JACOBIANstruct t;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber];
	extern STRU_WET_SURF ss[ssnumber]; //data structure used to store the properties on the structure wetted surface
	LOBATTOstruct bq;
	GLLQUADstruct gq;
	LOCAL_GSHAPEstruct lg;

	if (element_type == 0) { //hex element
		//initialize t.XS
		/*
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
		*/
		t.XS = new double**[NEL];
		for (i = 0; i < NEL; i++) {
			t.XS[i] = new double*[9];
			for (j = 0; j < 9; j++) {
				t.XS[i][j] = new double[NINT*NINT*NINT];
			}
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < 9; j++) {
				for (l = 0; l < NINT*NINT*NINT; l++) {
					t.XS[i][j][l] = 0.0;
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

		/*
		for (m = 0; m < NEL; m++) {
			std::cout << m << std::endl;
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						for (l = 0; l < 8; l++) {
							//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
							t.XS[m][0][0][i*NINT*NINT + j*NINT + k] = t.XS[m][0][0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
							t.XS[m][1][0][i*NINT*NINT + j*NINT + k] = t.XS[m][1][0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi
							t.XS[m][2][0][i*NINT*NINT + j*NINT + k] = t.XS[m][2][0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi
							t.XS[m][0][1][i*NINT*NINT + j*NINT + k] = t.XS[m][0][1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta
							t.XS[m][1][1][i*NINT*NINT + j*NINT + k] = t.XS[m][1][1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
							t.XS[m][2][1][i*NINT*NINT + j*NINT + k] = t.XS[m][2][1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta
							t.XS[m][0][2][i*NINT*NINT + j*NINT + k] = t.XS[m][0][2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta
							t.XS[m][1][2][i*NINT*NINT + j*NINT + k] = t.XS[m][1][2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta
							t.XS[m][2][2][i*NINT*NINT + j*NINT + k] = t.XS[m][2][2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
						}
						//The determinant of jacobian matrix for each node in the element
						t.JACOB[m][i*NINT*NINT + j*NINT + k] = t.XS[m][0][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][1][1][i*NINT*NINT + j*NINT + k] * t.XS[m][2][2][i*NINT*NINT + j*NINT + k] - t.XS[m][2][1][i*NINT*NINT + j*NINT + k] * t.XS[m][1][2][i*NINT*NINT + j*NINT + k])
							- t.XS[m][1][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][0][1][i*NINT*NINT + j*NINT + k] * t.XS[m][2][2][i*NINT*NINT + j*NINT + k] - t.XS[m][2][1][i*NINT*NINT + j*NINT + k] * t.XS[m][0][2][i*NINT*NINT + j*NINT + k])
							+ t.XS[m][2][0][i*NINT*NINT + j*NINT + k] * (t.XS[m][0][1][i*NINT*NINT + j*NINT + k] * t.XS[m][1][2][i*NINT*NINT + j*NINT + k] - t.XS[m][1][1][i*NINT*NINT + j*NINT + k] * t.XS[m][0][2][i*NINT*NINT + j*NINT + k]);
					}
				}
			}
		}
		*/
		bq = LOBATTO(N);
		gq = GLLQUAD(bq.Z, bq.WL, N, !FEM);
		lg = LOCAL_GSHAPE(gq.S, LNA, NINT);
		for (m = 0; m < NEL; m++) {
			//std::cout << m << std::endl;
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						for (l = 0; l < 8; l++) {
							//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
							/*
							t.XS[m][3 * 0 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
							t.XS[m][3 * 1 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi
							t.XS[m][3 * 2 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi
							t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta
							t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
							t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta
							t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta
							t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta
							t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
							*/
							t.XS[m][3 * 0 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dxi
							t.XS[m][3 * 1 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/deta 
							t.XS[m][3 * 2 + 0][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][0]; //dx/dzeta  
							t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dxi 
							t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/deta
							t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][1]; //dy/dzeta 
							t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dxi 
							t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/deta 
							t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][m] - 1][2]; //dz/dzeta
						}
						//The determinant of jacobian matrix for each node in the element
						t.JACOB[m][i*NINT*NINT + j*NINT + k] = t.XS[m][3 * 0 + 0][i*NINT*NINT + j*NINT + k] * (t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] - t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k])
							- t.XS[m][3 * 1 + 0][i*NINT*NINT + j*NINT + k] * (t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 2 + 2][i*NINT*NINT + j*NINT + k] - t.XS[m][3 * 2 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k])
							+ t.XS[m][3 * 2 + 0][i*NINT*NINT + j*NINT + k] * (t.XS[m][3 * 0 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 1 + 2][i*NINT*NINT + j*NINT + k] - t.XS[m][3 * 1 + 1][i*NINT*NINT + j*NINT + k] * t.XS[m][3 * 0 + 2][i*NINT*NINT + j*NINT + k]);
					}
				}
			}
		}
		/*
		for (z = 0; z < nrbsurfnumber; z++) {
			nr[z].xs = new double***[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].xs[i] = new double**[3];
				for (j = 0; j < 3; j++) {
					nr[z].xs[i][j] = new double*[3];
					for (k = 0; k < 3; k++) {
						nr[z].xs[i][j][k] = new double[NqINT*NqINT*NqINT];
					}
				}
			}
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 3; k++) {
						for (l = 0; l < NqINT*NqINT*NqINT; l++) {
							nr[z].xs[i][j][k][l] = 0.0;
						}
					}
				}
			}
			nr[z].JACOB = new double*[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].JACOB[i] = new double[NqINT*NqINT*NqINT];
			}
			for (m = 0; m < nr[z].NEL_nrb; m++) {
				//std::cout << m << std::endl;
				for (i = 0; i < NqINT; i++) {
					for (j = 0; j < NqINT; j++) {
						for (k = 0; k < NqINT; k++) {
							for (l = 0; l < 8; l++) {
								//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
								nr[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][0]; //dx/dxi
								nr[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][1]; //dy/dxi
								nr[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][2]; //dz/dxi
								nr[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][0]; //dx/deta
								nr[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][1]; //dy/deta
								nr[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][2]; //dz/deta
								nr[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][0]; //dx/dzeta
								nr[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][1]; //dy/dzeta
								nr[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][nr[z].NRBELE_ARR[m] - 1] - 1][2]; //dz/dzeta
							}
							//The determinant of jacobian matrix for each node in the element
							nr[z].JACOB[m][i*NqINT*NqINT + j*NqINT + k] = nr[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] * (nr[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - nr[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k])
								- nr[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] * (nr[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - nr[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k])
								+ nr[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] * (nr[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] - nr[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * nr[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k]);
						}
					}
				}
			}
		}
		*/

		//Derive the Jacobian determinant on boundary elements (wetted surface)
		/*
		for (z = 0; z < owsfnumber; z++) {
			ol[z].xs = new double***[ol[z].FSNEL];
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].xs[i] = new double**[3];
				for (j = 0; j < 3; j++) {
					ol[z].xs[i][j] = new double*[3];
					for (k = 0; k < 3; k++) {
						ol[z].xs[i][j][k] = new double[NqINT*NqINT*NqINT];
					}
				}
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 3; k++) {
						for (l = 0; l < NqINT*NqINT*NqINT; l++) {
							ol[z].xs[i][j][k][l] = 0.0;
						}
					}
				}
			}
			ol[z].JACOB = new double*[ol[z].FSNEL];
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].JACOB[i] = new double[NqINT*NqINT*NqINT];
			}
			for (m = 0; m < ol[z].FSNEL; m++) {
				//std::cout << m << std::endl;
				for (i = 0; i < NqINT; i++) {
					for (j = 0; j < NqINT; j++) {
						for (k = 0; k < NqINT; k++) {
							for (l = 0; l < 8; l++) {
								//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
								ol[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][0]; //dx/dxi
								ol[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][1]; //dy/dxi
								ol[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][2]; //dz/dxi
								ol[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][0]; //dx/deta
								ol[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][1]; //dy/deta
								ol[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][2]; //dz/deta
								ol[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][0]; //dx/dzeta
								ol[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][1]; //dy/dzeta
								ol[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][ol[z].GIDF[m] - 1] - 1][2]; //dz/dzeta
							}
							//The determinant of jacobian matrix for each node in the element
							ol[z].JACOB[m][i*NqINT*NqINT + j*NqINT + k] = ol[z].xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] * (ol[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - ol[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k])
								- ol[z].xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] * (ol[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - ol[z].xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k])
								+ ol[z].xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] * (ol[z].xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] - ol[z].xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * ol[z].xs[m][0][2][i*NqINT*NqINT + j*NqINT + k]);
						}
					}
				}
			}
		}
		*/

		if (mappingalgo == 2) {
			bq = LOBATTO(Nq);
			gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
			//the quadrature point is Gauss-Legendre point for FEM element and GLL for SEM element. 
			lg = LOCAL_GSHAPE(gq.S, LNA, NqINT);
			for (z = 0; z < owsfnumber; z++) {
				//GSHL_2D is the 2D shape function and its derivative with respect to 2D local coordinate 
				ol[z].GSHL_2D = new double***[3];
				for (i = 0; i < 3; i++) {
					ol[z].GSHL_2D[i] = new double**[4]; //4 points
					for (j = 0; j < 4; j++) {
						ol[z].GSHL_2D[i][j] = new double*[NqINT];
						for (k = 0; k < NqINT; k++) {
							ol[z].GSHL_2D[i][j][k] = new double[NqINT];
						}
					}
				}
				for (i = 0; i < 4; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
					for (j = 0; j < NqINT; j++) { // j k l are integration point where discrete value of shape function is derived
						for (k = 0; k < NqINT; k++) {
							ol[z].GSHL_2D[2][i][j][k] = (1.0 / 4.0)*(1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * gq.S[j])*(1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * gq.S[k]);
							//not derivative
							ol[z].GSHL_2D[0][i][j][k] = (1.0 / 4.0)*lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * (1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * gq.S[k]);
							//X direction derivative
							ol[z].GSHL_2D[1][i][j][k] = (1.0 / 4.0)*lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * (1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * gq.S[j]);
							//Z direction derivative
						}
					}
				}
				//If mapping algorithm==2, we could make a 2D Jacobian for fluid nodal force (to sent to structure) integration.  
				//Same number of quadrature point as the 3D integration Jacobian is used (Nq). 
				ol[z].Jacob_2D = new double*[ol[z].FSNEL];
				for (i = 0; i < ol[z].FSNEL; i++) {
					ol[z].Jacob_2D[i] = new double[NqINT*NqINT];
				}
				ol[z].xs_2D = new double***[ol[z].FSNEL];
				for (i = 0; i < ol[z].FSNEL; i++) {
					ol[z].xs_2D[i] = new double**[3];
					for (j = 0; j < 3; j++) {
						ol[z].xs_2D[i][j] = new double*[2];
						for (k = 0; k < 2; k++) {
							ol[z].xs_2D[i][j][k] = new double[NqINT*NqINT];
						}
					}
				}
				for (i = 0; i < ol[z].FSNEL; i++) {
					for (j = 0; j < 3; j++) {
						for (k = 0; k < 2; k++) {
							for (l = 0; l < NqINT*NqINT; l++) {
								ol[z].xs_2D[i][j][k][l] = 0.0;
							}
						}
					}
				}
				//We can determine the 2D Jacobian determinant for boundary condition here. 
				//Assuming linear geometric property mapping here. 
				for (m = 0; m < ol[z].FSNEL; m++) { //loop through each element on wetted surface
					for (i = 0; i < NqINT; i++) {
						for (j = 0; j < NqINT; j++) {
							for (l = 0; l < 4; l++) { //4 nodes on the linear element
								//Is the l in ol[z].GSHL_2D[0][l][i][j] and in ol[z].LNA_norm[l] must be consistent (the same surface). 
								ol[z].xs_2D[m][0][0][i * NqINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dnu
								ol[z].xs_2D[m][1][0][i * NqINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dnu
								ol[z].xs_2D[m][2][0][i * NqINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dnu
								ol[z].xs_2D[m][0][1][i * NqINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dmu
								ol[z].xs_2D[m][1][1][i * NqINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dmu
								ol[z].xs_2D[m][2][1][i * NqINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dmu
							}
							ol[z].Jacob_2D[m][i * NqINT + j] = pow(pow((ol[z].xs_2D[m][0][0][i * NqINT + j] * ol[z].xs_2D[m][1][1][i * NqINT + j] - ol[z].xs_2D[m][0][1][i * NqINT + j] * ol[z].xs_2D[m][1][0][i * NqINT + j]), 2) +
								pow((ol[z].xs_2D[m][0][0][i * NqINT + j] * ol[z].xs_2D[m][2][1][i * NqINT + j] - ol[z].xs_2D[m][0][1][i * NqINT + j] * ol[z].xs_2D[m][2][0][i * NqINT + j]), 2) +
								pow((ol[z].xs_2D[m][1][0][i * NqINT + j] * ol[z].xs_2D[m][2][1][i * NqINT + j] - ol[z].xs_2D[m][1][1][i * NqINT + j] * ol[z].xs_2D[m][2][0][i * NqINT + j]), 2), 0.5);
						}
					}
				}
				//check if the value is positive and check if i and j is exchanged, would the value of FPMASTER_2D be changed!!!
			}
		}
		if (mappingalgo == 5) {
			bq = LOBATTO(N);
			gq = GLLQUAD(bq.Z, bq.WL, N, 0);
			//the quadrature point is Gauss-Legendre point for whatever element type.  
			lg = LOCAL_GSHAPE(gq.S, LNA, N);
			for (z = 0; z < owsfnumber; z++) {
				//GSHL_2D is the 2D shape function and its derivative with respect to 2D local coordinate. 
				ol[z].GSHL_2D = new double***[3];
				for (i = 0; i < 3; i++) {
					ol[z].GSHL_2D[i] = new double**[4]; //4 points
					for (j = 0; j < 4; j++) {
						ol[z].GSHL_2D[i][j] = new double*[NINT];
						for (k = 0; k < NINT; k++) {
							ol[z].GSHL_2D[i][j][k] = new double[NINT];
						}
					}
				}
				for (i = 0; i < 4; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
					for (j = 0; j < NINT; j++) { // j k l are integration point where discrete value of shape function is derived
						for (k = 0; k < NINT; k++) {
							ol[z].GSHL_2D[2][i][j][k] = (1.0 / 4.0)*(1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * gq.S[j])*(1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * gq.S[k]);
							//not derivative
							ol[z].GSHL_2D[0][i][j][k] = (1.0 / 4.0)*lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * (1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * gq.S[k]);
							//X direction derivative
							ol[z].GSHL_2D[1][i][j][k] = (1.0 / 4.0)*lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[1]] * (1 + lg.MCOORD[ol[z].LNA_JB2D[i] - 1][ol[z].Jacob_face[0]] * gq.S[j]);
							//Z direction derivative
						}
					}
				}
				//If mapping algorithm==2, we could make a 2D Jacobian for fluid nodal force (to sent to structure) integration.  
				//Same number of quadrature point as the 3D integration Jacobian is used (Nq). 
				ol[z].Jacob_2D = new double*[ol[z].FSNEL];
				for (i = 0; i < ol[z].FSNEL; i++) {
					ol[z].Jacob_2D[i] = new double[NINT*NINT];
				}
				ol[z].xs_2D = new double***[ol[z].FSNEL];
				for (i = 0; i < ol[z].FSNEL; i++) {
					ol[z].xs_2D[i] = new double**[3];
					for (j = 0; j < 3; j++) {
						ol[z].xs_2D[i][j] = new double*[2];
						for (k = 0; k < 2; k++) {
							ol[z].xs_2D[i][j][k] = new double[NINT*NINT];
						}
					}
				}
				for (i = 0; i < ol[z].FSNEL; i++) {
					for (j = 0; j < 3; j++) {
						for (k = 0; k < 2; k++) {
							for (l = 0; l < NINT*NINT; l++) {
								ol[z].xs_2D[i][j][k][l] = 0.0;
							}
						}
					}
				}
				//We can determine the 2D Jacobian determinant for boundary condition here. 
				//Assuming linear geometric property mapping here. 
				for (m = 0; m < ol[z].FSNEL; m++) { //loop through each element on wetted surface
					for (i = 0; i < NINT; i++) {
						for (j = 0; j < NINT; j++) {
							for (l = 0; l < 4; l++) { //4 nodes on the linear element
								//Is the l in ol[z].GSHL_2D[0][l][i][j] and in ol[z].LNA_norm[l] must be consistent (the same surface). 
								ol[z].xs_2D[m][0][0][i * NINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dnu
								ol[z].xs_2D[m][1][0][i * NINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dnu
								ol[z].xs_2D[m][2][0][i * NINT + j] += ol[z].GSHL_2D[0][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dnu
								ol[z].xs_2D[m][0][1][i * NINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dmu
								ol[z].xs_2D[m][1][1][i * NINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dmu
								ol[z].xs_2D[m][2][1][i * NINT + j] += ol[z].GSHL_2D[1][l][i][j] * GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dmu
							}
							ol[z].Jacob_2D[m][i * NINT + j] = pow(pow((ol[z].xs_2D[m][0][0][i * NINT + j] * ol[z].xs_2D[m][1][1][i * NINT + j] - ol[z].xs_2D[m][0][1][i * NINT + j] * ol[z].xs_2D[m][1][0][i * NINT + j]), 2) +
								pow((ol[z].xs_2D[m][0][0][i * NINT + j] * ol[z].xs_2D[m][2][1][i * NINT + j] - ol[z].xs_2D[m][0][1][i * NINT + j] * ol[z].xs_2D[m][2][0][i * NINT + j]), 2) +
								pow((ol[z].xs_2D[m][1][0][i * NINT + j] * ol[z].xs_2D[m][2][1][i * NINT + j] - ol[z].xs_2D[m][1][1][i * NINT + j] * ol[z].xs_2D[m][2][0][i * NINT + j]), 2), 0.5);
						}
					}
				}
				//check if the value is positive and check if i and j is exchanged, would the value of FPMASTER_2D be changed!!!
			}
		}

		//For NRB boundary force integration, we still use the approach in algorithm 2 (insert one more quadrature point)
		bq = LOBATTO(Nq);
		gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
		lg = LOCAL_GSHAPE(gq.S, LNA, NqINT);
		for (z = 0; z < nrbsurfnumber; z++) {
			nr[z].GSHL_2D = new double***[3];
			for (i = 0; i < 3; i++) {
				nr[z].GSHL_2D[i] = new double**[4]; //4 points
				for (j = 0; j < 4; j++) {
					nr[z].GSHL_2D[i][j] = new double*[NqINT];
					for (k = 0; k < NqINT; k++) {
						nr[z].GSHL_2D[i][j][k] = new double[NqINT];
					}
				}
			}
			for (i = 0; i < 4; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
				for (j = 0; j < NqINT; j++) { // j k l are integration point where discrete value of shape function is derived
					for (k = 0; k < NqINT; k++) {
						nr[z].GSHL_2D[2][i][j][k] = (1.0 / 4.0)*(1 + lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[0]] * gq.S[j])*(1 + lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[1]] * gq.S[k]);
						//not derivative
						nr[z].GSHL_2D[0][i][j][k] = (1.0 / 4.0)*lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[0]] * (1 + lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[1]] * gq.S[k]);
						//X direction derivative
						nr[z].GSHL_2D[1][i][j][k] = (1.0 / 4.0)*lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[1]] * (1 + lg.MCOORD[nr[z].LNA_JB2D[i] - 1][nr[z].Jacob_face[0]] * gq.S[j]);
						//Z direction derivative
					}
				}
			}
			nr[z].Jacob_2D = new double*[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].Jacob_2D[i] = new double[NqINT*NqINT];
			}
			nr[z].xs_2D = new double***[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].xs_2D[i] = new double**[3];
				for (j = 0; j < 3; j++) {
					nr[z].xs_2D[i][j] = new double*[2];
					for (k = 0; k < 2; k++) {
						nr[z].xs_2D[i][j][k] = new double[NqINT*NqINT];
					}
				}
			}
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 2; k++) {
						for (l = 0; l < NqINT*NqINT; l++) {
							nr[z].xs_2D[i][j][k][l] = 0.0;
						}
					}
				}
			}
			for (m = 0; m < nr[z].NEL_nrb; m++) { //loop through each element on wetted surface
				for (i = 0; i < NqINT; i++) {
					for (j = 0; j < NqINT; j++) {
						for (l = 0; l < 4; l++) { //4 nodes on the linear element
							//Is the l in nr[z].GSHL_2D[0][l][i][j] and in nr[z].LNA_norm[l] must be consistent (the same surface). 
							nr[z].xs_2D[m][0][0][i * NqINT + j] += nr[z].GSHL_2D[0][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dnu
							nr[z].xs_2D[m][1][0][i * NqINT + j] += nr[z].GSHL_2D[0][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dnu
							nr[z].xs_2D[m][2][0][i * NqINT + j] += nr[z].GSHL_2D[0][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dnu
							nr[z].xs_2D[m][0][1][i * NqINT + j] += nr[z].GSHL_2D[1][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][0]; //dx/dmu
							nr[z].xs_2D[m][1][1][i * NqINT + j] += nr[z].GSHL_2D[1][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][1]; //dy/dmu
							nr[z].xs_2D[m][2][1][i * NqINT + j] += nr[z].GSHL_2D[1][l][i][j] * GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[l] - 1][m] - 1][2]; //dz/dmu
						}
						nr[z].Jacob_2D[m][i * NqINT + j] = pow(pow((nr[z].xs_2D[m][0][0][i * NqINT + j] * nr[z].xs_2D[m][1][1][i * NqINT + j] - nr[z].xs_2D[m][0][1][i * NqINT + j] * nr[z].xs_2D[m][1][0][i * NqINT + j]), 2) +
							pow((nr[z].xs_2D[m][0][0][i * NqINT + j] * nr[z].xs_2D[m][2][1][i * NqINT + j] - nr[z].xs_2D[m][0][1][i * NqINT + j] * nr[z].xs_2D[m][2][0][i * NqINT + j]), 2) +
							pow((nr[z].xs_2D[m][1][0][i * NqINT + j] * nr[z].xs_2D[m][2][1][i * NqINT + j] - nr[z].xs_2D[m][1][1][i * NqINT + j] * nr[z].xs_2D[m][2][0][i * NqINT + j]), 2), 0.5);
					}
				}
			}
		}
		//Finish the definition for hexahedral element 
	}
	if (element_type == 1) {
		t.JACOB_tet = new double[NEL]; //For first order tetrahedral element, the Jacobian determinant is the same throught the element
		double XS_tet[3][3];
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
			XS_tet[0][0] = x21; XS_tet[0][1] = x31; XS_tet[0][2] = x41;
			XS_tet[1][0] = y21; XS_tet[1][1] = y31; XS_tet[1][2] = y41;
			XS_tet[2][0] = z21; XS_tet[2][1] = z31; XS_tet[2][2] = z41;
			t.JACOB_tet[i] = XS_tet[0][0] * (XS_tet[1][1] * XS_tet[2][2] - XS_tet[2][1] * XS_tet[1][2])
				- XS_tet[1][0] * (XS_tet[0][1] * XS_tet[2][2] - XS_tet[2][1] * XS_tet[0][2])
				+ XS_tet[2][0] * (XS_tet[0][1] * XS_tet[1][2] - XS_tet[1][1] * XS_tet[0][2]);
		}
	}

	//Define the Jacobian matrix for structural boundary force integration
	bq = LOBATTO(hprefg);
	gq = GLLQUAD(bq.Z, bq.WL, hprefg, FEM);
	lg = LOCAL_GSHAPE(gq.S, LNA, hprefg + 1);
	for (z = 0; z < ssnumber; z++) {
		ss[z].GSHL_2D = new double***[3];
		for (i = 0; i < 3; i++) {
			ss[z].GSHL_2D[i] = new double**[4]; //4 points
			for (j = 0; j < 4; j++) {
				ss[z].GSHL_2D[i][j] = new double*[hprefg + 1];
				for (k = 0; k < hprefg + 1; k++) {
					ss[z].GSHL_2D[i][j][k] = new double[hprefg + 1];
				}
			}
		}
		for (i = 0; i < 4; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
			for (j = 0; j < hprefg + 1; j++) { // j k are integration point where discrete value of shape function is derived
				for (k = 0; k < hprefg + 1; k++) {
					ss[z].GSHL_2D[2][i][j][k] = (1.0 / 4.0)*(1 + lg.MCOORD_2D[i][0] * gq.S[j])*(1 + lg.MCOORD_2D[i][1] * gq.S[k]);
					//not derivative
					ss[z].GSHL_2D[0][i][j][k] = (1.0 / 4.0)*lg.MCOORD_2D[i][0] * (1 + lg.MCOORD_2D[i][1] * gq.S[k]);
					//xi direction derivative
					ss[z].GSHL_2D[1][i][j][k] = (1.0 / 4.0)*lg.MCOORD_2D[i][1] * (1 + lg.MCOORD_2D[i][0] * gq.S[j]);
					//eta direction derivatived
				}
			}
		}
		ss[z].Jacob_stru = new double*[ss[z].ELE_stru];
		for (i = 0; i < ss[z].ELE_stru; i++) {
			ss[z].Jacob_stru[i] = new double[(hprefg + 1)*(hprefg + 1)];
		}
		ss[z].xs_2D = new double***[ss[z].ELE_stru];
		for (i = 0; i < ss[z].ELE_stru; i++) {
			ss[z].xs_2D[i] = new double**[3];
			for (j = 0; j < 3; j++) {
				ss[z].xs_2D[i][j] = new double*[2];
				for (k = 0; k < 2; k++) {
					ss[z].xs_2D[i][j][k] = new double[(hprefg + 1)*(hprefg + 1)];
				}
			}
		}
		for (i = 0; i < ss[z].ELE_stru; i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 2; k++) {
					for (l = 0; l < (hprefg + 1)*(hprefg + 1); l++) {
						ss[z].xs_2D[i][j][k][l] = 0.0;
					}
				}
			}
		}
		for (m = 0; m < ss[z].ELE_stru; m++) { //loop through each element on wetted surface
			for (i = 0; i < (hprefg + 1); i++) {
				for (j = 0; j < (hprefg + 1); j++) {
					for (l = 0; l < 4; l++) { //4 nodes on the linear element
						//Is the l in ss[z].GSHL_2D[0][l][i][j] and in ss[z].LNA_norm[l] must be consistent (the same surface). 
						ss[z].xs_2D[m][0][0][i * (hprefg + 1) + j] += ss[z].GSHL_2D[0][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][0]; //dx/dnu
						ss[z].xs_2D[m][1][0][i * (hprefg + 1) + j] += ss[z].GSHL_2D[0][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][1]; //dy/dnu
						ss[z].xs_2D[m][2][0][i * (hprefg + 1) + j] += ss[z].GSHL_2D[0][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][2]; //dz/dnu
						ss[z].xs_2D[m][0][1][i * (hprefg + 1) + j] += ss[z].GSHL_2D[1][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][0]; //dx/dmu
						ss[z].xs_2D[m][1][1][i * (hprefg + 1) + j] += ss[z].GSHL_2D[1][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][1]; //dy/dmu
						ss[z].xs_2D[m][2][1][i * (hprefg + 1) + j] += ss[z].GSHL_2D[1][l][i][j] * ss[z].GCOORD_stru[ss[z].IEN_stru[l][m] - 1][2]; //dz/dmu
					}
					ss[z].Jacob_stru[m][i * (hprefg + 1) + j] = pow(pow((ss[z].xs_2D[m][0][0][i * (hprefg + 1) + j] * ss[z].xs_2D[m][1][1][i * (hprefg + 1) + j] - ss[z].xs_2D[m][0][1][i * (hprefg + 1) + j] * ss[z].xs_2D[m][1][0][i * (hprefg + 1) + j]), 2) +
						pow((ss[z].xs_2D[m][0][0][i * (hprefg + 1) + j] * ss[z].xs_2D[m][2][1][i * (hprefg + 1) + j] - ss[z].xs_2D[m][0][1][i * (hprefg + 1) + j] * ss[z].xs_2D[m][2][0][i * (hprefg + 1) + j]), 2) +
						pow((ss[z].xs_2D[m][1][0][i * (hprefg + 1) + j] * ss[z].xs_2D[m][2][1][i * (hprefg + 1) + j] - ss[z].xs_2D[m][1][1][i * (hprefg + 1) + j] * ss[z].xs_2D[m][2][0][i * (hprefg + 1) + j]), 2), 0.5);
				}
			}
		}
	}
	
	/*
	//The total value of the Jacobian determinant should be the 6 times of the total volume 
	double hd = 0.0; 
	for (i = 0; i < NEL; i++) {
		hd += t.JACOB_tet[i]; 
	}
	*/
	std::cout << " " << std::endl;
	return t;
}