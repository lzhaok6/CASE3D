//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>

struct JACOBIANstruct JACOBIAN(int NEL, double **GCOORD, int **IEN, int*** LNA) {
	int i, j, k, l, m, n;
	JACOBIANstruct t;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber]; 

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

	LOBATTOstruct bq;
	bq = LOBATTO(N);
	GLLQUADstruct gq;
	gq = GLLQUAD(bq.Z, bq.WL, N, !FEM);
	LOCAL_GSHAPEstruct lg;

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

	bq = LOBATTO(Nq);
	gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
	lg = LOCAL_GSHAPE(gq.S, LNA, NqINT);

	//Derive the Jacobian determinant on boundary elements (NRB)
	double**** xs; 
	xs = new double***[nr[0].NEL_nrb];
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		xs[i] = new double**[3];
		for (j = 0; j < 3; j++) {
			xs[i][j] = new double*[3];
			for (k = 0; k < 3; k++) {
				xs[i][j][k] = new double[NqINT*NqINT*NqINT];
			}
		}
	}
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < NqINT*NqINT*NqINT; l++) {
					xs[i][j][k][l] = 0.0;
				}
			}
		}
	}
	
	nr[0].JACOB = new double*[nr[0].NEL_nrb]; 
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		nr[0].JACOB[i] = new double[NqINT*NqINT*NqINT];
	}
	for (m = 0; m < nr[0].NEL_nrb; m++) {
		//std::cout << m << std::endl;
		for (i = 0; i < NqINT; i++) {
			for (j = 0; j < NqINT; j++) {
				for (k = 0; k < NqINT; k++) {
					for (l = 0; l < 8; l++) {
						//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
						xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] = xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/dxi
						xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] = xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/dxi
						xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] = xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/dxi
						xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] = xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/deta
						xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] = xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/deta
						xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] = xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/deta
						xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] = xs[m][0][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/dzeta
						xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] = xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/dzeta
						xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] = xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[nr[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/dzeta
					}
					//The determinant of jacobian matrix for each node in the element
					nr[0].JACOB[m][i*NqINT*NqINT + j*NqINT + k] = xs[m][0][0][i*NqINT*NqINT + j*NqINT + k] * (xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][1][2][i*NqINT*NqINT + j*NqINT + k])
						- xs[m][1][0][i*NqINT*NqINT + j*NqINT + k] * (xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][2][2][i*NqINT*NqINT + j*NqINT + k] - xs[m][2][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][0][2][i*NqINT*NqINT + j*NqINT + k])
						+ xs[m][2][0][i*NqINT*NqINT + j*NqINT + k] * (xs[m][0][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][1][2][i*NqINT*NqINT + j*NqINT + k] - xs[m][1][1][i*NqINT*NqINT + j*NqINT + k] * xs[m][0][2][i*NqINT*NqINT + j*NqINT + k]);
				}
			}
		}
	}

	//Derive the Jacobian determinant on boundary elements (wetted surface)
	double**** xs2;
	xs2 = new double***[ol[0].FSNEL];
	for (i = 0; i < ol[0].FSNEL; i++) {
		xs2[i] = new double**[3];
		for (j = 0; j < 3; j++) {
			xs2[i][j] = new double*[3];
			for (k = 0; k < 3; k++) {
				xs2[i][j][k] = new double[NqINT*NqINT*NqINT];
			}
		}
	}
	for (i = 0; i < ol[0].FSNEL; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < 3; k++) {
				for (l = 0; l < NqINT*NqINT*NqINT; l++) {
					xs2[i][j][k][l] = 0.0;
				}
			}
		}
	}
	ol[0].JACOB = new double*[ol[0].FSNEL];
	for (i = 0; i < ol[0].FSNEL; i++) {
		ol[0].JACOB[i] = new double[NqINT*NqINT*NqINT];
	}
	for (m = 0; m < ol[0].FSNEL; m++) {
		//std::cout << m << std::endl;
		for (i = 0; i < NqINT; i++) {
			for (j = 0; j < NqINT; j++) {
				for (k = 0; k < NqINT; k++) {
					for (l = 0; l < 8; l++) {
						//same with 2D code, the transpose of that in Mindmap (the determinant of two transpose matrix is the same)
						xs2[m][0][0][i*NqINT*NqINT + j*NqINT + k] = xs2[m][0][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/dxi
						xs2[m][1][0][i*NqINT*NqINT + j*NqINT + k] = xs2[m][1][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/dxi
						xs2[m][2][0][i*NqINT*NqINT + j*NqINT + k] = xs2[m][2][0][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/dxi
						xs2[m][0][1][i*NqINT*NqINT + j*NqINT + k] = xs2[m][0][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/deta
						xs2[m][1][1][i*NqINT*NqINT + j*NqINT + k] = xs2[m][1][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/deta
						xs2[m][2][1][i*NqINT*NqINT + j*NqINT + k] = xs2[m][2][1][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/deta
						xs2[m][0][2][i*NqINT*NqINT + j*NqINT + k] = xs2[m][0][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][0]; //dx/dzeta
						xs2[m][1][2][i*NqINT*NqINT + j*NqINT + k] = xs2[m][1][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][1]; //dy/dzeta
						xs2[m][2][2][i*NqINT*NqINT + j*NqINT + k] = xs2[m][2][2][i*NqINT*NqINT + j*NqINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[ol[0].IEN_gb[cn[l] - 1][m] - 1][2]; //dz/dzeta
					}
					//The determinant of jacobian matrix for each node in the element
					ol[0].JACOB[m][i*NqINT*NqINT + j*NqINT + k] = xs2[m][0][0][i*NqINT*NqINT + j*NqINT + k] * (xs2[m][1][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][2][2][i*NqINT*NqINT + j*NqINT + k] - xs2[m][2][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][1][2][i*NqINT*NqINT + j*NqINT + k])
						- xs2[m][1][0][i*NqINT*NqINT + j*NqINT + k] * (xs2[m][0][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][2][2][i*NqINT*NqINT + j*NqINT + k] - xs2[m][2][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][0][2][i*NqINT*NqINT + j*NqINT + k])
						+ xs2[m][2][0][i*NqINT*NqINT + j*NqINT + k] * (xs2[m][0][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][1][2][i*NqINT*NqINT + j*NqINT + k] - xs2[m][1][1][i*NqINT*NqINT + j*NqINT + k] * xs2[m][0][2][i*NqINT*NqINT + j*NqINT + k]);
				}
			}
		}
	}

	//If mapping algorithm==2, we could make a 2D Jacobian for fluid nodal force (to sent to structure) integration.  
	ol[0].Jacob_2D = new double*[ol[0].FSNEL];
	for (i = 0; i < ol[0].FSNEL; i++) {
		ol[0].Jacob_2D[i] = new double[NINT*NINT];
	}
	//Use second order GLL integration (one order higher than the element order which is one in algorithm 2 since the high-order element is splitted)
	int Nq2 = 2;
	int Nq2INT = Nq2 + 1;
	bq = LOBATTO(2);
	gq = GLLQUAD(bq.Z, bq.WL, Nq2, !FEM);
	lg = LOCAL_GSHAPE(gq.S, LNA, Nq2INT);
	t.XS_2D = new double***[ol[0].FSNEL];
	for (i = 0; i < ol[0].FSNEL; i++) {
		t.XS_2D[i] = new double**[2];
		for (j = 0; j < 2; j++) {
			t.XS_2D[i][j] = new double*[2];
			for (k = 0; k < 2; k++) {
				t.XS_2D[i][j][k] = new double[Nq2INT*Nq2INT];
			}
		}
	}
	//We can determine the 2D Jacobian determinant for boundary condition here. 
	//Assuming linear geometric property mapping here. 
	//why it is xi and deta? Could be be xi and eta for example? 
	for (m = 0; m < ol[0].FSNEL; m++) { //loop through each element on wetted surface
		for (i = 0; i < Nq2INT; i++) {
			for (j = 0; j < Nq2INT; j++) {
				for (l = 0; l < 4; l++) { //4 nodes on the linear element
					/*
					t.XS_2D[m][0][0][i*2 + j] = t.XS_2D[m][0][0][i*2 + j] + lg.GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][0]; //dx/dxi
					t.XS_2D[m][1][0][i*2 + j] = t.XS_2D[m][1][0][i*2 + j] + lg.GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][1]; //dy/dxi
					t.XS_2D[m][0][1][i*2 + j] = t.XS_2D[m][0][1][i*2 + j] + lg.GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][0]; //dx/deta
					t.XS_2D[m][1][1][i*2 + j] = t.XS_2D[m][1][1][i*2 + j] + lg.GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][1]; //dy/deta
					*/
					//We need to figure out which surface corresponds to the wetted surface here. 
					t.XS_2D[m][0][0][i * Nq2INT + j] = t.XS_2D[m][0][0][i * Nq2INT + j] + lg.GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][];
					t.XS_2D[m][1][0][i * Nq2INT + j] = t.XS_2D[m][1][0][i * Nq2INT + j] + lg.GSHL_2D[0][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][];
					t.XS_2D[m][0][1][i * Nq2INT + j] = t.XS_2D[m][0][1][i * Nq2INT + j] + lg.GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][];
					t.XS_2D[m][1][1][i * Nq2INT + j] = t.XS_2D[m][1][1][i * Nq2INT + j] + lg.GSHL_2D[1][l][i][j] * GCOORD[ol[0].IEN_gb[l][m] - 1][];
				}
				ol[0].Jacob_2D[m][i * Nq2INT + j] = t.XS_2D[m][0][0][i * Nq2INT + j] * t.XS_2D[m][1][1][i * Nq2INT + j] - t.XS_2D[m][1][0][i * Nq2INT + j] * t.XS_2D[m][0][1][i * Nq2INT + j];
			}
		}
	}

	std::cout << " " << std::endl;
	return t;
}