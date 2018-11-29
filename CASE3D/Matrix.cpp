//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <ctime>
#include <omp.h>

struct MATRIXstruct MATRIX(int NEL, int NNODE, double***SHL, double*W, int**IEN, int***LNA, double* JACOB_tet, double*** SHG_tet, double** GCOORD) {
	MATRIXstruct t;
	double QFUNC;  //CAPACITANCE MATRIX SUM VARIABLE (mass matrix)
	//double HFUNC;  //REACTANCE MATRIX SUM VARIABLE (stiffness matrix)
	int i, j, k, l, m, n, e, u, v;
	//double duration;
	std::clock_t start;

	int elenode3D = 0;
	if (element_type == 0) { //hex element
		elenode3D = NINT*NINT*NINT;
	}
	if (element_type == 1) { //tet element
		if (N == 1) {
			elenode3D = 4;
		}
		else {
			std::cout << "High-order tet element is not supported yet" << std::endl;
			system("PAUSE ");
		}
	}

	double** QMASTER;
	QMASTER = new double*[elenode3D]; //element matrix  
	for (i = 0; i < elenode3D; i++) {
		QMASTER[i] = new double[elenode3D];
	}
	t.Q = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		t.Q[i] = 0.0;
	}

	/*
	t.HMASTER = new double**[NEL];
	for (i = 0; i < NEL; i++) {
		t.HMASTER[i] = new double*[elenode3D];
		for (j = 0; j < elenode3D; j++) {
			t.HMASTER[i][j] = new double[elenode3D];
		}
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < elenode3D; j++) {
			for (k = 0; k < elenode3D; k++) {
				t.HMASTER[i][j][k] = 0.0;
			}
		}
	}
	*/
	t.HMASTER = new double*[NEL];
	for (i = 0; i < NEL; i++) {
		t.HMASTER[i] = new double[elenode3D*(elenode3D + 1) / 2];
	}
	for (i = 0; i < NEL; i++) {
		for (j = 0; j < elenode3D*(elenode3D + 1) / 2; j++) {
			t.HMASTER[i][j] = 0.0;
		}
	}

	double* W_new;
	W_new = new double[NINT*NINT*NINT];
	for (k = 0; k < NINT; k++) {       //k l m are for four points
		for (l = 0; l < NINT; l++) {
			for (m = 0; m < NINT; m++) {
				W_new[k*NINT*NINT + l*NINT + m] = W[k] * W[l] * W[m];
			}
		}
	}
	double* FUNC = new double[elenode3D];
	double* LMAXi = new double[NEL];
	double SUMH; 
	double HFUNC; 
	double* JACOB; 
	JACOB = new double[NINT*NINT*NINT];


	double** XS; //define XS to be a local variable which will be cleaned at the end of function
	double*** SHG;
	int*cn;
	LOBATTOstruct bq;
	GLLQUADstruct gq;
	LOCAL_GSHAPEstruct lg;
	if (element_type == 0) {
		cn = new int[8];
		cn[0] = LNA[0][0][0]; cn[1] = LNA[N][0][0];
		cn[2] = LNA[N][N][0]; cn[3] = LNA[0][N][0];
		cn[4] = LNA[0][0][N]; cn[5] = LNA[N][0][N];
		cn[6] = LNA[N][N][N]; cn[7] = LNA[0][N][N];
		SHG = new double**[4];
		for (i = 0; i < 4; i++) {
			SHG[i] = new double*[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				SHG[i][j] = new double[NINT*NINT*NINT];
			}
		}
		for (i = 0; i < 4; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				for (k = 0; k < NINT*NINT*NINT; k++) {
					SHG[i][j][k] = 0.0;
				}
			}
		}
		XS = new double*[9];
		for (i = 0; i < 9; i++) {
			XS[i] = new double[NINT*NINT*NINT];
		}
		for (i = 0; i < 9; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				XS[i][j] = 0.0;
			}
		}
		bq = LOBATTO(N);
		gq = GLLQUAD(bq.Z, bq.WL, N, !FEM);
		lg = LOCAL_GSHAPE(gq.S, LNA, NINT);
	}
	//Start looping through all the elements
	for (e = 0; e < NEL; e++) {
		//Derive the HMASTER first
		if (element_type == 0) {
			//derive the SHG first
			//derive XS first
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = 0;
						for (l = 0; l < 8; l++) {
							XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/dxi
							XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/deta 
							XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/dzeta  
							XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/dxi 
							XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/deta
							XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/dzeta 
							XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/dxi 
							XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/deta 
							XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/dzeta
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
							SHG[0][i][j*NINT*NINT + k*NINT + l] = ((XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dx
							SHG[1][i][j*NINT*NINT + k*NINT + l] = ((-XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 2][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 2][j*NINT*NINT + k*NINT + l] + XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 2][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dy
							SHG[2][i][j*NINT*NINT + k*NINT + l] = ((XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] - XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[0][i][j*NINT*NINT + k*NINT + l] + (-XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[2 * 3 + 1][j*NINT*NINT + k*NINT + l] + XS[2 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[1][i][j*NINT*NINT + k*NINT + l] + (XS[0 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[1 * 3 + 1][j*NINT*NINT + k*NINT + l] - XS[1 * 3 + 0][j*NINT*NINT + k*NINT + l] * XS[0 * 3 + 1][j*NINT*NINT + k*NINT + l])*SHL[2][i][j*NINT*NINT + k*NINT + l]) / JACOB[j*NINT*NINT + k*NINT + l]; //d(phi)/dz
							SHG[3][i][j*NINT*NINT + k*NINT + l] = SHL[3][i][j*NINT*NINT + k*NINT + l]; //phi=phi 
						}
					}
				}
			}
			for (int i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element 
				for (int j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
					if (i <= j) { //store only the upper triangle elements
						for (int k = 0; k < NINT*NINT*NINT; k++) {   //k l m are for four points
							//See https://stackoverflow.com/questions/9039189/make-efficient-the-copy-of-symmetric-matrix-in-c-sharp/9040526#9040526 for the storing algorithm
							t.HMASTER[e][i*elenode3D - i*(i + 1) / 2 + j] += W_new[k] * JACOB[k] * (SHG[2][i][k] * SHG[2][j][k] + SHG[1][i][k] * SHG[1][j][k] + SHG[0][i][k] * SHG[0][j][k]);
						}
					}
				}
			}
		}
		else if (element_type == 1) {
			//Currently only for linear tetrahedral element. 
			//Source: Pozrikidis, C Finite and Spectral Element Methods Using MATLAB Chapter 8
			double Vol = (1.0 / 6.0)*JACOB_tet[e];
			for (int i = 0; i < 4; i++) {     //totally NINT*NINT*NINT points for an element 
				for (int j = 0; j < 4; j++) { //the transpose of shape function
					if (i <= j) { //store only the upper triangle elements
						t.HMASTER[e][i*elenode3D - i*(i + 1) / 2 + j] = Vol * (SHG_tet[e][2][i] * SHG_tet[e][2][j] + SHG_tet[e][1][i] * SHG_tet[e][1][j] + SHG_tet[e][0][i] * SHG_tet[e][0][j]);
					}
				}
			}
		}
			//derive the QMASTER
		if (element_type == 0) {
			for (u = 0; u < NINT*NINT*NINT; u++) {
				for (v = 0; v < NINT*NINT*NINT; v++) {
					QMASTER[u][v] = 0.0;
				}
			}
			if (FEM == 1) {
				for (i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element
					for (j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
						for (k = 0; k < NINT; k++) {       //k l m are for four points
							for (l = 0; l < NINT; l++) {
								for (m = 0; m < NINT; m++) {
									//PERFORM GLL INTERGRATION FOR CAPACITANCE MATRIX
									QFUNC = JACOB[k*NINT*NINT + l*NINT + m] * (SHL[3][i][k*NINT*NINT + l*NINT + m] * SHL[3][j][k*NINT*NINT + l*NINT + m]);
									QMASTER[i][i] += W[k] * W[l] * W[m] * QFUNC;  //MULTIDIMENSIONAL GAUSSIAN QUADRATURE: INTEGRATE INNER INTERGRAL FIRST, THEN OUTER INTEGRAL
								}
							} //check point: whether the QMASTER matrix is diagonal, if not wrong
						}
					}
				}
				//assemble the local mass matrix 
				for (i = 0; i < NINT*NINT*NINT; i++) {
					t.Q[IEN[i][e] - 1] += QMASTER[i][i];
				}
			}
			else {
				for (i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element
					for (j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
						if (i == j) { //diagonal terms
							for (k = 0; k < NINT; k++) {       //k l m are for four points
								for (l = 0; l < NINT; l++) {
									for (m = 0; m < NINT; m++) {
										//PERFORM GLL INTERGRATION FOR CAPACITANCE MATRIX
										QFUNC = JACOB[k*NINT*NINT + l*NINT + m] * (SHL[3][i][k*NINT*NINT + l*NINT + m] * SHL[3][j][k*NINT*NINT + l*NINT + m]);
										QMASTER[i][j] = QMASTER[i][j] + W[k] * W[l] * W[m] * QFUNC;  //MULTIDIMENSIONAL GAUSSIAN QUADRATURE: INTEGRATE INNER INTERGRAL FIRST, THEN OUTER INTEGRAL
									}
								} //check point: whether the QMASTER matrix is diagonal, if not wrong
							}
							//ASSEMBLE GLOBAL CAPACITANCE MATRIX
							t.Q[IEN[i][e] - 1] = t.Q[IEN[i][e] - 1] + QMASTER[i][j];
						}
					}
				}
			}
		} //finish hexahedral mesh definition
		else if (element_type == 1) {
			//Currently only for linear tetrahedral element. 
			//Source: Pozrikidis, C Finite and Spectral Element Methods Using MATLAB Chapter 8
			double Vol = (1.0 / 6.0)*JACOB_tet[e];
			//Mass matrix
			for (u = 0; u < 4; u++) {
				for (v = 0; v < 4; v++) {
					QMASTER[u][v] = 0.0;
				}
			}
			for (i = 0; i < 4; i++) {     //totally NINT*NINT*NINT points for an element
				for (j = 0; j < 4; j++) { //the transpose of shape function
					if (i == j) { //diagonal terms
						QMASTER[i][j] = Vol*(1.0 / 4.0)*1.0;  //the mass matrix is directly diagonalized (Pozrikidis 8.2.38)
						t.Q[IEN[i][e] - 1] = t.Q[IEN[i][e] - 1] + QMASTER[i][j];
					}
				}
			}
		} //finish tetrahedral mesh definition

		//Derive the smallest eigenvalue to determine the criteria time step size
		for (i = 0; i < elenode3D; i++) { //loop through each line
			SUMH = 0.0;
			for (j = 0; j < elenode3D; j++) { //Column sum (Since HMASTER is symmetric, row sum and column sum would actually yield the same result)
				if (i < j) {
					HFUNC = abs(t.HMASTER[e][i*elenode3D - i*(i + 1) / 2 + j]);
					SUMH = SUMH + HFUNC;
				}
				else if (i > j) {
					HFUNC = abs(t.HMASTER[e][j*elenode3D - j*(j + 1) / 2 + i]);
					SUMH = SUMH + HFUNC;
				}
			}
			//add all element to diagonal element line by line
			FUNC[i] = (t.HMASTER[e][i*elenode3D - i*(i + 1) / 2 + i] + SUMH) / QMASTER[i][i]; //the HMASTER here doesn't have the absolute sign. 
		}
		//DERIVE THE MAX EIGEN VALUE IN THE ELEMENT e
		LMAXi[e] = FUNC[0];
		for (j = 1; j < elenode3D; j++) {
			if (FUNC[j] > LMAXi[e]) {
				LMAXi[e] = FUNC[j];
			}
		}
	} //finish looping all the elements
		//Obtain the largest eigenvalue of the elements
	t.LMAX = LMAXi[0];
	for (e = 1; e < NEL; e++) {
		if (LMAXi[e] > t.LMAX) {
			t.LMAX = LMAXi[e];
		}
	}
	//clean dynamic variables
	if (element_type == 0) {
		delete[] cn;
		for (i = 0; i < 9; i++) {
			delete[] XS[i];
		}
		delete[] XS;
		for (i = 0; i < 4; i++) {
			for (j = 0; j < NINT*NINT*NINT; j++) {
				delete[] SHG[i][j];
			}
			delete[] SHG[i];
		}
		delete[] SHG;
	}
	for (i = 0; i < elenode3D; i++) {
		delete[] QMASTER[i];
	}
	delete[] QMASTER;

	/*
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			std::cout << t.HMASTER[1][i][j] << " ";
		}
		std::cout << std::endl;
	}
	*/

	if (tensorfactorization == 1) {

		//If tensor product factorization is used, HMASTER is not used after the maximum eigenvalue is determined
		for (i = 0; i < NEL; i++) {
			delete[] t.HMASTER[i];
		}
		delete[] t.HMASTER;

		//tensors for the evaluation of stiffness terms
		t.gamma = new double***[NINT];
		for (i = 0; i < NINT; i++) {
			t.gamma[i] = new double**[NINT];
			for (j = 0; j < NINT; j++) {
				t.gamma[i][j] = new double*[NINT];
				for (k = 0; k < NINT; k++) {
					t.gamma[i][j][k] = new double[3];
				}
			}
		}
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < 3; l++) {
						t.gamma[i][j][k][l] = 0.0;
					}
				}
			}
		}
		/*
		t.gamman = new double*[NINT*NINT*NINT];
		for (i = 0; i < NINT*NINT*NINT; i++) {
			t.gamman[i] = new double[3];
		}
		*/
		/*
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < 3; j++) {
				t.gamman[i][j] = 0.0;
			}
		}
		*/
		for (i = 0; i < NINT*NINT*NINT * 3; i++) {
			t.gamman[i] = 0.0;
		}
		t.gamma_t = new double***[NINT];
		for (i = 0; i < NINT; i++) {
			t.gamma_t[i] = new double**[NINT];
			for (j = 0; j < NINT; j++) {
				t.gamma_t[i][j] = new double*[NINT];
				for (k = 0; k < NINT; k++) {
					t.gamma_t[i][j][k] = new double[3];
				}
			}
		}
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < 3; l++) {
						t.gamma_t[i][j][k][l] = 0.0;
					}
				}
			}
		}
		/*
		t.gamma_tn = new double*[NINT*NINT*NINT];
		for (i = 0; i < NINT*NINT*NINT; i++) {
			t.gamma_tn[i] = new double[3];
		}
		*/
		/*
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < 3; j++) {
				t.gamma_tn[i][j] = 0.0;
			}
		}
		*/
		for (i = 0; i < NINT*NINT*NINT * 3; i++) {
			t.gamma_tn[i] = 0.0;
		}

		//double*****t.G;
		t.G = new double****[3];
		for (i = 0; i < 3; i++) {
			t.G[i] = new double***[3];
			for (j = 0; j < 3; j++) {
				t.G[i][j] = new double**[NINT];
				for (k = 0; k < NINT; k++) {
					t.G[i][j][k] = new double*[NINT];
					for (l = 0; l < NINT; l++) {
						t.G[i][j][k][l] = new double[NINT];
					}
				}
			}
		}

		t.Gn = new double***[NEL];
		for (i = 0; i < NEL; i++) {
			t.Gn[i] = new double**[3];
			for (j = 0; j < 3; j++) {
				t.Gn[i][j] = new double*[3];
				for (k = 0; k < 3; k++) {
					t.Gn[i][j][k] = new double[NINT*NINT*NINT];
				}
			}
		}

		//need to initialize here!!!
		Eigen::MatrixXd JB(3, 3);
		Eigen::MatrixXd ga(3, 3);
		//Eigen::Matrix3d JB;
		//Eigen::Matrix3d ga;
		//std::clock_t start;
		//start = std::clock();

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
		int*cn = new int[8];
		cn[0] = LNA[0][0][0]; cn[1] = LNA[N][0][0];
		cn[2] = LNA[N][N][0]; cn[3] = LNA[0][N][0];
		cn[4] = LNA[0][0][N]; cn[5] = LNA[N][0][N];
		cn[6] = LNA[N][N][N]; cn[7] = LNA[0][N][N];
		LOBATTOstruct bq;
		GLLQUADstruct gq;
		LOCAL_GSHAPEstruct lg;
		bq = LOBATTO(N);
		gq = GLLQUAD(bq.Z, bq.WL, N, !FEM);
		lg = LOCAL_GSHAPE(gq.S, LNA, NINT);

		int ct;
		for (e = 0; e < NEL; e++) {
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = 0;
						XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = 0; XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = 0;
						for (l = 0; l < 8; l++) {
							XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/dxi
							XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/deta 
							XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][0]; //dx/dzeta  
							XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/dxi 
							XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/deta
							XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][1]; //dy/dzeta 
							XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[0][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/dxi 
							XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[1][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/deta 
							XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] = XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] + lg.GSHL[2][l][i][j][k] * GCOORD[IEN[cn[l] - 1][e] - 1][2]; //dz/dzeta
						}
						JACOB[i*NINT*NINT + j*NINT + k] = XS[3 * 0 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k])
							- XS[3 * 1 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 2 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 2 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k])
							+ XS[3 * 2 + 0][i*NINT*NINT + j*NINT + k] * (XS[3 * 0 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 1 + 2][i*NINT*NINT + j*NINT + k] - XS[3 * 1 + 1][i*NINT*NINT + j*NINT + k] * XS[3 * 0 + 2][i*NINT*NINT + j*NINT + k]);
					}
				}
			}
			//std::cout << e << std::endl; 
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						//
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								JB(l, m) = XS[l * 3 + m][i*NINT*NINT + j*NINT + k]; //Jacobian matrix
							}
						}
						//std::cout << JB << std::endl;
						//std::cout << JB.inverse() << std::endl;
						ga = (JB.inverse().transpose())*JB.inverse(); //3*3 matrix
						//std::cout << ga << std::endl;
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								t.G[l][m][i][j][k] = ga(l, m) * W_new[i*NINT*NINT + j*NINT + k] * JACOB[i*NINT*NINT + j*NINT + k];
								//std::cout << t.G[l][m][i][j][k] << std::endl;
							}
							//std::cout<<std::endl;
						}
					}
				}
			}

			//JB is the jacobian matrix and JACOB is the determinant of jacobian matrix
			ct = 0;
			//#pragma omp parallel for num_threads(6)
			for (int k = 0; k < NINT*NINT*NINT; k++) {
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						JB(l, m) = XS[l * 3 + m][k]; //Jacobian matrix
					}
				}
				ga = (JB.inverse().transpose())*JB.inverse(); //3*3 matrix
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						t.Gn[e][l][m][ct] = ga(l, m) * W_new[k] * JACOB[k];
					}
				}
				ct += 1;
			}
		}
		for (i = 0; i < 9; i++) {
			delete[] XS[i];
		}
		delete[] XS;
		delete[] cn; 
	}
	delete[] LMAXi;
	delete[] FUNC;
	delete[] JACOB; 
	delete[] W_new; 

	std::cout << " " << std::endl; 

	return t;
}