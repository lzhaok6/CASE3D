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

struct MATRIXstruct MATRIX(int NEL, int NNODE, double***SHL, double*W, int**IEN, int***LNA, double***XS, double****SHG, double**JACOB, double* JACOB_tet, double*** SHG_tet) {
	MATRIXstruct t;
	double QFUNC;  //CAPACITANCE MATRIX SUM VARIABLE (mass matrix)
	double HFUNC;  //REACTANCE MATRIX SUM VARIABLE (stiffness matrix)
	int i, j, k, l, m, n, e, u, v;
	double duration;
	std::clock_t start;

	t.QMASTER = new double**[NEL];
	if (element_type == 0) {
		for (i = 0; i < NEL; i++) {
			t.QMASTER[i] = new double*[NINT*NINT*NINT]; //element matrix  
			for (j = 0; j < NINT*NINT*NINT; j++) {
				t.QMASTER[i][j] = new double[NINT*NINT*NINT];
			}
		}
	}
	if (element_type == 1) {
		for (i = 0; i < NEL; i++) {
			t.QMASTER[i] = new double*[4]; //element matrix  
			for (j = 0; j < 4; j++) {
				t.QMASTER[i][j] = new double[4];
			}
		}
	}

	t.HMASTER = new double**[NEL]; //element matrix
	if (element_type == 0) {
		for (i = 0; i < NEL; i++) {
			t.HMASTER[i] = new double*[NINT*NINT*NINT];
			for (j = 0; j < NINT*NINT*NINT; j++) {
				t.HMASTER[i][j] = new double[NINT*NINT*NINT];
			}
		}
	}
	if (element_type == 1) {
		for (i = 0; i < NEL; i++) {
			t.HMASTER[i] = new double*[4];
			for (j = 0; j < 4; j++) {
				t.HMASTER[i][j] = new double[4];
			}
		}
	}

	t.Q = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		t.Q[i] = 0.0;
	}

	if (element_type == 0) {
		for (m = 0; m < NEL; m++) {
			for (i = 0; i < NINT*NINT*NINT; i++) {
				for (j = 0; j < NINT*NINT*NINT; j++) {
					t.HMASTER[m][i][j] = 0.0;
				}
			}
		}
	}
	if (element_type == 1) {
		for (m = 0; m < NEL; m++) {
			for (i = 0; i < 4; i++) {
				for (j = 0; j < 4; j++) {
					t.HMASTER[m][i][j] = 0.0;
				}
			}
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

	//----------------BEGIN MATRIX COMP.--------------------//
	if (element_type == 0) {
		if (FEM == 1) { //for FEM element, we apply row-sum technique to make the mass matrix diagonal
			for (e = 0; e < NEL; e++) {
				for (u = 0; u < NINT*NINT*NINT; u++) {
					for (v = 0; v < NINT*NINT*NINT; v++) {
						t.QMASTER[e][u][v] = 0.0;
					}
				}
				for (i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element
					for (j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
						for (k = 0; k < NINT; k++) {       //k l m are for four points
							for (l = 0; l < NINT; l++) {
								for (m = 0; m < NINT; m++) {
									//PERFORM GLL INTERGRATION FOR CAPACITANCE MATRIX
									QFUNC = JACOB[e][k*NINT*NINT + l*NINT + m] * (SHL[3][i][k*NINT*NINT + l*NINT + m] * SHL[3][j][k*NINT*NINT + l*NINT + m]);
									t.QMASTER[e][i][i] += W[k] * W[l] * W[m] * QFUNC;  //MULTIDIMENSIONAL GAUSSIAN QUADRATURE: INTEGRATE INNER INTERGRAL FIRST, THEN OUTER INTEGRAL
								}
							} //check point: whether the QMASTER matrix is diagonal, if not wrong
						}
					}
				}
				//assemble the local mass matrix 
				for (i = 0; i < NINT*NINT*NINT; i++) {
					t.Q[IEN[i][e] - 1] += t.QMASTER[e][i][i];
				}
			}
		}
		else {
			for (e = 0; e < NEL; e++) {
				//std::cout << e << std::endl;
				for (u = 0; u < NINT*NINT*NINT; u++) {
					for (v = 0; v < NINT*NINT*NINT; v++) {
						t.QMASTER[e][u][v] = 0.0;
					}
				}
				for (i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element
					for (j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
						if (i == j) { //diagonal terms
							for (k = 0; k < NINT; k++) {       //k l m are for four points
								for (l = 0; l < NINT; l++) {
									for (m = 0; m < NINT; m++) {
										//PERFORM GLL INTERGRATION FOR CAPACITANCE MATRIX
										QFUNC = JACOB[e][k*NINT*NINT + l*NINT + m] * (SHL[3][i][k*NINT*NINT + l*NINT + m] * SHL[3][j][k*NINT*NINT + l*NINT + m]);
										t.QMASTER[e][i][j] = t.QMASTER[e][i][j] + W[k] * W[l] * W[m] * QFUNC;  //MULTIDIMENSIONAL GAUSSIAN QUADRATURE: INTEGRATE INNER INTERGRAL FIRST, THEN OUTER INTEGRAL
									}
								} //check point: whether the QMASTER matrix is diagonal, if not wrong
							}
							//ASSEMBLE GLOBAL CAPACITANCE MATRIX
							t.Q[IEN[i][e] - 1] = t.Q[IEN[i][e] - 1] + t.QMASTER[e][i][j];
						}
					}
				}
				if (e % 1000 == 0) {
					std::cout << e << std::endl; //output which line is being read
				}
			}
		}

		//start = std::clock();
		/*
		for (e = 0; e < NEL; e++) {
			//std::cout << e << std::endl;
			for (i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element
				for (j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
					for (k = 0; k < NINT; k++) {       //k l m are for four points
						for (l = 0; l < NINT; l++) {
							for (m = 0; m < NINT; m++) {
								HFUNC = JACOB[e][k*NINT*NINT + l*NINT + m]*(SHG[e][2][i][k*NINT*NINT+l*NINT+m] * SHG[e][2][j][k*NINT*NINT + l*NINT + m] + SHG[e][1][i][k*NINT*NINT + l*NINT + m] * SHG[e][1][j][k*NINT*NINT + l*NINT + m] + SHG[e][0][i][k*NINT*NINT + l*NINT + m] * SHG[e][0][j][k*NINT*NINT + l*NINT + m]);
								//ASSEMBLE MASTER REACTANCE MATRIX
								t.HMASTER[e][i][j] = t.HMASTER[e][i][j] + W[k] * W[l] * W[m] * HFUNC;
								//QMASTER is diagonal because of the orthogonality of polynomial
								//HMASTER is not diagonal because the polynomial is differentiated
							}
						}
					}
				}
			}
			if (e % 1000 == 0) {
				std::cout << e << std::endl; //output which line is being read
			}
		}
		*/
		//about 25% faster than the algorithm above
		#pragma omp parallel for num_threads(6)
		for (int e = 0; e < NEL; e++) {
			//#pragma omp parallel for num_threads(6)
			for (int i = 0; i < NINT*NINT*NINT; i++) {     //totally NINT*NINT*NINT points for an element 
				//#pragma omp parallel for num_threads(6)
				for (int j = 0; j < NINT*NINT*NINT; j++) { //the transpose of shape function
					//#pragma omp parallel for num_threads(6)
					for (int k = 0; k < NINT*NINT*NINT; k++) {       //k l m are for four points
						//HFUNC = JACOB[e][k] * (SHG[e][2][i][k] * SHG[e][2][j][k] + SHG[e][1][i][k] * SHG[e][1][j][k] + SHG[e][0][i][k] * SHG[e][0][j][k]);
						//ASSEMBLE MASTER REACTANCE MATRIX
						//t.HMASTER[e][i][j] = t.HMASTER[e][i][j] + W_new[k] * HFUNC;
						t.HMASTER[e][i][j] = t.HMASTER[e][i][j] + W_new[k] * JACOB[e][k] * (SHG[e][2][i][k] * SHG[e][2][j][k] + SHG[e][1][i][k] * SHG[e][1][j][k] + SHG[e][0][i][k] * SHG[e][0][j][k]);
						//QMASTER is diagonal because of the orthogonality of polynomial 
						//HMASTER is not diagonal because the polynomial is differentiated 
					}
				}
			}
		}
		//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
		//std::cout << "total CPU time (ms): " << duration << std::endl;
		//std::cout << " " << std::endl;
	}

	if (element_type == 1) { //tetrahedral element mass (capacitance) and diffusion (reactance) matrix
		//Currently only for linear tetrahedral element. 
		//Source: Pozrikidis, C Finite and Spectral Element Methods Using MATLAB Chapter 8
		for (e = 0; e < NEL; e++) {
			double Vol = (1.0 / 6.0)*JACOB_tet[e];
			//Mass matrix
			for (u = 0; u < 4; u++) {
				for (v = 0; v < 4; v++) {
					t.QMASTER[e][u][v] = 0.0;
				}
			}
			for (i = 0; i < 4; i++) {     //totally NINT*NINT*NINT points for an element
				for (j = 0; j < 4; j++) { //the transpose of shape function
					if (i == j) { //diagonal terms
						t.QMASTER[e][i][j] = Vol*(1.0 / 4.0)*1.0;  //the mass matrix is directly diagonalized
						t.Q[IEN[i][e] - 1] = t.Q[IEN[i][e] - 1] + t.QMASTER[e][i][j];
					}
				}
			}
			//Diffusion matrix
			for (int i = 0; i < 4; i++) {     //totally NINT*NINT*NINT points for an element 
				for (int j = 0; j < 4; j++) { //the transpose of shape function
					t.HMASTER[e][i][j] = Vol * (SHG_tet[e][2][i] * SHG_tet[e][2][j] + SHG_tet[e][1][i] * SHG_tet[e][1][j] + SHG_tet[e][0][i] * SHG_tet[e][0][j]);
				}
			}
		}
	}
	
	if (tensorfactorization == 1) {
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
		//std::clock_t start;
		//start = std::clock();
		for (e = 0; e < NEL; e++) {
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						//
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								JB(l, m) = XS[e][l * 3 + m][i*NINT*NINT + j*NINT + k]; //Jacobian matrix
							}
						}
						//std::cout << JB << std::endl;
						//std::cout << JB.inverse() << std::endl;
						ga = (JB.inverse().transpose())*JB.inverse(); //3*3 matrix
																	  //std::cout << ga << std::endl;
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								t.G[l][m][i][j][k] = ga(l, m) * W_new[i*NINT*NINT + j*NINT + k] * JACOB[e][i*NINT*NINT + j*NINT + k];
								//std::cout << t.G[l][m][i][j][k] << std::endl;
							}
							//std::cout<<std::endl;
						}
					}
				}
			}
		}
		//JB is the jacobian matrix and JACOB is the determinant of jacobian matrix
		
		/*
		int ct;
		for (e = 0; e < NEL; e++) {
			//std::cout << e << std::endl;
			ct = 0;
			for (i = 0; i < NINT; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								JB(l, m) = XS[e][l][m][i*NINT*NINT + j*NINT + k]; //Jacobian matrix
							}
						}
						//std::cout << JB << std::endl;
						//std::cout << JB.inverse() << std::endl;
						ga = (JB.inverse().transpose())*JB.inverse(); //3*3 matrix
																	  //std::cout << ga << std::endl;
						for (l = 0; l < 3; l++) {
							for (m = 0; m < 3; m++) {
								t.Gn[e][l][m][ct] = ga(l, m) * W[i] * W[j] * W[k] * JACOB[e][i*NINT*NINT + j*NINT + k];
								//std::cout << t.G[l][m][i][j][k] << std::endl;
							}
							//std::cout<<std::endl;
						}
						ct += 1;
					}
				}
			}
			if (e % 1000 == 0) {
				std::cout << e << std::endl; //output which line is being read
			}
			//JB is the jacobian matrix and JACOB is the determinant of jacobian matrix
		}
		*/

		int ct;
		//start = std::clock();
		//#pragma omp parallel for num_threads(6)
		for (int e = 0; e < NEL; e++) {
			ct = 0;
			for (int k = 0; k < NINT*NINT*NINT; k++) {
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						JB(l, m) = XS[e][l * 3 + m][k]; //Jacobian matrix
					}
				}
				ga = (JB.inverse().transpose())*JB.inverse(); //3*3 matrix
				for (int l = 0; l < 3; l++) {
					for (int m = 0; m < 3; m++) {
						t.Gn[e][l][m][ct] = ga(l, m) * W_new[k] * JACOB[e][k];
					}
				}
				ct += 1;
			}
			/*
			if (e % 1000 == 0) {
				std::cout << e << std::endl; //output which line is being read
			}
			*/
			//JB is the jacobian matrix and JACOB is the determinant of jacobian matrix 
		}
	}

	//duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;
	//std::cout << "total CPU time (ms): " << duration << std::endl;
	std::cout << " " << std::endl;
	return t;
}