//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

//PURPOSE: DETERMINES THE LOCAL SHAPE FUNCTION OF ORDER N AND ITS DERIVATIVE AT POINTS BY TRANSFORMING LEGENDRE POLYNOMIALS INTO CARDINAL FUNCTIONS(I.E. LAGRANGE INTERPOLATES)
//The original subroutine assumes that the number of interpolation points and quadrature points is the same. 
//Jan. 8 2018: The above problem is fixed to allow different number of interpolation and quadrature points. 

struct LOCAL_SHAPEstruct LOCAL_SHAPE(int*** LNA, int NQUAD) {
	int i, j, k, l, m, n;
	LOCAL_SHAPEstruct t;

	t.SHL = new double**[4];
	for (i = 0; i < 4; i++) {
		t.SHL[i] = new double*[NINT*NINT*NINT];
		for (j = 0; j < NINT*NINT*NINT; j++) {
			t.SHL[i][j] = new double[(NQUAD + 1)*(NQUAD + 1)*(NQUAD + 1)];
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			for (k = 0; k < (NQUAD + 1)*(NQUAD + 1)*(NQUAD + 1); k++) {
				t.SHL[i][j][k] = 0.0;
			}
		}
	}


	t.SHOD = new double**[2];
	for (i = 0; i < 2;i++) {
		t.SHOD[i] = new double*[NINT];
		for (j = 0; j < NINT;j++) {
			t.SHOD[i][j] = new double[NQUAD + 1];
		}
	}
	for (i = 0; i < 2;i++) {
		for (j = 0; j < NINT;j++) {
			for (k = 0; k < NQUAD + 1;k++) {
				t.SHOD[i][j][k] = 0.0;
			}
		}
	}

	if (FEM == 1) {
		LOBATTOstruct b;
		b = LOBATTO(N);
		GLLQUADstruct g;
		g = GLLQUAD(b.Z, b.WL, N, 0); //Obtain Gauss-Legendre point (not Lobatto point)
		double femp[NINT];  //local Lagrange interpolation nodes
		for (i = 0; i < NINT; i++) {
			femp[i] = -1 + (2 / N)*i;
		}
		double nom; double denom;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				nom = 1.0; denom = 1.0;
				for (m = 0; m < NINT; m++) { //loop through every interpolation point 
					if (m != i) { //loop through nominator and denominator in basis function expression
						nom *= (g.S[j] - femp[m]);
						denom *= (femp[i] - femp[m]);
					}
				}
				t.SHOD[0][i][j] = nom / denom;
			}
		}
		//the calculation is performed by the code "Lagrange_polynomial.m" which is in the code folder
		t.SHOD[1][0][0] = -0.500000000000000;
		t.SHOD[1][0][1] = -0.500000000000000;
		t.SHOD[1][1][0] = 0.500000000000000;
		t.SHOD[1][1][1] = 0.500000000000000;
	}
	else {
		//----------------BEGAIN COMPUTATION----------------------------//
		//USE LEGENDRE POLYNOMIAL TO EVALUATE LOCAL SHAPE FUNCTION AT EACH NODE
		LEGENDREstruct hi;
		LEGENDREstruct hq;
		hi = LEGENDRE(N);
		hq = LEGENDRE(NQUAD);
		LOBATTOstruct bi;
		LOBATTOstruct bq;
		bi = LOBATTO(N);
		bq = LOBATTO(NQUAD);
		GLLQUADstruct gi;
		GLLQUADstruct gq;
		gi = GLLQUAD(bi.Z, bi.WL, N, 1);
		gq = GLLQUAD(bq.Z, bq.WL, NQUAD, 1);
		//DETERMINE 1D (one D) LOCAL SHAPE FUNCTIONS 
		for (i = 0; i < NINT; i++) { //interpolation point
			for (j = 0; j < NQUAD + 1; j++) { //quadrature point
				/*
				if (j != i) {
					t.SHOD[0][i][j] = -((1.0 - pow(S[j], 2))*h.LN[1][j]) / (N*(N + 1)*h.LN[0][i] * (S[j] - S[i]));
				}
				else {
					t.SHOD[0][i][j] = 1.0;
				}
				*/
				if (i == 0 && j == 0) {
					t.SHOD[0][i][j] = 1.0;
				}
				else if (i == NINT - 1 && j == NQUAD) {
					t.SHOD[0][i][j] = 1.0;
				}
				else if (gq.S[j] == gi.S[i]) {
					t.SHOD[0][i][j] = 1.0;
				}
				else {
					t.SHOD[0][i][j] = -((1.0 - pow(gq.S[j], 2))*hq.LN[1][j]) / (N*(N + 1)*hi.LN[0][i] * (gq.S[j] - gi.S[i]));
				}
			}
		}

		//USE LEGENDRE POLYNOMIAL CARDINAL FUNCTION TO EVALUATE LOCAL SHAPE FUNCTION DERIVATIVE AT EACH NODE
		//reference: Klenow thesis page 35
		for (l = 0; l < NINT; l++) { //interpolation point
			for (k = 0; k < NQUAD + 1; k++) { //quadrature point
				if (l == 0 && k == 0) {
					t.SHOD[1][l][k] = -double((N + 1)*N) / 4.0;
				}
				else if (l == NINT - 1 && k == NQUAD) {
					t.SHOD[1][l][k] = double((N + 1)*N) / 4.0;
				}
				else if (gq.S[k] == gi.S[l]) {
					t.SHOD[1][l][k] = 0.0;
				}
				else {
					t.SHOD[1][l][k] = (hq.LN[0][k] / hi.LN[0][l])*(1.0 / (gq.S[k] - gi.S[l]));
				}
			}
		}
	}

	//DETERMINE 3D LOCAL SHAPE FUNCTIONS
	//shape function discrete value is evaluated at integration points (used for Gauss-Lobatto integration later on)
	for (i = 0; i < NINT; i++) {   //Xi 
		for (j = 0; j < NINT; j++) {  //Eta
			for (k = 0; k < NINT; k++) { //Zeta
				//i j k are interpolation points and l m n are quadrature points 
				for (l = 0; l < NQUAD + 1; l++) { //the internal nodes on Xi 
					for (m = 0; m < NQUAD + 1; m++) { //the internal nodes on Eta 
						for (n = 0; n < NQUAD + 1; n++) { //the internal nodes on Zeta 
							t.SHL[3][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + n] = t.SHOD[0][i][l] * t.SHOD[0][j][m] * t.SHOD[0][k][n];
							//3D Nth ORDER SHAPE FUNCTION
							t.SHL[0][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + n] = t.SHOD[1][i][l] * t.SHOD[0][j][m] * t.SHOD[0][k][n];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Xi 
							t.SHL[1][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + n] = t.SHOD[0][i][l] * t.SHOD[1][j][m] * t.SHOD[0][k][n];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Eta 
							t.SHL[2][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + n] = t.SHOD[0][i][l] * t.SHOD[0][j][m] * t.SHOD[1][k][n];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Zeta
						}
					}
				}
			}
		}
	}

	//check point: if the shape function value at the corresponding point is 0 
	//std::cout << "shape function value on the corresponding point is: " << t.SHL[3][LNA[1][1][1] - 1][1][1][1] << std::endl;

	return t;
}
