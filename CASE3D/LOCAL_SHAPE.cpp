//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

//PURPOSE: DETERMINES THE LOCAL SHAPE FUNCTION OF ORDER n AND ITS DERIVATIVE AT POINTS BY TRANSFORMING LEGENDRE POLYNOMIALS INTO CARDINAL FUNCTIONS(I.E. LAGRANGE INTERPOLATES)
//The original subroutine assumes that the number of interpolation points and quadrature points is the same. 
//Jan. 8 2018: The above problem is fixed to allow different number of interpolation and quadrature points. 

struct LOCAL_SHAPEstruct LOCAL_SHAPE(int*** LNA, int n, int NQUAD) {
	int i, j, k, l, m, o;
	LOCAL_SHAPEstruct t;

	int nint = n + 1; 
	t.SHL = new double**[4];
	for (i = 0; i < 4; i++) {
		t.SHL[i] = new double*[nint*nint*nint];
		for (j = 0; j < nint*nint*nint; j++) {
			t.SHL[i][j] = new double[(NQUAD + 1)*(NQUAD + 1)*(NQUAD + 1)];
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < nint*nint*nint; j++) {
			for (k = 0; k < (NQUAD + 1)*(NQUAD + 1)*(NQUAD + 1); k++) {
				t.SHL[i][j][k] = 0.0;
			}
		}
	}

	t.SHOD = new double**[2];
	for (i = 0; i < 2;i++) {
		t.SHOD[i] = new double*[nint];
		for (j = 0; j < nint;j++) {
			t.SHOD[i][j] = new double[NQUAD + 1];
		}
	}
	for (i = 0; i < 2;i++) {
		for (j = 0; j < nint;j++) {
			for (k = 0; k < NQUAD + 1;k++) {
				t.SHOD[i][j][k] = 0.0;
			}
		}
	}

	if (FEM == 1) {
		LOBATTOstruct b;
		b = LOBATTO(n);
		GLLQUADstruct g;
		g = GLLQUAD(b.Z, b.WL, n, 0); //Obtain Gauss-Legendre point (not Lobatto point)
		double* femp; 
		femp = new double[nint];
		//double femp[nint];  //local Lagrange interpolation nodes
		for (i = 0; i < nint; i++) {
			femp[i] = -1 + (2 / n)*i;
		}
		double nom; double denom;
		for (i = 0; i < nint; i++) {
			for (j = 0; j < nint; j++) {
				nom = 1.0; denom = 1.0;
				for (m = 0; m < nint; m++) { //loop through every interpolation point 
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
		hi = LEGENDRE(n);
		hq = LEGENDRE(NQUAD);
		LOBATTOstruct bi;
		LOBATTOstruct bq;
		bi = LOBATTO(n);
		bq = LOBATTO(NQUAD);
		GLLQUADstruct gi;
		GLLQUADstruct gq;
		gi = GLLQUAD(bi.Z, bi.WL, n, 1);
		gq = GLLQUAD(bq.Z, bq.WL, NQUAD, 1);
		//DETERMINE 1D (one D) LOCAL SHAPE FUNCTIONS 
		for (i = 0; i < nint; i++) { //interpolation point
			for (j = 0; j < NQUAD + 1; j++) { //quadrature point
				/*
				if (j != i) {
					t.SHOD[0][i][j] = -((1.0 - pow(S[j], 2))*h.LN[1][j]) / (n*(n + 1)*h.LN[0][i] * (S[j] - S[i]));
				}
				else {
					t.SHOD[0][i][j] = 1.0;
				}
				*/
				if (i == 0 && j == 0) {
					t.SHOD[0][i][j] = 1.0;
				}
				else if (i == nint - 1 && j == NQUAD) {
					t.SHOD[0][i][j] = 1.0;
				}
				else if (gq.S[j] == gi.S[i]) {
					t.SHOD[0][i][j] = 1.0;
				}
				else {
					t.SHOD[0][i][j] = -((1.0 - pow(gq.S[j], 2))*hq.LN[1][j]) / (n*(n + 1)*hi.LN[0][i] * (gq.S[j] - gi.S[i]));
				}
			}
		}

		//USE LEGENDRE POLYNOMIAL CARDINAL FUNCTION TO EVALUATE LOCAL SHAPE FUNCTION DERIVATIVE AT EACH NODE
		//reference: Klenow thesis page 35
		for (l = 0; l < nint; l++) { //interpolation point
			for (k = 0; k < NQUAD + 1; k++) { //quadrature point
				if (l == 0 && k == 0) {
					t.SHOD[1][l][k] = -double((n + 1)*n) / 4.0;
				}
				else if (l == nint - 1 && k == NQUAD) {
					t.SHOD[1][l][k] = double((n + 1)*n) / 4.0;
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
	for (i = 0; i < nint; i++) {   //Xi 
		for (j = 0; j < nint; j++) {  //Eta
			for (k = 0; k < nint; k++) { //Zeta
				//i j k are interpolation points and l m o are quadrature points 
				for (l = 0; l < NQUAD + 1; l++) { //the internal nodes on Xi 
					for (m = 0; m < NQUAD + 1; m++) { //the internal nodes on Eta 
						for (o = 0; o < NQUAD + 1; o++) { //the internal nodes on Zeta 
							t.SHL[3][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + o] = t.SHOD[0][i][l] * t.SHOD[0][j][m] * t.SHOD[0][k][o];
							//3D Nth ORDER SHAPE FUNCTION
							t.SHL[0][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + o] = t.SHOD[1][i][l] * t.SHOD[0][j][m] * t.SHOD[0][k][o];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Xi 
							t.SHL[1][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + o] = t.SHOD[0][i][l] * t.SHOD[1][j][m] * t.SHOD[0][k][o];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Eta 
							t.SHL[2][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + o] = t.SHOD[0][i][l] * t.SHOD[0][j][m] * t.SHOD[1][k][o];
							//3D Nth ORDER SHAPE FUNCTION DERIVATIVE W/R TO Zeta
						}
					}
				}
			}
		}
	}

	std::cout << " " << std::endl; 
	//check point: if the shape function value at the corresponding point is 0 
	//std::cout << "shape function value on the corresponding point is: " << t.SHL[3][LNA[1][1][1] - 1][1][1][1] << std::endl;

	return t;
}
