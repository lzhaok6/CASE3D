//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//EIGENMAX estimates the maximum mesh Eigen value using Greschgorin's Theorem
double EIGENMAX(double** QMASTER, double*** HMASTER, int NEL) {
	double SUMH; //SUM VARIABLE
	double HFUNC; //REACTANCE MATRIX LOOP VARIABLE
	double FUNC[NINT*NINT*NINT]; //EIGENVALUE STORAGE ARRAY
	double LMAX; 
	double* LMAXi;
	int e, i, j;
	LMAXi = new double[NEL];
	//-------------BEGIN EIGENMAX COMP.----------------//
	for (e = 0; e < NEL;e++) {
		for (i = 0; i < NINT*NINT*NINT; i++) { //loop through each line
			SUMH = 0.0;
			for (j = 0; j < NINT*NINT*NINT; j++) { //column sum
				if (j != i) {
					HFUNC = abs(HMASTER[e][i][j]);
					SUMH = SUMH + HFUNC;
				}
			}
			//add all element to diagonal element line by line
			FUNC[i] = (HMASTER[e][i][i] + SUMH) / QMASTER[i][i]; //the HMASTER here doesn't have the absolute sign. 
		}

		//DERIVE THE MAX EIGEN VALUE
		LMAXi[e] = FUNC[0];
		for (j = 1; j < NINT*NINT*NINT; j++) {
			if (FUNC[j] > LMAXi[e]) {
				LMAXi[e] = FUNC[j];
			}
		}
	}

	LMAX = LMAXi[0];
	for (e = 1; e < NEL; e++) {
		if (LMAXi[e] > LMAX) {
			LMAX = LMAXi[e];
		}
	}

	std::cout << " " << std::endl;
	//check point: where the LMAX is a relatively valid value (compare to 2D problem)
	return(LMAX);
}