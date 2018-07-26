//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//EIGENMAX estimates the maximum mesh Eigen value using Greschgorin's Theorem
double EIGENMAX(double*** QMASTER, double*** HMASTER, int NEL) {
	double SUMH; //SUM VARIABLE
	double HFUNC; //REACTANCE MATRIX LOOP VARIABLE
	//double FUNC[NINT*NINT*NINT]; //EIGENVALUE STORAGE ARRAY
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
	double*FUNC; 
	FUNC = new double[elenode3D];
	double LMAX; 
	double* LMAXi;
	int e, i, j;
	LMAXi = new double[NEL];

	//-------------BEGIN EIGENMAX COMP.----------------//
	for (e = 0; e < NEL;e++) {
		for (i = 0; i < elenode3D; i++) { //loop through each line
			SUMH = 0.0;
			for (j = 0; j < elenode3D; j++) { //Column sum (Since HMASTER is symmetric, row sum and column sum would actually yield the same result)
				if (j != i) {
					HFUNC = abs(HMASTER[e][i][j]);
					SUMH = SUMH + HFUNC;
				}
			}
			//add all element to diagonal element line by line
			FUNC[i] = (HMASTER[e][i][i] + SUMH) / QMASTER[e][i][i]; //the HMASTER here doesn't have the absolute sign. 
		}

		//DERIVE THE MAX EIGEN VALUE IN THE ELEMENT e
		LMAXi[e] = FUNC[0];
		for (j = 1; j < elenode3D; j++) {
			if (FUNC[j] > LMAXi[e]) {
				LMAXi[e] = FUNC[j];
			}
		}
	}

	//Obtain the largest eigenvalue of the elements
	LMAX = LMAXi[0];
	for (e = 1; e < NEL; e++) {
		if (LMAXi[e] > LMAX) {
			LMAX = LMAXi[e];
		}
	}

	//std::cout << " " << std::endl;
	//check point: where the LMAX is a relatively valid value (compare to 2D problem)
	return(LMAX);
}