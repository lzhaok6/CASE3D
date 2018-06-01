//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>


//NRB determines the NRB local node numbering and the associated NRB arrays
void NRB(int NNODE, double **GCOORD, double* W, int*** LNA, int**IEN, int NEL, double***SHL, double***SHOD) {
	int i, j, k, l, m, e;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber];

	LOBATTOstruct bq;
	bq = LOBATTO(Nq); //one unit higher than the interpolation order
	GLLQUADstruct gq;
	gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
	LOCAL_SHAPEstruct ls;
	ls = LOCAL_SHAPE(LNA, N, Nq); //one order higher than the interpolation order

	//Construct the ADMASTER for the integration (insert one additional GLL quadrature point)
	//Can we use just one surface for NRBC? I think so  
	//We need to combine the physical group corresponding to NRBC first
	//Do we need to individual points on the NRBC? 
	//ADMASTER is for each element. 
	//We need NRBA to store all the nodes on NRB
	//We also need the total number of node in each NRB surface
	//We need to lump ADMASTER to get ADMASTERG 
	nr[0].ADMASTER = new double**[nr[0].NEL_nrb];
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		nr[0].ADMASTER[i] = new double*[NINT];
		for (j = 0; j < NINT; j++) {
			nr[0].ADMASTER[i][j] = new double[NINT];
		}
	}
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				nr[0].ADMASTER[i][j][k] = 0.0;
			}
		}
	}

	nr[0].ADMASTERG = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		nr[0].ADMASTERG[i] = 0.0;
	}

	double DUNC;
	//Check out Evernote "Some thoughts on boundary force integration" 
	for (e = 0; e < nr[0].NEL_nrb; e++) {
		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (i = 0; i < NqINT; i++) {
					for (j = 0; j < NqINT; j++) {
						for (k = 0; k < NqINT; k++) {
							DUNC += gq.W[i] * gq.W[j] * gq.W[k] * ls.SHL[3][nr[0].DP[m] - 1][i*NqINT*NqINT + j * NqINT + k] * ls.SHL[3][nr[0].DP[l] - 1][i*NqINT*NqINT + j * NqINT + k] * nr[0].JACOB[e][i*NqINT*NqINT + j*NqINT + k];
						}
					}
				}
				nr[0].ADMASTER[e][m][l] += DUNC;
			}
		}
	}

	for (e = 0; e < nr[0].NEL_nrb; e++) {
		for (j = 0; j < NINT*NINT; j++) {
			for (k = 0; k < NINT*NINT; k++) {
				nr[0].ADMASTERG[nr[0].IEN_gb[nr[0].DP[k] - 1][e] - 1] += nr[0].ADMASTER[e][j][k];
			}
		}
	}
	return;
}