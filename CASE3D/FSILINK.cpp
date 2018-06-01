//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>

//construct the link between structure and fluid module A C D separately
void FSILINK(double* W, int*** LNA, int**IEN, double***SHL, double**GCOORD, int NNODE, double***SHOD) {
	int i, j, k, l, m, n, o, e;
	int IC; //counter 
	//extern int owsfnumber;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[owsfnumber];

	//MpCCI model file is written in meshgeneration code. 
	//Need the node on wetted surface GIDN
	LOBATTOstruct bq;
	bq = LOBATTO(Nq); //one unit higher than the interpolation order
	GLLQUADstruct gq;
	gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
	LOCAL_SHAPEstruct ls;
	ls = LOCAL_SHAPE(LNA, N, Nq); //one order higher than the interpolation order
	
	double FUNC;
	//Check out Evernote "Some thoughts on boundary force integration" 
	for (e = 0; e < ol[0].FSNEL; e++) {
		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				FUNC = 0.0; //accumulator 
				for (i = 0; i < NqINT; i++) {
					for (j = 0; j < NqINT; j++) {
						for (k = 0; k < NqINT; k++) {
							FUNC += gq.W[i] * gq.W[j] * gq.W[k] * ls.SHL[3][ol[0].FP[m] - 1][i*NqINT*NqINT + j * NqINT + k] * ls.SHL[3][ol[0].FP[l] - 1][i*NqINT*NqINT + j * NqINT + k] * ol[0].JACOB[e][i*NqINT*NqINT + j*NqINT + k];
						}
					}
				}
				ol[0].FPMASTER[e][m][l] += FUNC;
			}
		}
	}

	//Derive the integration weight of 2D linear element (4 nodes) FPMASTER_2D. The size should be [4][NINT^2]
	//Initialize FPMASTER_2D first
	LOCAL_SHAPEstruct ls_ln; //ln means linear
	ls_ln = LOCAL_SHAPE(LNA, 1, Nq); //one order higher than the interpolation order
	//LNA is useless here. It only allows us to run the LOCAL_SHAPE function. We are just gonna use the SHOD from that function. 2
	for (e = 0; e < ol[0].FSNEL; e++) {
		for (m = 0; m < 2; m++) {
			for (n = 0; n < 2; n++) {
				for (l = 0; l < NINT; l++) {
					for (o = 0; o < NINT; o++) {
						FUNC = 0.0; //accumulator 
						for (i = 0; i < NqINT; i++) {
							for (j = 0; j < NqINT; j++) {
								FUNC += gq.W[i] * gq.W[j] * (ls_ln.SHOD[0][m][i] * ls_ln.SHOD[0][n][j]) * (ls.SHOD[0][l][i] * ls.SHOD[0][o][j]) * ol[0].Jacob_2D[e][i*NqINT*NqINT + j*NqINT];
								//FUNC += gq.W[i] * gq.W[j] * (linear 2D shape function defined on Nq order GLL quadrature points) * (Nth order 2D shape function defined on Nq order GLL quadrature points) * ol[0].Jacob_2D[e][i*NqINT*NqINT + j*NqINT];
							}
						}
						ol[0].FPMASTER_2D[e][ol[0].LNA_algo2[m][n] - 1][ol[0].LNA_2D[l][o] - 1] += FUNC;
					}
				}
			}
		}
	}

	std::cout << " " << std::endl;
	return;
}
