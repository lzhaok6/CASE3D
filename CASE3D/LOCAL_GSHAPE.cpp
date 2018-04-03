//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

//LOCAL_GSHAPE determines the local shape function and its derivatives of the geometry nodes 
//at the GLL quadrature points 

struct LOCAL_GSHAPEstruct LOCAL_GSHAPE(double* S, int*** LNA) {
	double MCOORD[8][3]; //MATRIX OF LOCAL ELEMENT COORDINATES
	int i, j, k, l, m;
	LOCAL_GSHAPEstruct t;
	t.GSHL = new double****[4];
	for (i = 0; i < 4; i++) {
		t.GSHL[i] = new double***[8];
		for (j = 0; j < 8; j++) { //NINT*NINT*NINT
			t.GSHL[i][j] = new double**[NINT];
			for (k = 0; k < NINT; k++) {
				t.GSHL[i][j][k] = new double*[NINT];
				for (l = 0; l < NINT;l++) {
					t.GSHL[i][j][k][l] = new double[NINT];
				}
			}
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < 8; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT; l++) {
					for (m = 0; m < NINT; m++) {
						t.GSHL[i][j][k][l][m] = 0.0;
					}
				}
			}
		}
	}

	//----------BEGIN COMPUT----------//
	//DEFINE LOCAL ELEMENT COORDINATES
	//LOCAL COORDINATE [-1,1]
	
	//X DIRECTION
	MCOORD[0][0] = -1.0; //good
	MCOORD[1][0] = 1.0; //good
	MCOORD[2][0] = 1.0; //good
	MCOORD[3][0] = -1.0; //good 
	MCOORD[4][0] = -1.0; //good
	MCOORD[5][0] = 1.0; //good 
	MCOORD[6][0] = 1.0; //good 
	MCOORD[7][0] = -1.0; //good
	//Y DIRECTION
	MCOORD[0][1] = -1.0;
	MCOORD[1][1] = -1.0;
	MCOORD[2][1] = 1.0;
	MCOORD[3][1] = 1.0;
	MCOORD[4][1] = -1.0;
	MCOORD[5][1] = -1.0;
	MCOORD[6][1] = 1.0;
	MCOORD[7][1] = 1.0;
	//Z DIRECTION
	MCOORD[0][2] = -1.0;
	MCOORD[1][2] = -1.0;
	MCOORD[2][2] = -1.0;
	MCOORD[3][2] = -1.0;
	MCOORD[4][2] = 1.0;
	MCOORD[5][2] = 1.0;
	MCOORD[6][2] = 1.0;
	MCOORD[7][2] = 1.0;
	
	/*
	//X DIRECTION
	MCOORD[LNA[0][0][0]-1][0] = -1.0; //good
	MCOORD[LNA[N][0][0]-1][0] = 1.0; //good
	MCOORD[LNA[N][N][0]-1][0] = 1.0; //good
	MCOORD[LNA[0][N][0]-1][0] = -1.0; //good 
	MCOORD[LNA[0][0][N]-1][0] = -1.0; //good
	MCOORD[LNA[N][0][N]-1][0] = 1.0; //good 
	MCOORD[LNA[N][N][N]-1][0] = 1.0; //good 
	MCOORD[LNA[0][N][N]-1][0] = -1.0; //good
	//Y DIRECTION
	MCOORD[LNA[0][0][0] - 1][1] = -1.0;
	MCOORD[LNA[N][0][0] - 1][1] = -1.0;
	MCOORD[LNA[N][N][0] - 1][1] = 1.0;
	MCOORD[LNA[0][N][0] - 1][1] = 1.0;
	MCOORD[LNA[0][0][N] - 1][1] = -1.0;
	MCOORD[LNA[N][0][N] - 1][1] = -1.0;
	MCOORD[LNA[N][N][N] - 1][1] = 1.0;
	MCOORD[LNA[0][N][N] - 1][1] = 1.0;
	//Z DIRECTION
	MCOORD[LNA[0][0][0] - 1][2] = -1.0;
	MCOORD[LNA[N][0][0] - 1][2] = -1.0;
	MCOORD[LNA[N][N][0] - 1][2] = -1.0;
	MCOORD[LNA[0][N][0] - 1][2] = -1.0;
	MCOORD[LNA[0][0][N] - 1][2] = 1.0;
	MCOORD[LNA[N][0][N] - 1][2] = 1.0;
	MCOORD[LNA[N][N][N] - 1][2] = 1.0;
	MCOORD[LNA[0][N][N] - 1][2] = 1.0;
	*/
	//EVALUATE LOCAL SHAPE FUNCTION AND LOCAL SHAPE FUNCTION DERIVATIVES AT QUAD POINTS 
	for (i = 0; i < 8; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
		for (j = 0; j < NINT; j++) { // j k l are integration point where discrete value of shape function is derived
			for (k = 0; k < NINT; k++) {	
				for (l = 0; l < NINT;l++) {
					t.GSHL[3][i][j][k][l] = (1.0 / 8.0)*(1 + MCOORD[i][0] * S[j])*(1 + MCOORD[i][1] * S[k])*(1 + MCOORD[i][2] * S[l]);
					//not derivative
					t.GSHL[0][i][j][k][l] = (1.0 / 8.0)*MCOORD[i][0] * (1 + MCOORD[i][1] * S[k])*(1 + MCOORD[i][2] * S[l]);
					//X direction derivative
					t.GSHL[1][i][j][k][l] = (1.0 / 8.0)*MCOORD[i][1] * (1 + MCOORD[i][0] * S[j])*(1 + MCOORD[i][2] * S[l]);
					//Y direction derivative
					t.GSHL[2][i][j][k][l] = (1.0 / 8.0)*MCOORD[i][2] * (1 + MCOORD[i][0] * S[j])*(1 + MCOORD[i][1] * S[k]);
					//Z direction derivative
				}
			}
		}
	}

	//check point: the shape function is zero at the corresponding point:
	//std::cout <<"Global shape function value: "<<t.GSHL[3][0][0][0][0] << std::endl; 
	//t.GSHL[3][0][0][0][0]=1
	std::cout << " " << std::endl;
	return t;
}
