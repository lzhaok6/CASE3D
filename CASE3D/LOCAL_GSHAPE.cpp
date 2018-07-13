//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

//LOCAL_GSHAPE determines the local shape function and its derivatives of the geometry nodes 
//at the GLL quadrature points 

struct LOCAL_GSHAPEstruct LOCAL_GSHAPE(double* S, int*** LNA, int NINT) {
	//double t.MCOORD[8][3]; //MATRIX OF LOCAL ELEMENT COORDINATES
	int i, j, k, l, m;
	LOCAL_GSHAPEstruct t;
	t.GSHL = new double****[4];
	for (i = 0; i < 4; i++) {
		t.GSHL[i] = new double***[8];
		for (j = 0; j < 8; j++) { //NINT*NINT*NINT
			t.GSHL[i][j] = new double**[NINT];
			for (k = 0; k < NINT; k++) {
				t.GSHL[i][j][k] = new double*[NINT];
				for (l = 0; l < NINT; l++) {
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

	/*
	t.GSHL_2D = new double***[3];
	for (i = 0; i < 3; i++) {
		t.GSHL_2D[i] = new double**[4]; //4 points
		for (j = 0; j < 4; j++) {
			t.GSHL_2D[i][j] = new double*[NINT];
			for (k = 0; k < NINT; k++) {
				t.GSHL_2D[i][j][k] = new double[NINT];
			}
		}
	}
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 4; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT; l++) {
					t.GSHL_2D[i][j][k][l] = 0.0;
				}
			}
		}
	}
	*/

	//----------BEGIN COMPUT----------//
	//DEFINE LOCAL ELEMENT COORDINATES
	//LOCAL COORDINATE [-1,1]

	//X DIRECTION
	t.MCOORD[0][0] = -1.0; //good
	t.MCOORD[1][0] = 1.0; //good
	t.MCOORD[2][0] = 1.0; //good
	t.MCOORD[3][0] = -1.0; //good 
	t.MCOORD[4][0] = -1.0; //good
	t.MCOORD[5][0] = 1.0; //good 
	t.MCOORD[6][0] = 1.0; //good 
	t.MCOORD[7][0] = -1.0; //good
	//Y DIRECTION
	t.MCOORD[0][1] = -1.0;
	t.MCOORD[1][1] = -1.0;
	t.MCOORD[2][1] = 1.0;
	t.MCOORD[3][1] = 1.0;
	t.MCOORD[4][1] = -1.0;
	t.MCOORD[5][1] = -1.0;
	t.MCOORD[6][1] = 1.0;
	t.MCOORD[7][1] = 1.0;
	//Z DIRECTION
	t.MCOORD[0][2] = -1.0;
	t.MCOORD[1][2] = -1.0;
	t.MCOORD[2][2] = -1.0;
	t.MCOORD[3][2] = -1.0;
	t.MCOORD[4][2] = 1.0;
	t.MCOORD[5][2] = 1.0;
	t.MCOORD[6][2] = 1.0;
	t.MCOORD[7][2] = 1.0;

	/*
	//X DIRECTION
	t.MCOORD[LNA[0][0][0]-1][0] = -1.0; //good
	t.MCOORD[LNA[N][0][0]-1][0] = 1.0; //good
	t.MCOORD[LNA[N][N][0]-1][0] = 1.0; //good
	t.MCOORD[LNA[0][N][0]-1][0] = -1.0; //good
	t.MCOORD[LNA[0][0][N]-1][0] = -1.0; //good
	t.MCOORD[LNA[N][0][N]-1][0] = 1.0; //good
	t.MCOORD[LNA[N][N][N]-1][0] = 1.0; //good
	t.MCOORD[LNA[0][N][N]-1][0] = -1.0; //good
	//Y DIRECTION
	t.MCOORD[LNA[0][0][0] - 1][1] = -1.0;
	t.MCOORD[LNA[N][0][0] - 1][1] = -1.0;
	t.MCOORD[LNA[N][N][0] - 1][1] = 1.0;
	t.MCOORD[LNA[0][N][0] - 1][1] = 1.0;
	t.MCOORD[LNA[0][0][N] - 1][1] = -1.0;
	t.MCOORD[LNA[N][0][N] - 1][1] = -1.0;
	t.MCOORD[LNA[N][N][N] - 1][1] = 1.0;
	t.MCOORD[LNA[0][N][N] - 1][1] = 1.0;
	//Z DIRECTION
	t.MCOORD[LNA[0][0][0] - 1][2] = -1.0;
	t.MCOORD[LNA[N][0][0] - 1][2] = -1.0;
	t.MCOORD[LNA[N][N][0] - 1][2] = -1.0;
	t.MCOORD[LNA[0][N][0] - 1][2] = -1.0;
	t.MCOORD[LNA[0][0][N] - 1][2] = 1.0;
	t.MCOORD[LNA[N][0][N] - 1][2] = 1.0;
	t.MCOORD[LNA[N][N][N] - 1][2] = 1.0;
	t.MCOORD[LNA[0][N][N] - 1][2] = 1.0;
	*/
	//EVALUATE LOCAL SHAPE FUNCTION AND LOCAL SHAPE FUNCTION DERIVATIVES AT QUAD POINTS 
	for (i = 0; i < 8; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
		for (j = 0; j < NINT; j++) { // j k l are integration point where discrete value of shape function is derived
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT; l++) {
					t.GSHL[3][i][j][k][l] = (1.0 / 8.0)*(1 + t.MCOORD[i][0] * S[j])*(1 + t.MCOORD[i][1] * S[k])*(1 + t.MCOORD[i][2] * S[l]);
					//not derivative
					t.GSHL[0][i][j][k][l] = (1.0 / 8.0)*t.MCOORD[i][0] * (1 + t.MCOORD[i][1] * S[k])*(1 + t.MCOORD[i][2] * S[l]);
					//xi direction derivative
					t.GSHL[1][i][j][k][l] = (1.0 / 8.0)*t.MCOORD[i][1] * (1 + t.MCOORD[i][0] * S[j])*(1 + t.MCOORD[i][2] * S[l]);
					//eta direction derivative
					t.GSHL[2][i][j][k][l] = (1.0 / 8.0)*t.MCOORD[i][2] * (1 + t.MCOORD[i][0] * S[j])*(1 + t.MCOORD[i][1] * S[k]);
					//zeta direction derivative
				}
			}
		}
	}
	//check point: the shape function is zero at the corresponding point:
	//std::cout <<"Global shape function value: "<<t.GSHL[3][0][0][0][0] << std::endl; 
	//t.GSHL[3][0][0][0][0]=1
	/*
	//We can move this part to the FSILINK code
	for (i = 0; i < 4; i++) { //ref: Pozrikidis IFSM P666 //i is the point where shape functions were derived (linear shape function for geometry discretization)
		for (j = 0; j < NINT; j++) { // j k l are integration point where discrete value of shape function is derived
			for (k = 0; k < NINT; k++) {
				t.GSHL_2D[2][i][j][k] = (1.0 / 4.0)*(1 + t.MCOORD[][0] * S[j])*(1 + t.MCOORD[][1] * S[k]);
				//not derivative
				t.GSHL_2D[0][i][j][k] = (1.0 / 4.0)*t.MCOORD[][0] * (1 + t.MCOORD[][1] * S[k]);
				//X direction derivative
				t.GSHL_2D[1][i][j][k] = (1.0 / 4.0)*t.MCOORD[][1] * (1 + t.MCOORD[][0] * S[j]);
				//Y direction derivative
			}
		}
	}
	*/
	std::cout << " " << std::endl;
	return t;
}
