//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

//ELE_GEN discretizes the 2D domain with NEL elements

using namespace std;

//-------------------OUTPUT VARIABLES-----------------------//
struct TD_ELE_GENstruct TD_ELE_GEN(int NEL, int XNODE, int NNODE, int XNEL, int YNEL, double XL, double YL) {
	TD_ELE_GENstruct t;
	//LOBATTOstruct b;
	TD_LOCAL_NODEstruct c;
	double MX;  //INTERPOLATION NODE COUNTER
	double DX;  //INTERPOLATION NODE COUNTER
	double MY;  //INTERPOLATION NODE COUNTER
	double DY;  //INTERPOLATION NODE COUNTER
	double XHE; double YHE;
	int ICX, ICY;
	int iix, iiy;
	int i, j;
	double *XE; double*YE;
	double **XIN; double **YIN;
	double* Z;
	XE = new double[XNEL + 1]; YE = new double[YNEL + 1];
	XIN = new double*[XNEL];
	for (i = 0; i < XNEL; i++) {
		XIN[i] = new double[NC + 1];
	}
	YIN = new double*[YNEL];
	for (i = 0; i < YNEL; i++) {
		YIN[i] = new double[NC + 1];
	}
	for (i = 0; i < XNEL; i++) {
		for (j = 0; j < NC + 1; j++) {
			XIN[i][j] = 0.0;
		}
	}
	for (i = 0; i < YNEL; i++) {
		for (j = 0; j < NC + 1; j++) {
			YIN[i][j] = 0.0;
		}
	}

	//-------------INITIALIZE POINTER VARIABLE---------------//
	//t.XCOORD = new double[NNODE];
	//t.YCOORD = new double[NNODE];
	t.GCOORD = new double*[NNODE];
	for (i = 0; i < NNODE; i++) {
		t.GCOORD[i] = new double[2];
	}

	/*
	for (i = 0; i < NNODE; i++) {
		t.XCOORD[i] = 0.0;
		t.YCOORD[i] = 0.0;
	}
	*/
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j <2; j++) {
			t.GCOORD[i][j] = 0.0;
		}
	}

	t.IEN = new int*[(NC + 1)*(NC + 1)];
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		t.IEN[i] = new int[NEL];
	}
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		for (j = 0; j < NEL; j++) {
			t.IEN[i][j] = 0;
		}
	}
	//--------------BEGIN ELE_GEN COMP-----------------------//
	//COMPUTE ELEMENT LENGTH
	XHE = XL / XNEL;
	YHE = YL / YNEL;

	//COMPUTE X-COORDINATES OF ELEMENT END NODES
	for (i = 0; i < XNEL + 1; i++) {
		XE[i] = (i)*XHE; //start from the left end node to right end node
	}

	//COMPUTE Y-COORDINATES OF ELEMENT END NODES
	for (i = 0; i < YNEL + 1; i++) {
		YE[i] = (i)*YHE;
	}

	Z = new double[(NC + 1) - 2];
	for (i = 0; i < (NC + 1) - 2; i++) {
		Z[i] = 0.0;
	}

	for (i = 0; i < (NC + 1) - 2; i++) {
		Z[i] = -1.0 + (i + 1)*(2.0 / ((NC + 1) - 1));
	}


	//COMPUTE ELEMENT INTERPOLATION NODE LOCATIONS (X-DIRECTION)
	for (i = 0; i < XNEL; i++) {
		MX = 0.5*(XE[i + 1] + XE[i]);
		DX = 0.5*(XE[i + 1] - XE[i]);

		XIN[i][0] = XE[i];

		if (NC>1) {

			//b=LOBATTO();

			for (j = 0; j < NC - 1; j++) {
				XIN[i][j + 1] = MX + Z[j] * DX;
			}
		}
		XIN[i][NC] = XE[i + 1];
	}

	//COMPUTE ELEMENT INTERPOLATION NODE LOCATIONS (Y-DIRECTION)
	for (i = 0; i <YNEL; i++) {
		MY = 0.5*(YE[i + 1] + YE[i]);
		DY = 0.5*(YE[i + 1] - YE[i]);

		YIN[i][0] = YE[i];

		if (NC>1) {
			for (j = 0; j < NC - 1; j++) {
				YIN[i][j + 1] = MY + Z[j] * DY;
			}
		}
		YIN[i][NC] = YE[i + 1];
	}

	//COMPUTE ELEMENT NODES ARRAY AND GLOBAL NODE LOCATIONS
	for (i = 0; i < YNEL; i++) {
		for (j = 0; j < XNEL; j++) {
			for (ICY = 0; ICY < NC + 1; ICY++) {
				for (ICX = 0; ICX < NC + 1; ICX++) {
					t.GCOORD[ICX + 1 + ICY*XNODE + NC*j + (NC*XNODE)*i - 1][0] = XIN[j][ICX];
					t.GCOORD[ICX + 1 + ICY*XNODE + NC*j + (NC*XNODE)*i - 1][1] = -YIN[i][ICY];
				}
			}
		}
	}

	/*
	for (i = 0; i < NNODE; i++) {
		t.XCOORD[i] = t.GCOORD[i][0];
		t.YCOORD[i] = t.GCOORD[i][1];
	}
	*/

	//COMPUTE LOCAL NODES ARRAY
	c = TD_LOCAL_NODE(NC);
	//LNA is defined from 0

	//COMPUTE ELEMENT NODES ARRAY AND GLOBAL NODE LOCATIONS
	for (iiy = 0; iiy < YNEL; iiy++) {
		for (iix = 0; iix < XNEL; iix++) {
			for (i = 0; i < NC + 1; i++) {
				for (j = 0; j < NC + 1; j++) {
					t.IEN[c.LNA[j][i] - 1][iiy*XNEL + (iix + 1) - 1] = 1 + XNODE*(NC + 1 - (i + 1)) + j + (NC*iix + NC*XNODE*iiy);
				}
			}
		}
	}
	//----------------------END ELE_GEN COMP------------------------------//

	for (i = 0; i < NC + 1; i++) {
		delete[] c.LNA[i];
	}
	delete[] c.LNA;

	delete[] XE;
	delete[] YE;
	for (i = 0; i < XNEL;i++) {
		delete[] XIN[i];
	}
	delete[] XIN;
	for (i = 0; i < YNEL;i++) {
		delete[] YIN[i];
	}
	delete[] YIN;
	delete[] Z;

	return t;
}