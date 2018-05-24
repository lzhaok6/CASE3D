//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>

void ELE_GEN(int NEL, double **GCOORD, int**IEN, int ***LNA ,double*Z) {
	int i, j, k, l, m;

    //Get information of wet surfaces
	//extern int owsfnumber;
	extern OWETSURF ol[5]; //defined in FSILINK 
	
	for (i = 0; i < owsfnumber; i++) {
		ol[i].FSNEL = 0;
	}

	//get FSNEL information 
	int count1 = 0;	int count4 = 0; int count2 = 0; int count5 = 0; int count3 = 0;
	int *GIDF1; int *GIDF2; int *GIDF3; int *GIDF4; 
	//get wet element No. on each wet surface
	GIDF1 = new int[SXNEL*SYNEL];
	GIDF2 = new int[SZNEL*SYNEL];
	GIDF3 = new int[SXNEL*SZNEL];
	GIDF4 = new int[SXNEL*SYNEL];
	for (i = 0; i < SXNEL*SYNEL;i++) {
		GIDF1[i] = 0;
	}
	for (i = 0; i < SZNEL*SYNEL;i++) {
		GIDF2[i] = 0;
	}
	for (i = 0; i < SXNEL*SZNEL; i++) {
		GIDF3[i] = 0;
	}
	for (i = 0; i < SXNEL*SYNEL; i++) {
		GIDF4[i] = 0;
	}

	//fill the GIDF array by node coordinate searching
	for (i = 0; i < NEL; i++) {
		//surface 1 
		if (abs(GCOORD[IEN[LNA[N][N][N] - 1][i] - 1][2] - (-SZ / 2)) < 1e-6 &&
			0.5*(GCOORD[IEN[LNA[0][0][N] - 1][i] - 1][1] + GCOORD[IEN[LNA[0][N][N] - 1][i] - 1][1]) > -SY  &&
			0.5*(GCOORD[IEN[LNA[0][0][N] - 1][i] - 1][0] + GCOORD[IEN[LNA[N][0][N] - 1][i] - 1][0]) > -SX
			&& 0.5*(GCOORD[IEN[LNA[0][0][N] - 1][i] - 1][0] + GCOORD[IEN[LNA[N][0][N] - 1][i] - 1][0]) < 0) {
			ol[0].FSNEL += 1;
			GIDF1[count1] = i + 1;
			count1 += 1;
		}
		//surface 4
		if (abs(GCOORD[IEN[LNA[N][N][0] - 1][i] - 1][2] - (SZ / 2)) < 1e-6 &&
			0.5*(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][1] + GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][1]) > -SY  &&
			0.5*(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][0] + GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][0]) > -SX
			&& 0.5*(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][0] + GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][0]) < 0) {
			ol[3].FSNEL += 1;
			GIDF4[count4] = i + 1;
			count4 += 1;
		}
		//surface 2
		if (abs(GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][0] - (-SX)) < 1e-6 &&
			0.5*(GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][1] + GCOORD[IEN[LNA[N][N][0] - 1][i] - 1][1]) > -SY  &&
			0.5*(GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][2] + GCOORD[IEN[LNA[N][0][N] - 1][i] - 1][2]) > -SZ / 2
			&& 0.5*(GCOORD[IEN[LNA[N][0][0] - 1][i] - 1][2] + GCOORD[IEN[LNA[N][0][N] - 1][i] - 1][2]) < SZ / 2) {
			ol[1].FSNEL += 1;
			GIDF2[count2] = i + 1;
			count2 += 1;
		}
		//surface 3
		if (abs(GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][1] - (-SY )) < 1e-6 &&
			0.5*(GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][0] + GCOORD[IEN[LNA[N][N][0] - 1][i] - 1][0]) > -SX &&
			0.5*(GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][0]+ GCOORD[IEN[LNA[N][N][0] - 1][i] - 1][0]) < 0 &&
			0.5*(GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][2]+ GCOORD[IEN[LNA[0][N][N] - 1][i] - 1][2]) > -SZ / 2
			&& 0.5*(GCOORD[IEN[LNA[0][N][0] - 1][i] - 1][2]+ GCOORD[IEN[LNA[0][N][N] - 1][i] - 1][2]) < SZ / 2) {
			ol[2].FSNEL += 1;
			GIDF3[count3] = i + 1;
			count3 += 1;
		}
	}
	if (AX > 0 && BZ > 0) {
		if (count1 != SXNEL*SYNEL || count2 != SZNEL*SYNEL || count3 != SXNEL*SZNEL || count4 != SXNEL*SYNEL) {
			std::cout << "the point recognition is inconsistent in ELE_GEN" << std::endl;
			system("PAUSE ");
		}
	}
	else if (AX == 0 && BZ == 0) {
		if (count3 != SXNEL*SZNEL) {
			std::cout << "the point recognition is inconsistent in ELE_GEN" << std::endl;
			system("PAUSE ");
		}
	}

	//declare and initialize GIDF
	for (i = 0; i < owsfnumber; i++) {
		ol[i].GIDF = new int[ol[i].FSNEL];
	}
	
	for (i = 0; i < owsfnumber; i++) {
		for (j = 0; j < ol[i].FSNEL; j++) {
			ol[i].GIDF[j] = 0;  
		}
	}

	//define GIDF for all wet surfaces 
	//surface 1 
	int count;
	count = 0;
	for (i = 0; i < SYNEL; i++) {
		for (j = 0; j < SXNEL; j++) {
			for (k = 0; k < ol[0].FSNEL;k++) {
				if (abs(GCOORD[IEN[LNA[0][N][N] - 1][GIDF1[k] - 1] - 1][0] - (-SX + j*XHE)) < 1e-5 &&
					abs(GCOORD[IEN[LNA[0][N][N] - 1][GIDF1[k] - 1] - 1][1] - ( - i*YHE)) < 1e-5) {
					ol[0].GIDF[count] = GIDF1[k];
					count += 1;
				}
			}
		}
	}
	if (count != ol[0].FSNEL) {
		std::cout << "ol[0].GIDF is wrong" << std::endl;
		system("PAUSE ");
	}
	//surface 4
	count = 0;
	for (i = 0; i < SYNEL; i++) {
		for (j = 0; j < SXNEL; j++) {
			for (k = 0; k < ol[3].FSNEL; k++) {
				if (abs(GCOORD[IEN[LNA[0][N][0] - 1][GIDF4[k] - 1] - 1][0] - (-SX + j*XHE)) < 1e-5 &&
					abs(GCOORD[IEN[LNA[0][N][0] - 1][GIDF4[k] - 1] - 1][1] - ( - i*YHE)) < 1e-5) {
					ol[3].GIDF[count] = GIDF4[k];
					count += 1;
				}
			}
		}
	}
	if (count != ol[3].FSNEL) {
		std::cout << "ol[3].GIDF is wrong" << std::endl;
		system("PAUSE ");
	}
	//surface 2
	count = 0;
	for (i = 0; i < SYNEL; i++) {
		for (j = 0; j < SZNEL; j++) {
			for (k = 0; k < ol[1].FSNEL; k++) {
				if (abs(GCOORD[IEN[LNA[N][N][0] - 1][GIDF2[k] - 1] - 1][2] - (-SZ / 2 + j*ZHE)) < 1e-5 &&
					abs(GCOORD[IEN[LNA[N][N][0] - 1][GIDF2[k] - 1] - 1][1] - ( - i*YHE)) < 1e-5) {
					ol[1].GIDF[count] = GIDF2[k];
					count += 1;
				}
			}
		}
	}
	if (count != ol[1].FSNEL) {
		std::cout << "ol[1].GIDF is wrong" << std::endl;
		system("PAUSE ");
	}

	count = 0;
	for (i = 0; i < SZNEL; i++) {
		for (j = 0; j < SXNEL; j++) {
			for (k = 0; k < ol[2].FSNEL; k++) {
				if (abs(GCOORD[IEN[LNA[0][N][0] - 1][GIDF3[k] - 1] - 1][0] - (-SX + j*XHE)) < 1e-5 &&
					abs(GCOORD[IEN[LNA[0][N][N] - 1][GIDF3[k] - 1] - 1][2] - (SZ / 2 - i*ZHE)) < 1e-5) {
					ol[2].GIDF[count] = GIDF3[k];
					count += 1;
				}
			}
		}
	}
	if (count != ol[2].FSNEL) {
		std::cout << "ol[2].GIDF is wrong" << std::endl;
		system("PAUSE ");
	}

	//get 2D element connectivity matrix IEN_2D using subroutine 2D_ELE_GEN.cpp 
	int XNEL1; int YNEL1; int NEL1; int XNEL2; int YNEL2; int NEL2; int XNEL3; int YNEL3; int NEL3; 
	if (mappingalgo == 1 || mappingalgo == 3 || mappingalgo == 4 || mappingalgo == 5) {
		//surface 1 
		XNEL1 = SXNEL / refine;
		YNEL1 = SYNEL / refine;
		//surface 2 
		XNEL2 = SZNEL / refine;
		YNEL2 = SYNEL / refine;
		//surface 3
		XNEL3 = SXNEL / refine;
		YNEL3 = SZNEL / refine;
	}
	else if (mappingalgo == 2) {
		//surface 1 
		XNEL1 = SXNEL*N;
		YNEL1 = SYNEL*N;
		//surface 2 
		XNEL2 = SZNEL*N;
		YNEL2 = SYNEL*N;
		//surface 3
		XNEL3 = SXNEL*N;
		YNEL3 = SZNEL*N;
	}

	NEL1 = XNEL1*YNEL1;
	NEL2 = XNEL2*YNEL2;
	NEL3 = XNEL3*YNEL3;
	if (AX > 0 && BZ > 0) {
		ol[0].FSNEL_fem = NEL1;
		ol[3].FSNEL_fem = NEL1;
		ol[1].FSNEL_fem = NEL2;
	}
	else {
		ol[0].FSNEL_fem = 0;
		ol[3].FSNEL_fem = 0;
		ol[1].FSNEL_fem = 0;
	}
	ol[2].FSNEL_fem = NEL3;
	int XNODE1 = XNEL1 + 1;
	int YNODE1 = YNEL1 + 1;
	int NODE1 = XNODE1*YNODE1;
	int XNODE2 = XNEL2 + 1;
	int YNODE2 = YNEL2 + 1;
	int NODE2 = XNODE2*YNODE2;
	int XNODE3 = XNEL3 + 1;
	int YNODE3 = YNEL3 + 1;
	int NODE3 = XNODE3*YNODE3;
	//IEN_2D is later used in interface_mapping.cpp
	//TD_ELE_GEN is also used in FSILINK subroutine used to generate the model file 
	for (i = 0; i < owsfnumber; i++) {
		ol[i].IEN_2D = new int*[(NC + 1)*(NC + 1)];
	}
	for (i = 0; i < (NC + 1)*(NC + 1); i++) {
		ol[0].IEN_2D[i] = new int[NEL1];
		ol[1].IEN_2D[i] = new int[NEL2];
		ol[2].IEN_2D[i] = new int[NEL3];
		ol[3].IEN_2D[i] = new int[NEL1];
	}

	for (j = 0; j < (NC + 1)*(NC + 1); j++) {
		for (k = 0; k < NEL1; k++) {
			ol[0].IEN_2D[j][k] = 0;
			ol[3].IEN_2D[j][k] = 0;
		}
		for (k = 0; k < NEL2; k++) {
			ol[1].IEN_2D[j][k] = 0;
		}
		for (k = 0; k < NEL3; k++) {
			ol[2].IEN_2D[j][k] = 0;
		}
	}

	int YNEL[5] = { YNEL1, YNEL2, YNEL3, YNEL1, YNEL2 };
	int XNEL[5] = { XNEL1, XNEL2, XNEL3, XNEL1, XNEL2 };
	int XNODE[5] = { XNODE1, XNODE2, XNODE3, XNODE1, XNODE2 };
	int nel[5] = { NEL1,NEL2,NEL3,NEL1,NEL2 };
	TD_LOCAL_NODEstruct c;
	c = TD_LOCAL_NODE(NC);
	for (m = 0; m < owsfnumber; m++) {
		for (k = 0; k < YNEL[m]; k++) {
			for (l = 0; l < XNEL[m]; l++) {
				for (i = 0; i < NC + 1; i++) {
					for (j = 0; j < NC + 1; j++) {
						ol[m].IEN_2D[c.LNA[j][i] - 1][k*XNEL[m] + (l + 1) - 1] = 1 + XNODE[m]*(NC + 1 - (i + 1)) + j + (NC*l + NC*XNODE[m]*k);
						if (c.LNA[j][i] > (NC + 1)*(NC + 1)) {
							std::cout << "LNA out of bound" << std::endl;
							system("PAUSE ");
						}
					}
				}
				if ((k*XNEL[m] + (l + 1)) > nel[m]) {
					std::cout << "IEN_2D out of bound" << std::endl;
					system("PAUSE ");
				}
			}
		}
	}
	
	//define the displacement direction of wet surfaces (x=0, y=1, z=2)
	ol[0].dir = 2;
	ol[1].dir = 0;
	ol[2].dir = 1;
	ol[3].dir = 2;

	//define the constant location of the wet surfaces 
	ol[0].location = -SZ / 2;
	ol[1].location = -SX;
	ol[2].location = -SY;
	ol[3].location = SZ / 2;

	ol[0].XNEL = SXNEL;
	ol[0].YNEL = SYNEL;
	ol[1].XNEL = SZNEL;
	ol[1].YNEL = SYNEL;
	ol[2].XNEL = SXNEL;
	ol[2].YNEL = SZNEL;
	ol[3].XNEL = SXNEL;
	ol[3].YNEL = SYNEL;

	delete[] GIDF1;
	delete[] GIDF2;
	delete[] GIDF3;
	delete[] GIDF4;
	
	//std::cout << " " << std::endl;
	return;
}

