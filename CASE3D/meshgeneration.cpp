#include "stdafx.h"
#include "Header.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>

/*
mesh generation explanation:
This code allows both self-generation of mesh and external mesh of the symmetric model of FSP with canopy.
It is controlled by the constant internalmesh in the header file. 
*/

struct meshgenerationstruct meshgeneration() {
	meshgenerationstruct t;
	LOBATTOstruct b;
	LOCAL_NODEstruct c;

	if (internalmesh == 1) {
		int iix, iiy, iiz, i, j, k, ICX, ICY, ICZ;
		int XNEL = round((AX + SX) / XHE);
		int YNEL = round((SY + DY) / YHE);
		int ZNEL = round((BZ * 2 + SZ) / ZHE);
		int SNEL = SXNEL*SYNEL*SZNEL;
		int NELP = XNEL*YNEL*ZNEL;

		int IC = 2;
		for (i = 0; i < XNEL; i++) {
			IC = IC - 1;
			for (j = 0; j < N + 1; j++) {
				IC = IC + 1;
			}
		}
		int XNODE = IC - 1;

		IC = 2;
		for (i = 0; i < YNEL; i++) {
			IC = IC - 1;
			for (j = 0; j < N + 1; j++) {
				IC = IC + 1;
			}
		}
		int YNODE = IC - 1;

		IC = 2;
		for (i = 0; i < ZNEL; i++) {
			IC = IC - 1;
			for (j = 0; j < N + 1; j++) {
				IC = IC + 1;
			}
		}
		int ZNODE = IC - 1;

		//total number of nodes
		int NNODEP = XNODE * YNODE * ZNODE;

		int **IENP = new int*[(N + 1)*(N + 1)*(N + 1)];
		for (i = 0; i < (N + 1)*(N + 1)*(N + 1); i++) {
			IENP[i] = new int[NELP];
		}
		for (i = 0; i < (N + 1)*(N + 1)*(N + 1); i++) {
			for (j = 0; j < NELP; j++) {
				IENP[i][j] = 0;
			}
		}

		c = LOCAL_NODE(N);

		for (iiy = 0; iiy < YNEL; iiy++) {
			for (iiz = 0; iiz < ZNEL; iiz++) {
				for (iix = 0; iix < XNEL; iix++) {
					for (j = 0; j < N + 1; j++) {
						for (k = 0; k < N + 1; k++) {
							for (i = 0; i < N + 1; i++) {
								IENP[c.LNA[i][j][k] - 1][iiy*XNEL*ZNEL + iiz*XNEL + iix] =
									(1 + i + XNODE*k + XNODE*ZNODE*(N + 1 - (j + 1)))
									+ N*iix + N*iiz*XNODE + N*iiy*XNODE*ZNODE;
								//the first line is for node number increased by nodes in that element 
								//the second line is for node number increased by previous elements 
							}
						}
					}
				}
			}
		}

		if (N > 1) {
			b = LOBATTO(N);
		}

		double*XE; double*YE; double*ZE;
		XE = new double[XNEL + 1]; YE = new double[YNEL + 1]; ZE = new double[ZNEL + 1];
		for (i = 0; i < XNEL + 1; i++) {
		}

		for (i = 0; i < XNEL + 1; i++) {
			XE[i] = i*XHE;
		}

		//COMPUTE Y-COORDINATES OF ELEMENT END NODES
		for (i = 0; i < YNEL + 1; i++) {
			YE[i] = i*YHE;
		}

		//COMPUTE Y - COORDINATES OF ELEMENT END NODES
		for (i = 0; i < ZNEL + 1; i++) {
			ZE[i] = i*ZHE;
		}

		double MX, MY, MZ; double dX, dY, dZ;
		//double AXIN[XNEL][N + 1]; double t.AYIN[YNEL][N + 1]; double AZIN[ZNEL][N + 1];
		double** AXIN; 
		//double** t.AYIN; 
		double** AZIN;
		AXIN = new double*[XNEL];
		for (i = 0; i < XNEL; i++) {
			AXIN[i] = new double[N + 1];
		}
		t.AYIN = new double*[YNEL];
		for (i = 0; i < YNEL; i++) {
			t.AYIN[i] = new double[N + 1];
		}
		AZIN = new double*[ZNEL];
		for (i = 0; i < ZNEL; i++) {
			AZIN[i] = new double[N + 1];
		}

		for (i = 0; i < XNEL; i++) {
			MX = 0.5*(XE[i + 1] + XE[i]);
			dX = 0.5*(XE[i + 1] - XE[i]);
			AXIN[i][0] = XE[i];
			if (N > 1) {
				for (j = 0; j < N - 1; j++) {
					AXIN[i][j + 1] = MX + b.Z[j] * dX;
				}
			}
			AXIN[i][N] = XE[i + 1];
		}

		//COMPUTE ELEMENT INTERPOLATION NODE LOCATIONS(Y - DIRECTION)
		for (i = 0; i < YNEL; i++) {
			MY = 0.5*(YE[i + 1] + YE[i]);
			dY = 0.5*(YE[i + 1] - YE[i]);
			t.AYIN[i][0] = YE[i];
			if (N > 1) {
				for (j = 0; j < N - 1; j++) {
					t.AYIN[i][j + 1] = MY + b.Z[j] * dY;
				}
			}
			t.AYIN[i][N] = YE[i + 1];
		}

		//COMPUTE ELEMENT INTERPOLATION NODE LOCATIONS(Z - DIRECTION)
		for (i = 0; i < ZNEL; i++) {
			MZ = 0.5*(ZE[i + 1] + ZE[i]);
			dZ = 0.5*(ZE[i + 1] - ZE[i]);
			AZIN[i][0] = ZE[i];
			if (N > 1) {
				for (j = 0; j < N - 1; j++) {
					AZIN[i][j + 1] = MZ + b.Z[j] * dZ;
				}
			}
			AZIN[i][N] = ZE[i + 1];
		}

		double **GCOORDP = new double*[NNODEP];
		for (i = 0; i < NNODEP; i++) {
			GCOORDP[i] = new double[3];
		}
		int ct = 0;
		for (j = 0; j < YNEL; j++) {
			for (k = 0; k < ZNEL; k++) {
				for (i = 0; i < XNEL; i++) {
					for (ICY = 0; ICY < N + 1; ICY++) {
						for (ICZ = 0; ICZ < N + 1; ICZ++) {
							for (ICX = 0; ICX < N + 1; ICX++) {
								GCOORDP[IENP[c.LNA[ICX][ICY][ICZ] - 1][i + k*XNEL + j*XNEL*ZNEL] - 1][0] = AXIN[i][ICX];
								GCOORDP[IENP[c.LNA[ICX][ICY][ICZ] - 1][i + k*XNEL + j*XNEL*ZNEL] - 1][1] = -t.AYIN[j][N - ICY];
								GCOORDP[IENP[c.LNA[ICX][ICY][ICZ] - 1][i + k*XNEL + j*XNEL*ZNEL] - 1][2] = AZIN[k][ICZ];
							}
						}
					}
				}
			}
		}

		//Set the flag of structural node to be -1 (negative)
		int*PN = new int[NNODEP];
		for (i = 0; i < NNODEP; i++) {
			PN[i] = i + 1; //for node
		}
		int*PE = new int[NELP];
		for (i = 0; i < NELP; i++) {
			PE[i] = i + 1; //for element
		}
		//make the nodes within the structural space to have -1 flag
		if (AX > 0 && BZ > 0) {  //multiple fuid modules
			for (i = 0; i < NELP; i++) {
				for (j = 0; j < NINT*NINT*NINT; j++) {
					/*
					if (GCOORDP[IENP[j][i] - 1][0] > AX && GCOORDP[IENP[j][i] - 1][1] > -SY &&
						GCOORDP[IENP[j][i] - 1][2] > BZ && GCOORDP[IENP[j][i] - 1][2] < BZ + SZ) {
						PN[IENP[j][i] - 1] = -1;
						PE[i] = -1;
					}
					*/
					
					if (GCOORDP[IENP[j][i] - 1][0] - AX > 1e-5 && GCOORDP[IENP[j][i] - 1][1] + SY > 1e-5 &&
						GCOORDP[IENP[j][i] - 1][2] - BZ > 1e-5 && GCOORDP[IENP[j][i] - 1][2] - BZ - SZ < -1e-5) {
						PN[IENP[j][i] - 1] = -1;
						PE[i] = -1;
					}
					
				}
			}
		}
		//for Bleich-Sandler configuration (only one fluid module under the structure)
		else if (AX == 0 && BZ == 0) {
			for (i = 0; i < NELP; i++) {
				for (j = 0; j < NINT*NINT*NINT; j++) {
					if (GCOORDP[IENP[j][i] - 1][0] >= AX && GCOORDP[IENP[j][i] - 1][1] > -SY &&
						GCOORDP[IENP[j][i] - 1][2] >= BZ && GCOORDP[IENP[j][i] - 1][2] <= BZ + SZ) {
						PN[IENP[j][i] - 1] = -1;
						PE[i] = -1;
					}
				}
			}
		}

		//renumber the positive element in array PN 
		ct = 0;
		for (i = 0; i < NNODEP; i++) {
			if (PN[i] > 0) {
				ct += 1;
				PN[i] = ct;
			}
		}
		t.NNODE = ct;
		t.NEL = NELP - SNEL;
		if (AX > 0 && BZ > 0) { //multiple fluid modules
			if ((NNODEP - ct) != (SXNEL*N)*(SYNEL*N)*(SZNEL*N - 1)) {
				std::cout << "new node number is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		else if (AX == 0 && BZ == 0) { //one fluid modules (BS configuration)
			if ((NNODEP - ct) != (SXNEL*N + 1)*(SYNEL*N)*(SZNEL*N + 1)) {
				std::cout << "new node number is wrong" << std::endl;
				system("PAUSE ");
			}
		}
		
		ct = 0;
		//renumber the positive element in array PE
		for (i = 0; i < NELP; i++) {
			if (PE[i] > 0) {
				ct += 1;
				PE[i] = ct;
			}
		}
		if ((NELP - ct) != SXNEL*SYNEL*SZNEL) {
			std::cout << "new element number is wrong" << std::endl;
			system("PAUSE ");
		}

		t.GCOORD = new double*[t.NNODE];
		for (i = 0; i < t.NNODE; i++) {
			t.GCOORD[i] = new double[3];
		}

		//renumber the GCOORD
		for (i = 0; i < NNODEP; i++) {
			if (PN[i] > 0) {
				for (j = 0; j < 3; j++) {
					t.GCOORD[PN[i] - 1][j] = GCOORDP[i][j];
				}
			}
		}

		//renumber the IEN 
		t.IEN = new int*[NINT*NINT*NINT];
		for (i = 0; i < NINT*NINT*NINT; i++) {
			t.IEN[i] = new int[t.NEL];
		}
		for (i = 0; i < NELP; i++) {
			if (PE[i] > 0) {
				for (j = 0; j < NINT*NINT*NINT; j++) {
					t.IEN[j][PE[i] - 1] = PN[IENP[j][i] - 1];
				}
			}
		}

		//offset the node location, make the origin identical with that of structure for X and Z (Y is not adjusted)
		double xoff = -(AX + SX);
		double zoff = -(2 * BZ + SZ) / 2;
		for (i = 0; i < t.NNODE; i++) {
			t.GCOORD[i][0] += xoff;
			t.GCOORD[i][2] += zoff;
		}

		//clean the dynamic array
		for (i = 0; i < (N + 1)*(N + 1)*(N + 1); i++) {
			delete[] IENP[i];
		}
		delete[] IENP;
		delete[] XE; delete[] YE; delete[] ZE;
		for (i = 0; i < XNEL; i++) {
			delete[] AXIN[i];
		}
		delete[] AXIN;
		/*
		for (i = 0; i < YNEL; i++) {
			delete[] t.AYIN[i];
		}
		delete[] t.AYIN;
		*/
		for (i = 0; i < ZNEL; i++) {
			delete[] AZIN[i];
		}
		delete[] AZIN;
		for (i = 0; i < NNODEP; i++) {
			delete[] GCOORDP[i];
		}
		delete[] GCOORDP;
		delete[] PN; delete[] PE;
	}

	else {   //external mesh from Gmsh
		int i, j;
		double ptholder[1 + NINT*NINT*NINT];
		for (i = 0; i < 1 + NINT*NINT*NINT; i++) {
			ptholder[i] = 0.0;
		}
		int row, col;
		double x;
		std::string lineA;
		std::string filename;
		std::cout << "reading the mesh file: " << std::endl;
		std::ifstream myfile("mesh_N=1_0.1_abaqus2_sym.txt");
		//std::ifstream myfile("Abaqus2mesh_debug.txt");
		if (!myfile) {
			std::cout << "can not open the mesh file" << std::endl;
			system("PAUSE ");
		}

		//int NEL = AXNEL*AYNEL*AZNEL + BXNEL*BYNEL*BZNEL + CXNEL*CYNEL*CZNEL + DXNEL*DYNEL*DZNEL + (AXNEL + BXNEL + CXNEL)*BZNEL*AYNEL;
		int NEL = round((AX + SX) / XHE)*round((DY + SY) / YHE)*round((SZ + 2 * BZ) / ZHE) - SXNEL*SYNEL*SZNEL;
		t.IEN = new int*[NINT*NINT*NINT];
		for (i = 0; i < NINT*NINT*NINT; i++) {
			t.IEN[i] = new int[NEL];
		}
		for (i = 0; i < NINT*NINT*NINT; i++) {
			for (j = 0; j < NEL; j++) {
				t.IEN[i][j] = 0;
			}
		}

		//time counter
		std::clock_t start;
		double duration;
		start = std::clock();
		/* Your algorithm here */

		int NEL_holder;
		t.NEL = 0;
		while (myfile.good()) {
			row = 0;
			while (std::getline(myfile, lineA)) { //read the file line by line
				std::istringstream streamA(lineA);
				col = 0;
				NEL_holder = 0;
				for (i = 0; i < 1 + NINT*NINT*NINT; i++) {
					ptholder[i] = 0.0;
				}
				while (streamA >> std::skipws >> x) {
					ptholder[col] = x;
					if (col >= 1 + NINT*NINT*NINT) {
						std::cout << " ptholder takes more values than it should" << std::endl;
						system("PAUSE ");
					}
					if (row == 0) {
						t.NNODE = ptholder[0];
						t.GCOORD = new double*[t.NNODE];
						for (i = 0; i < t.NNODE; i++) {
							t.GCOORD[i] = new double[3];
						}
						for (i = 0; i < t.NNODE; i++) {
							for (j = 0; j < 3; j++) {
								t.GCOORD[i][j] = 0.0;
							}
						}
						std::cout << "GCOORD is initialized" << std::endl;
					}
					if (row != 0 && row <= t.NNODE) {
						t.GCOORD[row - 1][0] = ptholder[0];
						t.GCOORD[row - 1][1] = ptholder[1];
						t.GCOORD[row - 1][2] = ptholder[2];
					}
					if (ptholder[NINT*NINT*NINT] != 0) {
						for (i = 1; i < 1 + NINT*NINT*NINT; i++) {
							t.IEN[i - 1][t.NEL] = ptholder[i];
						}
						NEL_holder = 1;
					}
					col++;
				}
				t.NEL += NEL_holder;
				row++;
			}
		}

		duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;

		std::cout << "printf: " << duration << '\n';

		if (t.NEL != NEL) {
			std::cout << "the prediction NEL is wrong" << std::endl;
			system("PAUSE ");
		}
	}

	return t;
}