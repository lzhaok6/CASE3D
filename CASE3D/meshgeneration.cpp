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
#include <iostream>
#include <vector>

/*
mesh generation explanation:
This code allows both self-generation of mesh and external mesh of the symmetric model of FSP with canopy.
It is controlled by the constant internalmesh in the header file. 
*/

struct meshgenerationstruct meshgeneration() {
	meshgenerationstruct t;
	LOBATTOstruct b;
	LOCAL_NODEstruct c;
	extern OWETSURF ol[owsfnumber];
	extern NRBSURF nr[nrbsurfnumber]; 
	c = LOCAL_NODE(N);

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

		//c = LOCAL_NODE(N);

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

	else {   //Import external mesh from Gmsh
		//We are assuming there is only one wetted surface. Thus, ol[0] is used throughout the loop below. 
		int i, j, k, l, m, n, e;
		std::string filename;
		std::cout << "reading the mesh file: " << std::endl;
		std::ifstream infile("frigate_N=1.msh");
		if (!infile) {
			std::cout << "can not open the mesh file" << std::endl;
			system("PAUSE ");
		}
		//read the file 
		int nodestart = 0;
		int nodeend = 0;
		int ele_line = 0;
		std::vector<std::vector<std::string>> output;
		int ct = -1;
		std::string csvLine;
		int physicalgroups = 1; //starting from the first physical group
		int endfile = 0;
		std::vector<int> phygrp_start; //the starting line number of the physical group
		int elestart = 0; //the flag denote the start of element connectivity definition
		
		//std::clock_t start;
		//start = std::clock();
		std::string csvElement;
		while (getline(infile, csvLine))
		{
			ct = ct + 1; //the current line number (starting from 0)
			std::istringstream csvStream(csvLine);
			std::vector<std::string> csvColumn;
			//std::string csvElement;
			while (getline(csvStream, csvElement, ' '))
			{
				csvColumn.push_back(csvElement);
			}
			output.push_back(csvColumn);
			if (csvColumn[0] == "$EndElements") {
				endfile = ct;
				break;
			}	
			if (csvColumn[0] == "$Nodes") {
				nodestart = ct + 2; //the node starts from the next line
			}
			if (csvColumn[0] == "$EndNodes") {
				nodeend = ct - 1; //the node starts from the next line
			}
			if (csvColumn[0] == "$Elements") {
				ele_line = ct + 2; //the line corresponding to the start of element connectivity definition
				phygrp_start.push_back(ele_line);
				elestart = 1;
			}
			//Get the information of physical groups
			if (elestart == 1 && ct >= ele_line) {
				if (csvColumn[3] == std::to_string(physicalgroups)) {
					//do nothing
				}
				else {
					physicalgroups += 1;
					//push the line number of the starting of the current physical group
					phygrp_start.push_back(ct);
				}
			}
			if (ct % 20000 == 0) {
				std::cout << ct << std::endl; //output which line is being read
			}
		}
		//double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC * 1000;

		//The physical groups except for the last one are the surface mesh (NINT*NINT). 
		//The last physical group is the volumn mesh (NINT*NINT*NINT) connectivity matrix
		t.NNODE = stoi(output[nodestart - 1][0]);
		t.NEL = endfile - phygrp_start[physicalgroups - 1];
		//Define connectivity matrix in the volumn mesh
		t.IEN = new int*[NINT*NINT*NINT]; 
		for (i = 0; i < NINT*NINT*NINT;i++) {
			t.IEN[i] = new int[t.NEL]; 
		}
		ct = 0;

		int nelstart = 5; //the column number where the element node definition starts
		if (N > 1) {
			nelstart = 6; 
		}
		for (i = phygrp_start[physicalgroups - 1]; i < endfile; i++) {
			for (j = nelstart; j < nelstart + NINT*NINT*NINT; j++) {
				t.IEN[j - nelstart][ct] = stoi(output[i][j]);
			}
			ct += 1;
		}

		//Define the rest of surface boundaries (for boundary condition definition)
		//std::vector<std::vector<std::vector<int>>> BCIEN;
		for (i = 0; i < physicalgroups - 1; i++) {
			std::vector<std::vector<int>>iens;
			for (j = phygrp_start[i]; j < phygrp_start[i + 1]; j++) {
				std::vector<int>local;
				for (k = 0; k < NINT*NINT; k++) {
					local.push_back(stoi(output[j][nelstart + k]));
					//BCIEN[i][j - phygrp_start[i]].push_back(stoi(output[j][5 + k]));
				}
				iens.push_back(local);
			}
			t.BCIEN.push_back(iens); 
		}

		//Define node coordinates
		t.GCOORD = new double*[t.NNODE];
		for (i = 0; i < t.NNODE; i++) {
			t.GCOORD[i] = new double[3];
		}

		for (i = nodestart; i < nodestart + t.NNODE; i++) {
			for (j = 0; j < 3; j++) {
				t.GCOORD[i - nodestart][j] = stod(output[i][1 + j]);
			}
		}
		//std::cout << " " << std::endl;
		
		//Since the current mesh is a little twisted (e.g, y and z directions are switched)
		double hd = 0.0; //place holder
		if (debug == 1) {
			//exchange y and z axis and revert the z coordinate
			for (i = 0; i < t.NNODE; i++) {
				hd = -t.GCOORD[i][1]; //z 
				t.GCOORD[i][1] = t.GCOORD[i][2];
				t.GCOORD[i][2] = hd;
				//drag the free surface to y=0 
				t.GCOORD[i][1] -= 6.0;
			}
		}

		int localnode[NINT*NINT];
		std::vector<std::vector<int>>global3d;
		//The following loop only attempts to extract the pattern of localnode[] from one corresponding global element since the assumption here is that all the elements in all physical groups
		//follow the same pattern. However, that is not rigorous enough since the corresponding local node for all 2D element in physical group do not necessary be the same. 
		for (i = 0; i < physicalgroups - 1; i++) {
			for (j = 0; j < phygrp_start[i + 1] - phygrp_start[i]; j++) { //phygrp_start[i + 1] - phygrp_start[i] is the element number in the ith physical group
				for (k = 0; k < t.NEL;k++) { 
					ct = 0;
					for (l = 0; l < NINT*NINT; l++) {
						for (m = 0; m < NINT*NINT*NINT; m++) {
							if (t.BCIEN[i][j][l] == t.IEN[m][k]) {
								localnode[l] = m + 1; //store the corresponding local node in BCIEN in global element connectivity matrix
								ct += 1;
							}
						}
					}
					if (ct==NINT*NINT) { //found the corresponding 3D element! 
						goto endloop; 
					}
				}
			}
		}
	endloop: 

		//=================extract the 2D local distribution of the local nodes LNA2D===================//
		//The section below trys to associate the relationship between localnode and wet surface local node distribution called nr[0].LNA_2D which is later used to separate high-order element to linear element
		//The local surface (left/right/front/behind/top/bottom) could be identified.
		//The assumption here is that all the wet surface corresponds to the same face in its corresponding local element (Check out Evernote: How to determine the wet elements on the fluid side). 
		int flag; 
		//int LNA_2D[NINT][NINT];

		//If the wet surface is the left face of the local element (i.e., i=0). 
		ct = 0; 
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[0][j][k]) {
						nr[0].LNA_2D[j][k] = l + 1;
						ol[0].LNA_2D[j][k] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = 0;
			goto endextraction; //If the face is found, jump to the end of process
		}
		//right face (i.e. i=N)
		ct = 0;
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[N][j][k]) {
						nr[0].LNA_2D[j][k] = l + 1;
						ol[0].LNA_2D[j][k] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = N;
			goto endextraction;
		}
		//back face (k=0)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[i][j][0]) {
						nr[0].LNA_2D[i][j] = l + 1;
						ol[0].LNA_2D[i][j] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = 0;
			goto endextraction;
		}
		//front face (k=N)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[i][j][N]) {
						nr[0].LNA_2D[i][j] = l + 1;
						ol[0].LNA_2D[i][j] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = N;
			goto endextraction;
		}
		//bottom face (j=0)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[i][0][k]) {
						nr[0].LNA_2D[i][k] = l + 1;
						ol[0].LNA_2D[i][k] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = 0;
			goto endextraction;
		}
		//top face (j=N)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[l] == c.LNA[i][N][k]) {
						nr[0].LNA_2D[i][k] = l + 1;
						ol[0].LNA_2D[i][k] = l + 1;
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			flag = N;
			goto endextraction;
		}

	endextraction:
		//================================end extraction===================================//
		//=============separate the high-order element into linear elements and obtain the normal direction================//
		int pys_num = 2; //the physical group number that corresponds to the wet surface (physical group 3)
		int elenum = phygrp_start[pys_num + 1] - phygrp_start[pys_num]; //number of high-order on the wet surface 
		int** IEN_py; //the linear element connectivity matrix of 2D physical groups (boundaries) separated from the original 2D high-order elements
		IEN_py = new int*[4];
		for (i = 0; i < 4; i++) {
			IEN_py[i] = new int[elenum];
		}
		if (mappingalgo == 2) {
			//Originally, we break the high-order element into linear elements. However, that would induce a lot of extra computation
			//Thus, we changed our plan to use the linear element contructed by the four corner nodes of the high-order element
			ct = 0;
			for (e = 0; e < elenum; e++) {
				if (flag == N) { //counter-clockwise
					IEN_py[0][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[0][0] - 1]; //oriente the nodes so that the normal direction is pointing out of the element
					ol[0].LNA_algo2[0][0] = 1;
					IEN_py[1][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[N][0] - 1];
					ol[0].LNA_algo2[1][0] = 2;
					IEN_py[2][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[N][N] - 1];
					ol[0].LNA_algo2[1][1] = 3;
					IEN_py[3][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[0][N] - 1];
					ol[0].LNA_algo2[0][1] = 4;
					ct += 1;
				}
				else if (flag == 0) { //clockwise
					IEN_py[0][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[0][0] - 1]; //oriente the nodes so that the normal direction is out of the element
					ol[0].LNA_algo2[0][0] = 1;
					IEN_py[1][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[0][N] - 1];
					ol[0].LNA_algo2[0][1] = 2;
					IEN_py[2][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[N][N] - 1];
					ol[0].LNA_algo2[1][1] = 3;
					IEN_py[3][ct] = t.BCIEN[pys_num][e][ol[0].LNA_2D[N][0] - 1];
					ol[0].LNA_algo2[1][0] = 4;
					ct += 1;
				}
			}
			if (ct != elenum) {
				std::cout << "Not all high-order elements are turned into wrapping linear IEN_py (connectivity)" << std::endl;
				system("PAUSE ");
			}
			
			std::vector<int>ele_num; //track the element number of the wetted surface elements
			for (e = 0; e < elenum; e++) {
				flag = 1;
				for (i = 0; i < 4; i++) {
					if (abs(t.GCOORD[IEN_py[i][e] - 1][1]) < 1e-1) {
						flag = 0;
					}
				}
				if (flag == 1) { //Then this element could be included the wetted surface element list
					ele_num.push_back(e);
				}
			}
			//ele_num is the elements for high-order element, we need to extract the corresponding linear elements from that. 

			//Store the total number of element on wetted surface in FSNEL
			ol[0].FSNEL = ele_num.size();

			ol[0].IEN_gb = new int*[NINT*NINT]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < NINT*NINT; i++) {
				ol[0].IEN_gb[i] = new int[ol[0].FSNEL];
			}
			for (i = 0; i < ol[0].FSNEL; i++) {
				for (j = 0; j < NINT*NINT; j++) {
					ol[0].IEN_gb[j][i] = t.BCIEN[pys_num][ele_num[i]][j];
				}
			}

			ol[0].IEN_algo2 = new int*[4]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < 4; i++) {
				ol[0].IEN_algo2[i] = new int[ol[0].FSNEL];
			}
			for (i = 0; i < ol[0].FSNEL; i++) {
				for (j = 0; j < 4; j++) {
					ol[0].IEN_algo2[j][i] = IEN_py[j][ele_num[i]];
				}
			}

			//Derive the IEN_2D to write the MpCCI model file (basically renumbering the node in IEN_algo2 to be recognized by MpCCI)
			ol[0].IEN_2D = new int*[4]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < 4; i++) {
				ol[0].IEN_2D[i] = new int[ol[0].FSNEL];
			}
			ct = 0; //count the node number assigned
			std::vector<int>dummy;
			for (i = 0; i < ol[0].FSNEL; i++) { //loop through each element
				for (j = 0; j < 4; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < 4; l++) {
							if (ol[0].IEN_algo2[l][k] == ol[0].IEN_algo2[j][i]) { //If this node has already been assigned, use the same numbering
								ol[0].IEN_2D[j][i] = ol[0].IEN_2D[l][k];
								flag = 0; //turn off the flag to assgin new number
							}
							else {
								//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
							}
						}
					}
					if (flag == 1) {
						ct += 1;
						dummy.push_back(ol[0].IEN_algo2[j][i]); //associate the local 2D node with the global node numbering 
						ol[0].IEN_2D[j][i] = ct;
					}
				}
			}
			ol[0].GIDNct = dummy.size();
			ol[0].GIDN = new int[ol[0].GIDNct];
			for (i = 0; i < ol[0].GIDNct; i++) {
				ol[0].GIDN[i] = dummy[i];
			}

			ol[0].norm = new double*[ol[0].FSNEL]; //store the normal direction of linear elements on the wetted surface
			for (i = 0; i < ol[0].FSNEL; i++) {
				ol[0].norm[i] = new double[3];
			}
			//Obtain the normal direction unit vector of the newly separated elements (numbering from 0 to elenum)
			//The normal direction calculation might be wrong. We need to validate it using a FSP mesh generated by BOLT.
			double ax, ay, az; double bx, by, bz;
			double n1, n2, n3; double absn;
			for (i = 0; i < ol[0].FSNEL; i++) {
				ax = t.GCOORD[ol[0].IEN_algo2[1][i] - 1][0] - t.GCOORD[ol[0].IEN_algo2[0][i] - 1][0];
				ay = t.GCOORD[ol[0].IEN_algo2[1][i] - 1][1] - t.GCOORD[ol[0].IEN_algo2[0][i] - 1][1];
				az = t.GCOORD[ol[0].IEN_algo2[1][i] - 1][2] - t.GCOORD[ol[0].IEN_algo2[0][i] - 1][2];
				bx = t.GCOORD[ol[0].IEN_algo2[2][i] - 1][0] - t.GCOORD[ol[0].IEN_algo2[1][i] - 1][0];
				by = t.GCOORD[ol[0].IEN_algo2[2][i] - 1][1] - t.GCOORD[ol[0].IEN_algo2[1][i] - 1][1];
				bz = t.GCOORD[ol[0].IEN_algo2[2][i] - 1][2] - t.GCOORD[ol[0].IEN_algo2[1][i] - 1][2];
				n1 = ay*bz - az*by; n2 = az*bx - ax*bz; n3 = ax*by - ay*bx;
				absn = sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
				ol[0].norm[i][0] = n1 / absn; ol[0].norm[i][1] = n2 / absn; ol[0].norm[i][2] = n3 / absn;
			}
			//check if the separated element has the same normal direction as the original mesh (N=1)

			//write the model file for MpCCI here
			std::ofstream myfile;
			myfile.open("model.txt");
			myfile << "EF wetsurface 3 1" << std::endl;
			myfile << "NODES " << ol[0].GIDNct << std::endl;
			for (i = 0; i < ol[0].GIDNct; i++) {
				myfile << i << " " << t.GCOORD[ol[0].GIDN[i] - 1][0] << " " << t.GCOORD[ol[0].GIDN[i] - 1][1] << " " << t.GCOORD[ol[0].GIDN[i] - 1][2] << " " << std::endl;
			}
			myfile << "ELEMENTS " << ol[0].FSNEL << std::endl;
			//output connectivity matrix
			for (i = 0; i < ol[0].FSNEL; i++) {
				myfile << i;
				for (j = 0; j < 4; j++) {
					myfile << " " << ol[0].IEN_2D[j][i] - 1; //node numbering starts from 0 in model file
				}
				myfile << std::endl;
			}
		}

		//==============Extract the connectivity matrix on NRB surface (glue all physical groups corresponding to the NRB together)==============//
		//All physical groups except for 3
		nr[0].NEL_nrb = (phygrp_start[physicalgroups - 1] - phygrp_start[0]) - elenum; //total element number. 
		//Please note that elenum here is for the physical group 3 which includes both free surface and wetted surface. 
		nr[0].IEN_gb = new int*[NINT*NINT]; //NRBELE_ARR
		for (i = 0; i < NINT*NINT; i++) {
			nr[0].IEN_gb[i] = new int[nr[0].NEL_nrb]; 
		}
		ct = 0; 
		for (i = 0; i < physicalgroups; i++) {
			if (i != pys_num) { //pys_num is the physical group number corresponding to the wetted surface
				for (e = 0; e < phygrp_start[i + 1] - phygrp_start[i]; e++) {
					for (j = 0; j < NINT*NINT; j++) {
						nr[0].IEN_gb[j][ct] = t.BCIEN[i][e][j];
					}
					ct += 1; 
				}
			}
		}
		if (ct!=nr[0].NEL_nrb) {
			std::cout << "The NRB element storage is wrong" << std::endl;
			system("PAUSE ");
		}

		//Obtain NRBA to store all the nodes on NRB surfaces by looping through all the NRB elements and filter out the points
		ct = 0; //count the node number assigned
		std::vector<int>dummy2;
		for (i = 0; i < nr[0].NEL_nrb; i++) { //loop through each element
			for (j = 0; j < NINT*NINT; j++) { //the nodes in current element
				flag = 1; //Initiate the flag to 1 
				for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
					for (l = 0; l < NINT*NINT; l++) {
						if (nr[0].IEN_gb[l][k] == nr[0].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
							flag = 0; //turn off the flag to assgin new number
						}
						else {
							//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
						}
					}
				}
				if (flag == 1) {
					ct += 1;
					dummy2.push_back(nr[0].IEN_gb[j][i]); //associate the local 2D node with the global node numbering 
				}
			}
		}
		nr[0].NRBNODE = dummy2.size();
		nr[0].NRBA = new int[nr[0].NRBNODE];
		for (i = 0; i < nr[0].NRBNODE; i++) {
			nr[0].NRBA[i] = dummy2[i];
		}
		//obtain the normal vector corresponding to every

	}

	
	return t;
}
