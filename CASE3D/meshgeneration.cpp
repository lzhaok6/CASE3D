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

	//We are assuming there is only one wetted surface. Thus, ol[0] is used throughout the loop below. 
	int i, j, k, l, m, n, e, z;
	std::string filename;
	std::cout << "reading the mesh file: " << std::endl;
	std::cout << "Have you configured the mesh file name correctly? If yes, hit Enter to proceed" << std::endl;
	system("PAUSE "); 
	std::ifstream infile("FSP_N=2.msh");
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
	for (i = 0; i < NINT*NINT*NINT; i++) {
		t.IEN[i] = new int[t.NEL];
	}
	ct = 0;

	/*
	int nelstart = 5; //the column number where the element node definition starts
	if (N > 1) {
		nelstart = 6;
	}
	*/
	//Find where does the node definition starts in physical groups
	//This method is based on the observation that the node starts after the first zero.
	ct = 0;
	while (stoi(output[phygrp_start[0]][ct]) != 0) {
		ct += 1;
	}
	int nelstart = ct + 1; 
	ct = 0;
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

	//Need to be modified, localnode needs to corresponds each physical group! we cannot make the assumption anymore
	//int localnode[NINT*NINT]; 
	int **localnode; 
	localnode = new int*[physicalgroups - 1];
	for (i = 0; i < physicalgroups - 1; i++) {
		localnode[i] = new int[NINT*NINT];
	}
	//The following loop only attempts to extract the pattern of localnode[] from one corresponding global element since the assumption here is that all the elements in all physical groups
	//follow the same pattern. However, that is not rigorous enough since the corresponding local node for all 2D element in physical group do not necessary be the same. 
	for (i = 0; i < physicalgroups - 1; i++) {
		for (j = 0; j < phygrp_start[i + 1] - phygrp_start[i]; j++) { //phygrp_start[i + 1] - phygrp_start[i] is the element number in the ith physical group
			for (k = 0; k < t.NEL; k++) {
				ct = 0;
				for (l = 0; l < NINT*NINT; l++) {
					for (m = 0; m < NINT*NINT*NINT; m++) {
						if (t.BCIEN[i][j][l] == t.IEN[m][k]) {
							localnode[i][l] = m + 1; //store the corresponding local node in BCIEN in global element connectivity matrix
							ct += 1;
						}
					}
				}
				if (ct == NINT*NINT) { //found the corresponding 3D element! 
					goto endloop;
				}
			}
		}
	endloop:;
	}


	//=================extract the 2D local distribution of the local nodes LNA2D===================//
	//The section below trys to associate the relationship between localnode and wet surface local node distribution called nr[0].LNA_2D which is later used to separate high-order element to linear element
	//The local surface (left/right/front/behind/top/bottom) could be identified.
	//The assumption here is that all the wet surface corresponds to the same face in its corresponding local element (Check out Evernote: How to determine the wet elements on the fluid side). 
	//int flag;
	//int LNA_norm[4];
	
	//The physical group number of wetted surface and NRB surfaces (Needs to be set up by the user)
	int wt_pys_num[4] = { 0,1,2,3 };  //the physical group number that corresponds to the wet surface (physical group 3)
	int nrb_pys_num[4] = { 4,5,6,7 };
	//If the wet surface is the left face of the local element (i.e., i=0). 
	std::cout << "Have you configured the wetted surface and NRB surface physical group number for this mesh? If so, hit Enter to proceed." << std::endl;
	//system("PAUSE ");
	
	int wet = 0; 
	int nrb = 0; 
	int wet_ct = 0;
	int nrb_ct = 0; 
	int wt_pys_size = sizeof(wt_pys_num) / sizeof(*wt_pys_num); //the number of wetted surfaces 
	int nrb_pys_size = sizeof(nrb_pys_num) / sizeof(*nrb_pys_num); //the number of nrb surfaces

	for (z = 0; z < physicalgroups - 1; z++) { //all physical groups except for connectivity
		//First determine if this is a wetted surface physical group or nrb physical group
		wet = 0; nrb = 0;
		for (i = 0; i < wt_pys_size; i++) {
			if (z == wt_pys_num[i]) {
				wet = 1; //this physical group corresponds to wetted surface
			}
		}
		for (i = 0; i < nrb_pys_size; i++) {
			if (z == nrb_pys_num[i]) {
				nrb = 1; //this physical group corresponds to wetted surface
			}
		}
		if (wet == 1 && nrb == 1) {
			std::cout << "a physical group cannot be both wetted surface and nrb!";
			system("PAUSE ");
		}
		ct = 0;
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[0][j][k]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[j][k] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[0][j][k];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[j][k];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[j][k] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[0][j][k];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[j][k];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_JB2D[0] = 1;  //Take a look at the header file for more information of LNA_JB2D
				nr[nrb_ct].LNA_JB2D[1] = 5;
				nr[nrb_ct].LNA_JB2D[2] = 8;
				nr[nrb_ct].LNA_JB2D[3] = 4;
				nrb_ct += 1; 
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_JB2D[0] = 1;
				ol[wet_ct].LNA_JB2D[1] = 5;
				ol[wet_ct].LNA_JB2D[2] = 8;
				ol[wet_ct].LNA_JB2D[3] = 4;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[0][1] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[1][0] = 4;
				ol[wet_ct].Jacob_face[0] = 1; //y coordinate
				ol[wet_ct].Jacob_face[1] = 2; //z coordinate
				wet_ct += 1;
			}
			goto endextraction; //If the face is found, jump to the end of process
		}
		//right face (i.e. i=N)
		ct = 0;
		for (j = 0; j < NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[N][j][k]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[j][k] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[N][j][k];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[j][k];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[j][k] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[N][j][k];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[j][k];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_JB2D[0] = 2;
				nr[nrb_ct].LNA_JB2D[1] = 3;
				nr[nrb_ct].LNA_JB2D[2] = 7;
				nr[nrb_ct].LNA_JB2D[3] = 6;
				nrb_ct += 1;
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_JB2D[0] = 2;
				ol[wet_ct].LNA_JB2D[1] = 3;
				ol[wet_ct].LNA_JB2D[2] = 7;
				ol[wet_ct].LNA_JB2D[3] = 6;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[1][0] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[0][1] = 4;
				ol[wet_ct].Jacob_face[0] = 1; //y coordinate
				ol[wet_ct].Jacob_face[1] = 2; //z coordinate
				wet_ct += 1;
			}
			goto endextraction;
		}
		//back face (k=0)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[i][j][0]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[i][j] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[i][j][0];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][j];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[i][j] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[i][j][0];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][j];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_JB2D[0] = 1;
				nr[nrb_ct].LNA_JB2D[1] = 4;
				nr[nrb_ct].LNA_JB2D[2] = 3;
				nr[nrb_ct].LNA_JB2D[3] = 2;
				nrb_ct += 1; 
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_JB2D[0] = 1;
				ol[wet_ct].LNA_JB2D[1] = 4;
				ol[wet_ct].LNA_JB2D[2] = 3;
				ol[wet_ct].LNA_JB2D[3] = 2;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[0][1] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[1][0] = 4;
				ol[wet_ct].Jacob_face[0] = 0; //x coordinate
				ol[wet_ct].Jacob_face[1] = 1; //y coordinate
				wet_ct += 1; 
			}
			goto endextraction;
		}
		//front face (k=N)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (j = 0; j < NINT; j++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[i][j][N]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[i][j] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[i][j][N];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][j];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[i][j] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[i][j][N];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][j];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_JB2D[0] = 5;
				nr[nrb_ct].LNA_JB2D[1] = 6;
				nr[nrb_ct].LNA_JB2D[2] = 7;
				nr[nrb_ct].LNA_JB2D[3] = 8;
				nrb_ct += 1; 
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_JB2D[0] = 5;
				ol[wet_ct].LNA_JB2D[1] = 6;
				ol[wet_ct].LNA_JB2D[2] = 7;
				ol[wet_ct].LNA_JB2D[3] = 8;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[1][0] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[0][1] = 4;
				ol[wet_ct].Jacob_face[0] = 0; //x coordinate
				ol[wet_ct].Jacob_face[1] = 1; //y coordinate
				wet_ct += 1;
			}
			goto endextraction;
		}
		//bottom face (j=0)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[i][0][k]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[i][k] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[i][0][k];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][k];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[i][k] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[i][0][k];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][k];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_JB2D[0] = 1;
				nr[nrb_ct].LNA_JB2D[1] = 2;
				nr[nrb_ct].LNA_JB2D[2] = 6;
				nr[nrb_ct].LNA_JB2D[3] = 5;
				nrb_ct += 1; 
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_JB2D[0] = 1;
				ol[wet_ct].LNA_JB2D[1] = 2;
				ol[wet_ct].LNA_JB2D[2] = 6;
				ol[wet_ct].LNA_JB2D[3] = 5;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[1][0] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[0][1] = 4;
				ol[wet_ct].Jacob_face[0] = 0; //x coordinate
				ol[wet_ct].Jacob_face[1] = 2; //z coordinate
				wet_ct += 1;
			}
			goto endextraction;
		}
		//top face (j=N)
		ct = 0;
		for (i = 0; i < NINT; i++) {
			for (k = 0; k < NINT; k++) {
				for (l = 0; l < NINT*NINT; l++) {
					if (localnode[z][l] == c.LNA[i][N][k]) {
						if (wet == 1) {
							ol[wet_ct].LNA_2D[i][k] = l + 1;
							ol[wet_ct].FP[ct] = c.LNA[i][N][k];
							ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][k];
						}
						else if (nrb == 1) {
							nr[nrb_ct].LNA_2D[i][k] = l + 1;
							nr[nrb_ct].DP[ct] = c.LNA[i][N][k];
							nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][k];
						}
						ct += 1;
					}
				}
			}
		}
		if (ct == NINT*NINT) {
			std::cout << " " << std::endl;
			if (nrb == 1) {
				nr[nrb_ct].LNA_norm[0] = nr[nrb_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				nr[nrb_ct].LNA_norm[1] = nr[nrb_ct].LNA_2D[0][N];
				nr[nrb_ct].LNA_norm[2] = nr[nrb_ct].LNA_2D[N][N];
				nr[nrb_ct].LNA_norm[3] = nr[nrb_ct].LNA_2D[N][0];
				nr[nrb_ct].LNA_JB2D[0] = 4;
				nr[nrb_ct].LNA_JB2D[1] = 8;
				nr[nrb_ct].LNA_JB2D[2] = 7;
				nr[nrb_ct].LNA_JB2D[3] = 3;
				nrb_ct += 1;
			}
			if (wet == 1) {
				ol[wet_ct].LNA_norm[0] = ol[wet_ct].LNA_2D[0][0]; //make sure the orientation is counter clockwise
				ol[wet_ct].LNA_norm[1] = ol[wet_ct].LNA_2D[0][N];
				ol[wet_ct].LNA_norm[2] = ol[wet_ct].LNA_2D[N][N];
				ol[wet_ct].LNA_norm[3] = ol[wet_ct].LNA_2D[N][0];
				ol[wet_ct].LNA_JB2D[0] = 4;
				ol[wet_ct].LNA_JB2D[1] = 8;
				ol[wet_ct].LNA_JB2D[2] = 7;
				ol[wet_ct].LNA_JB2D[3] = 3;
				ol[wet_ct].LNA_algo2[0][0] = 1;
				ol[wet_ct].LNA_algo2[0][1] = 2;
				ol[wet_ct].LNA_algo2[1][1] = 3;
				ol[wet_ct].LNA_algo2[1][0] = 4;
				ol[wet_ct].Jacob_face[0] = 0; //x coordinate
				ol[wet_ct].Jacob_face[1] = 2; //z coordinate
				wet_ct += 1; 
			}
			goto endextraction;
		}
	endextraction:;
	}
	if (nrb_ct != nrbsurfnumber || wet_ct != owsfnumber) {
		std::cout << "Not all wetted surface and NRB are captured" << std::endl; 
		system("PAUSE ");
	}
	//================================end extraction===================================//

	//=============separate the high-order element into linear elements and obtain the normal direction================//
	//int wt_pys_num = 2; //the physical group number that corresponds to the wet surface (physical group 3)
	if (mappingalgo == 2) {
		for (z = 0; z < wt_pys_size; z++) { //loop through each wetted surface
			int elenum = phygrp_start[wt_pys_num[z] + 1] - phygrp_start[wt_pys_num[z]]; //number of high-order on the wet surface s
			ol[z].IEN_py = new int*[4];
			for (i = 0; i < 4; i++) {
				ol[z].IEN_py[i] = new int[elenum];
			}
			//Originally, we break the high-order element into linear elements. However, that would induce a lot of extra computation
			//Thus, we changed our plan to use the linear element contructed by the four corner nodes of the high-order element
			ct = 0;
			for (e = 0; e < elenum; e++) {
				ol[z].IEN_py[0][ct] = t.BCIEN[wt_pys_num[z]][e][ol[z].LNA_norm[0] - 1];
				ol[z].IEN_py[1][ct] = t.BCIEN[wt_pys_num[z]][e][ol[z].LNA_norm[1] - 1];
				ol[z].IEN_py[2][ct] = t.BCIEN[wt_pys_num[z]][e][ol[z].LNA_norm[2] - 1];
				ol[z].IEN_py[3][ct] = t.BCIEN[wt_pys_num[z]][e][ol[z].LNA_norm[3] - 1];
				ct += 1;
			}
			if (ct != elenum) {
				std::cout << "Not all high-order elements are turned into wrapping linear ol[z].IEN_py (connectivity)" << std::endl;
				system("PAUSE ");
			}
			
			//track the element number of the wetted surface elements
			//If any node of the element is on the free surface, then this element is excluded from the wetted surface.
			std::cout << "The wetted surface element identification needs to be changed before proceeding with FSP" << std::endl;
			//system("PAUSE ");
			int flag;
			std::vector<int>ele_num; 
			for (e = 0; e < elenum; e++) {
				flag = 1;
				/*
				for (i = 0; i < 4; i++) {
					if (abs(t.GCOORD[ol[z].IEN_py[i][e] - 1][1]) < 1e-1) {
						flag = 0;
					}
				}
				*/
				if (flag == 1) { //Then this element could be included the wetted surface element list
					ele_num.push_back(e);
				}
			}
			//ele_num is the elements for high-order element, we need to extract the corresponding linear elements from that. 
			
			//Store the total number of element on wetted surface in FSNEL
			ol[z].FSNEL = ele_num.size();

			ol[z].IEN_gb = new int*[NINT*NINT]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < NINT*NINT; i++) {
				ol[z].IEN_gb[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < NINT*NINT; j++) {
					ol[z].IEN_gb[j][i] = t.BCIEN[wt_pys_num[z]][ele_num[i]][j];
				}
			}

			ol[z].IEN_algo2 = new int*[4]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < 4; i++) {
				ol[z].IEN_algo2[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < 4; j++) {
					ol[z].IEN_algo2[j][i] = ol[z].IEN_py[j][ele_num[i]];
				}
			}

			//Derive the IEN_2D to write the MpCCI model file (basically renumbering the node in IEN_algo2 to be recognized by MpCCI)
			ol[z].IEN_2D = new int*[4]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < 4; i++) {
				ol[z].IEN_2D[i] = new int[ol[z].FSNEL];
			}
			ct = 0; //count the node number assigned
			std::vector<int>dummy;
			for (i = 0; i < ol[z].FSNEL; i++) { //loop through each element
				for (j = 0; j < 4; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < 4; l++) {
							if (ol[z].IEN_algo2[l][k] == ol[z].IEN_algo2[j][i]) { //If this node has already been assigned, use the same numbering
								ol[z].IEN_2D[j][i] = ol[z].IEN_2D[l][k];
								flag = 0; //turn off the flag to assgin new number
							}
							else {
								//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
							}
						}
					}
					if (flag == 1) {
						ct += 1;
						dummy.push_back(ol[z].IEN_algo2[j][i]); //associate the local 2D node with the global node numbering 
						ol[z].IEN_2D[j][i] = ct;
					}
				}
			}
			ol[z].GIDNct = dummy.size();
			ol[z].GIDN = new int[ol[z].GIDNct];
			for (i = 0; i < ol[z].GIDNct; i++) {
				ol[z].GIDN[i] = dummy[i];
			}
			//Obtain GIDF
			ol[z].GIDF = new int[ol[z].FSNEL];
			for (i = 0; i < ol[z].FSNEL; i++) {
				flag = 0; 
				for (j = 0; j < t.NEL; j++) {
					ct = 0; 
					for (k = 0; k < NINT*NINT; k++) {
						if (t.IEN[ol[z].FP[k] - 1][j] == ol[z].IEN_gb[ol[z].FP_2D[k] - 1][i]) {
							ct += 1; 
						}
					}
					if (ct == NINT*NINT) { //find the corresponding 3D element
						ol[z].GIDF[i] = j + 1; 
						flag = 1; 
					}
				}
				if (flag == 0) {
					std::cout << "Not corresponding 3D element is found" << std::endl;
					system("PAUSE "); 
				}
			}

			ol[z].norm = new double*[ol[z].FSNEL]; //store the normal direction of linear elements on the wetted surface
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].norm[i] = new double[3];
			}
			//Obtain the normal direction unit vector of the newly separated elements (numbering from 0 to elenum)
			//The normal direction calculation might be wrong. We need to validate it using a FSP mesh generated by BOLT.
			double ax, ay, az; double bx, by, bz;
			double n1, n2, n3; double absn;
			for (i = 0; i < ol[z].FSNEL; i++) {
				ax = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][0] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[0] - 1][i] - 1][0];
				ay = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][1] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[0] - 1][i] - 1][1];
				az = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][2] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[0] - 1][i] - 1][2];
				bx = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[2] - 1][i] - 1][0] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][0];
				by = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[2] - 1][i] - 1][1] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][1];
				bz = t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[2] - 1][i] - 1][2] - t.GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[1] - 1][i] - 1][2];
				n1 = ay*bz - az*by; n2 = az*bx - ax*bz; n3 = ax*by - ay*bx;
				absn = sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
				ol[z].norm[i][0] = n1 / absn; ol[z].norm[i][1] = n2 / absn; ol[z].norm[i][2] = n3 / absn;
			}
			//check if the separated element has the same normal direction as the original mesh (N=1)
			//std::cout << " " << std::endl; 
		}
		//================================Write the model file for MpCCI here================================//
		std::ofstream myfile;
		myfile.open("model.txt");
		for (z = 0; z < wt_pys_size; z++) {
			std::string wetsurface_name; 
			wetsurface_name = "EF wetsurface" + std::to_string(z + 1) + " 3 1"; 
			myfile << wetsurface_name << std::endl;
			myfile << "NODES " << ol[z].GIDNct << std::endl;
			for (i = 0; i < ol[z].GIDNct; i++) {
				myfile << i << " " << t.GCOORD[ol[z].GIDN[i] - 1][0] << " " << t.GCOORD[ol[z].GIDN[i] - 1][1] << " " << t.GCOORD[ol[z].GIDN[i] - 1][2] << " " << std::endl;
			}
			myfile << "ELEMENTS " << ol[z].FSNEL << std::endl;
			//Output connectivity matrix
			for (i = 0; i < ol[z].FSNEL; i++) {
				myfile << i;
				for (j = 0; j < 4; j++) {
					myfile << " " << ol[z].IEN_2D[j][i] - 1; //node numbering starts from 0 in model file
				}
				myfile << std::endl;
			}
		}
	}

	//==============Extract the connectivity matrix on NRB surface (glue all physical groups corresponding to the NRB together)==============//
	//All physical groups except for 3
	for (z = 0; z < nrb_pys_size; z++) {
		nr[z].NEL_nrb = phygrp_start[nrb_pys_num[z] + 1] - phygrp_start[nrb_pys_num[z]];
		nr[z].IEN_gb = new int*[NINT*NINT]; //NRBELE_ARR
		for (i = 0; i < NINT*NINT; i++) {
			nr[z].IEN_gb[i] = new int[nr[z].NEL_nrb];
		}

		ct = 0;
		for (e = 0; e < nr[z].NEL_nrb; e++) {
			for (j = 0; j < NINT*NINT; j++) {
				nr[z].IEN_gb[j][ct] = t.BCIEN[nrb_pys_num[z]][e][j];
			}
			ct += 1;
		}
		if (ct != nr[z].NEL_nrb) {
			std::cout << "The NRB element storage is wrong" << std::endl;
			system("PAUSE ");
		}

		//Obtain NRBA to store all the nodes on NRB surfaces by looping through all the NRB elements and filter out the points
		ct = 0; //count the node number assigned
		std::vector<int>dummy2;
		int flag; 
		for (i = 0; i < nr[z].NEL_nrb; i++) { //loop through each element
			for (j = 0; j < NINT*NINT; j++) { //the nodes in current element
				flag = 1; //Initiate the flag to 1 
				for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
					for (l = 0; l < NINT*NINT; l++) {
						if (nr[z].IEN_gb[l][k] == nr[z].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
							flag = 0; //turn off the flag to assgin new number
						}
						else {
							//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
						}
					}
				}
				if (flag == 1) {
					ct += 1;
					dummy2.push_back(nr[z].IEN_gb[j][i]); //associate the local 2D node with the global node numbering 
				}
			}
		}
		nr[z].NRBNODE = dummy2.size();
		nr[z].NRBA = new int[nr[z].NRBNODE];
		for (i = 0; i < nr[z].NRBNODE; i++) {
			nr[z].NRBA[i] = dummy2[i];
		}

		//obtain NRBELE_ARR
		nr[z].NRBELE_ARR = new int[nr[z].NEL_nrb];
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			flag = 0;
			for (j = 0; j < t.NEL; j++) {
				ct = 0;
				for (k = 0; k < NINT*NINT; k++) {
					if (t.IEN[nr[z].DP[k] - 1][j] == nr[z].IEN_gb[nr[z].DP_2D[k] - 1][i]) {
						ct += 1;
					}
				}
				if (ct == NINT*NINT) { //find the corresponding 3D element
					nr[z].NRBELE_ARR[i] = j + 1;
					flag = 1;
				}
			}
			if (flag == 0) {
				std::cout << "Not corresponding 3D element is found" << std::endl;
				system("PAUSE ");
			}
		}
		//obtain the normal vector corresponding to every NRB element (high-order)
		nr[z].norm = new double*[nr[z].NEL_nrb]; //store the normal direction of linear elements on the wetted surface
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			nr[z].norm[i] = new double[3];
		}
		//Obtain the normal direction unit vector of the newly separated elements (numbering from 0 to elenum)
		//The normal direction calculation might be wrong. We need to validate it using a FSP mesh generated by BOLT.
		double ax, ay, az; double bx, by, bz;
		double n1, n2, n3; double absn;
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			ax = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][0] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[0] - 1][i] - 1][0];
			ay = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][1] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[0] - 1][i] - 1][1];
			az = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][2] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[0] - 1][i] - 1][2];
			bx = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[2] - 1][i] - 1][0] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][0];
			by = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[2] - 1][i] - 1][1] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][1];
			bz = t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[2] - 1][i] - 1][2] - t.GCOORD[nr[z].IEN_gb[nr[z].LNA_norm[1] - 1][i] - 1][2];
			n1 = ay*bz - az*by; n2 = az*bx - ax*bz; n3 = ax*by - ay*bx;
			absn = sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
			nr[z].norm[i][0] = n1 / absn; nr[z].norm[i][1] = n2 / absn; nr[z].norm[i][2] = n3 / absn;
		}
		//std::cout << " " << std::endl; 
	}

	std::cout << " " << std::endl;
	return t;
}
