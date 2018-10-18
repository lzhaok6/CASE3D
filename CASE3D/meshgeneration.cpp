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
This code takes Gmsh (.msh) format mesh file. 
*/

struct meshgenerationstruct meshgeneration() {
	meshgenerationstruct t;
	LOCAL_NODEstruct c;
	extern OWETSURF ol[owsfnumber];
	extern NRBSURF nr[nrbsurfnumber];
	c = LOCAL_NODE(N);

	//We are assuming there is only one wetted surface. Thus, ol[0] is used throughout the loop below. 
	int i, j, k, l, m, n, o, e, z;
	std::string filename;
	std::cout << "reading the mesh file: " << std::endl;
	std::cout << "Have you configured the mesh file name correctly? If yes, hit Enter to proceed" << std::endl;

	
	/*
	std::ifstream infile("C:/Users/lzhaok6/OneDrive/CASE_MESH/FSP_N=1.msh");
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
	int physicalgroups = 1; //starting from the first physical group
	int endfile = 0;
	std::vector<int> phygrp_start; //the starting line number of the physical group
	int elestart = 0; //the flag denote the start of element connectivity definition

	std::clock_t start;
	start = std::clock();
	std::string csvElement;
	std::vector<std::string> csvColumn;
	std::string csvLine;

	while (getline(infile, csvLine))
	{
		ct = ct + 1; //the current line number (starting from 0)
		std::istringstream csvStream(csvLine); //csvStream is a stream (turn the string csvLine to stream csvStream)
		csvColumn.clear();
		while (getline(csvStream, csvElement, ' ')) //Get line from stream (csvStream) into string (csvElement)
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
		//Get the information of physical groups (Here, we assume the physical group number is sequential from 1, 2, 3...)
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
	double duration = (std::clock() - start);

	int elenode3D = 0;
	int elenode2D = 0; 
	if (element_type == 0) { //hex element
		elenode3D = NINT*NINT*NINT;
		elenode2D = NINT*NINT;
	}
	if (element_type == 1) { //tet element
		if (N == 1) {
			elenode3D = 4;
			elenode2D = 3;
		}
		else {
			std::cout << "High-order tet element is not supported yet" << std::endl;
			system("PAUSE ");
		}
	}

	//The physical groups except for the last one are the surface mesh (elenode2D). 
	//The last physical group is the volumn mesh (elenode3D) connectivity matrix
	t.NNODE = stoi(output[nodestart - 1][0]);
	t.NEL = endfile - phygrp_start[physicalgroups - 1];
	//Define connectivity matrix in the volumn mesh
	t.IEN = new int*[elenode3D];
	for (i = 0; i < elenode3D; i++) {
		t.IEN[i] = new int[t.NEL];
	}
	int nelstart = 0;
	int tagnumber = stoi(output[phygrp_start[0]][2]); 
	nelstart = 2 + tagnumber + 1;
	//tagnumber is the number of tags (placed third in the line). After those tags the element connectivity starts
	//Refer to Evernote: Gmsh mesh file format for more information. 
	ct = 0;
	for (i = phygrp_start[physicalgroups - 1]; i < endfile; i++) {
		for (j = nelstart; j < nelstart + elenode3D; j++) {
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
			for (k = 0; k < elenode2D; k++) {
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
	*/
	int elenode3D = 0;
	int elenode2D = 0;
	if (element_type == 0) { //hex element
		elenode3D = NINT*NINT*NINT;
		elenode2D = NINT*NINT;
	}
	if (element_type == 1) { //tet element
		if (N == 1) {
			elenode3D = 4;
			elenode2D = 3;
		}
		else {
			std::cout << "High-order tet element is not supported yet" << std::endl;
			system("PAUSE ");
		}
	}
	std::cout << "reading the mesh file: " << std::endl;
	std::cout << "Have you configured the mesh file name correctly? If yes, hit Enter to proceed" << std::endl; 
	int ct = -1;

	FILE *fp = fopen("C:/Users/lzhaok6/OneDrive/CASE_MESH/FSP_N=2_mismatch.msh", "r");
	if (!fp) {
		printf("Cannot open the mesh file");
		system("PAUSE ");
	}

	if (N > 9) {
		std::cout << "the buf below may not be able to handle" << std::endl;
		system("PAUSE ");
	}

	char buf[1000];
	do {
		fgets(buf, 1000, fp);
	} while (!strstr(buf, "$Nodes"));
	fscanf(fp, "%d", &t.NNODE);
	t.GCOORD = new double*[t.NNODE];
	for (i = 0; i < t.NNODE; i++) {
		t.GCOORD[i] = new double[3];
	}
	double* VX = new double[t.NNODE]; //x coordinate
	double* VY = new double[t.NNODE]; //y coordinate
	double* VZ = new double[t.NNODE]; //z coordinate
	for (int i = 0; i < t.NNODE; i++) {
		fscanf(fp, "%*d %lf %lf %lf",
			VX + i, VY + i, VZ + i);
		if (i % 10000 == 0) {
			printf("%d %g %g %g\n", i, VX[i], VY[i], VZ[i]);
		}
	}
	//search until the $Elements is found (the start of element connecvitity definition)
	do {
		fgets(buf, 1000, fp);
	} while (!strstr(buf, "$Elements"));
	int NEL_tol; //total number of element including boundary 2D elements
	fscanf(fp, "%d", &NEL_tol);
	std::vector<std::vector<int>>iens;
	fgets(buf, 1000, fp);
	int* EToV; //element connectivity matrix for all physical groups (make it an 1D array to accelerate the memory access)
	EToV = new int[NEL_tol*elenode3D]; //more than enough
	for (i = 0; i < NEL_tol*elenode3D; i++) {
		EToV[i] = 0;
	}
	int* ele_type; //store the number of vertices for the corresponding element type
	ele_type = new int[200]; //more than enough
	for (i = 0; i < 200; i++) {
		ele_type[i] = 0;
	}
	//define element types (number start from 1)
	ele_type[2] = 3; //3-node triangle
	ele_type[3] = 2 * 2; //4-node quadrangle
	ele_type[4] = 2 * 2; //4-node tetrahedron
	ele_type[5] = 2 * 2 * 2; //8-node hexahedron
	ele_type[10] = 3 * 3; //2nd order 4-node quadrangle
	ele_type[12] = 3 * 3 * 3; //2nd order 8-node hexahedron
	ele_type[36] = 4 * 4; //3rd order 4-node quadrangle
	ele_type[92] = 4 * 4 * 4; //3rd order 8-node hexahedron

	int psy_curt; //store the current physical group number
	int ele_typ[20]; //maximum 20 physical groups
	int physicalgroups = 1; //number of physical groups
	int numele = 0;
	std::vector<int>numele_py; //total number of element in each physical group
	ct = 0;

	for (int e = 0; e < NEL_tol; e++) { //loop through all elements first
		if (e % 10000 == 0) {
			std::cout << "element: " << e << std::endl; //output which line is being read
		}
		int etype, Ntags;
		fscanf(fp, "%*d %d %d", &etype, &Ntags);
		//count the number of elements in the current physical group
		fscanf(fp, "%d", &psy_curt);
		if (psy_curt == physicalgroups) {
			numele += 1; //cumulatively store the number of element in the current physical group
						 //store the element type for each physical group
			ele_typ[psy_curt - 1] = etype;
		}
		else { //the current physical group is done, move on to the next physical group
			numele_py.push_back(numele);
			numele = 1;
			physicalgroups += 1; //add one more physical group
		}
		for (int t = 0; t < Ntags - 1; t++) {
			fscanf(fp, "%*d");
		}
		int Nverts = ele_type[etype];  //number of vertices of the element type in the current line
		if (Nverts == 0) {
			std::cout << "This element type is yet to be defined" << std::endl;
			system("PAUSE ");
		}
		for (int v = 0; v < Nverts; v++) { //store the element connectivity matrix for all physical groups
			fscanf(fp, "%d", EToV + ct + v);
		}
		ct += Nverts;
	}
	numele_py.push_back(numele); //number of element in each physical group

	//Populate the EToV matrix to BCIEN and IEN
	//Firstly, store the boundary element connectivity (BCIEN)
	int st = 0; //the starting point of the connectivity definition in EToV
	int nvert_py; //number of vertices in element for each physical group
	for (i = 0; i < physicalgroups - 1; i++) {
		nvert_py = ele_type[ele_typ[i]];
		std::vector<std::vector<int>>iens; //store the connectivity in the current physical group (2D vector)
		for (j = 0; j < numele_py[i]; j++) { //loop through the element in the current physical group
			std::vector<int>local;
			for (k = 0; k < nvert_py; k++) {
				local.push_back(EToV[st + j*nvert_py + k]);
			}
			iens.push_back(local);
		}
		t.BCIEN.push_back(iens);
		st += numele_py[i] * nvert_py;
	}
	//Secondly, store the global element connectivity (IEN)
	t.NEL = numele_py[physicalgroups - 1]; 
	t.IEN = new int*[elenode3D]; //element connectivity matrix
	for (i = 0; i < elenode3D; i++) { //Since some elements are boundary elements, they may only consume elenode2D in the first dimension. Thus, some space is wasted. 
		t.IEN[i] = new int[t.NEL];
	}
	for (j = 0; j < t.NEL; j++) { //loop through the element in the current physical group
		nvert_py = ele_type[ele_typ[physicalgroups - 1]];
		for (k = 0; k < nvert_py; k++) {
			t.IEN[k][j] = EToV[st + j*nvert_py + k];
		}
	}
	//populate GCOORD
	for (i = 0; i < t.NNODE; i++) {
		t.GCOORD[i][0] = VX[i];
		t.GCOORD[i][1] = VY[i];
		t.GCOORD[i][2] = VZ[i];
	}

	int* IEN_1D;
	if (element_type == 1) { //tet element
		IEN_1D = new int[t.NEL*elenode3D];
		for (i = 0; i < t.NEL; i++) {
			for (j = 0; j < elenode3D; j++) {
				IEN_1D[elenode3D*i + j] = t.IEN[j][i];
			}
		}
	}
	//=====================Changing the coordinate of GCOORD to SEM (shift the internal nodes)=====================//
	if (FEM == 0) {
		LOCAL_NODEstruct c_ln;
		c_ln = LOCAL_NODE(1);
		LOCAL_SHAPEstruct ls_ln; //ln means linear
		ls_ln = LOCAL_SHAPE(c_ln.LNA, 1, N, FEM); //Get the 3D linear shape function value on Nth order SEM nodes
		//t.SHL[3][LNA[i][j][k] - 1][l*(NQUAD + 1)*(NQUAD + 1) + m*(NQUAD + 1) + o]
		int LNA_ln[2][2][2];
		LNA_ln[0][0][0] = c.LNA[0][0][0]; LNA_ln[1][0][0] = c.LNA[N][0][0];
		LNA_ln[1][1][0] = c.LNA[N][N][0]; LNA_ln[0][1][0] = c.LNA[0][N][0];
		LNA_ln[0][0][1] = c.LNA[0][0][N]; LNA_ln[1][0][1] = c.LNA[N][0][N];
		LNA_ln[1][1][1] = c.LNA[N][N][N]; LNA_ln[0][1][1] = c.LNA[0][N][N];
		for (i = 0; i < t.NEL; i++) { //loop through all the elements
			ct = 0;
			for (m = 0; m < NINT; m++) { //m, n, o stands for internal nodes in one element (all but except for corner nodes)
				for (n = 0; n < NINT; n++) {
					for (o = 0; o < NINT; o++) {
						//jump over the corner points
						if (!((m == 0 && n == 0 && o == 0) || (m == N && n == 0 && o == 0) || (m == N && n == N && o == 0)
							|| (m == 0 && n == N && o == 0) || (m == 0 && n == 0 && o == N) || (m == N && n == 0 && o == N)
							|| (m == N && n == N && o == N) || (m == 0 && n == N && o == N))) {
							t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][0] = 0.0;
							t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][1] = 0.0;
							t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][2] = 0.0;
							for (j = 0; j < 2; j++) { //j, k, l stands for corner nodes
								for (k = 0; k < 2; k++) {
									for (l = 0; l < 2; l++) {
										t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][0] += t.GCOORD[t.IEN[LNA_ln[j][k][l] - 1][i] - 1][0] * ls_ln.SHL[3][c_ln.LNA[j][k][l] - 1][m*elenode2D + n*NINT + o];
										t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][1] += t.GCOORD[t.IEN[LNA_ln[j][k][l] - 1][i] - 1][1] * ls_ln.SHL[3][c_ln.LNA[j][k][l] - 1][m*elenode2D + n*NINT + o];
										t.GCOORD[t.IEN[c.LNA[m][n][o] - 1][i] - 1][2] += t.GCOORD[t.IEN[LNA_ln[j][k][l] - 1][i] - 1][2] * ls_ln.SHL[3][c_ln.LNA[j][k][l] - 1][m*elenode2D + n*NINT + o];
									}
								}
							}
							ct += 1;
						}
					}
				}
			}
			//std::cout << " " << std::endl;
		}
	}


	//Need to be modified, localnode needs to corresponds each physical group! we cannot make the assumption anymore
	//int localnode[elenode2D]; 
	int **localnode; 
	localnode = new int*[physicalgroups - 1];
	for (i = 0; i < physicalgroups - 1; i++) {
		localnode[i] = new int[elenode2D];
	}
	//The following loop only attempts to extract the pattern of localnode[] from one corresponding global element since the assumption here is that all the elements in all physical groups
	//follow the same pattern. However, that is not rigorous enough since the corresponding local node for all 2D element in physical group do not necessary be the same. 
	for (i = 0; i < physicalgroups - 1; i++) {
		for (j = 0; j < numele_py[i]; j++) { //phygrp_start[i + 1] - phygrp_start[i] is the element number in the ith physical group
			for (k = 0; k < t.NEL; k++) {
				ct = 0;
				for (l = 0; l < elenode2D; l++) {
					for (m = 0; m < elenode3D; m++) {
						if (t.BCIEN[i][j][l] == t.IEN[m][k]) {
							localnode[i][l] = m + 1; //store the corresponding local node in BCIEN in global element connectivity matrix
							ct += 1;
						}
					}
				}
				if (ct == elenode2D) { //found the corresponding 3D element! 
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
	//int wt_pys_num[4] = { 0,1,2,3 };  //the physical group number that corresponds to the wet surface (physical group 3)
	//int nrb_pys_num[4] = { 4,5,6,7 };
	//int wt_pys_num[1] = { 0 };
	//int nrb_pys_num[1] = { 1 };
	
	//If the wet surface is the left face of the local element (i.e., i=0). 
	std::cout << "Have you configured the wetted surface and NRB surface physical group number for this mesh? If so, hit Enter to proceed." << std::endl;
	//system("PAUSE ");
	
	int wet = 0; 
	int nrb = 0; 
	int wet_ct = 0;
	int nrb_ct = 0; 
	int wt_pys_size = sizeof(wt_pys_num) / sizeof(*wt_pys_num); //the number of wetted surfaces 
	int nrb_pys_size = sizeof(nrb_pys_num) / sizeof(*nrb_pys_num); //the number of nrb surfaces

	//The four nodes (hex element) and three nodes (tet element) arranged by LNA_norm trace the normal direction out of the fluid domain.
	//For linear tet element, this should just be the original sequency of the surface 2D triangular elements. 
	for (z = 0; z < physicalgroups - 1; z++) {
		nr[z].LNA_norm = new int[elenode2D];
		ol[z].LNA_norm = new int[elenode2D];
	}

	if (element_type == 1) {
		for (z = 0; z < physicalgroups - 1; z++) {
			nr[z].LNA_norm[0] = 1;
			ol[z].LNA_norm[0] = 1;
			nr[z].LNA_norm[1] = 2;
			ol[z].LNA_norm[1] = 2;
			nr[z].LNA_norm[2] = 3;
			ol[z].LNA_norm[2] = 3;
		}
	}
	if (element_type == 0) {
		int element_number = 0;
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
			//left face 
			ct = 0;
			for (j = 0; j < NINT; j++) {
				for (k = 0; k < NINT; k++) {
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[0][j][k]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[j][k] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[0][j][k];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[j][k];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[j][k] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[0][j][k];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[j][k];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 1; //y coordinate
					nr[nrb_ct].Jacob_face[1] = 2;
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
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[N][j][k]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[j][k] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[N][j][k];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[j][k];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[j][k] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[N][j][k];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[j][k];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 1; //y coordinate
					nr[nrb_ct].Jacob_face[1] = 2;
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
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[i][j][0]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[i][j] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[i][j][0];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][j];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[i][j] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[i][j][0];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][j];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 0; //x coordinate
					nr[nrb_ct].Jacob_face[1] = 1; //y coordinate
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
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[i][j][N]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[i][j] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[i][j][N];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][j];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[i][j] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[i][j][N];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][j];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 0; //x coordinate
					nr[nrb_ct].Jacob_face[1] = 1; //y coordinate
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
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[i][0][k]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[i][k] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[i][0][k];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][k];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[i][k] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[i][0][k];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][k];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 0; //x coordinate
					nr[nrb_ct].Jacob_face[1] = 2; //z coordinate
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
					for (l = 0; l < elenode2D; l++) {
						if (localnode[z][l] == c.LNA[i][N][k]) {
							if (wet == 1) {
								ol[wet_ct].LNA_2D[i][k] = l + 1;
								ol[wet_ct].FP_temp[ct] = c.LNA[i][N][k];
								ol[wet_ct].FP_2D[ct] = ol[wet_ct].LNA_2D[i][k];
							}
							else if (nrb == 1) {
								nr[nrb_ct].LNA_2D[i][k] = l + 1;
								nr[nrb_ct].DP_temp[ct] = c.LNA[i][N][k];
								nr[nrb_ct].DP_2D[ct] = nr[nrb_ct].LNA_2D[i][k];
							}
							ct += 1;
						}
					}
				}
			}
			if (ct == elenode2D) {
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
					nr[nrb_ct].Jacob_face[0] = 0; //x coordinate
					nr[nrb_ct].Jacob_face[1] = 2; //z coordinate
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
	}

	//================================end extraction===================================//

	//=============separate the high-order element into linear elements and obtain the normal direction================//
	//int wt_pys_num = 2; //the physical group number that corresponds to the wet surface (physical group 3)
	//if (mappingalgo == 2) {
		for (z = 0; z < wt_pys_size; z++) { //loop through each wetted surface

			int elenum = numele_py[wt_pys_num[z]]; //number of high-order on the wet surface s

			if (element_type == 0) { //hex element
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
			}

			//track the element number of the wetted surface elements
			//If any node of the element is on the free surface, then this element is excluded from the wetted surface.
			std::cout << "The wetted surface element identification needs to be changed before proceeding with FSP" << std::endl;
			//system("PAUSE ");
			int flag;
			std::vector<int>ele_num;
			ct = 0;
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
					ct += 1; 
				}
			}
			//ele_num is the elements for high-order element, we need to extract the corresponding linear elements from that. 

			//Store the total number of element on wetted surface in FSNEL
			//ol[z].FSNEL = ele_num.size();
			ol[z].FSNEL = ct;
			if (ct != ele_num.size()) {
				std::cout << "There is a problem with ol[z].FSNEL" << std::endl;
				system("PAUSE "); 
			}

			ol[z].IEN_gb = new int*[elenode2D]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			for (i = 0; i < elenode2D; i++) {
				ol[z].IEN_gb[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < elenode2D; j++) {
					ol[z].IEN_gb[j][i] = t.BCIEN[wt_pys_num[z]][ele_num[i]][j];
				}
			}

			ol[z].FP = new int*[ol[z].FSNEL];
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].FP[i] = new int[elenode2D];
			}
			//Define FP for hexahedral element
			if (element_type == 0) {
				for (i = 0; i < ol[z].FSNEL; i++) {
					for (j = 0; j < elenode2D; j++) {
						ol[z].FP[i][j] = ol[z].FP_temp[j];
					}
				}
			}

			//Define FP for tetrahedral element
			//Extract the local node pattern for tetrahedral element
			if (element_type == 1) {
				int node[3];
				ol[z].GIDF = new int[ol[z].FSNEL];
				//#pragma omp parallel for num_threads(6)
				for (l = 0; l < ol[z].FSNEL; l++) {
					if (l % 100 == 0) {
						std::cout << l << std::endl; //output which line is being read
					}
					node[0] = t.BCIEN[wt_pys_num[z]][l][0]; node[1] = t.BCIEN[wt_pys_num[z]][l][1]; node[2] = t.BCIEN[wt_pys_num[z]][l][2];
					for (i = 0; i < t.NEL; i++) { //loop through all the element
						ct = 0;
						for (k = 0; k < elenode2D; k++) { //loop through 2D local nodes
							for (j = 0; j < elenode3D; j++) { //loop through local 3D nodes
								//if (node[k] == t.IEN[j][i]) {
								if (node[k] == IEN_1D[i*elenode3D + j]) {
									ol[z].FP[l][k] = j + 1;
									ol[z].FP_2D[ct] = ct + 1;
									ct = ct + 1;
									if (ct == 3) { //found all local nodes in one global element (ready to move on to the next element)
										ol[z].GIDF[l] = i + 1;
										i = t.NEL; //Forcedly exit the i=0;i<t.NEL;i++ loop above 
									}
								}
							}
							if (ct == 0) { //if the first 2D local node does not have a corresponding node in the current gloabl element being searched, the search would move forward to the next gloabl element
								k = 3;
							}
						}
					}
				}
			}

			if (element_type == 0) { //hex element
				ol[z].IEN_algo2 = new int*[4]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
				for (i = 0; i < 4; i++) {
					ol[z].IEN_algo2[i] = new int[ol[z].FSNEL];
				}
				for (i = 0; i < ol[z].FSNEL; i++) {
					for (j = 0; j < 4; j++) {
						ol[z].IEN_algo2[j][i] = ol[z].IEN_py[j][ele_num[i]];
					}
				}
			}

			if (element_type == 0) {
				//Derive the IEN_2D to write the MpCCI model file (basically renumbering the node in IEN_algo2 sequentially to be recognized by MpCCI)
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
							ol[z].IEN_2D[j][i] = ct; //assign a new number to MpCCI element connectivity
						}

					}
				}
				//ol[z].GIDNct_MpCCI = dummy.size();
				ol[z].GIDNct_MpCCI = ct; 
				if (ct != dummy.size()) {
					std::cout << "There is a problem with ol[z].GIDNct_MpCCI" << std::endl;
					system("PAUSE "); 
				}
				ol[z].GIDN_MpCCI = new int[ol[z].GIDNct_MpCCI];
				for (i = 0; i < ol[z].GIDNct_MpCCI; i++) {
					ol[z].GIDN_MpCCI[i] = dummy[i];
				}
			}
			if (element_type == 1) {
				//Derive the IEN_2D to write the MpCCI model file (basically renumbering the node in IEN_algo2 sequentially to be recognized by MpCCI)
				ol[z].IEN_2D = new int*[3]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
				for (i = 0; i < 3; i++) {
					ol[z].IEN_2D[i] = new int[ol[z].FSNEL];
				}
				ct = 0; //count the node number assigned
				std::vector<int>dummy;
				for (i = 0; i < ol[z].FSNEL; i++) { //loop through each element
					for (j = 0; j < 3; j++) { //the nodes in current element
						flag = 1; //Initiate the flag to 1 
						for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
							for (l = 0; l < 3; l++) {
								if (ol[z].IEN_gb[l][k] == ol[z].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
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
							dummy.push_back(ol[z].IEN_gb[j][i]); //associate the local 2D node with the global node numbering 
							ol[z].IEN_2D[j][i] = ct; //assign a new number to MpCCI element connectivity
						}
					}
				}
				//ol[z].GIDNct_MpCCI = dummy.size();
				ol[z].GIDNct_MpCCI = ct; 
				if (ct != dummy.size()) {
					std::cout << "There is a problem with the ol[z].GIDNct_MpCCI" << std::endl;
					system("PAUSE "); 
				}
				ol[z].GIDN_MpCCI = new int[ol[z].GIDNct_MpCCI];
				for (i = 0; i < ol[z].GIDNct_MpCCI; i++) {
					ol[z].GIDN_MpCCI[i] = dummy[i];
				}
			}

			//Obtain GIDF (the global element number of each wetted surface)
			if (element_type == 0) {
				ol[z].GIDF = new int[ol[z].FSNEL];
				for (i = 0; i < ol[z].FSNEL; i++) {
					flag = 0;
					for (j = 0; j < t.NEL; j++) {
						ct = 0;
						for (k = 0; k < elenode2D; k++) {
							if (t.IEN[ol[z].FP[i][k] - 1][j] == ol[z].IEN_gb[ol[z].FP_2D[k] - 1][i]) {
								ct += 1;
							}
						}
						if (ct == elenode2D) { //find the corresponding 3D element
							ol[z].GIDF[i] = j + 1;
							flag = 1;
						}
					}
					if (flag == 0) {
						std::cout << "No corresponding 3D element is found" << std::endl;
						system("PAUSE ");
					}
				}
			}

			//Get the high-order global nodes on wetted surface (GIDN)
			std::vector<int>dummy3;
			ct = 0; 
			for (i = 0; i < ol[z].FSNEL; i++) { //loop through each element
				for (j = 0; j < elenode2D; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < elenode2D; l++) {
							if (ol[z].IEN_gb[l][k] == ol[z].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
								flag = 0; //turn off the flag to assgin new number
							}
							else {
								//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
							}
						}
					}
					if (flag == 1) {
						dummy3.push_back(ol[z].IEN_gb[j][i]); //associate the local 2D node with the global node numbering 
						ct += 1;
					}
				}
			}
			//ol[z].GIDNct = dummy3.size();
			ol[z].GIDNct = ct;
			if (ct != dummy3.size()) {
				std::cout << "There is a problem with GIDNct" << std::endl;
				system("PAUSE "); 
			}

			ol[z].GIDN = new int[ol[z].GIDNct];
			for (i = 0; i < ol[z].GIDNct; i++) {
				ol[z].GIDN[i] = dummy3[i];
			}

			ol[z].norm = new double*[ol[z].FSNEL]; //store the normal direction of linear elements on the wetted surface
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].norm[i] = new double[3];
			}
			//Obtain the normal direction unit vector of the newly separated elements (numbering from 0 to elenum)
			//The normal direction calculation might be wrong. We need to validate it using a FSP mesh generated by BOLT.
			double ax, ay, az; double bx, by, bz;
			double n1, n2, n3; double absn;
			ol[z].dimension = new double[ol[z].FSNEL]; //the dimension of the triangle formed by two vectors
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
				ol[z].dimension[i] = 0.5*pow((n1*n1 + n2*n2 + n3*n3), 0.5); //https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
			}
			//check if the separated element has the same normal direction as the original mesh (N=1)
			//std::cout << " " << std::endl; 
		}

		//================================Write the model file for MpCCI here================================//
		if (mappingalgo == 2) {
			std::ofstream myfile;
			myfile.open("model.txt");
			for (z = 0; z < wt_pys_size; z++) {
				std::string wetsurface_name;
				//wetsurface_name = "EF wetsurface" + std::to_string(z + 1) + " 3 1";
				wetsurface_name = "EF wetsurface" + std::to_string(z + 1) + " 3 1 " + std::to_string(element_type);
				myfile << wetsurface_name << std::endl;
				myfile << "NODES " << ol[z].GIDNct_MpCCI << std::endl;
				for (i = 0; i < ol[z].GIDNct_MpCCI; i++) {
					myfile << i << " " << t.GCOORD[ol[z].GIDN_MpCCI[i] - 1][0] << " " << t.GCOORD[ol[z].GIDN_MpCCI[i] - 1][1] << " " << t.GCOORD[ol[z].GIDN_MpCCI[i] - 1][2] << " " << std::endl;
				}
				myfile << "ELEMENTS " << ol[z].FSNEL << std::endl;
				//Output connectivity matrix
				for (i = 0; i < ol[z].FSNEL; i++) {
					myfile << i;
					if (element_type == 0) {
						for (j = 0; j < 4; j++) {
							myfile << " " << ol[z].IEN_2D[j][i] - 1; //node numbering starts from 0 in model file
						}
					}
					if (element_type == 1) {
						for (j = 0; j < 3; j++) {
							myfile << " " << ol[z].IEN_2D[j][i] - 1; //node numbering starts from 0 in model file
						}
					}
					myfile << std::endl;
				}
			}
		}
	//==============Extract the connectivity matrix on NRB surface (glue all physical groups corresponding to the NRB together)==============//
	//All physical groups except for 3
		for (z = 0; z < nrb_pys_size; z++) {
			nr[z].NEL_nrb = numele_py[nrb_pys_num[z]];
			nr[z].IEN_gb = new int*[elenode2D]; //NRBELE_ARR
			for (i = 0; i < elenode2D; i++) {
				nr[z].IEN_gb[i] = new int[nr[z].NEL_nrb];
			}

			//Extract the local node pattern for hexahedral and tetrahedral element
			nr[z].DP = new int*[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].DP[i] = new int[elenode2D];
			}
			if (element_type == 0) {
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					for (j = 0; j < elenode2D; j++) {
						nr[z].DP[i][j] = nr[z].DP_temp[j];
					}
				}
			}
			if (element_type == 1) {
				int node[3];
				nr[z].NRBELE_ARR = new int[nr[z].NEL_nrb];
				for (l = 0; l < nr[z].NEL_nrb; l++) {
					if (l % 100 == 0) {
						std::cout << l << std::endl; //output which line is being read
					}
					node[0] = t.BCIEN[nrb_pys_num[z]][l][0]; node[1] = t.BCIEN[nrb_pys_num[z]][l][1]; node[2] = t.BCIEN[nrb_pys_num[z]][l][2];
					for (i = 0; i < t.NEL; i++) { //loop through all the element
						ct = 0;
						for (k = 0; k < elenode2D; k++) { //loop through 2D local nodes
							for (j = 0; j < elenode3D; j++) { //loop through local nodes
								if (node[k] == IEN_1D[i*elenode3D + j]) { //found a corresponding node
									nr[z].DP[l][k] = j + 1;
									nr[z].DP_2D[ct] = ct + 1;
									ct = ct + 1;
									if (ct == 3) { //found all local nodes in one global element (ready to move on to the next element)
										nr[z].NRBELE_ARR[l] = i + 1;
										i = t.NEL; //Forcedly exit the i=0;i<t.NEL;i++ loop above 
									}
								}
							}
							if (ct == 0) {
								k = 3;
							}
						}
					}
				}
			}

			ct = 0;
			for (e = 0; e < nr[z].NEL_nrb; e++) {
				for (j = 0; j < elenode2D; j++) {
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
				for (j = 0; j < elenode2D; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < elenode2D; l++) {
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
			//nr[z].NRBNODE = dummy2.size();
			nr[z].NRBNODE = ct;
			if (ct != dummy2.size()) {
				std::cout << "There is a problem with nr[z].NRBNODE" << std::endl;
				system("PAUSE ");
			}

			nr[z].NRBA = new int[nr[z].NRBNODE];
			for (i = 0; i < nr[z].NRBNODE; i++) {
				nr[z].NRBA[i] = dummy2[i];
			}

			if (element_type == 0) {
				nr[z].NRBELE_ARR = new int[nr[z].NEL_nrb];
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					flag = 0;
					for (j = 0; j < t.NEL; j++) {
						ct = 0;
						for (k = 0; k < elenode2D; k++) {
							if (t.IEN[nr[z].DP[i][k] - 1][j] == nr[z].IEN_gb[nr[z].DP_2D[k] - 1][i]) {
								ct += 1;
							}
						}
						if (ct == elenode2D) { //find the corresponding 3D element
							nr[z].NRBELE_ARR[i] = j + 1;
							flag = 1;
						}
					}
					if (flag == 0) {
						std::cout << "Not corresponding 3D element is found" << std::endl;
						system("PAUSE ");
					}
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
			nr[z].dimension = new double[nr[z].NEL_nrb];
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
				nr[z].dimension[i] = 0.5*pow((n1*n1 + n2*n2 + n3*n3), 0.5); //https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
			}
			//std::cout << " " << std::endl; 
		}

	//For some reason, the normal side is flipped after using NekMesh to process the Gmsh generated all-tet mesh. 
	//If we don't use NekMesh, the surface element normal direction of the mesh directly from Gmsh is not correct (mixed).  
	
	/*
	if (element_type == 1) {
		for (z = 0; z < nrb_pys_size; z++) {
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				for (n = 0; n < 3; n++) {
					nr[z].norm[i][n] = -nr[z].norm[i][n];
				}
			}
		}
		for (z = 0; z < wt_pys_size; z++) {
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (n = 0; n < 3; n++) {
					ol[z].norm[i][n] = -ol[z].norm[i][n];
				}
			}
		}
	}
	*/

	/*
	//Check if the total dimension in Bleich-Sandler case is 1.0 for wet surface and NRB
	double hd = 0.0; 
	for (i = 0; i < nr[0].NEL_nrb; i++) {
		hd += nr[0].dimension[i];
	}
	hd = 0.0; 
	for (i = 0; i < ol[0].FSNEL; i++) {
		hd += ol[0].dimension[i];
	}
	*/
	if (Bleich == 1) {
		if (ol[0].norm[0][1] != 1 || nr[0].norm[0][1] != -1) {
			std::cout << "The normal direction of the wetted surface and NRB is wrong" << std::endl;
			system("PAUSE ");
		}
	}

	std::cout << " " << std::endl;
	return t;
}
