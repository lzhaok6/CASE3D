#include "stdafx.h"
#include "Header.h"
#include "data.h"
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
#include <algorithm>

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
	//std::string filename;
	std::cout << "reading the mesh file: " << std::endl;
	std::cout << "Have you configured the mesh file name correctly? If yes, hit Enter to proceed" << std::endl;

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
	//const char* filename = "C:/Users/lzhaok6/OneDrive/CASE_MESH/DDG_2ftftbasemesh_fs_150m_10m_24m.msh";
	const char* filename = "C:/Users/lzhaok6/OneDrive/CASE_MESH/FSP_N=1_0.25ft.msh";
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open the mesh file");
		system("PAUSE ");
	}
	std::string filenamestr(filename);
	std::string input_type; 
	if (filenamestr.find(".msh") != std::string::npos) { //Gmsh input file
		input_type = "Gmsh"; 
	}
	else if (filenamestr.find(".inp") != std::string::npos) {
		input_type = "Abaqus";
	}

	if (N > 9) {
		std::cout << "the buf below may not be able to handle" << std::endl;
		system("PAUSE ");
	}

	std::vector<int>numele_py; //total number of element in each physical group
	if (input_type == "Gmsh") {
		//Search for "$Nodes" (the start of node coordinate definition)
		char buf[1000];
		do {
			fgets(buf, 1000, fp); //read line by line
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
		//Search until the "$Elements" is found (the start of element connecvitity definition)
		do {
			fgets(buf, 1000, fp);
		} while (!strstr(buf, "$Elements"));
		int NEL_tol; //total number of element including boundary 2D elements
		fscanf(fp, "%d", &NEL_tol);
		std::vector<std::vector<int>>iens;
		fgets(buf, 1000, fp); //to change to a new line (does not seem relevant)
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
		ele_type[37] = 5 * 5; //4th order quadrature
		ele_type[38] = 6 * 6; //5th order quadrature
		ele_type[47] = 7 * 7; //6th order quadrature
		ele_type[93] = 5 * 5 * 5; //4th order hexahedron
		ele_type[94] = 6 * 6 * 6; //5th order hexahedron
		ele_type[95] = 7 * 7 * 7; //6th order hexahedron
		int psy_curt; //store the current physical group number
		int ele_typ[20]; //maximum 20 physical groups
		int physicalgroups = 1; //number of physical groups
		int numele = 0;
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
		//Finish the mesh reading
		//adjust the node location if necessary
		if (nodeadj == 1) {
			for (i = 0; i < t.NNODE; i++) {
				t.GCOORD[i][0] = VX[i];
				t.GCOORD[i][1] = VZ[i] + fs_offset;
				t.GCOORD[i][2] = -VY[i];
			}
		}
		//clean the dynamic allocated memory
		delete[] VX;
		delete[] VY;
		delete[] VZ;
		delete[] EToV; 
		delete[] ele_type; 
	}
	else if (input_type == "Abaqus") {
		char buf[1000];
		do {
			fgets(buf, 1000, fp); //read line by line
		} while (!strstr(buf, "*Node"));
		double a, b, c;
		int d, e, f, g;
		std::vector<double> VX;
		std::vector<double> VY;
		std::vector<double> VZ;
		ct = 1;
		while (fscanf(fp, "%*d, %lf, %lf, %lf", &a, &b, &c) != 0) {
			VX.push_back(a);
			VY.push_back(b);
			VZ.push_back(c);
			ct += 1;
			if (ct % 10000 == 0) {
				std::cout << "node: " << ct << std::endl; //output which line is being read
			}
		}
		t.NNODE = VX.size(); 
		//Find *Element to begin storing the element connectivity matrix
		do {
			fgets(buf, 1000, fp); //read line by line
		} while (!(buf[0] == '*'&& buf[1] == 'E' && buf[2] == 'l'));
		//start storing element connectivity matrix
		std::vector<int> EToV;
		ct = 1;
		if (element_type == 0) {
			std::cout << "the connectivity reader below does not support hex element yet" << std::endl;
			system(" PAUSE");
		}
		while (fscanf(fp, "%*d, %d, %d, %d, %d", &d, &e, &f, &g) != 0) {
			EToV.push_back(d);
			EToV.push_back(e);
			EToV.push_back(f);
			EToV.push_back(g);
			ct += 1;
			if (ct % 10000 == 0) {
				std::cout << "element: " << ct << std::endl; //output which line is being read
			}
		}
		//start searching for element sets (FSI_fluid or NRB) 
		std::vector<int> wt_py; //the physical group corresponds to wetted surface
		std::vector<int> nrb_py; //the physical group corresponds to nrb surface
		std::vector<std::vector<int>> ien_py; //the elements in each physical group
		std::vector<std::string> py_surface;
		std::string surf;
		ct = 0;
		//loop until the end of file
		while (!feof(fp)) {
			do {
				fgets(buf, 1000, fp); //read line by line and find *Elset and *Surface definition
				if (feof(fp)) {
					break;
				}
			} while (!((buf[0] == '*' && buf[2] == 'l' && buf[3] == 's') || (buf[0] == '*' && buf[2] == 'u' && buf[3] == 'r')));
			//search for the keyword FSI_fluid or NRB in the char arrary
			std::string currentline(buf); //convert the char array to string to search the key words
			std::vector<int>py_ele;
			if (currentline.find("FSI_fluid") != std::string::npos || currentline.find("NRB") != std::string::npos) { //found the definition element set for FSI_fluid
				//If this line is a surface definition, we store the surface face definition for each FSI_fluid pysical group
				if ((buf[0] == '*' && buf[2] == 'u' && buf[3] == 'r')) {
					std::cout << "Have you make the element set sequence in surface definition consistent with the element set definition?" << std::endl; 
					system(" PAUSE"); 
					//put surface definition
					int setnum;
					if (currentline.find("FSI_fluid") != std::string::npos) {
						setnum = wt_py.size();
					}
					else {
						setnum = nrb_py.size();
					}
					for (int i = 0; i < setnum; i++) {
						fgets(buf, 1000, fp);
						std::string haha(buf);
						py_surface.push_back(buf);
					}
					/*
					//fscanf does not work well with *FILE
					while (fscanf(fp, "%s", &surf) != 0) {
					py_surface.push_back(surf);
					std::cout << " " << std::endl;
					}
					*/
				}
				else {
					//store all the elements in this element set (BCIEN) 
					while (fscanf(fp, "%d,", &d) != 0) {
						py_ele.push_back(d); //store all the element number in the current physical group
					}
					if (py_ele.size() == 3 && py_ele[0] == 1) { //the elements are defined in a concise way (e.g., 1,100,1 meaning 1 to 100 with interval 1)
						for (int i = 0; i < py_ele[1] - 3; i++) {
							py_ele.push_back(py_ele[0] + 3 + i);
						}
						py_ele[1] = py_ele[0] + py_ele[2]; py_ele[2] = py_ele[1] + (py_ele[1] - py_ele[0]);
					}
					if (currentline.find("FSI_fluid") != std::string::npos) {
						wt_py.push_back(ct);
					}
					else {
						nrb_py.push_back(ct);
					}
					numele_py.push_back(py_ele.size()); 
					ct += 1;
					ien_py.push_back(py_ele);
				}
			}
		}
		//extract the element surface information and store the IEN and BCIEN matrices 
		std::vector<std::string>surface;
		for (int i = 0; i < py_surface.size(); i++) {
			ct = 0;
			std::stringstream ss(py_surface[i]);
			while (ss.good())
			{
				std::string substr;
				getline(ss, substr, ',');
				if (ct == 1) { //only store the second string
					substr.erase(std::remove(substr.begin(), substr.end(), '\n'), substr.end());
					substr.erase(std::remove(substr.begin(), substr.end(), ' '), substr.end());
					surface.push_back(substr);
				}
				ct += 1;
			}
		}
		std::cout << " " << std::endl;

		//define IEN
		t.NEL = EToV.size() / elenode3D;
		t.IEN = new int*[elenode3D];
		for (int i = 0; i < elenode3D; i++) {
			t.IEN[i] = new int[t.NEL];
		}
		for (int i = 0; i < t.NEL; i++) {
			for (int j = 0; j < elenode3D; j++) {
				t.IEN[j][i] = EToV[i*elenode3D + j];
			}
		}
		//define GIDF, NRBELE_ARR
		for (int z = 0; z < wt_py.size(); z++) {
			ol[0].FSNEL += ien_py[wt_py[z]].size();
		}
		for (int z = 0; z < nrb_py.size(); z++) {
			nr[0].NEL_nrb += ien_py[nrb_py[z]].size();
		}

		//define FP, FP_2D, DP, DP_2D
		for (int z = 0; z < owsfnumber; z++) {
			ol[z].FP = new int*[ol[z].FSNEL];
			for (int i = 0; i < ol[z].FSNEL; i++) {
				ol[z].FP[i] = new int[elenode2D];
			}
		}
		for (int z = 0; z < nrbsurfnumber; z++) {
			nr[z].DP = new int*[nr[z].NEL_nrb];
			for (int i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].DP[i] = new int[elenode2D];
			}
		}
		if (element_type == 0) {
			ct = 0;
			for (int z = 0; z < wt_py.size(); z++) {
				for (int i = 0; i < ien_py[wt_py[z]].size(); i++) {
					if (surface[wt_py[z]] == "S1") {
						ol[0].FP[ct][0] = 1; ol[0].FP[ct][1] = 2; ol[0].FP[ct][2] = 3; ol[0].FP[ct][3] = 4;
					}
					else if (surface[wt_py[z]] == "S2") {
						ol[0].FP[ct][0] = 5; ol[0].FP[ct][1] = 6; ol[0].FP[ct][2] = 7; ol[0].FP[ct][3] = 8;
					}
					else if (surface[wt_py[z]] == "S3") {
						ol[0].FP[ct][0] = 5; ol[0].FP[ct][1] = 1; ol[0].FP[ct][2] = 2; ol[0].FP[ct][3] = 6;
					}
					else if (surface[wt_py[z]] == "S4") {
						ol[0].FP[ct][0] = 6; ol[0].FP[ct][1] = 2; ol[0].FP[ct][2] = 3; ol[0].FP[ct][3] = 7;
					}
					else if (surface[wt_py[z]] == "S5") {
						ol[0].FP[ct][0] = 8; ol[0].FP[ct][1] = 7; ol[0].FP[ct][2] = 3; ol[0].FP[ct][3] = 4;
					}
					else if (surface[wt_py[z]] == "S6") {
						ol[0].FP[ct][0] = 5; ol[0].FP[ct][1] = 8; ol[0].FP[ct][2] = 4; ol[0].FP[ct][3] = 1;
					}
					ct += 1;
				}
			}
			ct = 0;
			for (int z = 0; z < nrb_py.size(); z++) {
				for (int i = 0; i < ien_py[nrb_py[z]].size(); i++) {
					if (surface[nrb_py[z]] == "S1") {
						nr[0].DP[ct][0] = 1; nr[0].DP[ct][1] = 2; nr[0].DP[ct][2] = 3; nr[0].DP[ct][3] = 4;
					}
					else if (surface[nrb_py[z]] == "S2") {
						nr[0].DP[ct][0] = 5; nr[0].DP[ct][1] = 6; nr[0].DP[ct][2] = 7; nr[0].DP[ct][3] = 8;
					}
					else if (surface[nrb_py[z]] == "S3") {
						nr[0].DP[ct][0] = 5; nr[0].DP[ct][1] = 1; nr[0].DP[ct][2] = 2; nr[0].DP[ct][3] = 6;
					}
					else if (surface[nrb_py[z]] == "S4") {
						nr[0].DP[ct][0] = 6; nr[0].DP[ct][1] = 2; nr[0].DP[ct][2] = 3; nr[0].DP[ct][3] = 7;
					}
					else if (surface[nrb_py[z]] == "S5") {
						nr[0].DP[ct][0] = 8; nr[0].DP[ct][1] = 7; nr[0].DP[ct][2] = 3; nr[0].DP[ct][3] = 4;
					}
					else if (surface[nrb_py[z]] == "S6") {
						nr[0].DP[ct][0] = 5; nr[0].DP[ct][1] = 8; nr[0].DP[ct][2] = 4; nr[0].DP[ct][3] = 1;
					}
					ct += 1;
				}
			}
		}
		else if (element_type == 1) {
			ct = 0; 
			for (int z = 0; z < wt_py.size(); z++) {
				for (int i = 0; i < ien_py[wt_py[z]].size(); i++) {
					if (surface[wt_py[z]] == "S1") {
						ol[0].FP[ct][0] = 1; ol[0].FP[ct][1] = 3; ol[0].FP[ct][2] = 2;
					}
					else if (surface[wt_py[z]] == "S2") {
						ol[0].FP[ct][0] = 4; ol[0].FP[ct][1] = 1; ol[0].FP[ct][2] = 2;
					}
					else if (surface[wt_py[z]] == "S3") {
						ol[0].FP[ct][0] = 4; ol[0].FP[ct][1] = 2; ol[0].FP[ct][2] = 3;
					}
					else if (surface[wt_py[z]] == "S4") {
						ol[0].FP[ct][0] = 1; ol[0].FP[ct][1] = 4; ol[0].FP[ct][2] = 3;
					}
					ct += 1; 
				}
				ol[0].FP_2D[0] = 1;
				ol[0].FP_2D[1] = 2;
				ol[0].FP_2D[2] = 3;
			}
			ct = 0; 
			for (int z = 0; z < nrb_py.size(); z++) {
				for (int i = 0; i < ien_py[nrb_py[z]].size(); i++) {
					if (surface[nrb_py[z]]== "S1") {
						nr[0].DP[ct][0] = 1; nr[0].DP[ct][1] = 3; nr[0].DP[ct][2] = 2;
					}
					else if (surface[nrb_py[z]] == "S2") {
						nr[0].DP[ct][0] = 4; nr[0].DP[ct][1] = 1; nr[0].DP[ct][2] = 2;
					}
					else if (surface[nrb_py[z]] == "S3") {
						nr[0].DP[ct][0] = 4; nr[0].DP[ct][1] = 2; nr[0].DP[ct][2] = 3;
					}
					else if (surface[nrb_py[z]] == "S4") {
						nr[0].DP[ct][0] = 1; nr[0].DP[ct][1] = 4; nr[0].DP[ct][2] = 3;
					}
					ct += 1;
				}
				nr[0].DP_2D[0] = 1;
				nr[0].DP_2D[1] = 2;
				nr[0].DP_2D[2] = 3;
			}
		}

		//define BCIEN
		std::vector<std::vector<int>>iens_wt; //store the connectivity in the current physical group (2D vector)
		ct = 0; 
		for (int z = 0; z < wt_py.size(); z++) {
			for (j = 0; j < numele_py[wt_py[z]]; j++) { //loop through the element in the current physical group
				std::vector<int>local;
				for (k = 0; k < elenode2D; k++) {
					local.push_back(t.IEN[ol[0].FP[ct][k] - 1][ien_py[wt_py[z]][j] - 1]);
				}
				iens_wt.push_back(local);
				ct += 1; 
			}
		}
		t.BCIEN.push_back(iens_wt);
		std::vector<std::vector<int>>iens_nr; //store the connectivity in the current physical group (2D vector)
		ct = 0; 
		for (int z = 0; z < nrb_py.size(); z++) {
			for (j = 0; j < numele_py[nrb_py[z]]; j++) { //loop through the element in the current physical group
				std::vector<int>local;
				for (k = 0; k < elenode2D; k++) {
					local.push_back(t.IEN[nr[0].DP[ct][k] - 1][ien_py[nrb_py[z]][j] - 1]);
				}
				iens_nr.push_back(local);
				ct += 1; 
			}
		}
		t.BCIEN.push_back(iens_nr);
		t.GCOORD = new double*[t.NNODE];
		for (i = 0; i < t.NNODE; i++) {
			t.GCOORD[i] = new double[3];
		}
		for (int i = 0; i < t.NNODE; i++) {
			t.GCOORD[i][0] = VX[i];
			t.GCOORD[i][1] = VY[i];
			t.GCOORD[i][2] = VZ[i];
		}
		//Finish the mesh reading
		//adjust the node location if necessary
		if (nodeadj == 1) {
			for (int i = 0; i < t.NNODE; i++) {
				t.GCOORD[i][0] = VX[i];
				t.GCOORD[i][1] = VZ[i] + fs_offset;
				t.GCOORD[i][2] = -VY[i];
			}
		}
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
	if (element_type == 0) {
		localnode = new int*[owsfnumber + nrbsurfnumber];
		for (i = 0; i < owsfnumber + nrbsurfnumber; i++) {
			localnode[i] = new int[elenode2D];
		}
		//The following loop only attempts to extract the pattern of localnode[] from one corresponding global element since the assumption here is that all the elements in all physical groups
		//follow the same pattern. However, that is not rigorous enough since the corresponding local node for all 2D element in physical group do not necessary be the same. 
		for (i = 0; i < owsfnumber + nrbsurfnumber; i++) {
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
					if (ct == elenode2D) { //found the corresponding 3D element! (only one element per physical group is sampled)
						goto endloop;
					}
				}
			}
		endloop:;
		}
	}

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

	for (z = 0; z < owsfnumber; z++) {
		ol[z].LNA_norm = new int[elenode2D];
	}
	for (z = 0; z < nrbsurfnumber; z++) {
		nr[z].LNA_norm = new int[elenode2D];
	}

	if (element_type == 1) {
		for (z = 0; z < owsfnumber; z++) {
			ol[z].LNA_norm[0] = 1;
			ol[z].LNA_norm[1] = 2;
			ol[z].LNA_norm[2] = 3;
		}
		for (z = 0; z < nrbsurfnumber; z++) {
			nr[z].LNA_norm[0] = 1;
			nr[z].LNA_norm[1] = 2;
			nr[z].LNA_norm[2] = 3;
		}
	}

	if (element_type == 0) {
		int element_number = 0;
		for (z = 0; z < owsfnumber + nrbsurfnumber; z++) { //all physical groups except for connectivity
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
	for (z = 0; z < wt_pys_size; z++) { //loop through each wetted surface
		int elenum; 
		if (input_type == "Gmsh") {
			elenum = numele_py[wt_pys_num[z]]; //number of high-order on the wet surface s
		}
		else if (input_type == "Abaqus") {
			elenum = ol[0].FSNEL; 
			if (element_type == 0) {
				std::cout << "The Abaqus hex mesh is currently buggy, don't use it now. See the comment below for details" << std::endl;
				system("PAUSE "); 
				//IEN_py is not fully assigned. 
				//localnode is wrong as well. 
			}
		}

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

		if (input_type == "Gmsh") {
			ol[z].FP = new int*[ol[z].FSNEL];
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].FP[i] = new int[elenode2D];
			}
			//Define FP for hexahedral element (assume all the hexahedral element in the same wetted surface has the same FP definition)
			if (element_type == 0) {
				for (i = 0; i < ol[z].FSNEL; i++) {
					for (j = 0; j < elenode2D; j++) {
						ol[z].FP[i][j] = ol[z].FP_temp[j];
					}
				}
			}
		}

		//Define FP for tetrahedral element
		//Extract the local node pattern for tetrahedral element
		if (element_type == 1 && input_type == "Gmsh" && (mappingalgo == 4 || mappingalgo == 5)) {
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
			if (debug_IEN == 0) {
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
			int half; 
			for (i = 0; i < ol[z].FSNEL; i++) { //loop through each element
				if (i % 10000 == 0) {
					std::cout << i << std::endl;
				}
				for (j = 0; j < 3; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					half = round(i / 2.0); //search the second half first
					//search the second half first 
					for (k = half; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < 3; l++) {
							if (ol[z].IEN_gb[l][k] == ol[z].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
								ol[z].IEN_2D[j][i] = ol[z].IEN_2D[l][k];
								flag = 0; //turn off the flag used to assgin new number
								goto endloop2;
							}
						}
					}
					//search the second half
					for (k = 0; k < half; k++) {
						for (l = 0; l < 3; l++) {
							if (ol[z].IEN_gb[l][k] == ol[z].IEN_gb[j][i]) { //If this node has already been assigned, use the same numbering
								ol[z].IEN_2D[j][i] = ol[z].IEN_2D[l][k];
								flag = 0; //turn off the flag used to assgin new number
								goto endloop2;
							}
						}
					}
				endloop2:;
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

		//define GIDN (the counterpart of NRBA on non-reflecting boundary)
		std::vector<int>GIDN_total;
		for (i = 0; i < ol[z].FSNEL; i++) { //loop through each element
			for (j = 0; j < elenode2D; j++) { //the nodes in current element
				GIDN_total.push_back(ol[z].IEN_gb[j][i]);
			}
		}
		std::sort(GIDN_total.begin(), GIDN_total.end());
		//add the non repeated number into a new arrary
		std::vector<int>GIDN_unrepeat;
		GIDN_unrepeat.push_back(GIDN_total[0]);
		for (i = 1; i < GIDN_total.size(); i++) {
			flag = 1; //the flag used to store a new number
					  //search back if this number is alreay defined
			for (k = i - 1; k >= 0; k--) {
				if (GIDN_total[k] == GIDN_total[i]) { //this number is already defined
					flag = 0;
				}
				else {
					break;
				}
			}
			if (flag == 1) {
				GIDN_unrepeat.push_back(GIDN_total[i]);
			}
		}
		//store GIDN_total into GIDN
		ol[z].GIDNct = GIDN_unrepeat.size();
		ol[z].GIDN = new int[ol[z].GIDNct];
		for (i = 0; i < ol[z].GIDNct; i++) {
			ol[z].GIDN[i] = GIDN_unrepeat[i];
		}

		//Obtain GIDF (the global element numbering of each local wetted surface)
		if (element_type == 0 && (mappingalgo == 4 || mappingalgo == 5)) {
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
		std::cout << " " << std::endl; 
	}

	//================================Write the model file for MpCCI here================================//
	if (mappingalgo == 2) {
		nodesperelem = new int*[owsfnumber];
		for (z = 0; z < owsfnumber; z++) {
			nodesperelem[z] = new int[ol[z].FSNEL];
		}
		for (z = 0; z < owsfnumber; z++) {
			for (i = 0; i < ol[z].FSNEL; i++) {
				if (element_type == 0) {
					nodesperelem[z][i] = 4;
				}
				else if (element_type == 1) {
					nodesperelem[z][i] = 3;
				}
			}
		}

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
		if (input_type == "Gmsh") {
			nr[z].NEL_nrb = numele_py[nrb_pys_num[z]];
		}
		nr[z].IEN_gb = new int*[elenode2D]; //NRBELE_ARR
		for (i = 0; i < elenode2D; i++) {
			nr[z].IEN_gb[i] = new int[nr[z].NEL_nrb];
		}

		//Extract the local node pattern for hexahedral and tetrahedral element
		if (input_type == "Gmsh") {
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
		}

		if (element_type == 1 && input_type == "Gmsh") {
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

		int flag;

		std::vector<int>NRBA_total;
		for (i = 0; i < nr[z].NEL_nrb; i++) { //loop through each element
			for (j = 0; j < elenode2D; j++) { //the nodes in current element
				NRBA_total.push_back(nr[z].IEN_gb[j][i]);
			}
		}
		std::sort(NRBA_total.begin(), NRBA_total.end());
		//add the non repeated number into a new arrary
		std::vector<int>NRBA_unrepeat;
		NRBA_unrepeat.push_back(NRBA_total[0]);
		for (i = 1; i < NRBA_total.size(); i++) {
			flag = 1; //the flag used to store a new number
			//search back if this number is alreay defined
			for (k = i - 1; k >= 0; k--) {
				if (NRBA_total[k] == NRBA_total[i]) { //this number is already defined
					flag = 0;
				}
				else {
					break;
				}
			}
			if (flag == 1) {
				NRBA_unrepeat.push_back(NRBA_total[i]);
			}
		}
		//store NRBA_total into NRBA
		nr[z].NRBNODE = NRBA_unrepeat.size();
		nr[z].NRBA = new int[nr[z].NRBNODE];
		for (i = 0; i < nr[z].NRBNODE; i++) {
			nr[z].NRBA[i] = NRBA_unrepeat[i];
		}

		if (improvednrb == 1) {
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
						std::cout << "No corresponding 3D element is found" << std::endl;
						system("PAUSE ");
					}
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
		std::cout << " " << std::endl; 
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

	/*
	std::string name3 = "normal side check.txt";
	std::ofstream myfile2;
	myfile2.open(name3);

	for (z = 0; z < owsfnumber; z++) {
		myfile2 << "wetted surface " << z << std::endl;
		for (i = 0; i < ol[z].FSNEL; i++) {
			myfile2 << ol[z].norm[i][0] << " " << ol[z].norm[i][1] << " " << ol[z].norm[i][2] << std::endl;
		}
	} 
	for (z = 0; z < nrbsurfnumber; z++) {
		myfile2 << "NRB surface " << z << std::endl;
		for (i = 0; i < nr[z].NEL_nrb; i++) {
			myfile2 << nr[z].norm[i][0] << " " << nr[z].norm[i][1] << " " << nr[z].norm[i][2] << std::endl;
		}
	}
	*/

	//Clean the dynamic allocated arrary
	if (element_type == 1) {
		delete[] IEN_1D; 
	}

	std::cout << " " << std::endl;
	return t;
}
