//#include "stdafx.h"
#include "header.h"
#include "data.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <sstream>
#include <cstring>

//This function serves to associate the fluid interface mesh with structure interface mesh. The relationship is used for interface mapping later on
//Please note that each the sequence of the structural wetted surface has to be the same with the fluid wetted surface. For example, if the front structure wetted surface (ss[0])
//is the first surface, then the fluid wetted surface has to be the first surface ol[0] as well. 
//int** nodesperelem; 
void Neighborhood_search(double** GCOORD, int***LNA, int**IEN_flu, int NEL_flu) {
	extern OWETSURF ol[owsfnumber];
	extern STRU_WET_SURF ss[ssnumber]; //data structure used to store the properties on the structure wetted surface
	//extern int** nodesperelem;
	//extern int test; 
	int i, j, k, l, m, n, o, h, e, z, q, ii;

	std::string orphannodfile = "orphannode.txt";
	std::ofstream orphannodehd;
	orphannodehd.open(orphannodfile);

	int elenode2D;
	int elenode2D_ln; 
	if (element_type == 0) { //hex element
		elenode2D = NINT*NINT; 
		elenode2D_ln = 4; 
	}
	if (element_type == 1) { //tet element
		elenode2D = 3;
		elenode2D_ln = 3;
	}

	//Get the corner point of the fluid elements
	if (element_type == 0) { //hex element
		for (z = 0; z < owsfnumber; z++) {
			ol[z].IEN_flu_3D = new int[ol[z].FSNEL * 8];
			for (i = 0; i < ol[z].FSNEL; i++) {
				ol[z].IEN_flu_3D[i * 8 + 0] = IEN_flu[LNA[0][0][0] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 1] = IEN_flu[LNA[N][0][0] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 2] = IEN_flu[LNA[N][N][0] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 3] = IEN_flu[LNA[0][N][0] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 4] = IEN_flu[LNA[0][0][N] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 5] = IEN_flu[LNA[N][0][N] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 6] = IEN_flu[LNA[N][N][N] - 1][ol[z].GIDF[i] - 1];
				ol[z].IEN_flu_3D[i * 8 + 7] = IEN_flu[LNA[0][N][N] - 1][ol[z].GIDF[i] - 1];
			}
		}
		for (z = 0; z < owsfnumber; z++) {
			ol[z].IEN_flu_2D = new int*[4];
			for (i = 0; i < 4; i++) {
				ol[z].IEN_flu_2D[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < 4; j++) {
					//ol[z].IEN_flu_2D[j][i] = IEN_flu[ol[z].LNA_norm[j] - 1][ol[z].GIDF[i] - 1];
					ol[z].IEN_flu_2D[j][i] = ol[z].IEN_gb[ol[z].LNA_norm[i * 4 + j] - 1][i];
				}
			}
		}
		
	}
	if (element_type == 1) {
		ol[z].IEN_flu_3D = new int[ol[z].FSNEL * 4];
		for (i = 0; i < ol[z].FSNEL; i++) {
			ol[z].IEN_flu_3D[i * 4 + 0] = IEN_flu[0][ol[z].GIDF[i] - 1];
			ol[z].IEN_flu_3D[i * 4 + 1] = IEN_flu[1][ol[z].GIDF[i] - 1];
			ol[z].IEN_flu_3D[i * 4 + 2] = IEN_flu[2][ol[z].GIDF[i] - 1];
			ol[z].IEN_flu_3D[i * 4 + 3] = IEN_flu[3][ol[z].GIDF[i] - 1];
		}
		for (z = 0; z < owsfnumber; z++) {
			ol[z].IEN_flu_2D = new int*[3];
			for (i = 0; i < 3; i++) {
				ol[z].IEN_flu_2D[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < 3; j++) {
					ol[z].IEN_flu_2D[j][i] = IEN_flu[ol[z].FP[i][j] - 1][ol[z].GIDF[i] - 1];
				}
			}
		}
	}

	//First bring in the structural wetted surface mesh into the code
	//std::ifstream infile_algo5("C:/Users/lzhaok6/Desktop/DDG_datacheck.inp"); //The Abaqus input file
	std::ifstream infile_algo5("C:/Users/lzhaok6/Desktop/FSP_canopy_0.15_abaqus_MpCCI_explicit_sym_0.3048wl.inp");
	if (!infile_algo5) {
		std::cout << "can not open the structure input file" << std::endl;
		system("PAUSE ");
	}
	//read the file 
	int nodestart = 0;
	int nodeend = 0;
	int eleend = 0;
	int ele_line = 0;
	std::vector<std::vector<std::string>> output;
	int ct = -1;
	int endfile = 0;
	int elestart = 0; //the flag denote the start of element connectivity definition
	std::string csvElement;
	std::vector<std::string> csvColumn;
	std::string csvLine;
	std::vector<int> sidesets_start;
	std::vector<std::string> sidesets_name;
	std::vector<int> surface;
	int elementct = 1; 
	while (getline(infile_algo5, csvLine))
	{
		ct = ct + 1; //the current line number (starting from 0)
		std::istringstream csvStream(csvLine); //csvStream is a stream (turn the string csvLine to stream csvStream)
		csvColumn.clear();
		while (getline(csvStream, csvElement, ',')) //Get line from stream (csvStream) into string (csvElement)
		{
			csvColumn.push_back(csvElement);
		}
		output.push_back(csvColumn);
		if (csvColumn[0] == "*NODE" || csvColumn[0] == "*Node") {
			nodestart = ct + 1; //the node starts from the next line
		}
		if (csvColumn[0] == "*ELEMENT" || csvColumn[0] == "*Element") {
			if (elementct == 1) {
				nodeend = ct - 1;
				elestart = ct + 1; //the line where element definition starts
			}
			elementct += 1; 
		}
		if (csvColumn[0] == "*ELSET" || csvColumn[0] == "*Elset") {
			sidesets_start.push_back(ct + 1); //the line where sidesets (physical group) definition starts
			sidesets_name.push_back(csvColumn[1]);
			if (sidesets_start.size() == 1) {
				eleend = ct - 1; //the place where the first element set starts is where the element connectivity ends 
			}
		}
		if (csvColumn[0] == "*Surface") {
			surface.push_back(ct);
		}
		if (ct % 20000 == 0) {
			std::cout << ct << std::endl; //output which line is being read
		}
	}
	int ele_start_num = stoi(output[elestart][0]); //In case the first element is not numbered 1
	int node_start_num = stoi(output[nodestart][0]); //In case the first node is not numbered 1

	int NNODE_stru = nodeend - nodestart + 1; //total node number in structure model
	//int NEL_stru = eleend - elestart + 1; //total element number in structure model
	
	//Define connectivity matrix in the volumn mesh
	//First we need to know the largest Numbering of the element definition. (Often times in Abaqus input file, the element is not numbered sequentially. Thus, we need to know the largest element numbering to create the matrix to store the raw connectivities)
	std::vector<int> eleno; 
	for (i = elestart; i < eleend + 1; i++) {
		if (isdigit(output[i][0].back())) {
			eleno.push_back(stoi(output[i][0]));
		}
	}
	std::sort(eleno.begin(), eleno.end());
	int NEL_stru = eleno.back(); 

	int **IEN = new int*[4]; //4-node shell element or 3-node triangular element
	for (i = 0; i < 4; i++) {
		IEN[i] = new int[NEL_stru];
	}
	for (i = 0; i < 4; i++) {
		for (j = 0; j < NEL_stru; j++) {
			IEN[i][j] = 0; 
		}
	}
	std::vector<int> element_num; 
	for (i = elestart; i < eleend + 1; i++) {
		//check if this line is indeed the definition of connectivity
		if (isdigit(output[i][0].back())) { //check if the last character of the first string is a digit
			for (j = 0; j < output[i].size() - 1; j++) {
				//IEN[j][i - elestart] = stoi(output[i][j + 1]) - (node_start_num - 1); //change the string to integer (numbering from 1)
				IEN[j][stoi(output[i][0]) - ele_start_num] = stoi(output[i][j + 1]) - (node_start_num - 1); //change the string to integer (numbering from 1)
			}
			element_num.push_back(stoi(output[i][0]));
		}
	}
	std::sort(element_num.begin(), element_num.end());
	if (element_num[0] != ele_start_num) {
		std::cout << "ele_start_num is not the correct numbering of the first element in connectivity matrix definition" << std::endl;
		system(" PAUSE");
	}

	std::cout << "Remember to check the offset (under the *instance)" << std::endl;
	system("PAUSE ");
	for (z = 0; z < ssnumber; z++) {
		ss[z].GCOORD_stru = new double*[NNODE_stru];
		for (i = 0; i < NNODE_stru; i++) {
			ss[z].GCOORD_stru[i] = new double[3];
		}
		for (i = nodestart; i < nodeend + 1; i++) {
			//ss[z].GCOORD_stru[i - nodestart][0] = 144.220001 - stod(output[i][1]);
			//ss[z].GCOORD_stru[i - nodestart][1] = stod(output[i][2]) - 6.0;
			//ss[z].GCOORD_stru[i - nodestart][2] = -stod(output[i][3]);
			ss[z].GCOORD_stru[i - nodestart][0] = stod(output[i][1]);
			ss[z].GCOORD_stru[i - nodestart][1] = stod(output[i][2]) - 0.3048;
			ss[z].GCOORD_stru[i - nodestart][2] = stod(output[i][3]);
		}
	}

	//Find the element set name corresponding to the structural wetted surface
	//one surface could have several element sets
	std::vector<std::string> surface_name_1D;
	std::vector<std::vector<std::string>> surface_name;
	std::vector<int> eleset_orientation_1D;
	std::vector<std::vector<int>> eleset_orientation; //If the surface is positive, the value is 1; If the value is negative, the value is -1
	
	//store the element set information of each surface (name and orientation)
	for (i = 0; i < surface.size(); i++) {
		ct = 1;
		surface_name_1D.clear(); 
		eleset_orientation_1D.clear();
		while (output[surface[i] + ct][0].at(0) != '*') {
			surface_name_1D.push_back(output[surface[i] + ct][0]); //http://www.cplusplus.com/reference/string/string/substr/
			if (output[surface[i] + ct][1] == " SNEG") { //negative surface
				eleset_orientation_1D.push_back(-1);
			}
			else if (output[surface[i] + ct][1] == " SPOS") { //positive surface 
				eleset_orientation_1D.push_back(1);
			}
			else {
				std::cout << "The surface side is not recognized" << std::endl;
				system("PAUSE ");
			}
			ct += 1;
		}
		surface_name.push_back(surface_name_1D);
		eleset_orientation.push_back(eleset_orientation_1D); 
	}
	//Store the element number on the wetted surface
	//Currently, the code assumes one surface only has one element set. We need to accomodate the case with mutiple element sets corresponds to one surface. 
	std::vector<std::vector<int>> ele_num;
	std::vector<int> ele_num_clm; 
	std::vector<int> ele_num_sideset; 
	std::vector<std::vector<int>> surface_orientation;
	std::vector<int> surface_orientation_1D;
	ct = 0; 
	int eleset_linenum = 0;  //how many lines for the element set definition
	for (j = 0; j < surface.size(); j++) {
		//Read the nodes in all the subsets corresponding to this surface 
		ele_num_clm.clear();
		surface_orientation_1D.clear(); 
		ct = 0;
		for (n = 0; n < surface_name[j].size(); n++) { //One surface could corresponds to several sidesets
			for (i = 0; i < sidesets_start.size(); i++) { //Loop through side sets to find the one corresponding to the surface being searched. 
				if (sidesets_name[i].substr(7, std::string::npos) == surface_name[j][n]) { //If the name of this sideset is the same with one of the wetted surfaces
					//find the how many lines does the current sideset definition have
					eleset_linenum = 0;
					while (output[sidesets_start[i] + eleset_linenum][0].at(0) != '*') {
						eleset_linenum += 1;
					}
					if (eleset_linenum == 1) { //The element number is expressed as initial,end,internal way (only one line is used for definition) 
						for (k = stoi(output[sidesets_start[i]][0]); k < stoi(output[sidesets_start[i]][1]) + 1; k += stoi(output[sidesets_start[i]][2])) {
							ele_num_clm.push_back(k);
							surface_orientation_1D.push_back(eleset_orientation[j][n]);
							ct += 1;
						}
					}
					else {
						for (k = sidesets_start[i]; k < sidesets_start[i] + eleset_linenum; k++) {
							for (l = 0; l < output[k].size(); l++) {
								ele_num_clm.push_back(stoi(output[k][l]));
								surface_orientation_1D.push_back(eleset_orientation[j][n]);
								ct += 1;
							}
						}
					}
				}
			}
		}
		ele_num.push_back(ele_num_clm);
		ele_num_sideset.push_back(ct);
		surface_orientation.push_back(surface_orientation_1D);
	}
	for (i = 0; i < surface.size(); i++) {
		if (ele_num_sideset[i] != ele_num[i].size()) {
			std::cout << "There is a problem with the counting of structural elements in each sideset" << std::endl;
			system("PAUSE "); 
		}
	}
	
	//Store the element connectivity on the wetted surface
	nodesperelem = new int*[ssnumber];
	//get the number of structural wetted elements on each surface
	for (z = 0; z < ssnumber; z++) {
		//Define nodesperelem (used by MpCCI functions in data.c and adapter.c)
		nodesperelem[z] = new int[ele_num_sideset[z]];
		ss[z].elenode = new int[ele_num_sideset[z]]; 
		ss[z].ELE_stru = ele_num_sideset[z];
		ss[z].IEN_stru = new int[4 * ss[z].ELE_stru];
		ss[z].IEN_stru_norm = new int*[4];
		for (i = 0; i < 4; i++) {
			ss[z].IEN_stru_norm[i] = new int[ss[z].ELE_stru];
		}
		ct = 0;
		for (j = 0; j < ele_num_sideset[z]; j++) { //loop through all the element in the wetted 
			ss[z].IEN_stru[ct * 4 + 0] = IEN[0][ele_num[z][j] - 1 - (ele_start_num - 1)];
			ss[z].IEN_stru[ct * 4 + 1] = IEN[1][ele_num[z][j] - 1 - (ele_start_num - 1)];
			ss[z].IEN_stru[ct * 4 + 2] = IEN[2][ele_num[z][j] - 1 - (ele_start_num - 1)];
			ss[z].IEN_stru[ct * 4 + 3] = IEN[3][ele_num[z][j] - 1 - (ele_start_num - 1)];
			ct += 1;
		}
		if (ct != ss[z].ELE_stru) {
			std::cout << "The wetted surface elements are not fully defined" << std::endl;
			system("PAUSE ");
		}
		ct = 0;
		for (j = 0; j < ele_num_sideset[z]; j++) { //loop through all the element in the wetted surface
			if (IEN[3][ele_num[z][j] - 1 - (ele_start_num - 1)] == 0) { //If the last element is zero, it means this element is triangular
				ss[z].elenode[j] = 3; //triangular element
				nodesperelem[z][j] = 3;
			}
			else {
				ss[z].elenode[j] = 4; //quad element
				nodesperelem[z][j] = 4;
			}
			for (k = 0; k < ss[z].elenode[j]; k++) {
				if (surface_orientation[z][j] == 1) { //out of the structure, revert the node sequence make it into the structure
					ss[z].IEN_stru_norm[ss[z].elenode[j] - k - 1][ct] = IEN[k][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[3][ct] = IEN[0][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[2][ct] = IEN[1][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[1][ct] = IEN[2][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[0][ct] = IEN[3][ele_num[z][j] - 1 - (ele_start_num - 1)];
				}
				else if (surface_orientation[z][j] == -1) { //into the structure, keep the original node sequence
					ss[z].IEN_stru_norm[k][ct] = IEN[k][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[0][ct] = IEN[0][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[1][ct] = IEN[1][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[2][ct] = IEN[2][ele_num[z][j] - 1 - (ele_start_num - 1)];
					//ss[z].IEN_stru_norm[3][ct] = IEN[3][ele_num[z][j] - 1 - (ele_start_num - 1)];
				}
			}
			ct += 1;
		}
		if (ct != ss[z].ELE_stru) {
			std::cout << "The wetted surface elements are not fully defined" << std::endl;
			system("PAUSE ");
		}
	}
	
	
	//Create the pressure on structure gauss nodes (the pressure to be received) 
	//Used in interface_mapping.cpp 
	for (z = 0; z < ssnumber; z++) {
		ss[z].P_gs = new double[ss[z].ELE_stru*(hprefg + 1)*(hprefg + 1)];
		for (i = 0; i < ss[z].ELE_stru*(hprefg + 1)*(hprefg + 1); i++) {
			ss[z].P_gs[i] = 0.0;
		}
	}
	
	//Determine the normal direction of the structure wetted surface elements for the calculation of the fluid force (pointing into the structure)
	for (z = 0; z < ssnumber; z++) {
		ss[z].norm_stru = new double*[ss[z].ELE_stru];
		for (i = 0; i < ss[z].ELE_stru; i++) {
			ss[z].norm_stru[i] = new double[3];
		}
		for (i = 0; i < ss[z].ELE_stru; i++) {
			for (j = 0; j < 3; j++) {
				ss[z].norm_stru[i][j] = 0.0;
			}
		}
	}
	if (mappingalgo == 4 || mappingalgo == 5) {
		double x1, x2, x3; double y1, y2, y3; double z1, z2, z3;
		double ax, ay, az; double bx, by, bz;
		double n1, n2, n3; double absn;
		for (z = 0; z < ssnumber; z++) {
			for (i = 0; i < ss[z].ELE_stru; i++) {
				x1 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[0][i] - 1][0]; x2 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[1][i] - 1][0]; x3 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[2][i] - 1][0];
				y1 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[0][i] - 1][1]; y2 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[1][i] - 1][1]; y3 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[2][i] - 1][1];
				z1 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[0][i] - 1][2]; z2 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[1][i] - 1][2]; z3 = ss[z].GCOORD_stru[ss[z].IEN_stru_norm[2][i] - 1][2];
				ax = x2 - x1; ay = y2 - y1; az = z2 - z1;
				bx = x3 - x2; by = y3 - y2; bz = z3 - z2;
				n1 = ay*bz - az*by; n2 = az*bx - ax*bz; n3 = ax*by - ay*bx;
				absn = sqrt(pow(n1, 2) + pow(n2, 2) + pow(n3, 2));
				ss[z].norm_stru[i][0] = n1 / absn; ss[z].norm_stru[i][1] = n2 / absn; ss[z].norm_stru[i][2] = n3 / absn;
				//std::cout << " " << std::endl;
			}
		}
	}

	//Define the local node connectivity on structural wetted surface mesh
	TD_LOCAL_NODEstruct td;
	td = TD_LOCAL_NODE(1);
	for (z = 0; z < ssnumber; z++) {
		ss[z].LNA_stru = new int[2 * 2];
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				ss[z].LNA_stru[i * 2 + j] = td.LNA[i][j];
			}
		}
	}
	//local connectivity of the quadrature nodes
	td = TD_LOCAL_NODE(hprefg);
	for (z = 0; z < ssnumber; z++) {
		ss[z].LNA_gs = new int*[hprefg + 1];
		for (i = 0; i < hprefg + 1; i++) {
			ss[z].LNA_gs[i] = new int[hprefg + 1];
		}
		for (i = 0; i < hprefg + 1; i++) {
			for (j = 0; j < hprefg + 1; j++) {
				ss[z].LNA_gs[i][j] = td.LNA[i][j];
			}
		}
	}

	//Define the quadrature node on linear structure element 
	LOBATTOstruct bb;
	bb = LOBATTO(hprefg);
	GLLQUADstruct ff;
	ff = GLLQUAD(bb.Z, bb.WL, hprefg, FEM);
	for (z = 0; z < ssnumber; z++) {
		for (i = 0; i < hprefg + 1; i++) {
			ss[z].W_stru[i] = ff.W[i];
		}
	}
	//Define the linear shape function on structure element (the value is on the gauss node)
	double* basep; //interpolation node local coordinate
	double* gsp; //quadrature node local coordinate
	double nomx, nomy;
	double denomx, denomy;
	basep = new double[2]; 
	for (i = 0; i < 2; i++) {
		basep[i] = 0.0;
	}
	for (i = 0; i < 2; i++) {
		basep[i] = -1.0 + i*(2.0 / (2 - 1));
	}
	gsp = new double[hprefg + 1]; 
	for (i = 0; i < hprefg + 1; i++) {
		gsp[i] = ff.S[i];
	}
	for (z = 0; z < ssnumber; z++) {
		ss[z].phi_stru = new double**[4];
		for (i = 0; i < 4; i++) { //for all the points in that element
			ss[z].phi_stru[i] = new double*[hprefg + 1];
			for (j = 0; j < hprefg + 1; j++) {
				ss[z].phi_stru[i][j] = new double[hprefg + 1];
			}
		}
		for (i = 0; i < 4; i++) {
			for (j = 0; j < hprefg + 1; j++) {
				for (k = 0; k < hprefg + 1; k++) {
					ss[z].phi_stru[i][j][k] = 0.0;
				}
			}
		}
		for (h = 0; h < 2; h++) {
			for (k = 0; k < 2; k++) {
				for (i = 0; i < hprefg + 1; i++) {
					for (j = 0; j < hprefg + 1; j++) {
						nomx = 1.0; nomy = 1.0; //multiplier initialization
						denomx = 1.0; denomy = 1.0; //multiplier initialization
						for (e = 0; e < 2; e++) { //loop through nominator and denominator in basis function expression
							if (e != h) {
								nomx *= (gsp[i] - basep[e]);
								denomx *= (basep[h] - basep[e]);
							}
							if (e != k) {
								nomy *= (gsp[j] - basep[e]);
								denomy *= (basep[k] - basep[e]);
							}
						}
						ss[z].phi_stru[ss[z].LNA_stru[h * 2 + k] - 1][i][j] = (nomx / denomx)*(nomy / denomy); //tensor product
					}
				}
			}
		}
	}

	TD_LOCAL_NODEstruct ct_gs;
	//Interpolate the ss.GCOORD_stru to obtain the GCOORD_stru_gs (the coordinate of gauss nodes)
	for (z = 0; z < ssnumber; z++) {
		ss[z].GCOORD_stru_gs = new double[ss[z].ELE_stru*(hprefg + 1)*(hprefg + 1) * 3];
		//Define int**IEN_stru_gs for the later nodal force integration
		ss[z].IEN_stru_gs = new int*[(hprefg + 1)*(hprefg + 1)];
		for (i = 0; i < (hprefg + 1)*(hprefg + 1); i++) {
			ss[z].IEN_stru_gs[i] = new int[ss[z].ELE_stru];
		}
		ct = 1;
		//Since we will use Gauss-Legendre nodes, different elements don't share nodes. As a result, the nodes could be numbered sequentially
		for (i = 0; i < ss[z].ELE_stru; i++) {
			for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {
				ss[z].IEN_stru_gs[j][i] = ct;
				ct += 1;
			}
		}
		//Derive the coordinate of the gauss nodes 
		//TD_LOCAL_NODEstruct ct_gs;
		ct_gs = TD_LOCAL_NODE(hprefg);
		TD_LOCAL_NODEstruct ct_ln;
		ct_ln = TD_LOCAL_NODE(1);
		LOCAL_NODEstruct lna;
		lna = LOCAL_NODE(1);
		LOCAL_SHAPEstruct ls_ln; //ln means linear
		ls_ln = LOCAL_SHAPE(lna.LNA, 1, hprefg, 1); //Get the 2D linear shape function value on Nth order Gauss-Legendre nodes
		for (i = 0; i < ss[z].ELE_stru; i++) { //loop through all the elements
			//quad element
			if (ss[z].elenode[i] == 4) {
				for (m = 0; m < hprefg + 1; m++) { //m, n, o stands for internal nodes in one element (all but except for corner nodes)
					for (n = 0; n < hprefg + 1; n++) {
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 0] = 0.0;
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 1] = 0.0;
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 2] = 0.0;
						for (j = 0; j < 2; j++) { //j, k, l stands for corner nodes
							for (k = 0; k < 2; k++) {
								ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 0] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + ct_ln.LNA[j][k] - 1] - 1][0] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
								ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 1] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + ct_ln.LNA[j][k] - 1] - 1][1] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
								ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1) + 2] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + ct_ln.LNA[j][k] - 1] - 1][2] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
							}
						}
					}
				}
			}
			//triangular element (we assume the quad node of triangular element is less than the quad element)
			else {
				//Derive the 3-point quadrature node location on the interface triangle elements 
				double xi[3]; double eta[3]; double phi[3];
				xi[0] = 1.0 / 6.0; xi[1] = 2.0 / 3.0; xi[2] = 1.0 / 6.0;  //From Ragab FEM course note
				eta[0] = 1.0 / 6.0; eta[1] = 1.0 / 6.0; eta[2] = 2.0 / 3.0;
				for (m = 0; m < 3; m++) { //m stands for the gauss point
					ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 0] = 0.0;
					ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 1] = 0.0;
					ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 2] = 0.0;
					phi[0] = 1 - xi[m] - eta[m]; phi[1] = xi[m]; phi[2] = eta[m];
					for (j = 0; j < 3; j++) { //j stands for corner nodes
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 0] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + j] - 1][0] * phi[j];
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 1] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + j] - 1][1] * phi[j];
						ss[z].GCOORD_stru_gs[3 * (ss[z].IEN_stru_gs[m][i] - 1) + 2] += ss[z].GCOORD_stru[ss[z].IEN_stru[i * 4 + j] - 1][2] * phi[j];
					}
				}
			}
		}
	}

	for (z = 0; z < ssnumber; z++) {
		ss[z].gs_num = ss[z].ELE_stru * (hprefg + 1) * (hprefg + 1); //the total number of gauss points on the structural surface
		ss[z].gs_flu_global = new double*[ss[z].gs_num]; //used to store the global coordinate of the projected gauss point on fluid surface
		for (i = 0; i < ss[z].gs_num; i++) {
			ss[z].gs_flu_global[i] = new double[3];
		}
		ss[z].gs_flu_local = new double*[ss[z].gs_num];
		for (i = 0; i < ss[z].gs_num; i++) {
			ss[z].gs_flu_local[i] = new double[3];
		}
		ss[z].gs_flu = new int[ss[z].gs_num];
		ss[z].orphan_flag_gs = new int[ss[z].gs_num];
	}

	for (z = 0; z < owsfnumber; z++) {
		ol[z].flu_local = new double*[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
		ol[z].flu_stru_global = new double*[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
		ol[z].flu_stru = new int[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
		for (i = 0; i < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); i++) {
			ol[z].flu_local[i] = new double[2];
			ol[z].flu_stru_global[i] = new double[3];
		}
		ol[z].orphan_flag_flu = new int[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
	}

	for (z = 0; z < ssnumber; z++) {
		ss[z].FP_flu = new int*[ss[z].ELE_stru * (hprefg + 1) * (hprefg + 1)];
		for (i = 0; i < ss[z].ELE_stru * (hprefg + 1) * (hprefg + 1); i++) {
			ss[z].FP_flu[i] = new int[elenode2D];
		}
	}

	int* FP_searched = new int[elenode2D];
	for (i = 0; i < elenode2D; i++) {
		FP_searched[i] = 0.0;
	}

	//Generate the MpCCI model file (read in the readfile function in data.c)
	if (mappingalgo == 4 || mappingalgo == 5) {
		int flag;
		for (z = 0; z < owsfnumber; z++) {
			ss[z].IEN_stru_MpCCI = new int[4 * ss[z].ELE_stru]; //Connecvitity matrix of wetted surface (after removing the free surface elements)
			ct = 0; //count the node number assigned
			std::vector<int>dummy;
			for (i = 0; i < ss[z].ELE_stru; i++) { //loop through each element
				for (j = 0; j < ss[z].elenode[i]; j++) { //the nodes in current element
					flag = 1; //Initiate the flag to 1 
					for (k = 0; k < i; k++) { //see if the number has already been assigned by the nodes in previous elements
						for (l = 0; l < 4; l++) {
							if (ss[z].IEN_stru[k * 4 + l] == ss[z].IEN_stru[i * 4 + j]) { //If this node has already been assigned, use the same numbering
								ss[z].IEN_stru_MpCCI[i * 4 + j] = ss[z].IEN_stru_MpCCI[k * 4 + l];
								flag = 0; //turn off the flag to assgin new number
							}
							else {
								//If the number has not assigned yet the flag is still 1, thus a new number could be assigned. 
							}
						}
					}
					if (flag == 1) {
						ct += 1;
						dummy.push_back(ss[z].IEN_stru[i * 4 + j]); //associate the local 2D node with the global node numbering 
						ss[z].IEN_stru_MpCCI[i * 4 + j] = ct; //assign a new number to MpCCI element connectivity
					}
				}
			}
			//ss.Node_stru = dummy.size();
			ss[z].Node_stru = ct;
			if (ct != dummy.size()) {
				std::cout << "There is a problem with ss.Node_stru" << std::endl;
				system("PAUSE ");
			}
			ss[z].Node_glob = new int[ss[z].Node_stru];
			for (i = 0; i < ss[z].Node_stru; i++) {
				ss[z].Node_glob[i] = dummy[i];
			}
		}
	}

	//start writing the file
	if (mappingalgo == 4 || mappingalgo == 5) {
		std::ofstream myfile_algo5;
		myfile_algo5.open("model.txt");
		for (z = 0; z < owsfnumber; z++) {
			std::string wetsurface_name;
			wetsurface_name = "EF wetsurface" + std::to_string(z + 1) + " 3 1 " + std::to_string(0); //gonna be quad element no matter what element type is used by the fluid
			myfile_algo5 << wetsurface_name << std::endl;
			myfile_algo5 << "NODES " << ss[z].Node_stru << std::endl;
			for (i = 0; i < ss[z].Node_stru; i++) {
				myfile_algo5 << i << " " << ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][0] << " " << ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][1] << " " << ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][2] << " " << std::endl;
			}
			myfile_algo5 << "ELEMENTS " << ss[z].ELE_stru << std::endl;
			//Output connectivity matrix
			for (i = 0; i < ss[z].ELE_stru; i++) {
				myfile_algo5 << i;
				for (j = 0; j < ss[z].elenode[i]; j++) {
					myfile_algo5 << " " << ss[z].IEN_stru_MpCCI[i * 4 + j] - 1; //node numbering starts from 0 in model file
				}
				myfile_algo5 << std::endl;
			}
		}
	}

	int gs_nearest; //used to store the nearest fluid node on the fluid FSI interface
	int flu_nearest = 0; //used to store the nearest structure node on the structural wetted surface
	//Container to store the searching result (one gauss point corresponds to a fluid element)
	int** FSNEL_stru_gs; //[gauss point][structural wetted surface element]

	//Eigen::Matrix3d fdot; 
	double fdot00, fdot01, fdot02, fdot10, fdot11, fdot12, fdot20, fdot21, fdot22;
	double f0, f1, f2; 
	double dx[3]; 
	//Eigen::Matrix3d A; 
	double A00, A01, A02, A10, A11, A12, A20, A21, A22; 
	double b0, b1, b2; 
	//Eigen::Vector3d f;
	//Eigen::Vector3d b; //RHS coefficient
	//Eigen::Vector3d dx;
	Eigen::Vector3d x; 
	//Eigen::Matrix3d AA;
	double AA00, AA01, AA02, AA10, AA11, AA12, AA20, AA21, AA22;
	//Eigen::Vector3d BB;
	double BB0, BB1, BB2;
	//Eigen::Vector3d local; //store the local coordinate
	double local[3]; 
	//Eigen::Matrix2d AA_f;
	double AA_f00, AA_f01, AA_f10, AA_f11; 
	//Eigen::Vector2d BB_f;
	double BB_f0, BB_f1; 
	//Eigen::Vector2d local_f; //store the local coordinate
	double local_f[2]; 
	double search_range = 1; //set the searching range to 1m (only the element with all points in the range will be searched)
	double range[2];
	range[1] = search_range; 
	double distance[2]; //the distance from the strutural gauss node to its projected point
	distance[0] = search_range; 
	int flag = 1; 
	double tol = 1e-6; 
	double diff = 1.0; 
	double xu = 0.0; double yu = 0.0; double zu = 0.0; //the cooridnate of the projected point
	double xu_searched = 0.0; double yu_searched = 0.0; double zu_searched = 0.0; //The projected node closest to the structure gauss node
	double xu_k = 0.0; double yu_k = 0.0; double zu_k = 0.0; //the projected point coordinate at step k 
	double xu_k_1 = 0.0; double yu_k_1 = 0.0; double zu_k_1 = 0.0; //the projected point coordinate at step k + 1
	double x1, x2, x3, x4, x5, x6, x7, x8;
	double y1, y2, y3, y4, y5, y6, y7, y8;
	double z1, z2, z3, z4, z5, z6, z7, z8;
	double xn; double yn; double zn; 
	int searched_ele; //the fluid element found after node projection
	int orphan; //the flag indicating if the node is an orphan. 
	int inelement; 
	//Node projection from structure to fluid (Project the structural Gauss points to fluid elements)
	//Structural wetted elements do not share gauss points
	//We assume the structural wetted surface is imported as a whole surface and the fluid wetted surface could have several groups
	
	int localnode[2];
	localnode[0] = 3; localnode[1] = (hprefg + 1)*(hprefg + 1); 
	std::vector<int> localnode_seq_1D;
	std::vector<std::vector<int>> localnode_seq;
	std::vector<int> LNA_seq_1D;
	std::vector<std::vector<int>> LNA_seq;
	localnode_seq_1D.push_back(0); localnode_seq_1D.push_back(1); localnode_seq_1D.push_back(2);
	localnode_seq.push_back(localnode_seq_1D);
	localnode_seq_1D.clear(); 
	for (j = 0; j < hprefg + 1; j++) {
		for (q = 0; q < hprefg + 1; q++) {
			localnode_seq_1D.push_back(j*(hprefg + 1) + q);
		}
	}
	localnode_seq.push_back(localnode_seq_1D);
	LNA_seq_1D.push_back(0); LNA_seq_1D.push_back(1); LNA_seq_1D.push_back(2);
	LNA_seq.push_back(LNA_seq_1D); 
	LNA_seq_1D.clear();
	for (j = 0; j < hprefg + 1; j++) {
		for (q = 0; q < hprefg + 1; q++) {
			LNA_seq_1D.push_back(ct_gs.LNA[j][q] - 1);
		}
	}
	LNA_seq.push_back(LNA_seq_1D);
	int ele_id; 

	/*
	double xs = 6.2405066666666702;
	double ys = -4.2162466866666666;
	double zs = -0.65991649066666680;
	double dist = 0.0;
	double smallestdist = 1000.0; 
	int nodefound; 
	e = 0; 
	for (k = 0; k < ol[0].FSNEL; k++) {
		for (l = 0; l < elenode2D_ln; l++) {	
			dist = pow(pow(GCOORD[ol[e].IEN_flu_2D[l][k] - 1][0] - xs, 2) + pow(GCOORD[ol[e].IEN_flu_2D[l][k] - 1][1] - ys, 2) + pow(GCOORD[ol[e].IEN_flu_2D[l][k] - 1][2] - zs, 2), 0.5);
			if (dist < smallestdist) {
				smallestdist = dist; 
				nodefound = ol[e].IEN_flu_2D[l][k]; 
			}
		}
	}
	*/


	//Create fluid boundary coordinate array for fast access
	for (z = 0; z < owsfnumber; z++) {
		ol[z].GCOORD_fs = new double[ol[z].FSNEL*elenode2D_ln * 3];
	}
	for (z = 0; z < owsfnumber; z++) {
		for (i = 0; i < ol[z].FSNEL; i++) {
			for (l = 0; l < elenode2D_ln; l++) {
				for (j = 0; j < 3; j++) {
					ol[z].GCOORD_fs[i*elenode2D_ln * 3 + l * 3 + j] = GCOORD[ol[z].IEN_flu_2D[l][i] - 1][j]; 
				}
			}
		}
	}

	//There are two types of orphans. The first type does not have corresponding element (the projected node does not belong to any element). However, the node is still in the searching range. Thus, the value on the closest point is assigned. 
	//The second type is out of the searching range, no value should be assigned on this point. 
	int flag2; //the flag to determine if the node is within the searching range. (loop through all the fluid interface nodes and find that the whole structure element is within at least one fluid node's searching range)
	//the difference between flag and flag2 is that flag determines if the current should be projected to the element being searched and flag2 determins if the current node does not fall within the searching range of any node on target surface. 
	int ite_num; 
	int orphan_ct = 0; //count how many total orphan node there is
	for (e = 0; e < ssnumber; e++) {
		for (i = 0; i < ss[e].ELE_stru; i++) {
			if (i % 1000 == 0) {
				std::cout << i << " structure elements scanned out of: " << ss[e].ELE_stru << std::endl;
			}
			if (ss[e].elenode[i] == 3) { //quad element
				ele_id = 0; 
			}
			else {
				ele_id = 1;
			}
			for (j = 0; j < localnode[ele_id]; j++) { //j stands for gauss points
				flag2 = 0; //If it falls in one of the fluid nodes' searching range, turn on this flag. 
				orphan = 1; //let's first assume the node is an orphan. If the corresponding element is found, the flag is turned to 0 
				//Start searching for fluid element for gauss node IEN_gs[j][i]
				//Calculate the distance from the current the structural gauss point to the fluid element in the search range
				distance[0] = search_range;
				range[1] = search_range; //initiate the fluid node distance to be equal to the searching range for every gauss node
				double xs = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 0];
				double ys = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 1];
				double zs = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 2];
				for (k = 0; k < ol[e].FSNEL; k++) { //looping through all fluid wetted elements
					flag = 1;
					//First determine if the fluid element is inside the searching range (If not, we would not project the structural gauss node)
					for (l = 0; l < elenode2D_ln; l++) {
						range[0] = pow(pow(ol[e].GCOORD_fs[k*elenode2D_ln * 3 + l * 3 + 0] - xs, 2) + pow(ol[e].GCOORD_fs[k*elenode2D_ln * 3 + l * 3 + 1] - ys, 2) + pow(ol[e].GCOORD_fs[k*elenode2D_ln * 3 + l * 3 + 2] - zs, 2), 0.5);
						if (range[0] < range[1]) {
							flag2 = 1; 
							range[1] = range[0]; //range[1] is used to store the shortest distance so far
							gs_nearest = ol[e].IEN_flu_2D[l][k];
						}
						//At the same time, store the node with the shorted distance to the gauss node under searching. If the node is orphan, we use the value on this fluid node for the structure gauss node
						if (range[0] > search_range) { //One of the corner nodes in the element is out of the searching range. Jump through this element
							flag = 0;
						}
					}
					if (flag == 1) { //the element should be searched (within the searching range)
						//Find the global coordinate of the orthogonally projected point on the current element
						//Find the coordinate of the orthogonally projected fluid interpolation point on the current structure element
						//Find the mathematica code in /Users/zhaokuanlu/OneDrive/post_processing/Node_projection_method2.nb
						//Also check out the theory in /Users/zhaokuanlu/OneDrive/post_processing/Node_projection_method2.pdf
						x1 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 0 * 3 + 0]; x2 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 1 * 3 + 0]; x3 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 2 * 3 + 0];
						y1 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 0 * 3 + 1]; y2 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 1 * 3 + 1]; y3 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 2 * 3 + 1];
						z1 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 0 * 3 + 2]; z2 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 1 * 3 + 2]; z3 = ol[e].GCOORD_fs[k*elenode2D_ln * 3 + 2 * 3 + 2];
						xn = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 0]; yn = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 1]; zn = ss[e].GCOORD_stru_gs[3 * (ss[e].IEN_stru_gs[LNA_seq[ele_id][j]][i] - 1) + 2];
						b0 = -(-(-x1 + x3)*xn - (-y1 + y3)*yn - (-z1 + z3)*zn);
						b1 = -(-(-x1 + x2)*xn - (-y1 + y2)*yn - (-z1 + z2)*zn);
						b2 = -(-(x2*y1 - x3*y1 - x1*y2 + x3 *y2 + x1 *y3 - x2 *y3) *z1 -
							y1*(-x2 *z1 + x3 *z1 + x1 *z2 - x3 *z2 - x1 *z3 + x2 *z3) -
							x1*(y2 *z1 - y3 *z1 - y1 *z2 + y3 *z2 + y1 *z3 - y2 *z3));
						A00 = -x1 + x3;
						A01 = -y1 + y3;
						A02 = -z1 + z3;
						A10 = -x1 + x2;
						A11 = -y1 + y2;
						A12 = -z1 + z2;
						A20 = y2 *z1 - y3 *z1 - y1 *z2 + y3 *z2 + y1 *z3 - y2 *z3;
						A21 = -x2 *z1 + x3 *z1 + x1 *z2 - x3 *z2 - x1 *z3 + x2 *z3;
						A22 = x2 *y1 - x3 *y1 - x1 *y2 + x3 *y2 + x1 *y3 - x2 *y3;
						//x = A.inverse()*b;
						//xu = x(0); yu = x(1); zu = x(2);
						//We calculated the equation above using matlab instead of the matrix calculation in Eigen
						xu = (b2*(A01*A12 - A02*A11)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b1*(A01*A22 - A02*A21)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) + (b0*(A11*A22 - A12*A21)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						yu = (b1*(A00*A22 - A02*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b2*(A00*A12 - A02*A10)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b0*(A10*A22 - A12*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						zu = (b2*(A00*A11 - A01*A10)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b1*(A00*A21 - A01*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) + (b0*(A10*A21 - A11*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						//We need to determine at this stage if the projected node is within the searched element, if not we need to pass this element. 
						double xi = 0.0;
						double eta = 0.0;
						double zeta = 0.0;  //converged local coordinate
						//the local coordinate in tetrahedral and hexahedral fluid elements is derived differently. 
						if (element_type == 1) { //tetrahedral element
							//Need to change the name of the matrix and vectors. also need to change the if(){} below 
							//determine the local coordinate in the tetrahedral element
							//There is no need to perform iteration like hex element
							x1 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 0] - 1][0]; x2 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 1] - 1][0]; x3 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 2] - 1][0]; x4 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 3] - 1][0];
							y1 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 0] - 1][1]; y2 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 1] - 1][1]; y3 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 2] - 1][1]; y4 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 3] - 1][1];
							z1 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 0] - 1][2]; z2 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 1] - 1][2]; z3 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 2] - 1][2]; z4 = GCOORD[ol[e].IEN_flu_3D[k * 4 + 3] - 1][2];
							AA00 = x2 - x1; AA00 = x3 - x1; AA02 = x4 - x1;
							AA10 = y2 - y1; AA11 = y3 - y1; AA12 = y4 - y1;
							AA20 = z2 - z1; AA21 = z3 - z1; AA22 = z4 - z1;
							BB0 = xu - x1; BB1 = yu - y1; BB2 = zu - z1;
							//local = AA.inverse()*BB;
							local[0] = (BB2*(AA01*AA12 - AA02*AA11)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) - (BB1*(AA01*AA22 - AA02*AA21)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) + (BB0*(AA11*AA22 - AA12*AA21)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20);
							local[1] = (BB1*(AA00*AA22 - AA02*AA20)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) - (BB2*(AA00*AA12 - AA02*AA10)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) - (BB0*(AA10*AA22 - AA12*AA20)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20);
							local[2] = (BB2*(AA00*AA11 - AA01*AA10)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) - (BB1*(AA00*AA21 - AA01*AA20)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20) + (BB0*(AA10*AA21 - AA11*AA20)) / (AA00*AA11*AA22 - AA00*AA12*AA21 - AA01*AA10*AA22 + AA01*AA12*AA20 + AA02*AA10*AA21 - AA02*AA11*AA20);
							//determine if the point is inside the searched element
							if (local[0] > -1e-5 && (local[0] - 1) < 1e-5 && local[1] > -1e-5 && (local[1] - 1) < 1e-5 && local[2] > -1e-5 && (local[2] - 1) < 1e-5
								&& (local[0] + local[1] + local[2] - 1) < 1e-5) {
								//The point is indeed inside the searched element and we can store the local coordinate in that element
								xi = local[0]; //xi
								eta = local[1]; //eta
								zeta = local[2]; //zeta
								orphan = 0;
								inelement = 1;
							}
							else { //The node is an orphan (the orphan flag is not turned off)
								//Search for the nearest point on the fluid FSI interface
								//orphan = 1;
								inelement = 0;
							}
						}
						if (element_type == 0) { //hex element
							double xi_k, eta_k, zeta_k;
							double xi_k_1, eta_k_1, zeta_k_1;
							//Determine the local coordinate in the hexahedral element to see if the point is within the element being searched
							//Since the equation system is nonlinear, we use the Newton iteration and set the initial value to be 0
							//The derivative (fdot) is obtained by Mathematica: /Users/zhaokuanlu/OneDrive/post_processing/Hex_global_to_local_coordinate_Newton_iteration.nb
							//See the file for theory behind this: /Users/zhaokuanlu/OneDrive/post_processing/Hex_global_to_local_coordinate_Newton_iteration.pdf
							xi_k = 0.1; eta_k = 0.1; zeta_k = 0.1;
							x1 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 0] - 1][0]; x2 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 1] - 1][0]; x3 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 2] - 1][0]; x4 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 3] - 1][0]; x5 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 4] - 1][0]; x6 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 5] - 1][0]; x7 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 6] - 1][0]; x8 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 7] - 1][0];
							y1 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 0] - 1][1]; y2 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 1] - 1][1]; y3 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 2] - 1][1]; y4 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 3] - 1][1]; y5 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 4] - 1][1]; y6 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 5] - 1][1]; y7 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 6] - 1][1]; y8 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 7] - 1][1];
							z1 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 0] - 1][2]; z2 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 1] - 1][2]; z3 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 2] - 1][2]; z4 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 3] - 1][2]; z5 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 4] - 1][2]; z6 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 5] - 1][2]; z7 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 6] - 1][2]; z8 = GCOORD[ol[e].IEN_flu_3D[k * 8 + 7] - 1][2];
							diff = 1;
							ite_num = 0; 
							while (diff > tol) {
								fdot00 = -(1.0 / 8.0) *(1.0 - eta_k) *x1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 + zeta_k);
								fdot01 = -(1.0 / 8.0) *x1*(1.0 - xi_k) *(1.0 - zeta_k) + 1.0 / 8.0 *x4*(1.0 - xi_k) *(1.0 - zeta_k) -
									1.0 / 8.0 *x2*(1.0 + xi_k) *(1.0 - zeta_k) + 1.0 / 8.0 *x3*(1.0 + xi_k) *(1.0 - zeta_k) -
									1.0 / 8.0 *x5*(1.0 - xi_k) *(1.0 + zeta_k) + 1.0 / 8.0 *x8*(1.0 - xi_k) *(1.0 + zeta_k) -
									1.0 / 8.0 *x6*(1.0 + xi_k) *(1.0 + zeta_k) + 1.0 / 8.0 *x7*(1.0 + xi_k) *(1.0 + zeta_k);
								fdot02 = -(1.0 / 8.0) *(1.0 - eta_k) *x1*(1.0 - xi_k) - 1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - xi_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 - xi_k) + 1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 - xi_k) -
									1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 + xi_k) - 1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 + xi_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + xi_k) + 1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + xi_k);
								fdot10 = -(1.0 / 8.0) *(1.0 - eta_k) *y1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *y3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *y4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *y5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *y7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *y8*(1.0 + zeta_k);
								fdot11 = -(1.0 / 8.0) *(1.0 - xi_k) *y1*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *y3*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *y4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - xi_k) *y5*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *y7*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *y8*(1.0 + zeta_k);
								fdot12 = -(1.0 / 8.0) *(1.0 - eta_k) *(1.0 - xi_k) *y1 - 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y2 -
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y3 - 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y4 +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y5 + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y6 +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y7 + 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y8;
								fdot20 = -(1.0 / 8.0) *(1.0 - eta_k) *z1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *z3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *z4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *z5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *z7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *z8*(1.0 + zeta_k);
								fdot21 = -(1.0 / 8.0) *(1.0 - xi_k) *z1*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *z3*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *z4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - xi_k) *z5*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *z7*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *z8*(1.0 + zeta_k);
								fdot22 = -(1.0 / 8.0) *(1.0 - eta_k) *(1.0 - xi_k) *z1 - 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z2 -
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z3 - 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z4 +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z5 + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z6 +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z7 + 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z8;
								f0 = -xu + 1.0 / 8.0 *(1.0 - eta_k) *x1*(1.0 - xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 + xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 + xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 - xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 - xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + xi_k) *(1.0 + zeta_k);
								f1 = -yu + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y1*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y3*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y4*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y5*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y7*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y8*(1.0 + zeta_k);
								f2 = -zu + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z1*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z3*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z4*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z5*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z7*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z8*(1.0 + zeta_k);
								//dx = -fdot.inverse()*f; //delta x^k (checked)
								dx[0] = (f1*(fdot01*fdot22 - fdot02*fdot21)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) - (f2*(fdot01*fdot12 - fdot02*fdot11)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) - (f0*(fdot11*fdot22 - fdot12*fdot21)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20);
								dx[1] = (f2*(fdot00*fdot12 - fdot02*fdot10)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) - (f1*(fdot00*fdot22 - fdot02*fdot20)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) + (f0*(fdot10*fdot22 - fdot12*fdot20)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20);
								dx[2] = (f1*(fdot00*fdot21 - fdot01*fdot20)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) - (f2*(fdot00*fdot11 - fdot01*fdot10)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20) - (f0*(fdot10*fdot21 - fdot11*fdot20)) / (fdot00*fdot11*fdot22 - fdot00*fdot12*fdot21 - fdot01*fdot10*fdot22 + fdot01*fdot12*fdot20 + fdot02*fdot10*fdot21 - fdot02*fdot11*fdot20);
								diff = pow(pow(dx[0], 2) + pow(dx[1], 2) + pow(dx[2], 2), 0.5);
								xi_k_1 = xi_k + dx[0]; eta_k_1 = eta_k + dx[1]; zeta_k_1 = zeta_k + dx[2];
								xi_k = xi_k_1; eta_k = eta_k_1; zeta_k = zeta_k_1; //update the value
								ite_num += 1;
								if (ite_num > 1000) {
									//std::cout << "Maximum iteration reached (1000)." << std::endl;
									//The convergence is not achieved because the point (xu, yu, zu) is too far away from the element. 
									//orphannodehd << x1 << ", " << x2 << ", " << x3 << ", " << x4 << ", " << x5 << ", " << x6 << ", " << x7 << ", " << x8 << std::endl;
									//orphannodehd << y1 << ", " << y2 << ", " << y3 << ", " << y4 << ", " << y5 << ", " << y6 << ", " << y7 << ", " << y8 << std::endl;
									//orphannodehd << z1 << ", " << z2 << ", " << z3 << ", " << z4 << ", " << z5 << ", " << z6 << ", " << z7 << ", " << z8 << std::endl;
									//orphannodehd << xu << ", " << yu << ", " << zu << std::endl;
									//system("PAUSE ");
									xi_k = 2; eta_k = 2; zeta_k = 2; //Enfore a large enough value so that it would not pass the orphan test outside the iteration.
									break;
								}
							}
							//std::cout << "Iteration number " << ite_num << std::endl; 
							//orphannodehd << x1 << ", " << x2 << ", " << x3 << ", " << x4 << ", " << x5 << ", " << x6 << ", " << x7 << ", " << x8 << std::endl;
							//orphannodehd << y1 << ", " << y2 << ", " << y3 << ", " << y4 << ", " << y5 << ", " << y6 << ", " << y7 << ", " << y8 << std::endl;
							//orphannodehd << z1 << ", " << z2 << ", " << z3 << ", " << z4 << ", " << z5 << ", " << z6 << ", " << z7 << ", " << z8 << std::endl;
							//orphannodehd << xu << ", " << yu << ", " << zu << std::endl;
							xi = xi_k; eta = eta_k; zeta = zeta_k;
							if (xi + 1.1 > -1e-5 && xi - 1.1 < 1e-5 && eta + 1.1 > -1e-5 && eta - 1.1 < 1e-5 && zeta + 1.1 > -1e-5 && zeta - 1.1 < 1e-5) {
								//since the interface element might be warped, we need to leave a margin. (See 2018 Fall report 15). 
								//The point is indeed inside the searched element. We can turn off the orphan flag
								//orphannodehd << x1 << ", " << x2 << ", " << x3 << ", " << x4 << ", " << x5 << ", " << x6 << ", " << x7 << ", " << x8 << std::endl;
								//orphannodehd << y1 << ", " << y2 << ", " << y3 << ", " << y4 << ", " << y5 << ", " << y6 << ", " << y7 << ", " << y8 << std::endl;
								//orphannodehd << z1 << ", " << z2 << ", " << z3 << ", " << z4 << ", " << z5 << ", " << z6 << ", " << z7 << ", " << z8 << std::endl;
								//orphannodehd << xu << ", " << yu << ", " << zu << std::endl;
								orphan = 0;
								inelement = 1;
							}
							else { //The node is an orphan (the orphan flag is not turned off)
								inelement = 0;
							}
						}
						//End determining if the projected node is in the searched element and its local coordinate

						//If the projected node is indeed within the searched element, calculate the orthogonal distance 
						//Calculate the distance between the structural gauss nodes and the projected node
						if (inelement == 1) { //if the node is within the current fluid element
							distance[1] = pow(pow(xn - xu, 2) + pow(yn - yu, 2) + pow(zn - zu, 2), 0.5);
							//if (distance[1] < distance[0]) {
							if (distance[1] - distance[0] < -1e-5) { //the newly found projection has the shorted distance, need to update the projected node 
								//store the local coordinate in the fluid mesh of the projected structure gauss node
								ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][0] = xi; //xi
								ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][1] = eta; //eta
								ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][2] = zeta; //zeta
								//orphannodehd << xi << " " << eta << " " << zeta << std::endl; 
								//store the global coordinate in the fluid mesh of the projected structure gauss node
								xu_searched = xu; yu_searched = yu; zu_searched = zu;
								distance[0] = distance[1]; //update the shortest orthogonal distance so far
								searched_ele = ol[e].GIDF[k];
								for (ii = 0; ii < elenode2D; ii++) {
									FP_searched[ii] = ol[e].FP[k][ii];
								}
							}
						}
					}
				} //finish searching all the fluid elements in the current group
				if (orphan == 1 && flag2 == 1) { //If the node is still an orphan after looping through all the fluid elements, we need to search the fluid mesh to determine which is the closest node
					//the orphan node and the other nodes should be treated separately
					//Instead of storing the corresponding element in gs_flu, the stored value becomes the fluid mesh "node" number 
					//At the same time, we keep track of which gauss node is orphan
					ss[e].gs_flu[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = gs_nearest;
					//ss.orphan_flag_gs.push_back(i*(hprefg + 1)*(hprefg + 1) + j);
					ss[e].orphan_flag_gs[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = 1;
					orphan_ct += 1; 
				}
				else if (orphan == 1 && flag2 == 0) { //The node is an orphan and also it is not within the searching range. Thus, no value should be given to that quadrature point
					//we set ss[e].gs_flu to 0 to show that this node should not receive any value
					ss[e].gs_flu[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = 0;
					ss[e].orphan_flag_gs[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = 1;
					orphan_ct += 1;
				}
				else {
					//the corresponding fluid element and the projected node coordinate should be found at this stage
					//store the corresponding element 
					ss[e].gs_flu[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = searched_ele;
					ss[e].orphan_flag_gs[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]] = 0;
					for (ii = 0; ii < elenode2D; ii++) {
						ss[e].FP_flu[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][ii] = FP_searched[ii];
					}
				}
				//Check if the searched node is on the element surface
				if (element_type == 1) {
					double xi_fd = ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][0];
					double eta_fd = ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][1];
					double zeta_fd = ss[e].gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][2];
					if (abs(xi_fd) < 1e-5 || abs(eta_fd) < 1e-5 || abs(zeta_fd) < 1e-5 || abs(xi_fd + eta_fd + zeta_fd - 1) < 1e-5) {
						//the node is on one of the surfaces
					}
					else {
						std::cout << "The local node is not on one of the surfaces" << std::endl;
						system("PAUSE ");
					}
				}
				//store the projected node location
				ss[e].gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][0] = xu_searched;
				ss[e].gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][1] = yu_searched;
				ss[e].gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + localnode_seq[ele_id][j]][2] = zu_searched;
			}
		}
	}

	for (z = 0; z < owsfnumber; z++) {
		delete[] ol[z].GCOORD_fs;
	}

	//Derive the Gauss-Legendre points on fluid elements (for displacement integration)
	//We can first use equal amount of interpolation node for GL nodes. This would guarantee the accuracy of the integrated nodal boundary force.
	//The same connectivity matrix could be used (ol[z].IEN_gb[j][i])
	//Derive the coordinate of the gauss nodes 
	if (element_type == 0) {
		TD_LOCAL_NODEstruct ct_ln;
		ct_ln = TD_LOCAL_NODE(1);
		LOCAL_NODEstruct lna; 
		lna = LOCAL_NODE(1);
		LOCAL_SHAPEstruct ls_ln; //ln means linear
		ls_ln = LOCAL_SHAPE(lna.LNA, 1, hprefg_flu, 1); //Get the 2D linear shape function value on Nth order Gauss-Legendre nodes
		for (z = 0; z < owsfnumber; z++) {
			ol[z].GCOORD_flu_gs = new double*[ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1)];
			for (i = 0; i < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); i++) {
				ol[z].GCOORD_flu_gs[i] = new double[3];
			}
			int LNA_2D[2][2];
			//LNA_2D[0][0] = ol[z].LNA_2D[0][0]; LNA_2D[1][0] = ol[z].LNA_2D[N][0]; LNA_2D[0][1] = ol[z].LNA_2D[0][N]; LNA_2D[1][1] = ol[z].LNA_2D[N][N];
			for (i = 0; i < ol[z].FSNEL; i++) { //loop through all the elements
				LNA_2D[0][0] = ol[z].LNA_2D[i*NINT*NINT + 0 * NINT + 0]; LNA_2D[1][0] = ol[z].LNA_2D[i*NINT*NINT + N * NINT + 0]; LNA_2D[0][1] = ol[z].LNA_2D[i*NINT*NINT + 0 * NINT + N]; LNA_2D[1][1] = ol[z].LNA_2D[i*NINT*NINT + N * NINT + N];
				for (m = 0; m < (hprefg_flu + 1); m++) { //m, n, o stands for internal nodes in one element (all but except for corner nodes)
					for (n = 0; n < (hprefg_flu + 1); n++) {
						ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][0] = 0.0; //m and n equivalent to LNA[m][n]
						ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][1] = 0.0;
						ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][2] = 0.0;
						for (j = 0; j < 2; j++) { //j, k stands for corner nodes
							for (k = 0; k < 2; k++) {
								//ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][0] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][0] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg_flu + 1) + n];
								//ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][1] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][1] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg_flu + 1) + n];
								//ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][2] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][2] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg_flu + 1) + n];
								ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][0] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][0] * ls_ln.SHL_2D[2][LNA_2D[j][k] - 1][m*(hprefg_flu + 1) + n];
								ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][1] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][1] * ls_ln.SHL_2D[2][LNA_2D[j][k] - 1][m*(hprefg_flu + 1) + n];
								ol[z].GCOORD_flu_gs[i*(hprefg_flu + 1)*(hprefg_flu + 1) + m*(hprefg_flu + 1) + n][2] += GCOORD[ol[z].IEN_gb[LNA_2D[j][k] - 1][i] - 1][2] * ls_ln.SHL_2D[2][LNA_2D[j][k] - 1][m*(hprefg_flu + 1) + n];
							}
						}
					}
				}
			}
		}
	}
	if (element_type == 1) {
		//Derive the 3-point quadrature node location on the interface triangle elements 
		double xi[3]; double eta[3]; double phi[3];
		xi[0] = 1.0 / 6.0; xi[1] = 2.0 / 3.0; xi[2] = 1.0 / 6.0;  //From Ragab FEM course note
		eta[0] = 1.0 / 6.0; eta[1] = 1.0 / 6.0; eta[2] = 2.0 / 3.0;
		for (z = 0; z < owsfnumber; z++) {
			ol[z].GCOORD_flu_gs = new double*[ol[z].FSNEL * 3];
			for (i = 0; i < ol[z].FSNEL * 3; i++) {
				ol[z].GCOORD_flu_gs[i] = new double[3];
			}
			for (i = 0; i < ol[z].FSNEL; i++) { //loop through all the elements
				for (m = 0; m < 3; m++) { //m stands for the gauss point to be determined
					ol[z].GCOORD_flu_gs[i * 3 + m][0] = 0.0;
					ol[z].GCOORD_flu_gs[i * 3 + m][1] = 0.0;
					ol[z].GCOORD_flu_gs[i * 3 + m][2] = 0.0;
					phi[0] = 1 - xi[m] - eta[m]; phi[1] = xi[m]; phi[2] = eta[m];
					for (j = 0; j < 3; j++) { //j stands for corner nodes
						ol[z].GCOORD_flu_gs[i * 3 + m][0] += GCOORD[ol[z].IEN_gb[j][i] - 1][0] * phi[j];
						ol[z].GCOORD_flu_gs[i * 3 + m][1] += GCOORD[ol[z].IEN_gb[j][i] - 1][1] * phi[j];
						ol[z].GCOORD_flu_gs[i * 3 + m][2] += GCOORD[ol[z].IEN_gb[j][i] - 1][2] * phi[j];
					}
				}
			}
		}
	}

	//Create a coordinate array of structure boundary nodes for fast access
	for (z = 0; z < ssnumber; z++) {
		ss[z].GCOORD_stru_fs = new double*[NNODE_stru];
		for (k = 0; k < NNODE_stru; k++) {
			ss[z].GCOORD_stru_fs[k] = new double[3];
		}
		ct = 0;
		for (k = 0; k < ss[z].ELE_stru; k++) { //looping through structure wetted elements
			for (l = 0; l < ss[z].elenode[k]; l++) {
				for (i = 0; i < 3; i++) {
					ss[z].GCOORD_stru_fs[ct][i] = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + l] - 1][i];
				}
				ct += 1;
			}
		}
	}
	
	//Node projection from fluid to structure (Project the fluid interpolation points to structural elements
	//Eigen::Matrix2d fdot_f;
	double fdot_f00, fdot_f01, fdot_f10, fdot_f11; 
	//Eigen::Vector2d f_f;
	double f_f0, f_f1; 
	//Eigen::Vector2d dx_f;
	double dx_f[2]; 
	diff = 1.0;
	range[1] = search_range;
	distance[0] = search_range;
	int elenode2D_gs = 0; 
	if (element_type == 0) {
		elenode2D_gs = (hprefg_flu + 1)*(hprefg_flu + 1);
	}
	if (element_type == 1) {
		elenode2D_gs = 3;
	}
	orphan_ct = 0; 
	for (z = 0; z < owsfnumber; z++) {
		for (i = 0; i < ol[z].FSNEL; i++) { //loop through fluid elements
			if (i % 1000 == 0) {
				std::cout << i << " fluid elements scanned out of: " << ol[z].FSNEL << std::endl;
			}
			for (j = 0; j < elenode2D_gs; j++) { //loop through each (gauss) node on the fluid FSI 
				flag2 = 0;
				orphan = 1; //let's first assume the node is an orphan. If the corresponding element is found, the flag is turned to 0 and stays the same for the current node
				distance[0] = search_range;
				range[1] = search_range;
				//Start searching for structural element for fluid interpolation node ol[z].IEN_gb[j][i]
				//Calculate the distance from the current fluid interpolation node to the structural element in the search range
				//First determine if the structure element is inside the searching range
				double xf = ol[z].GCOORD_flu_gs[i*elenode2D_gs + j][0];
				double yf = ol[z].GCOORD_flu_gs[i*elenode2D_gs + j][1];
				double zf = ol[z].GCOORD_flu_gs[i*elenode2D_gs + j][2];
				ct = 0;
				for (k = 0; k < ss[z].ELE_stru; k++) { //looping through structure wetted elements
					flag = 1;
					for (l = 0; l < ss[z].elenode[k]; l++) {
						range[0] = pow(pow(xf - ss[z].GCOORD_stru_fs[ct][0], 2) + pow(yf - ss[z].GCOORD_stru_fs[ct][1], 2) + pow(zf - ss[z].GCOORD_stru_fs[ct][2], 2), 0.5);
						if (range[0] < range[1]) {
							flag2 = 1;
							range[1] = range[0]; //range[1] is used to store the shortest distance so far
							if (debug_algo5 == 0) {
								flu_nearest = ss[z].IEN_stru_MpCCI[k * 4 + l];
							}
						}
						//At the same time, store the node with the shorted distance to the gauss node under searching. If the node is orphan, we use the value on this fluid node for the structure gauss node
						if (range[0] > search_range) { //One of the corner nodes in the element is out of the searching range. Jump through this element
							flag = 0;
						}
						ct += 1;
					}
					if (flag == 1) { //the element should be searched (within the searching range)
						//Find the coordinate of the orthogonally projected fluid interpolation point on the current structure 4-node shell element
						//Find the mathematica code in /Users/zhaokuanlu/OneDrive/post_processing/Node_projection_method2.nb
						//Also check out the theory in /Users/zhaokuanlu/OneDrive/post_processing/Node_projection_method2.pdf
						x1 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 0] - 1][0]; x2 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 1] - 1][0]; x3 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 2] - 1][0];
						y1 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 0] - 1][1]; y2 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 1] - 1][1]; y3 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 2] - 1][1];
						z1 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 0] - 1][2]; z2 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 1] - 1][2]; z3 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 2] - 1][2];
						xn = xf; yn = yf; zn = zf; //the fluid node to be projected
						b0 = -(-(-x1 + x3)*xn - (-y1 + y3)*yn - (-z1 + z3)*zn);
						b1 = -(-(-x1 + x2)*xn - (-y1 + y2)*yn - (-z1 + z2)*zn);
						b2 = -(-(x2*y1 - x3*y1 - x1*y2 + x3 *y2 + x1 *y3 - x2 *y3) *z1 -
							y1*(-x2 *z1 + x3 *z1 + x1 *z2 - x3 *z2 - x1 *z3 + x2 *z3) -
							x1*(y2 *z1 - y3 *z1 - y1 *z2 + y3 *z2 + y1 *z3 - y2 *z3));
						A00 = -x1 + x3; 
						A01 = -y1 + y3; 
						A02 = -z1 + z3;
						A10 = -x1 + x2; 
						A11 = -y1 + y2; 
						A12 = -z1 + z2;
						A20 = y2 *z1 - y3 *z1 - y1 *z2 + y3 *z2 + y1 *z3 - y2 *z3;
						A21 = -x2 *z1 + x3 *z1 + x1 *z2 - x3 *z2 - x1 *z3 + x2 *z3;
						A22 = x2 *y1 - x3 *y1 - x1 *y2 + x3 *y2 + x1 *y3 - x2 *y3;
						//x = A.inverse()*b;
						//xu = x(0); yu = x(1); zu = x(2);
						xu = (b2*(A01*A12 - A02*A11)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b1*(A01*A22 - A02*A21)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) + (b0*(A11*A22 - A12*A21)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						yu = (b1*(A00*A22 - A02*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b2*(A00*A12 - A02*A10)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b0*(A10*A22 - A12*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						zu = (b2*(A00*A11 - A01*A10)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) - (b1*(A00*A21 - A01*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20) + (b0*(A10*A21 - A11*A20)) / (A00*A11*A22 - A00*A12*A21 - A01*A10*A22 + A01*A12*A20 + A02*A10*A21 - A02*A11*A20);
						//We need to determine at this stage if the projected node is within the searched element, if not we need to pass this element. 
						double xi_k, eta_k;
						double xi_k_1, eta_k_1;
						double xi, eta; //converged local coordinate
						double xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4; //pseudo structural node coordinate 
						double xup, yup; //pseudo projected point coordinate
						//The derivative is obtained by Mathematica: /Users/zhaokuanlu/OneDrive/post_processing/Quad_global_to_local_coordinate_Newton_iteration.nb
						//See the: /Users/zhaokuanlu/OneDrive/post_processing/Quad_global_to_local_coordinate_Newton_iteration.pdf for more information
						//determine the local coordinate in the hexahedral element
						//Since the equation system is nonlinear, we use the Newton iteration and set the initial value to be 0 
						xi_k = 0.5; eta_k = 0.5; //initialize the iteration
						if (ss[z].elenode[k] == 4) { //if the structure element is 4 node, then we define the fourth node, otherwise, no need to define it. 
							x4 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 3] - 1][0]; y4 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 3] - 1][1]; z4 = ss[z].GCOORD_stru[ss[z].IEN_stru[k * 4 + 3] - 1][2];
						}
						if (x1 == x2 && x1 == x3 && x2 == x3) { //the element is parallel to yz plane (same x coordinate), use y(xi,eta) and z(xi,eta) to determine xi and eta
							xp1 = y1; xp2 = y2; xp3 = y3; xp4 = y4;
							yp1 = z1; yp2 = z2; yp3 = z3; yp4 = z4;
							xup = yu; yup = zu;
						}
						else if (y1 == y2 && y1 == y3 && y2 == y3) { //the element is parallel to xz plane, use x(xi,eta) and z(xi,eta) to determine xi and eta
							xp1 = x1; xp2 = x2; xp3 = x3; xp4 = x4;
							yp1 = z1; yp2 = z2; yp3 = z3; yp4 = z4;
							xup = xu; yup = zu;
						}
						else if (z1 == z2 && z1 == z3 && z2 == z3) { //the element is parallel to xy plane, use x(xi,eta) and y(xi,eta) to determine xi and eta
							xp1 = x1; xp2 = x2; xp3 = x3; xp4 = x4;
							yp1 = y1; yp2 = y2; yp3 = y3; yp4 = y4;
							xup = xu; yup = yu;
						}
						else { //the element is a 3D surface, use any two dimensions (we choose to use x and y here)
							xp1 = x1; xp2 = x2; xp3 = x3; xp4 = x4;
							yp1 = y1; yp2 = y2; yp3 = y3; yp4 = y4;
							xup = xu; yup = yu;
						}
						if (ss[z].elenode[k] == 4) { //the structure element is a quad element
							diff = 1.0;
							ite_num = 0;
							while (diff > tol) {
								fdot_f00 = -(1.0 / 4.0) *(1.0 - eta_k) *xp1 + 1.0 / 4.0 *(1.0 - eta_k) *xp2 + 1.0 / 4.0 *(1.0 + eta_k) *xp3 -
									1.0 / 4.0 *(1.0 + eta_k) *xp4;
								fdot_f01 = -(1.0 / 4.0) *xp1*(1.0 - xi_k) + 1.0 / 4.0 *xp4*(1.0 - xi_k) - 1.0 / 4.0 *xp2*(1.0 + xi_k) +
									1.0 / 4.0 *xp3*(1.0 + xi_k);
								fdot_f10 = -(1.0 / 4.0) *(1.0 - eta_k) *yp1 + 1.0 / 4.0 *(1.0 - eta_k) *yp2 + 1.0 / 4.0 *(1.0 + eta_k) *yp3 -
									1.0 / 4.0 *(1.0 + eta_k) *yp4;
								fdot_f11 = -(1.0 / 4.0) *(1.0 - xi_k) *yp1 - 1.0 / 4.0 *(1.0 + xi_k) *yp2 + 1.0 / 4.0 *(1.0 + xi_k) *yp3 +
									1.0 / 4.0 *(1.0 - xi_k) *yp4;
								f_f0 = -xup + 1.0 / 4.0 *(1.0 - eta_k) *xp1*(1.0 - xi_k) + 1.0 / 4.0 *(1.0 + eta_k) *xp4*(1.0 - xi_k) +
									1.0 / 4.0 *(1.0 - eta_k) *xp2*(1.0 + xi_k) + 1.0 / 4.0 *(1.0 + eta_k) *xp3*(1.0 + xi_k);
								f_f1 = -yup + 1.0 / 4.0 *(1.0 - eta_k) *(1.0 - xi_k) *yp1 + 1.0 / 4.0 *(1.0 - eta_k) *(1.0 + xi_k) *yp2 +
									1.0 / 4.0 *(1.0 + eta_k) *(1.0 + xi_k) *yp3 + 1.0 / 4.0 *(1.0 + eta_k) *(1.0 - xi_k) *yp4;
								dx_f[0] = (f_f1*fdot_f01) / (fdot_f00*fdot_f11 - fdot_f01*fdot_f10) - (f_f0*fdot_f11) / (fdot_f00*fdot_f11 - fdot_f01*fdot_f10);
								dx_f[1] = (f_f0*fdot_f10) / (fdot_f00*fdot_f11 - fdot_f01*fdot_f10) - (f_f1*fdot_f00) / (fdot_f00*fdot_f11 - fdot_f01*fdot_f10);
								//dx_f = -fdot_f.inverse()*f_f; //delta x^k
								diff = pow(pow(dx_f[0], 2) + pow(dx_f[1], 2), 0.5);
								xi_k_1 = xi_k + dx_f[0]; eta_k_1 = eta_k + dx_f[1];
								xi_k = xi_k_1; eta_k = eta_k_1; //update the value
								ite_num += 1;
								if (ite_num > 1000) {
									//std::cout << "Maximum iteration reached (1000)." << std::endl;
									xi_k = 2; eta_k = 2; //Enfore a large enough value so that it would not pass the orphan test outside the iteration.
									break;
								}
							}
							//std::cout << ite_num << std::endl; 
							xi = xi_k; eta = eta_k; //the local coordinate obtained! 
							if (xi + 1 > -1e-5 && xi - 1 < 1e-5 && eta + 1 > -1e-5 && eta - 1 < 1e-5) {
								//The point is indeed inside the searched element and we can store the local coordinate in that element
								//orphannodehd << x1 << ", " << x2 << ", " << x3 << ", " << x4 << std::endl;
								//orphannodehd << y1 << ", " << y2 << ", " << y3 << ", " << y4 << std::endl;
								//orphannodehd << z1 << ", " << z2 << ", " << z3 << ", " << z4 << std::endl;
								//orphannodehd << xu << ", " << yu << ", " << zu << std::endl;
								//orphannodehd << xi << ", " << eta << std::endl;
								orphan = 0; //turn off the orphan flag
								inelement = 1;
							}
							else { //The node is an orphan
								   //orphan = 1;
								inelement = 0;
							}
							//End determining if the projected node is in the searched element and its local coordinate
						}
						else { //the structure element is a triangular element (no need to iterate to find the local coordinate)
							//see "post_processing/Local_coordinate_triangular_element.pdf" for the derivation
							AA_f00 = xp2 - xp1; AA_f01 = xp3 - xp1;
							AA_f10 = yp2 - yp1; AA_f11 = yp3 - yp1;
							BB_f0 = xup - xp1; BB_f1 = yup - yp1;
							//local_f = AA_f.inverse()*BB_f;
							local_f[0] = (AA_f11*BB_f0) / (AA_f00*AA_f11 - AA_f01*AA_f10) - (AA_f01*BB_f1) / (AA_f00*AA_f11 - AA_f01*AA_f10);
							local_f[1] = (AA_f00*BB_f1) / (AA_f00*AA_f11 - AA_f01*AA_f10) - (AA_f10*BB_f0) / (AA_f00*AA_f11 - AA_f01*AA_f10); 
							xi = local_f[0]; eta = local_f[1];
							if (xi > -1e-5 && (xi - 1) < 1e-5 && eta > -1e-5 && (eta - 1) < 1e-5 && (xi + eta - 1) < 1e-5) {
								//The point is indeed inside the searched element and we can store the local coordinate in that element
								orphan = 0; //turn off the orphan flag
								inelement = 1;
							}
							else { //The node is an orphan
								//orphan = 1;
								inelement = 0;
							}
						}
						//If the projected node is indeed within the searched element, calculate the orthogonal distance 
						//Calculate the distance between the fluid interpolation nodes and the projected node
						if (inelement == 1) { //if the node is not an orphan
							distance[1] = pow(pow(xn - xu, 2) + pow(yn - yu, 2) + pow(zn - zu, 2), 0.5);
							//if (distance[1] < distance[0]) {
							if (distance[1] - distance[0] < -1e-5) {
								xu_searched = xu; yu_searched = yu; zu_searched = zu;
								distance[0] = distance[1];
								ol[z].flu_local[i*elenode2D_gs + j][0] = xi; //xi
								ol[z].flu_local[i*elenode2D_gs + j][1] = eta; //eta
								searched_ele = k + 1;
							}
						}
					}
				}//finish the searching through the structural wetted elements 
				if (orphan == 1 && flag2 == 1) { //If the node is still an orphan, we need to search the structure mesh to determine which is the closest node
					//the orphan node and the other nodes should be treated separately
					//Instead of storing the corresponding element in flu_stru, the stored value becomes the structure "node" number 
					//At the same time, we keep track of which fluid node is orphan
					ol[z].flu_stru[i*elenode2D_gs + j] = flu_nearest;
					ol[z].orphan_flag_flu[i*elenode2D_gs + j] = 1;
					orphan_ct += 1; 
				}
				else if (orphan == 1 && flag2 == 0) {
					ol[z].flu_stru[i*elenode2D_gs + j] = 0;
					ol[z].orphan_flag_flu[i*elenode2D_gs + j] = 1;
					orphan_ct += 1;
				}
				else {
					//store the corresponding element 
					ol[z].flu_stru[i*elenode2D_gs + j] = searched_ele;
					ol[z].orphan_flag_flu[i*elenode2D_gs + j] = 0;
				}
				//store the projected node location
				ol[z].flu_stru_global[i*elenode2D_gs + j][0] = xu_searched;
				ol[z].flu_stru_global[i*elenode2D_gs + j][1] = yu_searched;
				ol[z].flu_stru_global[i*elenode2D_gs + j][2] = zu_searched;
			}
			//Finish with all the nodes in the current fluid element
		}
		//Finish with all the nodes in the current physical group
	}
	//Finish all the fluid nodes
	std::cout << " " << std::endl;

	//clean the dynamic variables
	for (i = 0; i < 4; i++) {
		delete[] IEN[i];
	}
	delete[] IEN; 
	for (z = 0; z < owsfnumber; z++) {
		delete[] ol[z].IEN_flu_3D;
	}
	if (element_type == 0) { //hex element
		for (z = 0; z < owsfnumber; z++) {
			for (i = 0; i < 4; i++) {
				delete[] ol[z].IEN_flu_2D[i];
			}
			delete[] ol[z].IEN_flu_2D;
		}
	}
	if (element_type == 1) { //hex element
		for (z = 0; z < owsfnumber; z++) {
			for (i = 0; i < 3; i++) {
				delete[] ol[z].IEN_flu_2D[i];
			}
			delete[] ol[z].IEN_flu_2D;
		}
	}
	for (z = 0; z < ssnumber; z++) {
		for (i = 0; i < 4; i++) {
			delete[] ss[z].IEN_stru_norm[i];
		}
		delete[] ss[z].IEN_stru_norm;
	}
	for (z = 0; z < ssnumber; z++) {
		delete[] ss[z].GCOORD_stru_gs;
	}
	if (element_type == 0) {
		for (z = 0; z < owsfnumber; z++) {
			for (i = 0; i < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); i++) {
				delete[] ol[z].GCOORD_flu_gs[i];
			}
			delete[] ol[z].GCOORD_flu_gs;
		}
	}
	if (element_type == 1) {
		for (z = 0; z < owsfnumber; z++) {
			for (i = 0; i < ol[z].FSNEL * 3; i++) {
				delete[] ol[z].GCOORD_flu_gs[i];
			}
			delete[] ol[z].GCOORD_flu_gs;
		}
	}
	for (z = 0; z < ssnumber; z++) {
		for (k = 0; k < NNODE_stru; k++) {
			delete[] ss[z].GCOORD_stru_fs[k];
		}
		delete[] ss[z].GCOORD_stru_fs;
	}
	for (z = 0; z < ssnumber; z++) {
		for (i = 0; i < ss[z].gs_num; i++) {
			delete[] ss[z].gs_flu_global[i];
		}
		delete[] ss[z].gs_flu_global;
	}
	for (z = 0; z < owsfnumber; z++) {
		for (i = 0; i < ol[z].FSNEL*(hprefg_flu + 1)*(hprefg_flu + 1); i++) {
			delete[] ol[z].flu_stru_global[i];
		}
		delete[] ol[z].flu_stru_global;
	}
	for (z = 0; z < ssnumber; z++) {
		for (i = 0; i < ss[z].ELE_stru * (hprefg + 1) * (hprefg + 1); i++) {
			delete[] ss[z].FP_flu[i];
		}
		delete[] ss[z].FP_flu;
	}
	std::cout << std::endl; 

	return;
}