//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <vector>
#include <sstream>
#include <cstring>

//This function serves to associate the fluid interface mesh with structure interface mesh. The relationship is used for interface mapping later on
void Neighborhood_search(double** GCOORD, int***LNA, int**IEN_flu, int NEL_flu) {
	extern OWETSURF ol[owsfnumber];
	extern stru_wet_surf ss; //data structure used to store the properties on the structure wetted surface
	int i, j, k, l, m, n, o, z;

	/*
	Eigen::Matrix3d A;
	A(0, 0) = 1; A(0, 1) = 2; A(0, 2) = 3;
	A(1, 0) = 4; A(1, 1) = 6; A(1, 2) = 2; 
	A(2, 0) = 1; A(2, 1) = 9; A(2, 2) = 0; 
	std::cout << A.inverse() << std::endl; 
	*/

	int elenode2D;
	if (element_type == 0) { //hex element
		elenode2D = NINT*NINT;
	}
	if (element_type == 1) { //tet element
		elenode2D = 3;
	}

	//Get the corner point of the fluid elements 
	int** IEN_flu_3D; 
	//int** IEN_flu_2D; 
	if (element_type == 0) { //hex element
		IEN_flu_3D = new int*[8];
		for (i = 0; i < 8; i++) {
			IEN_flu_3D[i] = new int[NEL_flu];
		}
		for (i = 0; i < NEL_flu; i++) {
			IEN_flu_3D[0][i] = IEN_flu[LNA[0][0][0] - 1][i];
			IEN_flu_3D[1][i] = IEN_flu[LNA[N][0][0] - 1][i];
			IEN_flu_3D[2][i] = IEN_flu[LNA[N][N][0] - 1][i];
			IEN_flu_3D[3][i] = IEN_flu[LNA[0][N][0] - 1][i];
			IEN_flu_3D[4][i] = IEN_flu[LNA[0][0][N] - 1][i];
			IEN_flu_3D[5][i] = IEN_flu[LNA[N][0][N] - 1][i];
			IEN_flu_3D[6][i] = IEN_flu[LNA[N][N][N] - 1][i];
			IEN_flu_3D[7][i] = IEN_flu[LNA[0][N][N] - 1][i];
		}

		for (z = 0; z < owsfnumber; z++) {
			ol[z].IEN_flu_2D = new int*[4];
			for (i = 0; i < 4; i++) {
				ol[z].IEN_flu_2D[i] = new int[ol[z].FSNEL];
			}
			for (i = 0; i < ol[z].FSNEL; i++) {
				for (j = 0; j < 4; j++) {
					ol[z].IEN_flu_2D[j][i] = IEN_flu[ol[z].LNA_norm[j] - 1][ol[z].GIDF[i] - 1];
				}
			}
		}
		
	}
	if (element_type == 1) {
		IEN_flu_3D = new int*[4];
		for (i = 0; i < 4; i++) {
			IEN_flu_3D[i] = new int[NEL_flu];
		}
		for (i = 0; i < NEL_flu; i++) {
			IEN_flu_3D[0][i] = IEN_flu[0][i];
			IEN_flu_3D[1][i] = IEN_flu[1][i];
			IEN_flu_3D[2][i] = IEN_flu[2][i];
			IEN_flu_3D[3][i] = IEN_flu[3][i];
		}
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

	//First bring in the structural wetted surface mesh into the code
	std::ifstream infile("C:/Users/lzhaok6/Desktop/FSP_canopy_0.3_abaqus_MpCCI_explicit_sym_0.3048wl.inp"); //The Abaqus input file
	if (!infile) {
		std::cout << "can not open the mesh file" << std::endl;
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
	while (getline(infile, csvLine))
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
			nodeend = ct - 1;
			elestart = ct + 1; //the line where element definition starts
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
	int NEL_stru = eleend - elestart + 1; //total element number in structure model
	//Define connectivity matrix in the volumn mesh
	int **IEN = new int*[4]; //4-node shell element
	for (i = 0; i < 4; i++) {
		IEN[i] = new int[NEL_stru];
	}
	for (i = elestart; i < eleend + 1; i++) {
		for (j = 0; j < 4; j++) {
			IEN[j][i - elestart] = stoi(output[i][j + 1]) - (node_start_num - 1); //change the string to integer
		}
	}
	double **GCOORD_stru; //the coordinate of structure nodes
	GCOORD_stru = new double*[NNODE_stru];
	for (i = 0; i < NNODE_stru; i++) {
		GCOORD_stru[i] = new double[3];
	}

	std::cout << "Remember to check the offset" << std::endl;
	system("PAUSE ");
	double xoffset = 0.0;
	double yoffset = -0.3048;
	double zoffset = 0.0;
	for (i = nodestart; i < nodeend + 1; i++) {
		GCOORD_stru[i - nodestart][0] = stod(output[i][1]) + xoffset;
		GCOORD_stru[i - nodestart][1] = stod(output[i][2]) + yoffset;
		GCOORD_stru[i - nodestart][2] = stod(output[i][3]) + zoffset;
	}
	 
	//Find the element set name corresponding to the structural wetted surface
	std::vector<std::string> surface_name; 
	std::vector<int> surface_orientation; //If the surface is positive, the value is 1; If the value is negative, the value is -1
	for (i = 0; i < surface.size(); i++) {
		surface_name.push_back(output[surface[i] + 1][0]); //http://www.cplusplus.com/reference/string/string/substr/
		if (output[surface[i] + 1][1] == " SNEG") { //negative surface
			//surface_orientation.push_back(-1);
			surface_orientation.push_back(1);
		}
		else if (output[surface[i] + 1][1] == " SPOS") { //positive surface 
			surface_orientation.push_back(1);
		}
		else {
			std::cout << "The surface side is not recognized" << std::endl;
			system("PAUSE ");
		}
	}
	//store the element number on the wetted surface
	std::vector<std::vector<int>> ele_num;
	std::vector<int> ele_num_clm; 
	ct = 0; 
	for (i = 0; i < sidesets_start.size(); i++) { //loop through side sets to find the structural wetted surface element set
		for (j = 0; j < surface.size(); j++) { //If the name of this sideset is the same with one of the wetted surfaces
			if (sidesets_name[i].substr(7, std::string::npos) == surface_name[j]) {
				ele_num_clm.clear(); 
				//read all the nodes in this sideset
				if (i < sidesets_start.size() - 1) {
					for (k = sidesets_start[i]; k < sidesets_start[i + 1] - 1; k++) {
						for (l = 0; l < output[k].size(); l++) {
							ele_num_clm.push_back(stoi(output[k][l]));
						}
					}
				}
				else { //for the last sideset, there is no sidesets_start[i + 1], we assume the surface definition follows the last side set
					for (k = sidesets_start[i]; k < surface[0]; k++) {
						for (l = 0; l < output[k].size(); l++) {
							ele_num_clm.push_back(stoi(output[k][l]));
						}
					}
				}
				ele_num.push_back(ele_num_clm);
				std::cout << " " << std::endl;
			}
		}
	}
	//Store the element connectivity 
	//get the total number of structural wetted surface elements
	int ELE_stru = 0; //total number of structural wetted surface elements
	for (i = 0; i < surface.size(); i++) {
		ELE_stru += ele_num[i].size();
	}
	int** IEN_stru; 
	IEN_stru = new int*[4];
	for (i = 0; i < 4; i++) {
		IEN_stru[i] = new int[ELE_stru];
	}
	ct = 0; 
	for (i = 0; i < surface.size(); i++) { //loop through structural wetted surfaces
		for (j = 0; j < ele_num[i].size(); j++) { //loop through all the element in the wetted 
			if (surface_orientation[i] == 1) { //positive surface, keep the original node sequence
				IEN_stru[0][ct] = IEN[0][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[1][ct] = IEN[1][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[2][ct] = IEN[2][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[3][ct] = IEN[3][ele_num[i][j] - 1 - (ele_start_num - 1)];
				ct += 1;
			}
			else { //revert the node sequence
				IEN_stru[0][ct] = IEN[3][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[1][ct] = IEN[2][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[2][ct] = IEN[1][ele_num[i][j] - 1 - (ele_start_num - 1)];
				IEN_stru[3][ct] = IEN[0][ele_num[i][j] - 1 - (ele_start_num - 1)];
				ct += 1;
			}
		}
	}
	if (ct != ELE_stru) {
		std::cout << "The wetted surface elements are not fully defined" << std::endl;
		system("PAUSE ");
	}
	//Interpolate the GCOORD_stru to obtain the GCOORD_stru_gs
	double **GCOORD_stru_gs; //the coordinate of structure wetted surface gauss points
	GCOORD_stru_gs = new double*[ELE_stru*(hprefg + 1)*(hprefg + 1)];
	for (i = 0; i < ELE_stru*(hprefg + 1)*(hprefg + 1); i++) {
		GCOORD_stru_gs[i] = new double[3];
	}
	//Define int**IEN_stru_gs for the later nodal force integration
	int** IEN_stru_gs; 
	IEN_stru_gs = new int*[(hprefg + 1)*(hprefg + 1)];
	for (i = 0; i < (hprefg + 1)*(hprefg + 1); i++) {
		IEN_stru_gs[i] = new int[ELE_stru];
	}

	ct = 1; 
	//Since we will use Gauss-Legendre nodes, different elements don't share nodes. As a result, the nodes could be numbered sequentially
	for (i = 0; i < ELE_stru; i++) {
		for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {
			IEN_stru_gs[j][i] = ct;
			ct += 1;
		}
	}

	//Derive the coordinate of the gauss nodes 
	TD_LOCAL_NODEstruct ct_gs;
	ct_gs = TD_LOCAL_NODE(hprefg);
	TD_LOCAL_NODEstruct ct_ln;
	ct_ln = TD_LOCAL_NODE(1);
	LOCAL_SHAPEstruct ls_ln; //ln means linear
	ls_ln = LOCAL_SHAPE(LNA, 1, hprefg, 1); //Get the 2D linear shape function value on Nth order Gauss-Legendre nodes
	for (i = 0; i < ELE_stru; i++) { //loop through all the elements
		for (m = 0; m < hprefg + 1; m++) { //m, n, o stands for internal nodes in one element (all but except for corner nodes)
			for (n = 0; n < hprefg + 1; n++) {
				GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][0] = 0.0;
				GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][1] = 0.0;
				GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][2] = 0.0;
				for (j = 0; j < 2; j++) { //j, k, l stands for corner nodes
					for (k = 0; k < 2; k++) {
						GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][0] += GCOORD_stru[IEN_stru[ct_ln.LNA[j][k] - 1][i] - 1][0] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
						GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][1] += GCOORD_stru[IEN_stru[ct_ln.LNA[j][k] - 1][i] - 1][1] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
						GCOORD_stru_gs[IEN_stru_gs[ct_gs.LNA[m][n] - 1][i] - 1][2] += GCOORD_stru[IEN_stru[ct_ln.LNA[j][k] - 1][i] - 1][2] * ls_ln.SHL_2D[2][ct_ln.LNA[j][k] - 1][m*(hprefg + 1) + n];
					}
				}
			}
		}
	}

	int gs_num = ELE_stru * (hprefg + 1) * (hprefg + 1); //the total number of gauss points on the structural surface

	ss.gs_flu_global = new double*[gs_num]; //used to store the global coordinate of the projected gauss point on fluid surface
	for (i = 0; i < gs_num; i++) {
		ss.gs_flu_global[i] = new double[3];
	}
	ss.gs_flu_local = new double*[gs_num]; 
	for (i = 0; i < gs_num; i++) {
		ss.gs_flu_local[i] = new double[3];
	}
	ss.gs_flu = new int[gs_num];

	for (z = 0; z < owsfnumber; z++) {
		ol[z].flu_local = new double*[ol[z].FSNEL*elenode2D];
		ol[z].flu_stru_global = new double*[ol[z].FSNEL*elenode2D];
		ol[z].flu_stru = new int[ol[z].FSNEL*elenode2D];
		for (i = 0; i < ol[z].FSNEL*elenode2D; i++) {
			ol[z].flu_local[i] = new double[2];
			ol[z].flu_stru_global[i] = new double[3];
		}
	}

	int gs_nearest; //used to store the nearest fluid node on the fluid FSI interface
	int flu_nearest; //used to store the nearest structure node on the structural wetted surface
	//Container to store the searching result (one gauss point corresponds to a fluid element)
	int** FSNEL_stru_gs; //[gauss point][structural wetted surface element]

	Eigen::Matrix3d fdot; 
	Eigen::Vector3d f;
	Eigen::Vector3d dx;
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
	for (i = 0; i < ELE_stru; i++) {
		for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {
			orphan = 1; //let's first assume the node is an orphan. If the corresponding element is found, the flag is turned to 0 
			//Start searching for fluid element for gauss node IEN_gs[j][i]
			//Calculate the distance from the current the structural gauss point to the fluid element in the search range
			distance[0] = search_range;
			range[1] = search_range; //initiate the fluid node distance to be equal to the searching range for every gauss node
			for (z = 0; z < owsfnumber; z++) { //search the fluid wetted element group by group
				for (k = 0; k < ol[z].FSNEL; k++) { //looping through fluid wetted elements
					flag = 1;
					//First determine if the fluid element is inside the searching range
					for (l = 0; l < elenode2D; l++) {
						range[0] = pow(pow(GCOORD[ol[z].IEN_flu_2D[l][k] - 1][0] - GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][0], 2) + pow(GCOORD[ol[z].IEN_flu_2D[l][k] - 1][1] - GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][1], 2) + pow(GCOORD[ol[z].IEN_flu_2D[l][k] - 1][2] - GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][2], 2), 0.5);
						if (range[0] < range[1]) {
							range[1] = range[0]; //range[1] is used to store the shortest distance so far
							gs_nearest = ol[z].IEN_flu_2D[l][k];
						}
						//At the same time, store the node with the shorted distance to the gauss node under searching. If the node is orphan, we use the value on this fluid node for the structure gauss node
						if (range[0] > search_range) { //One of the corner nodes in the element is out of the searching range. Jump through this element
							flag = 0;
						}
					}
					if (flag == 1) { //the element should be searched (within the searching range)
						//Find the global coordinate of the orthogonally projected point on the current element
						//initiate the coordinate of the projected point to be the same with the structural gauss point
						xu_k = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][0]; yu_k = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][1]; zu_k = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][2];
						x1 = GCOORD[ol[z].IEN_flu_2D[0][k] - 1][0]; x2 = GCOORD[ol[z].IEN_flu_2D[1][k] - 1][0]; x3 = GCOORD[ol[z].IEN_flu_2D[2][k] - 1][0];
						y1 = GCOORD[ol[z].IEN_flu_2D[0][k] - 1][1]; y2 = GCOORD[ol[z].IEN_flu_2D[1][k] - 1][1]; y3 = GCOORD[ol[z].IEN_flu_2D[2][k] - 1][1];
						z1 = GCOORD[ol[z].IEN_flu_2D[0][k] - 1][2]; z2 = GCOORD[ol[z].IEN_flu_2D[1][k] - 1][2]; z3 = GCOORD[ol[z].IEN_flu_2D[2][k] - 1][2];
						xn = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][0]; yn = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][1]; zn = GCOORD_stru_gs[IEN_stru_gs[j][i] - 1][2];
						diff = 1; 
						while (diff > tol) {
							fdot(0, 0) = -x1 - xn + 2 * xu_k; fdot(0, 1) = -y1 - yn + 2 * yu_k; fdot(0, 2) = -z1 - zn + 2 * zu_k; //1st point
							fdot(1, 0) = -x2 - xn + 2 * xu_k; fdot(1, 1) = -y2 - yn + 2 * yu_k; fdot(1, 2) = -z2 - zn + 2 * zu_k; //2nd point
							fdot(2, 0) = -x3 - xn + 2 * xu_k; fdot(2, 1) = -y3 - yn + 2 * yu_k; fdot(2, 2) = -z3 - zn + 2 * zu_k; //3rd point
							f(0) = (xu_k - x1)*(xu_k - xn) + (yu_k - y1)*(yu_k - yn) + (zu_k - z1)*(zu_k - zn);
							f(1) = (xu_k - x2)*(xu_k - xn) + (yu_k - y2)*(yu_k - yn) + (zu_k - z2)*(zu_k - zn);
							f(2) = (xu_k - x3)*(xu_k - xn) + (yu_k - y3)*(yu_k - yn) + (zu_k - z3)*(zu_k - zn);
							//std::cout << fdot.inverse() << std::endl; 
							dx = -fdot.inverse()*f; //delta x^k 
							//diff = pow(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2), 0.5);
							diff = pow(pow(dx(0), 2) + pow(dx(1), 2) + pow(dx(2), 2), 0.5);
							//xu_k_1 = xu_k + dx(0, 0); yu_k_1 = yu_k + dx(1, 0); zu_k_1 = zu_k + dx(2, 0);
							xu_k_1 = xu_k + dx(0); yu_k_1 = yu_k + dx(1); zu_k_1 = zu_k + dx(2);
							xu_k = xu_k_1; yu_k = yu_k_1; zu_k = zu_k_1; //update the value
						}
						xu = xu_k; yu = yu_k; zu = zu_k;
						//Found the projected point after Newton iteration 

						//We need to determine at this stage if the projected node is within the searched element, if not we need to pass this element. 
						if (element_type == 1) { //tetrahedral element
							Eigen::Matrix3d b(3, 1);
							Eigen::Matrix3d A(3, 3);
							Eigen::Matrix3d local(3, 1); //store the local coordinate 
							//determine the local coordinate in the tetrahedral element
							A(0, 0) = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][0] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][0];
							A(0, 1) = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][0] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][0];
							A(0, 2) = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][0] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][0];
							A(1, 0) = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][1] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][1];
							A(1, 1) = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][1] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][1];
							A(1, 2) = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][1] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][1];
							A(2, 0) = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][2] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][2];
							A(2, 1) = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][2] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][2];
							A(2, 2) = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][2] - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][2];
							b(0, 0) = xu - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][0];
							b(1, 0) = yu - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][1];
							b(2, 0) = zu - GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][2];
							local = A.inverse()*b;
							//determine if the point is inside the searched element
							if (local(0, 0) >= 0 && local(0, 0) <= 1 && local(1, 0) >= 0 && local(1, 0) <= 1 && local(2, 0) >= 0 && local(2, 0) <= 1
								&& local(0, 0) + local(1, 0) + local(2, 0) <= 1) {
								//The point is indeed inside the searched element and we can store the local coordinate in that element
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][0] = local(0, 0); //xi
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][1] = local(1, 0); //eta
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][2] = local(2, 0); //zeta
								orphan = 0;
							}
							else { //The node is an orphan
								   //Search for the nearest point on the fluid FSI interface
								orphan = 1;
							}
						}
						double xi, eta, zeta; //converged local coordinate
						if (element_type == 0) { //hex element
							double xi_k, eta_k, zeta_k;
							double xi_k_1, eta_k_1, zeta_k_1;
							//The derivative is obtained by Mathematica: /Users/zhaokuanlu/OneDrive/post_processing/Hex_global_to_local_coordinate_Newton_iteration.nb
							//determine the local coordinate in the hexahedral element
							//Since the equation system is nonlinear, we use the Newton iteration and set the initial value to be 0 
							//xi_k = 0.5; eta_k = 0.5; zeta_k = 0.5;
							xi_k = 0.1; eta_k = 0.1; zeta_k = 0.1;
							x1 = GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][0]; x2 = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][0]; x3 = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][0]; x4 = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][0]; x5 = GCOORD[IEN_flu_3D[4][ol[z].GIDF[k] - 1] - 1][0]; x6 = GCOORD[IEN_flu_3D[5][ol[z].GIDF[k] - 1] - 1][0]; x7 = GCOORD[IEN_flu_3D[6][ol[z].GIDF[k] - 1] - 1][0]; x8 = GCOORD[IEN_flu_3D[7][ol[z].GIDF[k] - 1] - 1][0];
							y1 = GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][1]; y2 = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][1]; y3 = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][1]; y4 = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][1]; y5 = GCOORD[IEN_flu_3D[4][ol[z].GIDF[k] - 1] - 1][1]; y6 = GCOORD[IEN_flu_3D[5][ol[z].GIDF[k] - 1] - 1][1]; y7 = GCOORD[IEN_flu_3D[6][ol[z].GIDF[k] - 1] - 1][1]; y8 = GCOORD[IEN_flu_3D[7][ol[z].GIDF[k] - 1] - 1][1];
							z1 = GCOORD[IEN_flu_3D[0][ol[z].GIDF[k] - 1] - 1][2]; z2 = GCOORD[IEN_flu_3D[1][ol[z].GIDF[k] - 1] - 1][2]; z3 = GCOORD[IEN_flu_3D[2][ol[z].GIDF[k] - 1] - 1][2]; z4 = GCOORD[IEN_flu_3D[3][ol[z].GIDF[k] - 1] - 1][2]; z5 = GCOORD[IEN_flu_3D[4][ol[z].GIDF[k] - 1] - 1][2]; z6 = GCOORD[IEN_flu_3D[5][ol[z].GIDF[k] - 1] - 1][2]; z7 = GCOORD[IEN_flu_3D[6][ol[z].GIDF[k] - 1] - 1][2]; z8 = GCOORD[IEN_flu_3D[7][ol[z].GIDF[k] - 1] - 1][2];
							diff = 1;
							while (diff > tol) {
								fdot(0, 0) = -(1.0 / 8.0) *(1.0 - eta_k) *x1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 + zeta_k);
								fdot(0, 1) = -(1.0 / 8.0) *x1*(1.0 - xi_k) *(1.0 - zeta_k) + 1.0 / 8.0 *x4*(1.0 - xi_k) *(1.0 - zeta_k) -
									1.0 / 8.0 *x2*(1.0 + xi_k) *(1.0 - zeta_k) + 1.0 / 8.0 *x3*(1.0 + xi_k) *(1.0 - zeta_k) -
									1.0 / 8.0 *x5*(1.0 - xi_k) *(1.0 + zeta_k) + 1.0 / 8.0 *x8*(1.0 - xi_k) *(1.0 + zeta_k) -
									1.0 / 8.0 *x6*(1.0 + xi_k) *(1.0 + zeta_k) + 1.0 / 8.0 *x7*(1.0 + xi_k) *(1.0 + zeta_k);
								fdot(0, 2) = -(1.0 / 8.0) *(1.0 - eta_k) *x1*(1.0 - xi_k) - 1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - xi_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 - xi_k) + 1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 - xi_k) -
									1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 + xi_k) - 1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 + xi_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + xi_k) + 1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + xi_k);
								fdot(1, 0) = -(1.0 / 8.0) *(1.0 - eta_k) *y1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *y3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *y4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *y5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *y7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *y8*(1.0 + zeta_k);
								fdot(1, 1) = -(1.0 / 8.0) *(1.0 - xi_k) *y1*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *y3*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *y4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - xi_k) *y5*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *y7*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *y8*(1.0 + zeta_k);
								fdot(1, 2) = -(1.0 / 8.0) *(1.0 - eta_k) *(1.0 - xi_k) *y1 - 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y2 -
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y3 - 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y4 +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y5 + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y6 +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y7 + 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y8;
								fdot(2, 0) = -(1.0 / 8.0) *(1.0 - eta_k) *z1*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *z3*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *z4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - eta_k) *z5*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - eta_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *z7*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + eta_k) *z8*(1.0 + zeta_k);
								fdot(2, 1) = -(1.0 / 8.0) *(1.0 - xi_k) *z1*(1.0 - zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *z3*(1.0 - zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *z4*(1.0 - zeta_k) -
									1.0 / 8.0 *(1.0 - xi_k) *z5*(1.0 + zeta_k) - 1.0 / 8.0 *(1.0 + xi_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + xi_k) *z7*(1.0 + zeta_k) + 1.0 / 8.0 *(1.0 - xi_k) *z8*(1.0 + zeta_k);
								fdot(2, 2) = -(1.0 / 8.0) *(1.0 - eta_k) *(1.0 - xi_k) *z1 - 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z2 -
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z3 - 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z4 +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z5 + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z6 +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z7 + 1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z8;
								f(0) = -xu + 1.0 / 8.0 *(1.0 - eta_k) *x1*(1.0 - xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x4*(1.0 - xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x2*(1.0 + xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x3*(1.0 + xi_k) *(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x5*(1.0 - xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x8*(1.0 - xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *x6*(1.0 + xi_k) *(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *x7*(1.0 + xi_k) *(1.0 + zeta_k);
								f(1) = -yu + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y1*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y3*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y4*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *y5*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *y6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *y7*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *y8*(1.0 + zeta_k);
								f(2) = -zu + 1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z1*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z2*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z3*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z4*(1.0 - zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 - xi_k) *z5*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 - eta_k) *(1.0 + xi_k) *z6*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 + xi_k) *z7*(1.0 + zeta_k) +
									1.0 / 8.0 *(1.0 + eta_k) *(1.0 - xi_k) *z8*(1.0 + zeta_k);
								dx = -fdot.inverse()*f; //delta x^k
								diff = pow(pow(dx(0), 2) + pow(dx(1), 2) + pow(dx(2), 2), 0.5);
								xi_k_1 = xi_k + dx(0); eta_k_1 = eta_k + dx(1); zeta_k_1 = zeta_k + dx(2);
								xi_k = xi_k_1; eta_k = eta_k_1; zeta_k = zeta_k_1; //update the value
							}
							xi = xi_k; eta = eta_k; zeta = zeta_k;
							if (xi + 1 > -1e-5 && xi - 1 < 1e-5 && eta + 1 > -1e-5 && eta - 1 < 1e-5 && zeta + 1 > -1e-5 && zeta - 1 < 1e-5) {
								//The point is indeed inside the searched element. We can turn off the orphan flag
								orphan = 0;
								inelement = 1;
							}
							else { //The node is not in the current searched element (might be an orphan)
								inelement = 0; 
							}
						}
						//End determining if the projected node is in the searched element and its local coordinate

						//If the projected node is indeed within the searched element, calculate the orthogonal distance 
						//Calculate the distance between the structural gauss nodes and the projected node
						if (inelement == 1) { //if the node is within the current fluid element
							distance[1] = pow(pow(xn - xu, 2) + pow(yn - yu, 2) + pow(zn - zu, 2), 0.5);
							if (distance[1] < distance[0]) {
								//store the local coordinate in the fluid mesh of the projected structure gauss node
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][0] = xi; //xi
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][1] = eta; //eta
								ss.gs_flu_local[i*(hprefg + 1)*(hprefg + 1) + j][2] = zeta; //zeta
								//store the global coordinate in the fluid mesh of the projected structure gauss node
								xu_searched = xu; yu_searched = yu; zu_searched = zu;
								distance[0] = distance[1]; //update the shortest orthogonal distance so far
								searched_ele = ol[z].GIDF[k];
							}
						}
					}
					//std::cout << " " << std::endl;
				} //finish searching all the fluid elements in the current group
			} //finish the searching through the fluid interface elements (several groups)
			if (orphan == 1) { //If the node is still an orphan after looping through all the fluid elements, we need to search the fluid mesh to determine which is the closest node
				//the orphan node and the other nodes should be treated separately
				//Instead of storing the corresponding element in gs_flu, the stored value becomes the fluid mesh "node" number 
				//At the same time, we keep track of which gauss node is orphan
				ss.gs_flu[i*(hprefg + 1)*(hprefg + 1) + j] = gs_nearest;
				ss.orphan_num_gs.push_back(i*(hprefg + 1)*(hprefg + 1) + j);
			}
			else {
				//the corresponding fluid element and the projected node coordinate should be found at this stage
				//store the corresponding element 
				ss.gs_flu[i*(hprefg + 1)*(hprefg + 1) + j] = searched_ele;
			}
			//store the projected node location
			ss.gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][0] = xu_searched;
			ss.gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][1] = yu_searched;
			ss.gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][2] = zu_searched;
		}
	}

	//Node projection from fluid to structure (Project the fluid interpolation points to structural elements)
	Eigen::Matrix2d fdot_f;
	Eigen::Vector2d f_f;
	Eigen::Vector2d dx_f;
	diff = 1.0;
	range[1] = search_range;
	distance[0] = search_range;
	for (z = 0; z < owsfnumber; z++) {
		for (i = 0; i < ol[z].FSNEL; i++) { //loop through fluid elements
			for (j = 0; j < elenode2D; j++) { //loop through each node on the fluid FSI 
				orphan = 1; //let's first assume the node is an orphan. If the corresponding element is found, the flag is turned to 0 
				distance[0] = search_range;
				range[1] = search_range;
				//Start searching for structural element for fluid interpolation node ol[z].IEN_gb[j][i]
				//Calculate the distance from the current fluid interpolation node to the structural element in the search range
				//First determine if the structure element is inside the searching range
				for (k = 0; k < ELE_stru; k++) { //looping through structure wetted elements
					flag = 1;
					for (l = 0; l < 4; l++) {
						range[0] = pow(pow(GCOORD[ol[z].IEN_gb[j][i] - 1][0] - GCOORD_stru[IEN_stru[l][k] - 1][0], 2) + pow(GCOORD[ol[z].IEN_gb[j][i] - 1][1] - GCOORD_stru[IEN_stru[l][k] - 1][1], 2) + pow(GCOORD[ol[z].IEN_gb[j][i] - 1][2] - GCOORD_stru[IEN_stru[l][k] - 1][2], 2), 0.5);
						if (range[0] < range[1]) {
							range[1] = range[0]; //range[1] is used to store the shortest distance so far
							flu_nearest = IEN_stru[l][k];
						}
						//At the same time, store the node with the shorted distance to the gauss node under searching. If the node is orphan, we use the value on this fluid node for the structure gauss node
						if (range[0] > search_range) { //One of the corner nodes in the element is out of the searching range. Jump through this element
							flag = 0;
						}
					}
					if (flag == 1) { //the element should be searched (within the searching range)
						//Find the coordinate of the orthogonally projected fluid interpolation point on the current structure element
						//initiate the coordinate of the projected point to be the same with the fluid interpolation point
						//xu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][0]; yu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][1]; zu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][2];
						xu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][0] + 0.1; yu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][1] + 1; zu_k = GCOORD[ol[z].IEN_gb[j][i] - 1][2] + 1;
						x1 = GCOORD_stru[IEN_stru[0][k] - 1][0]; x2 = GCOORD_stru[IEN_stru[1][k] - 1][0]; x3 = GCOORD_stru[IEN_stru[2][k] - 1][0];
						y1 = GCOORD_stru[IEN_stru[0][k] - 1][1]; y2 = GCOORD_stru[IEN_stru[1][k] - 1][1]; y3 = GCOORD_stru[IEN_stru[2][k] - 1][1];
						z1 = GCOORD_stru[IEN_stru[0][k] - 1][2]; z2 = GCOORD_stru[IEN_stru[1][k] - 1][2]; z3 = GCOORD_stru[IEN_stru[2][k] - 1][2];
						xn = GCOORD[ol[z].IEN_gb[j][i] - 1][0]; yn = GCOORD[ol[z].IEN_gb[j][i] - 1][1]; zn = GCOORD[ol[z].IEN_gb[j][i] - 1][2];
						diff = 1.0;
						while (diff > tol) {
							fdot(0, 0) = -x1 - xn + 2 * xu_k; fdot(0, 1) = -y1 - yn + 2 * yu_k; fdot(0, 2) = -z1 - zn + 2 * zu_k; //1st point
							fdot(1, 0) = -x2 - xn + 2 * xu_k; fdot(1, 1) = -y2 - yn + 2 * yu_k; fdot(1, 2) = -z2 - zn + 2 * zu_k; //2nd point
							fdot(2, 0) = -x3 - xn + 2 * xu_k; fdot(2, 1) = -y3 - yn + 2 * yu_k; fdot(2, 2) = -z3 - zn + 2 * zu_k; //3rd point
							f(0) = (xu_k - x1)*(xu_k - xn) + (yu_k - y1)*(yu_k - yn) + (zu_k - z1)*(zu_k - zn);
							f(1) = (xu_k - x2)*(xu_k - xn) + (yu_k - y2)*(yu_k - yn) + (zu_k - z2)*(zu_k - zn);
							f(2) = (xu_k - x3)*(xu_k - xn) + (yu_k - y3)*(yu_k - yn) + (zu_k - z3)*(zu_k - zn);
							dx = -fdot.inverse()*f; //delta x^k
							diff = pow(pow(dx(0), 2) + pow(dx(1), 2) + pow(dx(2), 2), 0.5);
							xu_k_1 = xu_k + dx(0); yu_k_1 = yu_k + dx(1); zu_k_1 = zu_k + dx(2);
							xu_k = xu_k_1; yu_k = yu_k_1; zu_k = zu_k_1; //update the value
						}
						xu = xu_k; yu = yu_k; zu = zu_k;
						//Found the projected point after Newton iteration 

						//We need to determine at this stage if the projected node is within the searched element, if not we need to pass this element. 
						double xi_k, eta_k;
						double xi_k_1, eta_k_1;
						double xi, eta; //converged local coordinate
						double xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4; //pseudo structural node coordinate 
						double xup, yup; //pseudo projected point coordinate
						//The derivative is obtained by Mathematica: /Users/zhaokuanlu/OneDrive/post_processing/.nb
						//determine the local coordinate in the hexahedral element
						//Since the equation system is nonlinear, we use the Newton iteration and set the initial value to be 0 
						xi_k = 0.5; eta_k = 0.5;
						x1 = GCOORD_stru[IEN_stru[0][k] - 1][0]; x2 = GCOORD_stru[IEN_stru[1][k] - 1][0]; x3 = GCOORD_stru[IEN_stru[2][k] - 1][0]; x4 = GCOORD_stru[IEN_stru[3][k] - 1][0];
						y1 = GCOORD_stru[IEN_stru[0][k] - 1][1]; y2 = GCOORD_stru[IEN_stru[1][k] - 1][1]; y3 = GCOORD_stru[IEN_stru[2][k] - 1][1]; y4 = GCOORD_stru[IEN_stru[3][k] - 1][1];
						z1 = GCOORD_stru[IEN_stru[0][k] - 1][2]; z2 = GCOORD_stru[IEN_stru[1][k] - 1][2]; z3 = GCOORD_stru[IEN_stru[2][k] - 1][2]; z4 = GCOORD_stru[IEN_stru[3][k] - 1][2];
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
						else { //the element is a 3D surface, use any two dimension (we choose to use x and y here)
							xp1 = x1; xp2 = x2; xp3 = x3; xp4 = x4;
							yp1 = y1; yp2 = y2; yp3 = y3; yp4 = y4;
							xup = xu; yup = yu; 
						}
						diff = 1.0;
						while (diff > tol) {
							fdot_f(0, 0) = -(1.0 / 4.0) *(1.0 - eta_k) *xp1 + 1.0 / 4.0 *(1.0 - eta_k) *xp2 + 1.0 / 4.0 *(1.0 + eta_k) *xp3 -
								1.0 / 4.0 *(1.0 + eta_k) *xp4;
							fdot_f(0, 1) = -(1.0 / 4.0) *xp1*(1.0 - xi_k) + 1.0 / 4.0 *xp4*(1.0 - xi_k) - 1.0 / 4.0 *xp2*(1.0 + xi_k) +
								1.0 / 4.0 *xp3*(1.0 + xi_k);
							fdot_f(1, 0) = -(1.0 / 4.0) *(1.0 - eta_k) *yp1 + 1.0 / 4.0 *(1.0 - eta_k) *yp2 + 1.0 / 4.0 *(1.0 + eta_k) *yp3 -
								1.0 / 4.0 *(1.0 + eta_k) *yp4;
							fdot_f(1, 1) = -(1.0 / 4.0) *(1.0 - xi_k) *yp1 - 1.0 / 4.0 *(1.0 + xi_k) *yp2 + 1.0 / 4.0 *(1.0 + xi_k) *yp3 +
								1.0 / 4.0 *(1.0 - xi_k) *yp4;
							f_f(0) = -xup + 1.0 / 4.0 *(1.0 - eta_k) *xp1*(1.0 - xi_k) + 1.0 / 4.0 *(1.0 + eta_k) *xp4*(1.0 - xi_k) +
								1.0 / 4.0 *(1.0 - eta_k) *xp2*(1.0 + xi_k) + 1.0 / 4.0 *(1.0 + eta_k) *xp3*(1.0 + xi_k);
							f_f(1) = -yup + 1.0 / 4.0 *(1.0 - eta_k) *(1.0 - xi_k) *yp1 + 1.0 / 4.0 *(1.0 - eta_k) *(1.0 + xi_k) *yp2 +
								1.0 / 4.0 *(1.0 + eta_k) *(1.0 + xi_k) *yp3 + 1.0 / 4.0 *(1.0 + eta_k) *(1.0 - xi_k) *yp4;
							dx_f = -fdot_f.inverse()*f_f; //delta x^k
							diff = pow(pow(dx_f(0), 2) + pow(dx_f(1), 2), 0.5);
							xi_k_1 = xi_k + dx_f(0); eta_k_1 = eta_k + dx_f(1);
							xi_k = xi_k_1; eta_k = eta_k_1; //update the value
						}
						xi = xi_k; eta = eta_k;
						if (xi + 1 > -1e-5 && xi - 1 < 1e-5 && eta + 1 > -1e-5 && eta - 1 < 1e-5) {
							//The point is indeed inside the searched element and we can store the local coordinate in that element
							orphan = 0; //turn off the orphan flag
							inelement = 1;
						}
						else { //The node is an orphan
							//orphan = 1;
							inelement = 0;
						}
						//End determining if the projected node is in the searched element and its local coordinate

						//If the projected node is indeed within the searched element, calculate the orthogonal distance 
						//Calculate the distance between the fluid interpolation nodes and the projected node
						if (orphan == 0) { //if the node is not an orphan
							distance[1] = pow(pow(xn - xu, 2) + pow(yn - yu, 2) + pow(zn - zu, 2), 0.5);
							if (distance[1] < distance[0]) {
								xu_searched = xu; yu_searched = yu; zu_searched = zu;
								distance[0] = distance[1];
								ol[z].flu_local[i*elenode2D + j][0] = xi; //xi
								ol[z].flu_local[i*elenode2D + j][1] = eta; //eta
								searched_ele = k + 1;
							}
						}
					}
				}//finish the searching through the structural wetted elements 
				if (orphan == 1) { //If the node is still an orphan, we need to search the structure mesh to determine which is the closest node
					//the orphan node and the other nodes should be treated separately
					//Instead of storing the corresponding element in flu_stru, the stored value becomes the structure "node" number 
					//At the same time, we keep track of which fluid node is orphan
					ol[z].flu_stru[i*elenode2D + j] = flu_nearest;
					ol[z].orphan_num_flu.push_back(i*elenode2D + j);
				}
				//the corresponding structural element and the projected node coordinate should be found at this stage
				//store the corresponding element 
				ol[z].flu_stru[i*elenode2D + j] = searched_ele;
				//store the projected node location
				ol[z].flu_stru_global[i*elenode2D + j][0] = xu_searched;
				ol[z].flu_stru_global[i*elenode2D + j][1] = yu_searched;
				ol[z].flu_stru_global[i*elenode2D + j][2] = zu_searched;
			}
			//Finish with all the nodes in the current fluid element
		}
		//Finish with all the nodes in the current physical group
	}
	//Finish all the fluid nodes
	std::cout << " " << std::endl; 

	return;
}