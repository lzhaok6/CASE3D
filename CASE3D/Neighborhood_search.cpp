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

//This function serves to associate the fluid interface mesh with structure interface mesh. The relationship is used for interface mapping later on
void Neighborhood_search() {
	
	int i, j, k, l, m, n; 
	//First bring in the structural wetted surface mesh into the code
	
	//Define the connectivity matrix on structure wetted surface (imported)
	int** IEN_stru; 
	int** IEN_gs; //Gauss points in the wetted element
	//count how many gauss points are there
	int** IEN_gs_flu; //the corresponding fluid element number of the gauss point
	int ELE_stru; 
	int gs_num = ELE_stru * 4; //the total number of gauss points 
	int*gs_flu = new int[gs_num]; //used to store the corresponding fluid element
	double**gs_flu_global = new double*[gs_num]; //used to store the global coordinate of the projected point
	for (i = 0; i < gs_num; i++) {
		gs_flu_global[i] = new double[3];
	}
	double**gs_flu_local; //used to store the local coordinate of the projected point
	for (i = 0; i < gs_num; i++) {
		gs_flu_local[i] = new double[3];
	}
	double **GCOORD_stru; //the coordinate of structure wetted surface points
	//Container to store the searching result (one gauss point corresponds to a fluid element)
	int** FSNEL_stru_gs; //[gauss point][structural wetted surface element]

	//Define the connecvitiy matrix no fluid FSI surface (Should be available from meshgeneration.cpp)
	int** IEN_flu; //the corner nodes of the fluid interface element 
	int FSNEL; //the wetted surface element number on fluid surface
	double** GCOORD; //the coordinate in the fluid domain
	int elenode2D;

	Eigen::MatrixXd fdot(3, 3);
	Eigen::MatrixXd f(3, 1);
	Eigen::MatrixXd dx(3, 1);
	double search_range = 1; //set the searching range to 1m (only the element with all points in the range will be searched)
	double range = 0.0;
	double distance[2]; //the distance from the strutural gauss node to its projected point
	distance[0] = search_range; 
	int flag = 1; 
	double tol = 1e-6; 
	double diff = 0.0; 
	double xu = 0.0; double yu = 0.0; double zu = 0.0; //the cooridnate of the projected point
	double xu_searched = 0.0; double yu_searched = 0.0; double zu_searched = 0.0; //The projected node closest to the structure gauss node
	double xu_k = 0.0; double yu_k = 0.0; double zu_k = 0.0; //the projected point coordinate at step k 
	double xu_k_1 = 0.0; double yu_k_1 = 0.0; double zu_k_1 = 0.0; //the projected point coordinate at step k + 1
	double x1; double y1; double z1; double x2; double y2; double z2; double x3; double y3; double z3; 
	double xn; double yn; double zn; 
	int searched_ele; //the fluid element found after node projection
	//Node projection from structure to fluid (Project the structural Gauss points to fluid elements)
	//Structural wetted elements do not share gauss points
	for (i = 0; i < ELE_stru; i++) {
		for (j = 0; j < (hprefg + 1)*(hprefg + 1); j++) {


			//Start searching for fluid element for gauss node IEN_gs[j][i]
			//Calculate the distance from the current the structural gauss point to the fluid element in the search range
			//First determine if the element is inside the searching range
			for (k = 0; k < FSNEL; k++) { //looping through fluid wetted elements
				flag = 1;
				for (l = 0; l < elenode2D; l++) {
					range = pow(pow(GCOORD[IEN_flu[l][k] - 1][0] - GCOORD_stru[IEN_gs[j][i] - 1][0], 2) + pow(GCOORD[IEN_flu[l][k] - 1][1] - GCOORD_stru[IEN_gs[j][i] - 1][1], 2) + pow(GCOORD[IEN_flu[l][k] - 1][2] - GCOORD_stru[IEN_gs[j][i] - 1][2], 2), 0.5);
					if (range > 1) { //One of the corner nodes in the element is out of the searching range. Jump through this element
						flag = 0;
					}
				}
				if (flag == 1) { //the element should be searched (within the searching range)
					//Find the coordinate of the orthogonally projected point on the current element
					//initiate the coordinate of the projected point to be the same with the structural gauss point
					xu_k = GCOORD_stru[IEN_gs[j][i] - 1][0]; yu_k = GCOORD_stru[IEN_gs[j][i] - 1][1]; zu_k = GCOORD_stru[IEN_gs[j][i] - 1][2];
					x1 = GCOORD[IEN_flu[0][k] - 1][0]; x2 = GCOORD[IEN_flu[1][k] - 1][0]; x3 = GCOORD[IEN_flu[2][k] - 1][0];
					y1 = GCOORD[IEN_flu[0][k] - 1][1]; y2 = GCOORD[IEN_flu[1][k] - 1][1]; y3 = GCOORD[IEN_flu[2][k] - 1][1];
					z1 = GCOORD[IEN_flu[0][k] - 1][2]; z2 = GCOORD[IEN_flu[1][k] - 1][2]; z3 = GCOORD[IEN_flu[2][k] - 1][2];
					xn = GCOORD_stru[IEN_gs[j][i] - 1][0]; yn = GCOORD_stru[IEN_gs[j][i] - 1][1]; zn = GCOORD_stru[IEN_gs[j][i] - 1][2];
					while (diff > tol) {
						fdot(0, 0) = -x1 - xn + 2 * xu_k; fdot(0, 1) = -y1 - yn + 2 * yu_k; fdot(0, 2) = -z1 - zn + 2 * zu_k; //1st point
						fdot(1, 0) = -x2 - xn + 2 * xu_k; fdot(1, 1) = -y2 - yn + 2 * yu_k; fdot(1, 2) = -z2 - zn + 2 * zu_k; //2nd point
						fdot(2, 0) = -x3 - xn + 2 * xu_k; fdot(2, 1) = -y3 - yn + 2 * yu_k; fdot(2, 2) = -z3 - zn + 2 * zu_k; //3rd point
						f(0, 0) = (xu_k - x1)*(xu_k - xn) + (yu_k - y1)*(yu_k - yn) + (zu_k - z1)*(zu_k - zn);
						f(1, 0) = (xu_k - x2)*(xu_k - xn) + (yu_k - y2)*(yu_k - yn) + (zu_k - z2)*(zu_k - zn);
						f(2, 0) = (xu_k - x3)*(xu_k - xn) + (yu_k - y3)*(yu_k - yn) + (zu_k - z3)*(zu_k - zn);
						dx = -fdot.inverse()*f; //delta x^k
						diff = pow(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2), 0.5);
						xu_k_1 = xu_k + dx(0, 0); yu_k_1 = yu_k + dx(1, 0); zu_k_1 = zu_k + dx(2, 0);
						xu_k = xu_k_1; yu_k = yu_k_1; zu_k = zu_k_1; //update the value
					}
					xu = xu_k; yu = yu_k; zu = zu_k;
					//Find the projected point after Newton iteration 
					//Calculate the distance between the structural gauss nodes and the projected node
					distance[1] = pow(pow(xn - xu, 2) + pow(yn - yu, 2) + pow(zn - zu, 2), 0.5);
					if (distance[1] < distance[0]) {
						xu_searched = xu; yu_searched = yu; zu_searched = zu;
						//We need to determine at this stage if the projected node is within the searched element, if not we need to pass this element. 
						



						distance[0] = distance[1];
						searched_ele = i + 1;
					}
				}
			}//finish the searching through the fluid interface elements 
			//the corresponding fluid element and the projected node coordinate should be found at this stage
			//store the corresponding element 
			gs_flu[i*(hprefg + 1)*(hprefg + 1) + j] = searched_ele;
			//store the projected node location
			gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][0] = xu_searched;
			gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][1] = yu_searched;
			gs_flu_global[i*(hprefg + 1)*(hprefg + 1) + j][2] = zu_searched;


		}
	}
	
	//Next, we determine if the point is inside the element just found. If not, it will be considered as an orphan node and the closest node would be found
	if (element_type == 1) { //tetrahedral element
		Eigen::MatrixXd b(3, 1);
		Eigen::MatrixXd A(3, 3);
		Eigen::MatrixXd local(3, 1); //store the local coordinate 
		//determine the local coordinate in the tetrahedral element
		for (i = 0; i < gs_num; i++) { //for every structure gauss node
			A(0, 0) = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][0] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][0];
			A(0, 1) = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][0] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][0];
			A(0, 2) = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][0] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][0];
			A(1, 0) = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][1] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][1];
			A(1, 1) = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][1] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][1];
			A(1, 2) = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][1] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][1];
			A(2, 0) = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][2] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][2];
			A(2, 1) = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][2] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][2];
			A(2, 2) = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][2] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][2];
			b(0, 0) = gs_flu_global[i][0] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][0];
			b(1, 0) = gs_flu_global[i][1] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][1];
			b(2, 0) = gs_flu_global[i][2] - GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][2];
			local = A.inverse()*b;
			//determine if the point is inside the searched element
			if (local(0, 0) >= 0 && local(0, 0) <= 1 && local(1, 0) >= 0 && local(1, 0) <= 1 && local(2, 0) >= 0 && local(2, 0) <= 1
				&& local(0, 0) + local(1, 0) + local(2, 0) <= 1) {
				//The point is indeed inside the searched element and we can store the local coordinate in that element
				gs_flu_local[i][0] = local(0, 0); //xi
				gs_flu_local[i][1] = local(1, 0); //eta
				gs_flu_local[i][2] = local(2, 0); //zeta
			}
			else { //The node is an orphan
			//Search for the nearest point on the fluid FSI interface

			}
		}
	}
	if (element_type == 0) { //hex element
		double x1, x2, x3, x4, x5, x6, x7, x8;
		double y1, y2, y3, y4, y5, y6, y7, y8;
		double z1, z2, z3, z4, z5, z6, z7, z8;
		double x, y, z;
		double xi_k, eta_k, zeta_k;
		double xi_k_1, eta_k_1, zeta_k_1;
		double xi, eta, zeta; //converged local coordinate
		//The derivative is obtained by Mathematica: /Users/zhaokuanlu/OneDrive/post_processing/Hex_global_to_local_coordinate_Newton_iteration.nb
		//determine the local coordinate in the hexahedral element
		//Since the equation system is nonlinear, we use the Newton iteration and set the initial value to be 0 
		for (i = 0; i < gs_num; i++) { //for every structure gauss node
			xi_k = 0.0; eta_k = 0.0; zeta_k = 0.0; 
			x1 = GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][0]; x2 = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][0]; x3 = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][0]; x4 = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][0]; x5 = GCOORD[IEN_flu[4][gs_flu[i] - 1] - 1][0]; x6 = GCOORD[IEN_flu[5][gs_flu[i] - 1] - 1][0]; x7 = GCOORD[IEN_flu[6][gs_flu[i] - 1] - 1][0]; x8 = GCOORD[IEN_flu[7][gs_flu[i] - 1] - 1][0];
			y1 = GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][1]; y2 = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][1]; y3 = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][1]; y4 = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][1]; y5 = GCOORD[IEN_flu[4][gs_flu[i] - 1] - 1][1]; y6 = GCOORD[IEN_flu[5][gs_flu[i] - 1] - 1][1]; y7 = GCOORD[IEN_flu[6][gs_flu[i] - 1] - 1][1]; y8 = GCOORD[IEN_flu[7][gs_flu[i] - 1] - 1][1];
			z1 = GCOORD[IEN_flu[0][gs_flu[i] - 1] - 1][2]; z2 = GCOORD[IEN_flu[1][gs_flu[i] - 1] - 1][2]; z3 = GCOORD[IEN_flu[2][gs_flu[i] - 1] - 1][2]; z4 = GCOORD[IEN_flu[3][gs_flu[i] - 1] - 1][2]; z5 = GCOORD[IEN_flu[4][gs_flu[i] - 1] - 1][2]; z6 = GCOORD[IEN_flu[5][gs_flu[i] - 1] - 1][2]; z7 = GCOORD[IEN_flu[6][gs_flu[i] - 1] - 1][2]; z8 = GCOORD[IEN_flu[7][gs_flu[i] - 1] - 1][2];
			x = gs_flu_global[i][0]; y = gs_flu_global[i][1]; z = gs_flu_global[i][2];
			while (diff > tol) {
				fdot(0, 0) = -(1.0 / 8.0)* (1 - eta_k)* x1*(1 - zeta_k) + 1.0 / 8.0* (1 - eta_k)* x2*(1 - zeta_k) - 1.0 / 8.0 *(1 + eta_k)* x3*(1 - zeta_k) + 1.0 / 8.0 *(1 + eta_k)* x4*(1 - zeta_k) -
					1.0 / 8.0* (1 - eta_k)* x5*(1 + zeta_k) + 1.0 / 8.0* (1 - eta_k)* x6*(1 + zeta_k) - 1.0 / 8.0* (1 + eta_k)* x7*(1 + zeta_k) + 1.0 / 8.0* (1 + eta_k)* x8*(1 + zeta_k);
				fdot(0, 1) = -(1.0 / 8.0)* x1*(1 - xi_k)* (1 - zeta_k) + 1.0 / 8.0* x3*(1 - xi_k)* (1 - zeta_k) - 1.0 / 8.0* x2*(1 + xi_k)* (1 - zeta_k) + 1.0 / 8.0* x4*(1 + xi_k)* (1 - zeta_k) -
					1.0 / 8.0* x5*(1 - xi_k)* (1 + zeta_k) + 1.0 / 8.0* x7*(1 - xi_k)* (1 + zeta_k) - 1.0 / 8.0 *x6*(1 + xi_k)* (1 + zeta_k) + 1.0 / 8.0* x8*(1 + xi_k)* (1 + zeta_k);
				fdot(0, 2) = -(1.0 / 8.0)* (1 - eta_k)* x1*(1 - xi_k) - 1.0 / 8.0* (1 + eta_k) *x3*(1 - xi_k) + 1.0 / 8.0* (1 - eta_k)* x5*(1 - xi_k) + 1.0 / 8.0 *(1 + eta_k)* x7*(1 - xi_k) -
					1.0 / 8.0 *(1 - eta_k)* x2*(1 + xi_k) - 1.0 / 8.0 *(1 + eta_k)* x4*(1 + xi_k) + 1.0 / 8.0* (1 - eta_k) *x6*(1 + xi_k) + 1.0 / 8.0* (1 + eta_k)* x8*(1 + xi_k);
				fdot(1, 0) = -(1.0 / 8.0) *(1 - eta_k)* y1*(1 - zeta_k) + 1.0 / 8.0* (1 - eta_k)* y2*(1 - zeta_k) - 1.0 / 8.0* (1 + eta_k)* y3*(1 - zeta_k) + 1.0 / 8.0 *(1 + eta_k)* y4*(1 - zeta_k) -
					1.0 / 8.0 *(1 - eta_k)* y5*(1 + zeta_k) + 1.0 / 8.0* (1 - eta_k) *y6*(1 + zeta_k) - 1.0 / 8.0* (1 + eta_k)* y7*(1 + zeta_k) + 1.0 / 8.0* (1 + eta_k) *y8*(1 + zeta_k);
				fdot(1, 1) = -(1.0 / 8.0) *(1 - xi_k) *y1*(1 - zeta_k) - 1.0 / 8.0* (1 + xi_k)* y2*(1 - zeta_k) + 1.0 / 8.0 *(1 - xi_k)* y3*(1 - zeta_k) + 1.0 / 8.0* (1 + xi_k) *y4*(1 - zeta_k) -
					1.0 / 8.0* (1 - xi_k)* y5*(1 + zeta_k) - 1.0 / 8.0*(1 + xi_k)* y6*(1 + zeta_k) + 1.0 / 8.0* (1 - xi_k)* y7*(1 + zeta_k) + 1.0 / 8.0*(1 + xi_k)* y8*(1 + zeta_k);
				fdot(1, 2) = -(1.0 / 8.0)* (1 - eta_k)* (1 - xi_k)* y1 - 1.0 / 8.0* (1 - eta_k)* (1 + xi_k)* y2 - 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* y3 - 1.0 / 8.0* (1 + eta_k)* (1 + xi_k)* y4 +
					1.0 / 8.0* (1 - eta_k)* (1 - xi_k)* y5 + 1.0 / 8.0* (1 - eta_k)* (1 + xi_k)* y6 + 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* y7 + 1.0 / 8.0* (1 + eta_k)* (1 + xi_k)*y8;
				fdot(2, 0) = -(1.0 / 8.0) *(1 - eta_k)* z1*(1 - zeta_k) + 1.0 / 8.0* (1 - eta_k)* z2*(1 - zeta_k) - 1.0 / 8.0* (1 + eta_k) *z3*(1 - zeta_k) + 1.0 / 8.0* (1 + eta_k)* z4*(1 - zeta_k) -
					1.0 / 8.0* (1 - eta_k)* z5*(1 + zeta_k) + 1.0 / 8.0* (1 - eta_k)* z6*(1 + zeta_k) - 1.0 / 8.0* (1 + eta_k)* z7*(1 + zeta_k) + 1.0 / 8.0* (1 + eta_k)* z8*(1 + zeta_k);
				fdot(2, 1) = -(1.0 / 8.0)* (1 - xi_k)* z1*(1 - zeta_k) - 1.0 / 8.0* (1 + xi_k)* z2*(1 - zeta_k) + 1.0 / 8.0* (1 - xi_k) *z3*(1 - zeta_k) + 1.0 / 8.0* (1 + xi_k)* z4*(1 - zeta_k) -
					1.0 / 8.0* (1 - xi_k) *z5*(1 + zeta_k) - 1.0 / 8.0* (1 + xi_k)* z6*(1 + zeta_k) + 1.0 / 8.0* (1 - xi_k) *z7*(1 + zeta_k) + 1.0 / 8.0* (1 + xi_k)* z8*(1 + zeta_k);
				fdot(2, 2) = -(1.0 / 8.0) *(1 - eta_k)* (1 - xi_k)* z1 - 1.0 / 8.0 *(1 - eta_k)* (1 + xi_k)* z2 - 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* z3 - 1.0 / 8.0* (1 + eta_k) *(1 + xi_k) *z4 +
					1.0 / 8.0 *(1 - eta_k)* (1 - xi_k)* z5 + 1.0 / 8.0* (1 - eta_k)* (1 + xi_k) *z6 + 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* z7 + 1.0 / 8.0* (1 + eta_k)* (1 + xi_k) *z8;
				f(0, 0) = -x + 1.0 / 8.0 * (1 - eta_k)*x1*(1 - xi_k)*(1 - zeta_k) + 1.0 / 8.0*(1 + eta_k)*x3*(1 - xi_k)*(1 - zeta_k) +
					1.0 / 8.0* (1 - eta_k)* x2*(1 + xi_k)* (1 - zeta_k) + 1.0 / 8.0 *(1 + eta_k)* x4*(1 + xi_k) *(1 - zeta_k) +
					1.0 / 8.0*(1 - eta_k)*x5*(1 - xi_k) *(1 + zeta_k) + 1.0 / 8.0*(1 + eta_k)*x7*(1 - xi_k)*(1 + zeta_k) +
					1.0 / 8.0*(1 - eta_k)*x6*(1 + xi_k)*(1 + zeta_k) + 1.0 / 8.0*(1 + eta_k)*x8*(1 + xi_k)*(1 + zeta_k);
				f(1, 0) = -y + 1.0 / 8.0* (1 - eta_k)* (1 - xi_k)* y1*(1 - zeta_k) +
					1.0 / 8.0 *(1 - eta_k) *(1 + xi_k)* y2*(1 - zeta_k) + 1.0 / 8.0 *(1 + eta_k)* (1 - xi_k)* y3*(1 - zeta_k) +
					1.0 / 8.0* (1 + eta_k)* (1 + xi_k)* y4*(1 - zeta_k) + 1.0 / 8.0 *(1 - eta_k)* (1 - xi_k)* y5*(1 + zeta_k) +
					1.0 / 8.0* (1 - eta_k)* (1 + xi_k)* y6*(1 + zeta_k) + 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* y7*(1 + zeta_k) +
					1.0 / 8.0* (1 + eta_k)* (1 + xi_k)* y8*(1 + zeta_k);
				f(2, 0) = -z + 1.0 / 8.0 *(1 - eta_k)* (1 - xi_k) *z1*(1 - zeta_k) +
					1.0 / 8.0 *(1 - eta_k) *(1 + xi_k)* z2*(1 - zeta_k) + 1.0 / 8.0* (1 + eta_k) *(1 - xi_k)* z3*(1 - zeta_k) +
					1.0 / 8.0 *(1 + eta_k) *(1 + xi_k)* z4*(1 - zeta_k) + 1.0 / 8.0* (1 - eta_k)* (1 - xi_k) *z5*(1 + zeta_k) +
					1.0 / 8.0* (1 - eta_k)* (1 + xi_k)* z6*(1 + zeta_k) + 1.0 / 8.0* (1 + eta_k)* (1 - xi_k)* z7*(1 + zeta_k) +
					1.0 / 8.0 *(1 + eta_k)* (1 + xi_k) *z8*(1 + zeta_k);
				dx = -fdot.inverse()*f; //delta x^k
				diff = pow(pow(dx(0, 0), 2) + pow(dx(1, 0), 2) + pow(dx(2, 0), 2), 0.5);
				xi_k_1 = xi_k + dx(0, 0); eta_k_1 = eta_k + dx(1, 0); zeta_k_1 = zeta_k + dx(2, 0);
				xi_k = xi_k_1; eta_k = eta_k_1; zeta_k = zeta_k_1; //update the value
			}
			xi = xi_k; eta = eta_k; zeta = zeta_k; 
			if (xi >= -1 && xi <= 1 && eta >= -1 && eta <= 1 && zeta >= -1 && zeta <= 1) {
				//The point is indeed inside the searched element and we can store the local coordinate in that element
				gs_flu_local[i][0] = xi; //xi
				gs_flu_local[i][1] = eta; //eta
				gs_flu_local[i][2] = zeta; //zeta
			}
			else { //The node is an orphan

			}
		}
	}



	return;
}