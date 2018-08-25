//#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "header.h"
#include <fstream>
#include "data.h"
//set up two cases, the first case is map from user mesh to MpCCI mesh; the second case is map from MpCCI mesh to user mesh. 
//input all the force and displacement 
//currently it is in 3 directions set the other two direction to be zero. (For example on surface A, nodisp has only component in x direction)

//Mapping algorithm explanation
/*
There are totally 4 algorithms:
The first one is the h/p refinement based on Farhat's energy conservation algorithm (1998). The code maps the force from FEM/SEM mesh to the base fluid mesh which is identical with the structural mesh.
The second one breaks the higher order element into sub-elements and let the MpCCI do the interpolation. 
The third one is the algorithm 1 in an average sense. The force is averaged on h/p refined mesh and then given to the base fluid mesh. The force is still conservative.
The fourth one is the Sprague's pressure averaging approach (also used by Klenow), which can also be comprehended as the consistent mapping algorithm in Farhat's paper (1998)
The difference from the 8-4 version:
The algorithm 1 and 3 are combined into algorithm 1 which is capable of dealing with h & p refinement. 
The algorithm 2 stays the same. 
The algorithm 4 is now the algorithm 3. 
A new algorithm 4 is added. 
*/

struct interface_mappingstruct interface_mapping(int fluid2structure, double ** GCOORD, double* WP, int** IEN, int***LNA) {
	interface_mappingstruct t;
	int i, j, k, l, m, n, h, q, ii, jj, kk, z;
	int u, v;
	extern OWETSURF ol[owsfnumber];
	extern STRU_WET_SURF ss[ssnumber]; //data structure used to store the properties on the structure wetted surface
	double accu = 0.0;
	int ct = 0;
	
	//This subroutine is temporarily changed to test the BS fluid domain (without side modules)
	switch (fluid2structure) //sem points to fem points
	{
	case 1: 
		if (mappingalgo == 2) {
			//In this code, we use the 2nd order GLL integration instead of 1st order GL integration in FSP case since we have to include a 2D Jacobian determinate to accommodate the geometric mapping
			ct = 0;
			for (z = 0; z < owsfnumber; z++) {
				if (ol[z].FSNEL > 0) {
					for (l = 0; l < ol[z].FSNEL; l++) { //loop through each element
						for (i = 0; i < 2; i++) {
							for (j = 0; j < 2; j++) {
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[i][j] - 1][l] - 1) + 0] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[i][j] - 1][l] - 1) + 1] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[i][j] - 1][l] - 1) + 2] = 0.0;
							}
						}
					}
					ct += 1;
				}
			}

			ct = 0;
			//The problem remains to be fixed is that the pressure value is at the interpolation point instead of the quadrature point
			for (z = 0; z < owsfnumber; z++) {
				if (ol[z].FSNEL > 0) {
					for (l = 0; l < ol[z].FSNEL; l++) {
						for (i = 0; i < 2; i++) {
							for (j = 0; j < 2; j++) {
								for (u = 0; u < NINT; u++) { //u v stands for y x
									for (v = 0; v < NINT; v++) {
										for (n = 0; n < 3; n++) {
											wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[i][j] - 1][l] - 1) + n] //force in z direction
												+= ol[z].norm[l][n] * WP[ol[z].IEN_gb[ol[z].LNA_2D[u][v] - 1][l] - 1] * ol[z].FPMASTER_2D[l][ol[z].LNA_algo2[i][j] - 1][ol[z].LNA_2D[u][v] - 1];
										}
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
		}

		else if (mappingalgo == 4) {
			if (element_type == 1) { //tetrahedral element
				ct = 0;
				for (z = 0; z < owsfnumber; z++) {
					if (ol[z].FSNEL > 0) {
						for (l = 0; l < ol[z].FSNEL; l++) { //loop through each element
							for (i = 0; i < 3; i++) {
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[i][l] - 1) + 0] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[i][l] - 1) + 1] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[i][l] - 1) + 2] = 0.0;
							}
						}
						ct += 1;
					}
				}
				ct = 0;
				//The problem remains to be fixed is that the pressure value is at the interpolation point instead of the quadrature point
				for (z = 0; z < owsfnumber; z++) {
					if (ol[z].FSNEL > 0) {
						for (l = 0; l < ol[z].FSNEL; l++) {
							for (i = 0; i < 3; i++) {
								for (j = 0; j < 3; j++) {
									for (n = 0; n < 3; n++) {
										wsflist[ct]->nodeforce[3 * (ol[z].IEN_2D[i][l] - 1) + n] //force in z direction
											+= ol[z].norm[l][n] * (1.0 / 3.0) * WP[ol[z].IEN_gb[j][l] - 1] * (1.0 / 3.0) * ol[z].dimension[l];
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
		}

		else if (mappingalgo == 5) {
			//Interpolate the pressure from fluid interpolation node to structure gauss node ss[z].P_gs[ss[z].LNA_gs[v][u] - 1][l]
			LOBATTOstruct b;
			b = LOBATTO(N);
			GLLQUADstruct f;
			f = GLLQUAD(b.Z, b.WL, N, !FEM);

			double basep[NINT];
			if (FEM == 0) { //SEM interpolation points
				for (i = 0; i < NINT; i++) {
					basep[i] = f.S[i]; //GLL point 
				}
			}
			else { //FEM interpolation points
				basep[0] = -1;
				basep[1] = 1;
			}
			double lcx, lcy, lcz;
			double nomx, nomy, nomz; double denomx, denomy, denomz; 
			double phig; 

			
			if (debug_algo5 == 1) {
				for (z = 0; z < ssnumber; z++) {
					for (l = 0; l < ss[z].ELE_stru*(hprefg + 1)*(hprefg + 1); l++) {
						for (ii = 0; ii < NINT; ii++) {
							for (jj = 0; jj < NINT; jj++) {
								for (kk = 0; kk < NINT; kk++) {
									WP[IEN[LNA[ii][jj][kk] - 1][ss[z].gs_flu[l] - 1] - 1] = 1.0;
								}
							}
						}
					}
				}
			}

			if (element_type == 0) {
				for (z = 0; z < ssnumber; z++) {
					for (l = 0; l < ss[z].ELE_stru; l++) {
						for (h = 0; h < hprefg + 1; h++) {
							for (q = 0; q < hprefg + 1; q++) {
								ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] = 0.0; //initialize the value 
								if (ss[z].orphan_flag_gs[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] == 0) { //the node is not an orphan
									//local coordinate of the projected structure gauss point
									lcx = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][0];
									lcy = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][1];
									lcz = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][2];
									for (ii = 0; ii < NINT; ii++) {
										for (jj = 0; jj < NINT; jj++) {
											for (kk = 0; kk < NINT; kk++) {
												nomx = 1.0; nomy = 1.0; nomz = 1.0; //multiplier initialization
												denomx = 1.0; denomy = 1.0; denomz = 1.0; //multiplier initialization
												for (m = 0; m < NINT; m++) {
													if (m != ii) {
														nomx *= (lcx - basep[m]);
														denomx *= (basep[ii] - basep[m]);
													}
													if (m != jj) {
														nomy *= (lcy - basep[m]);
														denomy *= (basep[jj] - basep[m]);
													}
													if (m != kk) {
														nomz *= (lcz - basep[m]);
														denomz *= (basep[kk] - basep[m]);
													}
												}
												phig = (nomx / denomx)*(nomy / denomy)*(nomz / denomz);
												ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] += WP[IEN[LNA[ii][jj][kk] - 1][ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1] - 1] * phig;
												//P_gs is created in Neighborhood_search.cpp 
											}
										}
									}
								}
								else { //the node is an orphan, directly give the value on the nearest node 
									ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] = WP[ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1];
								}
							}
						}
					}
				}
			}
			if (element_type == 1) {
				for (z = 0; z < ssnumber; z++) {
					for (l = 0; l < ss[z].ELE_stru; l++) {
						for (h = 0; h < hprefg + 1; h++) {
							for (q = 0; q < hprefg + 1; q++) {
								ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] = 0.0; //initialize the value 
								if (ss[z].orphan_flag_gs[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] == 0) { //the node is not an orphan
									//local coordinate of the projected structure gauss point
									lcx = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][0];
									lcy = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][1];
									lcz = ss[z].gs_flu_local[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q][2];
									ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] += WP[IEN[0][ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1] - 1] * (1 - lcx - lcy - lcz) +
										WP[IEN[1][ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1] - 1] * lcx +
										WP[IEN[2][ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1] - 1] * lcy +
										WP[IEN[3][ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1] - 1] * lcz;
								}
								else { //the node is an orphan, directly give the value on the nearest node 
									ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] = WP[ss[z].gs_flu[l*(hprefg + 1)*(hprefg + 1) + h*(hprefg + 1) + q] - 1];
								}
							}
						}
					}
				}
			}
			/*
			double hd = 0.0;
			if (debug_algo5 == 1) {
				for (l = 0; l < ss[z].ELE_stru; l++) {
					for (h = 0; h < hprefg + 1; h++) {
						for (q = 0; q < hprefg + 1; q++) {
							hd += ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l];
						}
					}
				}
			} //wrong 
			*/
			/*
			if (debug_algo5 == 1) {
				for (l = 0; l < ss[z].ELE_stru; l++) {
					for (h = 0; h < hprefg + 1; h++) {
						for (q = 0; q < hprefg + 1; q++) {
							ss[z].P_gs[ss[z].LNA_gs[h][q] - 1][l] = 0.0;
						}
					}
				}
				for (h = 0; h < hprefg + 1; h++) {
					for (q = 0; q < hprefg + 1; q++) {
						ss[z].P_gs[ss[z].LNA_gs[0][0] - 1][0] = 1.0;
						ss[z].P_gs[ss[z].LNA_gs[1][0] - 1][0] = 2.0;
						ss[z].P_gs[ss[z].LNA_gs[1][1] - 1][0] = 3.0;
						ss[z].P_gs[ss[z].LNA_gs[0][1] - 1][0] = 4.0;
					}
				}
			}
			*/
			for (z = 0; z < ssnumber; z++) {
				for (l = 0; l < ss[z].ELE_stru; l++) { //loop through each element first
					for (i = 0; i < 4; i++) { //i,j stands for fem points y m 
						wsflist[z]->nodeforce[3 * (ss[z].IEN_stru_MpCCI[i][l] - 1) + 0] = 0.0;
						wsflist[z]->nodeforce[3 * (ss[z].IEN_stru_MpCCI[i][l] - 1) + 1] = 0.0;
						wsflist[z]->nodeforce[3 * (ss[z].IEN_stru_MpCCI[i][l] - 1) + 2] = 0.0;
					}
				}
				//phi_femg: the linear shape function value at base mesh Gauss-Legendre nodes
				for (l = 0; l < ss[z].ELE_stru; l++) {
					for (i = 0; i < 2; i++) {
						for (j = 0; j < 2; j++) {
							for (u = 0; u < hprefg + 1; u++) { //u v stands for y x
								for (v = 0; v < hprefg + 1; v++) {
									for (n = 0; n < 3; n++) {
										wsflist[z]->nodeforce[3 * (ss[z].IEN_stru_MpCCI[ss[z].LNA_stru[i][j] - 1][l] - 1) + n] //force in m direction
											+= ss[z].norm_stru[l][n] * ss[z].W_stru[u] * ss[z].W_stru[v] * ss[z].P_gs[ss[z].LNA_gs[v][u] - 1][l] * ss[z].phi_stru[ss[z].LNA_stru[i][j] - 1][v][u] * ss[z].Jacob_stru[l][u*(hprefg + 1) + v];
										//ss[z].phi_stru is the linear shape function value on gauss node (Gauss-Legendre)
									}
								}
							}
						}
					}
				}
			}
			double BF_val[ssnumber];
			for (z = 0; z < ssnumber; z++) {
				BF_val[z] = 0.0; 
				for (j = 0; j < wsflist[z]->nnodes; j++) {
					for (k = 0; k < 3; k++) {
						BF_val[z] += wsflist[z]->nodeforce[j * 3 + k];
					}
				}
			}
			std::cout << " " << std::endl;
		}
		break;

	case 0: //map disp from coupling mesh(fem mesh) to user defined mesh(sem mesh)
		//The displacement is mapped from the linear element to the corresponding high-order element
		double hd = 0.0; 
		if (mappingalgo == 2) {
			hd = 0.0;
			for (z = 0; z < owsfnumber; z++) {
				for (j = 0; j < wsflist[z]->nnodes; j++) {
					for (k = 0; k < 3; k++) {
						hd += wsflist[z]->nodecoord[j * 3 + k];
					}
				}
			}
			std::cout << " " << std::endl;

			double DISPTEMP = 0.0;

			ct = 0;
			for (z = 0; z < owsfnumber; z++) {
				if (ol[z].FSNEL > 0) {
					for (l = 0; l < ol[z].FSNEL; l++) {  //loop through each element first
						for (i = 0; i < NINT; i++) { //loop through each point on element surface 
							for (j = 0; j < NINT; j++) {
								ol[z].DISP[ol[z].IEN_gb[ol[z].LNA_2D[i][j] - 1][l] - 1][0] = 0.0;
								ol[z].DISP[ol[z].IEN_gb[ol[z].LNA_2D[i][j] - 1][l] - 1][1] = 0.0;
								ol[z].DISP[ol[z].IEN_gb[ol[z].LNA_2D[i][j] - 1][l] - 1][2] = 0.0;
								for (u = 0; u < 2; u++) {
									for (v = 0; v < 2; v++) {
										for (n = 0; n < 3; n++) {
											//DISPTEMP = wsflist[ct]->nodecoord[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[u][v] - 1][l] - 1) + n] - GCOORD[ol[z].IEN_gb[ol[z].LNA_2D[u][v] - 1][l] - 1][n];
											DISPTEMP = wsflist[ct]->nodecoord[3 * (ol[z].IEN_2D[ol[z].LNA_algo2[u][v] - 1][l] - 1) + n] - GCOORD[ol[z].IEN_gb[ol[z].LNA_norm[ol[z].LNA_algo2[u][v] - 1] - 1][l] - 1][n];
											ol[z].DISP[ol[z].IEN_gb[ol[z].LNA_2D[i][j] - 1][l] - 1][n] += DISPTEMP * ol[z].phi_fem[ol[z].LNA_algo2[u][v] - 1][i][j];
										}
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
			hd = 0.0;
			for (z = 0; z < owsfnumber; z++) {
				for (i = 0; i < ol[z].GIDNct; i++) {
					for (n = 0; n < 3; n++) {
						hd += ol[z].DISP[ol[z].GIDN[i] - 1][n];
					}
				}
			}
		}

		else if (mappingalgo == 4) {
			/*
			if (debug2 == 1) {
				for (z = 0; z < owsfnumber; z++) {
					for (l = 0; l < ol[z].FSNEL; l++) {  //loop through each element first
						for (i = 0; i < 3; i++) { //loop through each point on element surface 
							wsflist[z]->nodecoord[3 * (ol[z].IEN_2D[i][l] - 1) + 1] = GCOORD[ol[z].IEN_gb[i][l] - 1][1] + 1;
						}
					}
				}
			}
			*/
			if (element_type == 1) { //tetrahedral element
				double DISPTEMP = 0.0;
				ct = 0;
				for (z = 0; z < owsfnumber; z++) {
					if (ol[z].FSNEL > 0) {
						for (l = 0; l < ol[z].FSNEL; l++) {  //loop through each element first
							for (i = 0; i < 3; i++) { //loop through each point on element surface 
								ol[z].DISP[ol[z].IEN_gb[i][l] - 1][0] = 0.0;
								ol[z].DISP[ol[z].IEN_gb[i][l] - 1][1] = 0.0;
								ol[z].DISP[ol[z].IEN_gb[i][l] - 1][2] = 0.0;
								for (n = 0; n < 3; n++) {
									DISPTEMP = wsflist[ct]->nodecoord[3 * (ol[z].IEN_2D[i][l] - 1) + n] - GCOORD[ol[z].IEN_gb[i][l] - 1][n];
									ol[z].DISP[ol[z].IEN_gb[i][l] - 1][n] = DISPTEMP;
								}
							}
						}
						ct += 1;
					}
				}
			}
		}

		else if (mappingalgo == 5) {
			//Interpolate the structural displacement to fluid interpolation node
			double basep[2]; 
			basep[0] = -1;
			basep[1] = 1;
			double lcx, lcy;
			double nomx, nomy; double denomx, denomy;
			double phig;
			double DISPTEMP = 0.0;
			
			if (debug_algo5 == 1) {
				for (z = 0; z < ssnumber; z++) {
					for (i = 0; i < ss[z].Node_stru; i++) {
						for (n = 0; n < 3; n++) {
							wsflist[z]->nodecoord[3 * i + n] = ss[z].GCOORD_stru[ss[z].Node_glob[i] - 1][n] + 1;
						}
					}
				}
			}
			
			int elenode2D_gs;
			if (element_type == 0) { //hex element
				elenode2D_gs = NINT*NINT;
			}
			if (element_type == 1) { //tet element
				elenode2D_gs = 3;
			}

			for (z = 0; z < owsfnumber; z++) {
				for (l = 0; l < ol[z].FSNEL; l++) {  //loop through each element first
					for (i = 0; i < elenode2D_gs; i++) { //loop through each point on element surface
						ol[z].DISP_gs[l*elenode2D_gs + i][0] = 0.0;
						ol[z].DISP_gs[l*elenode2D_gs + i][1] = 0.0;
						ol[z].DISP_gs[l*elenode2D_gs + i][2] = 0.0;
						if (ol[z].orphan_flag_flu[l*elenode2D_gs + i] == 0) { //the node is not an orphan
							lcx = ol[z].flu_local[l*elenode2D_gs + i][0];
							lcy = ol[z].flu_local[l*elenode2D_gs + i][1];
							for (u = 0; u < 2; u++) { //u,v stands for fem points 
								for (v = 0; v < 2; v++) {
									nomx = 1.0; nomy = 1.0; //multiplier initialization
									denomx = 1.0; denomy = 1.0; //multiplier initialization
									for (m = 0; m < 2; m++) {
										if (m != u) {
											nomx *= (lcx - basep[m]);
											denomx *= (basep[u] - basep[m]);
										}
										if (m != v) {
											nomy *= (lcy - basep[m]);
											denomy *= (basep[v] - basep[m]);
										}
									}
									phig = (nomx / denomx)*(nomy / denomy);
									for (n = 0; n < 3; n++) {
										DISPTEMP = wsflist[z]->nodecoord[3 * (ss[z].IEN_stru_MpCCI[ss[z].LNA_stru[u][v] - 1][ol[z].flu_stru[l*elenode2D_gs + i] - 1] - 1) + n] - ss[z].GCOORD_stru[ss[z].IEN_stru[ss[z].LNA_stru[u][v] - 1][ol[z].flu_stru[l*elenode2D_gs + i] - 1] - 1][n];
										ol[z].DISP_gs[l*elenode2D_gs + i][n] += DISPTEMP * phig;
									}
								}
							}
						}
						else { //the node is an orphan
							for (n = 0; n < 3; n++) {
								ol[z].DISP_gs[l*elenode2D_gs + i][n] = wsflist[z]->nodecoord[3 * (ol[z].flu_stru[l*elenode2D_gs + i] - 1) + n] - ss[z].GCOORD_stru[ss[z].Node_glob[ol[z].flu_stru[l*elenode2D_gs + i] - 1] - 1][n];
								//flu_stru for orphan node is actually the node number instead of the element number 
							}
						}
					}
				}
			}

			hd = 0.0;
			for (z = 0; z < owsfnumber; z++) {
				for (i = 0; i < ol[z].FSNEL*elenode2D_gs; i++) {
					for (n = 0; n < 3; n++) {
						hd += ol[z].DISP_gs[i][n];
					}
				}
			}
			std::cout << " " << std::endl;
		}

		
		break;
	}
	t.energy_rec = 0;
	t.energy_sent = 0; 

	return t;
}

