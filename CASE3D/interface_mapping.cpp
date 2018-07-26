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

struct interface_mappingstruct interface_mapping(int fluid2structure, double ** GCOORD) {
	interface_mappingstruct t;
	int i, j, k, l, n, z;
	int u, v;
	extern OWETSURF ol[owsfnumber];
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
												+= ol[z].norm[l][n] * ol[z].WP[ol[z].IEN_gb[ol[z].LNA_2D[u][v] - 1][l] - 1] * ol[z].FPMASTER_2D[l][ol[z].LNA_algo2[i][j] - 1][ol[z].LNA_2D[u][v] - 1];
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
											+= ol[z].norm[l][n] * (1.0 / 3.0) * ol[z].WP[ol[z].IEN_gb[j][l] - 1] * (1.0 / 3.0) * ol[z].dimension[l];
									}
								}
							}
						}
					}
					ct += 1;
				}
			}

		}

		double BF_val[4];
		ct = 0;
		for (i = 0; i < owsfnumber; i++) {
			if (ol[i].FSNEL > 0) {
				BF_val[i] = 0.0;
				for (j = 0; j < wsflist[ct]->nnodes; j++) {
					for (k = 0; k < 3; k++) {
						BF_val[i] += wsflist[ct]->nodeforce[j * 3 + k];
					}
				}
				/*
				if (abs(BF_val[i] - ol[i].OBF_val) > 1) {
					std::cout << "the force mapping is inconsistent on wet surface: " << i << std::endl;
					//system("PAUSE ");
				}
				*/
				ct += 1;
			}
		}
		std::cout << " " << std::endl;

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
								for (u = 0; u < 2; u++) { //u,v stands for fem points 
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
			if (debug2 == 1) {
				for (z = 0; z < owsfnumber; z++) {
					for (l = 0; l < ol[z].FSNEL; l++) {  //loop through each element first
						for (i = 0; i < 3; i++) { //loop through each point on element surface 
							wsflist[z]->nodecoord[3 * (ol[z].IEN_2D[i][l] - 1) + 1] = GCOORD[ol[z].IEN_gb[i][l] - 1][1] + 1;
						}
					}
				}
			}
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

		std::cout << " " << std::endl;
		break;
	}
	t.energy_rec = 0;
	t.energy_sent = 0; 

	return t;
}

