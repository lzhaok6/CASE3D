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

double BF_val[4];
struct interface_mappingstruct interface_mapping(int fluid2structure, int**IEN_3D, int***LNA_3D, int**LNA_2D, int** LNA_base, int**LNA_basealgo5, double *Z, int TIME, double ** GCOORD, double ***phi_fem, double* W, double*** phi_femg, double*** phi_sem2) {
	interface_mappingstruct t;
	int i, j, k, l, m, n;
	int u, v;
	extern OWETSURF ol[owsfnumber];
	double accu = 0.0;
	int ct = 0;
	
	//This subroutine is temporarily changed to test the BS fluid domain (without side modules)
	switch (fluid2structure) //sem points to fem points
	{
	case 1: 
		ct = 0;
		if (mappingalgo == 2) {
			//In this code, we use the 2nd order GLL integration instead of 1st order GL integration in FSP case since we have to include a 2D Jacobian determinate to accommodate the geometric mapping
			LOBATTOstruct b;
			b = LOBATTO(2);
			GLLQUADstruct f;
			f = GLLQUAD(b.Z, b.WL, 2, 1); //linear integration (Gauss-Legendre points) 
			ct = 0;
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL; l++) { //loop through each element
						for (i = 0; i < 2; i++) {
							for (j = 0; j < 2; j++) {
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 0] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 1] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 2] = 0.0;
							}
						}
					}
					ct += 1;
				}
			}

			ct = 0;
			//The problem remains to be fixed is that the pressure value is at the interpolation point instead of the quadrature point
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL; l++) {
						for (i = 0; i < 2; i++) {
							for (j = 0; j < 2; j++) {
								for (u = 0; u < 3; u++) { //u v stands for y x
									for (v = 0; v < 3; v++) {
										for (n = 0; n < 3; n++) {
											wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + n] //force in z direction
												+= ol[k].norm[l][n] * f.W[u] * f.W[v] * ol[k].WP[ol[k].IEN_gb[LNA_2D[u][v] - 1][l] - 1] * phi_sem2[LNA_2D[i][j] - 1][v][u] * ol[k].Jacob_2D[l][u * 3 + v];
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
		
		std::cout << " " << std::endl;
		break;

	case 0: //map disp from coupling mesh(fem mesh) to user defined mesh(sem mesh)
		double DISPTEMP;
		
		if (mappingalgo == 2) {
			ct = 0;
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL_fem > 0) {
					for (l = 0; l < ol[m].FSNEL_fem; l++) {  //loop through each element
						for (u = 0; u < 2; u++) {
							for (v = 0; v < 2; v++) {
								for (n = 0; n < 3; n++) {
									ol[m].DISP[ol[m].IEN_gb[LNA_2D[u][v] - 1][l] - 1][n]
										= wsflist[ct]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[u][v] - 1][l] - 1) + n] - GCOORD[ol[m].IEN_gb[LNA_2D[u][v] - 1][l] - 1][n];
								}
							}
						}
					}
					ct += 1;
				}
			}
		}

		break;
	}

	return t;
}

