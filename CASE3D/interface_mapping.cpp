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
struct interface_mappingstruct interface_mapping(int fluid2structure, int**IEN_3D, int***LNA_3D, int**LNA_2D, int** LNA_base, int**LNA_basealgo5, double *Z, int TIME, double ** GCOORD, double ***phi_fem, double* W, double*** phi_femg, double*** phi_fem2) {
	interface_mappingstruct t;
	int i, j, k, l, m, h;
	int u, v;
	//extern int owsfnumber;
	extern OWETSURF ol[owsfnumber];
	double accu = 0.0;
	int ct = 0;
	t.energy_rec = 0.0;
	t.energy_sent = 0.0;
	
	//This subroutine is temporarily changed to test the BS fluid domain (without side modules)
	switch (fluid2structure) //sem points to fem points
	{
	case 1: 
		ct = 0;
		if (mappingalgo == 1) {
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL / refine / refine; l++) { //loop through each element first
						for (i = 0; i < NCINT; i++) { //i,j stands for fem points y z 
							for (j = 0; j < NCINT; j++) {
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 0] = 0.0; //error prone
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 1] = 0.0;
								wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + 2] = 0.0;
							}
						}
					}
					ct += 1;
				}
			}
			ct = 0;
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {
						for (i = 0; i < NCINT; i++) {
							for (j = 0; j < NCINT; j++) {
								for (u = 0; u < hpref + 1; u++) { //u v stands for y x
									for (v = 0; v < hpref + 1; v++) {
										wsflist[ct]->nodeforce[3 * (ol[m].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[m].dir] //force in z direction
											+= ol[m].BF3[LNA_base[v][u] - 1][l] * phi_fem[LNA_2D[i][j] - 1][v][u];
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
		}
		
		else if (mappingalgo == 2) {  
			
			if (nodeforcemap2 == 1) {
				LOBATTOstruct b;
				b = LOBATTO(1);
				GLLQUADstruct f;
				f = GLLQUAD(b.Z, b.WL, 1, 0); //linear integration (Gauss-Legendre points) 
				ct = 0;
				for (k = 0; k < owsfnumber; k++) {
					if (ol[k].FSNEL > 0) {
						for (l = 0; l < ol[k].FSNEL_fem; l++) { //loop through each element
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
				for (k = 0; k < owsfnumber; k++) {
					if (ol[k].FSNEL > 0) {
						for (l = 0; l < ol[k].FSNEL_fem; l++) {
							for (i = 0; i < 2; i++) {
								for (j = 0; j < 2; j++) {
									for (u = 0; u < 2; u++) { //u v stands for y x
										for (v = 0; v < 2; v++) {
											wsflist[ct]->nodeforce[3 * (ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[k].dir] //force in z direction
												+= ol[k].NORM_BF * f.W[u] * f.W[v] * ol[k].WP[ol[k].IEN_gb[LNA_2D[u][v] - 1][l] - 1] * phi_fem2[LNA_2D[i][j] - 1][v][u] * (ol[k].XYHE_gb[0][l] / 2.0)*(ol[k].XYHE_gb[1][l] / 2.0);
										}
									}
								}
							}
						}
						ct += 1;
					}
				}
				//std::cout << "" << std::endl;
			}
			else {
				ct = 0;
				for (k = 0; k < owsfnumber; k++) {
					if (ol[k].FSNEL > 0) {
						for (l = 0; l < ol[k].FSNEL_fem; l++) { //loop through each element
							for (i = 0; i < 2; i++) {
								for (j = 0; j < 2; j++) {
									wsflist[ct]->nodepressure[ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1] = 0.0;
									wsflist[ct]->nodepressure[ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1] = 0.0;
									wsflist[ct]->nodepressure[ol[k].IEN_2D[LNA_2D[i][j] - 1][l] - 1] = 0.0;
								}
							}
						}
						ct += 1;
					}
				}
				ct = 0;
				for (k = 0; k < owsfnumber; k++) {
					if (ol[k].FSNEL > 0) {
						for (l = 0; l < ol[k].FSNEL_fem; l++) {
							for (u = 0; u < 2; u++) {
								for (v = 0; v < 2; v++) {
									wsflist[ct]->nodepressure[ol[k].IEN_2D[LNA_2D[u][v] - 1][l] - 1] 
										= ol[k].WP[ol[k].IEN_gb[LNA_2D[u][v] - 1][l] - 1];
								}
							}
						}
						ct += 1;
					}
				}
			}
			//ol[m].DISP[ol[m].IEN_gb[LNA_2D[u][v] - 1][l] - 1][ol[m].dir][1]
				//= wsflist[ct]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[u][v] - 1][l] - 1) + ol[m].dir] - ol[m].location;
			//std::cout << " " << std::endl;
		}
		
		else if (mappingalgo == 3) {
			//Average force mapping for h/p refinement
			ct = 0;
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL / refine / refine; l++) { //loop through each element first
						for (i = 0; i < NCINT; i++) { //i,j stands for fem points y z 
							for (j = 0; j < NCINT; j++) {
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
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {
						for (i = 0; i < NCINT; i++) {
							for (j = 0; j < NCINT; j++) {
								accu = 0.0;
								for (u = 0; u < hpref + 1; u++) { //u v stands for y x
									for (v = 0; v < hpref + 1; v++) {
										accu += ol[m].BF3[LNA_base[v][u] - 1][l];
									}
								}
								wsflist[ct]->nodeforce[3 * (ol[m].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[m].dir] += accu / NCINT / NCINT;
							}
						}
					}
					ct += 1;
				}
			}
		}
		
		else if (mappingalgo == 4) {
			ct = 0;
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL / refine / refine; l++) { //loop through each element first
						for (i = 0; i < NCINT; i++) { //i,j stands for fem points y z 
							for (j = 0; j < NCINT; j++) {
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
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {
						for (i = 0; i < NCINT; i++) {
							for (j = 0; j < NCINT; j++) {
								accu = 0.0;
								for (u = 0; u < hpref + 1; u++) { //u v stands for y x
									for (v = 0; v < hpref + 1; v++) {
										accu += ol[m].BP[LNA_base[v][u] - 1][l]; //total pressure on base fluid element
									}
								}
								wsflist[ct]->nodeforce[3 * (ol[m].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[m].dir] += (accu / (hpref + 1) / (hpref + 1))*(XHE*refine*YHE*refine / NCINT / NCINT)*ol[m].NORM;
							}
						}
					}
					ct += 1;
				}
			}
			//In the loop above, we assume the mesh size is the same in all directions. Thus, we used XHE*YHE to represent elemennt dimension on all wetted surfaces. 
			//However, that's problematic if the mesh size is different in x, y and z direction. 
			//std::cout << " " << std::endl;
		}

		else if (mappingalgo == 5) {
			ct = 0;
			for (k = 0; k < owsfnumber; k++) {
				if (ol[k].FSNEL > 0) {
					for (l = 0; l < ol[k].FSNEL / refine / refine; l++) { //loop through each element first
						for (i = 0; i < NCINT; i++) { //i,j stands for fem points y z 
							for (j = 0; j < NCINT; j++) {
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
			//BPG: the pressure on Gauss-Legendre nodes
			//phi_femg: the linear shape function value at base mesh Gauss-Legendre nodes
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {
						for (i = 0; i < NCINT; i++) {
							for (j = 0; j < NCINT; j++) {
								for (u = 0; u < hprefg + 1; u++) { //u v stands for y x
									for (v = 0; v < hprefg + 1; v++) {
										wsflist[ct]->nodeforce[3 * (ol[m].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[m].dir] //force in z direction
											+= ol[m].NORM_BF * W[u] * W[v] * ol[m].BPG[LNA_basealgo5[v][u] - 1][l] * phi_femg[LNA_2D[i][j] - 1][v][u] * (XHE * refine / 2.0)*(YHE * refine / 2.0);
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
		}
		
		//check if the force mapping is conservative 
		//double BF_val[4];
		ct = 0;
		for (i = 0; i < 4; i++) {
			if (ol[i].FSNEL > 0) {
				BF_val[i] = 0.0;
				for (j = 0; j < wsflist[ct]->nnodes; j++) {
					BF_val[i] += wsflist[ct]->nodeforce[j * 3 + ol[i].dir];
				}
				if (abs(BF_val[i] - ol[i].OBF_val) > 1) {
					std::cout << "the force mapping is inconsistent on wet surface: " << i << std::endl;
					//system("PAUSE ");
				}
				ct += 1;
			}
		}
		std::cout << " " << std::endl;
		break;

	case 0: //map disp from coupling mesh(fem mesh) to user defined mesh(sem mesh)
		double DISPTEMP;
		if (mappingalgo == 1 || mappingalgo == 5) {
			
			if (debug == 1) {
				for (m = 0; m < owsfnumber; m++) {
					for (l = 0; l < ol[m].FSNEL; l++) {
						for (i = 0; i < NCINT; i++) {
							for (j = 0; j < NCINT; j++) {
								wsflist[m]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[i][j] - 1][l] - 1) + ol[m].dir] = ol[m].location + ((double)rand() / (RAND_MAX));
							}
						}
					}
				}
			}
			
			ct = 0;
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {  //loop through each element first
						for (i = 0; i < hpref + 1; i++) { //loop through each point on element surface 
							for (j = 0; j < hpref + 1; j++) {
								ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][0][1] = 0.0;
								ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][1][1] = 0.0;
								ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][2][1] = 0.0;
								for (u = 0; u < NCINT; u++) { //u,v stands for fem points 
									for (v = 0; v < NCINT; v++) {
										DISPTEMP = wsflist[ct]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[u][v] - 1][l] - 1) + ol[m].dir] - ol[m].location;
										//ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][ol[m].dir]
											//+= DISPTEMP * phi_fem_alias[LNA_2D[u][v] - 1][j + i*(hpref + 1)][m];
										//phi_fem_alias is equivalent to phi_fem[LNA_2D[u][v] - 1][i][j]
										ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][ol[m].dir][1]
											+= DISPTEMP * phi_fem[LNA_2D[u][v] - 1][i][j];
									}
								}
							}
						}
					}
					ct += 1;
				}
			}
		}
		
		else if (mappingalgo == 2) {
			ct = 0;
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL_fem > 0) {
					for (l = 0; l < ol[m].FSNEL_fem; l++) {  //loop through each element
						for (u = 0; u < 2; u++) {
							for (v = 0; v < 2; v++) {
								ol[m].DISP[ol[m].IEN_gb[LNA_2D[u][v] - 1][l] - 1][ol[m].dir][1]
									= wsflist[ct]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[u][v] - 1][l] - 1) + ol[m].dir] - ol[m].location;
							}
						}
					}
					ct += 1;
				}
			}
		}
		
		else if (mappingalgo == 3 || mappingalgo == 4) {
			//for debug purpose
			//wsflist[3]->nodecoord[3 * (ol[3].IEN_2D[LNA_2D[0][1] - 1][0] - 1) + ol[3].dir] = ol[3].location + 1.0;
			/*
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[0][0] - 1][0] - 1) + ol[2].dir] = ol[2].location + 1.0;
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[0][1] - 1][0] - 1) + ol[2].dir] = ol[2].location + 2.0;
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[0][0] - 1][1] - 1) + ol[2].dir] = ol[2].location + 3.0;
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[0][1] - 1][1] - 1) + ol[2].dir] = ol[2].location + 4.0;
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[1][0] - 1][1] - 1) + ol[2].dir] = ol[2].location + 1.0;
			wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[LNA_2D[1][1] - 1][1] - 1) + ol[2].dir] = ol[2].location + 2.0;
			*/
			/*
			for (i = 0; i < ol[2].FSNEL / refine / refine;i++) {
				for (j = 0; j < 4; j++) {
					wsflist[0]->nodecoord[3 * (ol[2].IEN_2D[j][i] - 1) + ol[2].dir] = ol[2].location + 2.0;
				}
			}
			*/
			for (m = 0; m < owsfnumber; m++) {
				for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {  //loop through each element first
					for (i = 0; i < hpref + 1; i++) { //loop through each point on element surface 
						for (j = 0; j < hpref + 1; j++) {
							ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][ol[m].dir][1] = 0.0;
						}
					}
				}
			}

			ct = 0;
			for (m = 0; m < owsfnumber; m++) {
				if (ol[m].FSNEL > 0) {
					for (l = 0; l < ol[m].FSNEL / refine / refine; l++) {  //loop through each element first
						for (i = 0; i < hpref + 1; i++) { //loop through each point on element surface 
							for (j = 0; j < hpref + 1; j++) {
								accu = 0.0;
								for (u = 0; u < NCINT; u++) { //u,v stands for fem points 
									for (v = 0; v < NCINT; v++) {
										accu += wsflist[ct]->nodecoord[3 * (ol[m].IEN_2D[LNA_2D[u][v] - 1][l] - 1) + ol[m].dir] - ol[m].location;
									}
								}
								ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][ol[m].dir][1] = accu / NCINT / NCINT;
								//ol[m].DISP[ol[m].IEN_base[LNA_base[i][j] - 1][l] - 1][ol[m].dir] += (accu / NCINT / NCINT) * ol[m].NW[LNA_base[i][j] - 1][l];
							}
						}
					}
					ct += 1;
				}
			}
		}

		//Observe the energy input and output (d(E)=sum(d(disp)*F)).
		//Energy released (positive)/gained(negative) by fluid (energy_sent) and energy gained(positive)/released(negative) by structure (energy_rec). 
		t.energy_sent = 0.0;
		for (m = 0; m < owsfnumber; m++) {
			for (l = 0; l < ol[m].FSNEL; l++) {
				for (i = 0; i < NINT; i++) {
					for (j = 0; j < NINT; j++) {
						DISPTEMP = (ol[m].DISP[IEN_3D[ol[m].FP[i*NINT + j] - 1][ol[m].GIDF[l] - 1] - 1][ol[m].dir][1] -
							ol[m].DISP[IEN_3D[ol[m].FP[i*NINT + j] - 1][ol[m].GIDF[l] - 1] - 1][ol[m].dir][0]); //DISP is defined on NNODE 
							//DISPTEMP = ol[m].DISP[IEN_3D[ol[m].FP[i*NINT + j] - 1][ol[m].GIDF[l] - 1] - 1][ol[m].dir][1]; //DISP is defined on NNODE 
						t.energy_sent += ol[m].BF1[ol[m].FP[i*NINT + j] - 1][l] * DISPTEMP * ol[m].NORM;
					}
				}
			}
		}

		//total received energy by structural nodes
		t.energy_rec = 0.0;
		for (m = 0; m < owsfnumber; m++) {
			for (i = 0; i < ol[m].GIDNct_st; i++) {
				DISPTEMP = wsflist[m]->nodecoord[3 * i + ol[m].dir] - ol[m].nodecoord_mpcci[3 * i + ol[m].dir];
				//DISPTEMP = wsflist[m]->nodecoord[3 * i + ol[m].dir]- ol[m].location;
				t.energy_rec += wsflist[m]->nodeforce[3 * i + ol[m].dir] * DISPTEMP * ol[m].NORM;
			}
		}
		std::cout << " " << std::endl;
		break;
	}

	return t;
}

