//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>


//NRB determines the NRB local node numbering and the associated NRB arrays
struct NRBstruct NRB(int NNODE, double **GCOORD, int*** LNA) {
	NRBstruct t;
	int i, j, k, l, m, n, e, z, o;
	extern OWETSURF ol[owsfnumber]; //defined in FSILINK 
	extern NRBSURF nr[nrbsurfnumber];
	
	t.ADMASTERG = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		t.ADMASTERG[i] = 0.0;
	}

	if (element_type == 0) {
		LOBATTOstruct bq;
		bq = LOBATTO(Nq); //one unit higher than the interpolation order
		GLLQUADstruct gq;
		gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
		LOCAL_SHAPEstruct ls;
		ls = LOCAL_SHAPE(LNA, N, Nq, FEM); //one order higher than the interpolation order
		//Construct the ADMASTER for the integration (insert one additional GLL quadrature point)
		//Can we use just one surface for NRBC? I think so  
		//We need to combine the physical group corresponding to NRBC first
		//Do we need to individual points on the NRBC? 
		//ADMASTER is for each element. 
		//We need NRBA to store all the nodes on NRB
		//We also need the total number of node in each NRB surface
		//We need to lump ADMASTER to get ADMASTERG 

		for (z = 0; z < nrbsurfnumber; z++) {
			nr[z].ADMASTER = new double**[nr[z].NEL_nrb];
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				nr[z].ADMASTER[i] = new double*[NINT*NINT];
				for (j = 0; j < NINT*NINT; j++) {
					nr[z].ADMASTER[i][j] = new double[NINT*NINT];
				}
			}
			for (i = 0; i < nr[z].NEL_nrb; i++) {
				for (j = 0; j < NINT*NINT; j++) {
					for (k = 0; k < NINT*NINT; k++) {
						nr[z].ADMASTER[i][j][k] = 0.0;
					}
				}
			}

			double DUNC;
			//Check out Evernote "Some thoughts on boundary force integration" 
			for (e = 0; e < nr[z].NEL_nrb; e++) {
				for (m = 0; m < NINT; m++) {
					for (n = 0; n < NINT; n++) {
						for (l = 0; l < NINT; l++) {
							for (o = 0; o < NINT; o++) {
								DUNC = 0.0; //accumulator 
								for (i = 0; i < NqINT; i++) {
									for (j = 0; j < NqINT; j++) {
										DUNC += gq.W[i] * gq.W[j] * (ls.SHOD[0][m][i] * ls.SHOD[0][n][j]) * (ls.SHOD[0][l][i] * ls.SHOD[0][o][j]) * nr[z].Jacob_2D[e][i*NqINT + j];
										//FUNC += gq.W[i] * gq.W[j] * (ls_ln.SHOD[0][m][i] * ls_ln.SHOD[0][n][j]) * (ls.SHOD[0][l][i] * ls.SHOD[0][o][j]) * (XHE / 2.0)*(YHE / 2.0);
									}
								}
								//ol[z].FPMASTER[e][ol[z].LNA_2D[m][n] - 1][ol[z].LNA_2D[l][o] - 1] += FUNC;
								nr[z].ADMASTER[e][m*NINT + n][l*NINT + o] += DUNC;
							}
						}
					}
				}
			}

			for (e = 0; e < nr[z].NEL_nrb; e++) {
				for (j = 0; j < NINT*NINT; j++) {
					for (k = 0; k < NINT*NINT; k++) {
						t.ADMASTERG[nr[z].IEN_gb[nr[z].DP_2D[k] - 1][e] - 1] += nr[z].ADMASTER[e][j][k];
					}
				}
			}

		}
	}
	if (element_type == 1) {

		if (mappingalgo == 4) {
			for (z = 0; z < nrbsurfnumber; z++) {
				nr[z].ADMASTER = new double**[nr[z].NEL_nrb];
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					nr[z].ADMASTER[i] = new double*[3];
					for (j = 0; j < 3; j++) {
						nr[z].ADMASTER[i][j] = new double[3];
					}
				}
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					for (j = 0; j < 3; j++) {
						for (k = 0; k < 3; k++) {
							nr[z].ADMASTER[i][j][k] = 0.0;
						}
					}
				}
				//Check out the Goodnote "Bounary force integration" or Evernote "Tetrahedral average integration"
				for (e = 0; e < nr[z].NEL_nrb; e++) {
					for (m = 0; m < 3; m++) {
						for (n = 0; n < 3; n++) {
							nr[z].ADMASTER[e][m][n] = (1.0 / 9.0) * nr[z].dimension[e];
						}
					}
				}
				for (e = 0; e < nr[z].NEL_nrb; e++) {
					for (j = 0; j < 3; j++) {
						for (k = 0; k < 3; k++) {
							t.ADMASTERG[nr[z].IEN_gb[k][e] - 1] += nr[z].ADMASTER[e][j][k];
						}
					}
				}
			}
		}

		if (mappingalgo == 5) {
			double xi[3]; double eta[3]; double phi[3][3];
			xi[0] = 1.0 / 6.0; xi[1] = 2.0 / 3.0; xi[2] = 1.0 / 6.0;
			eta[0] = 1.0 / 6.0; eta[1] = 1.0 / 6.0; eta[2] = 2.0 / 3.0;
			double w[3]; //integration weight
			w[0] = 1.0 / 3.0; w[1] = 1.0 / 3.0; w[2] = 1.0 / 3.0;
			for (z = 0; z < owsfnumber; z++) {
				nr[z].ADMASTER = new double**[nr[z].NEL_nrb];
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					nr[z].ADMASTER[i] = new double*[3];
					for (j = 0; j < 3; j++) {
						nr[z].ADMASTER[i][j] = new double[3];
					}
				}
				for (i = 0; i < nr[z].NEL_nrb; i++) {
					for (j = 0; j < 3; j++) {
						for (k = 0; k < 3; k++) {
							nr[z].ADMASTER[i][j][k] = 0.0;
						}
					}
				}
				for (n = 0; n < 3; n++) {
					phi[0][n] = 1 - xi[n] - eta[n];
					phi[1][n] = xi[n];
					phi[2][n] = eta[n];
				}
				for (e = 0; e < nr[z].NEL_nrb; e++) {
					for (m = 0; m < 3; m++) {
						for (n = 0; n < 3; n++) {
							for (o = 0; o < 3; o++) { //o is the 3-point quadrature point
								nr[z].ADMASTER[e][m][n] += (1.0 / 2.0) * w[o] * phi[m][o] * phi[n][o] * (2 * ol[z].dimension[e]);
							}
						}
					}
				}
			}
			for (e = 0; e < nr[z].NEL_nrb; e++) {
				for (j = 0; j < 3; j++) {
					for (k = 0; k < 3; k++) {
						t.ADMASTERG[nr[z].IEN_gb[k][e] - 1] += nr[z].ADMASTER[e][j][k];
					}
				}
			}
		}

	}

	//Obtain NRB node numbering for all nodes from ol[z].NRBA
	//The goal is to remove the repetitive nodes on the intersection of different NRB surfaces
	int* holder;
	holder = new int[NNODE];
	for (i = 0; i < NNODE; i++) {
		holder[i] = 0.0;
	}
	for (z = 0; z < nrbsurfnumber; z++) {
		for (j = 0; j < nr[z].NRBNODE; j++) {
			holder[nr[z].NRBA[j] - 1] += 1;  //If the final is >1, it means the node shows up in different surfaces (i.e., on the intersection)
		}
	}
	//Count the total number of node on NRBC
	t.NNODE_nrb = 0;
	for (i = 0; i < NNODE; i++) {
		if (holder[i] != 0) {
			t.NNODE_nrb += 1;
		}
	}
	int ct = 0;
	//int* NRBA_t; //node numbering for all NRB nodes
	t.NRBA_t = new int[t.NNODE_nrb];
	for (i = 0; i < NNODE; i++) {
		if (holder[i] != 0) { //The nodes showed up in ol[z].NRBA, regardless how many times, are included in NRBA_t. This way, the same node is not counted twice.
			t.NRBA_t[ct] = i + 1;
			ct += 1; 
		}
	}
	if (ct != t.NNODE_nrb) {
		std::cout << "The NRBA_t is wrong" << std::endl;
		system("PAUSE ");
	}
	std::cout << " " << std::endl;
	//The node counting routine is verified by FSP code

	std::cout << " " << std::endl;
	return t;
}