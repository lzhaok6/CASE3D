//#include "stdafx.h"
#include "header.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>


//NRB determines the NRB local node numbering and the associated NRB arrays
struct NRBstruct NRB(int NNODE, double **GCOORD, double* W, int*** LNA, int**IEN, int NEL, double***SHL, double***SHOD) {
	NRBstruct t;
	//int* NRBA;
	int i, j, k, l, m;
	int count;
	double DUNC;
	//int NRBELE;
	//int NRBELE1 = 0; int NRBELE2 = 0; int NRBELE3 = 0; int NRBELE4 = 0;
	//int *NRBELE_ARR; //the NRB element numbering
	extern int owsfnumber;
	extern OWETSURF ol[4]; //defined in FSILINK 
	t.NRBNODE_loc = new int[owsfnumber];
	for (i = 0; i < owsfnumber + 1; i++) {
		t.NRBNODE_loc[i] = 0;
	}

	int XNNODE = N*round((AX + SX) / XHE) + 1; int XNEL = round((AX + SX) / XHE);
	int ZNNODE = N*round((2 * BZ + SZ) / ZHE) + 1; int ZNEL = round((2 * BZ + SZ) / ZHE);
	int YNNODE = N*round((SY + DY) / YHE) + 1; int YNEL = round((SY + DY) / YHE);


	if (WAVE == 1) { //plane wave
		t.NRBNODE = XNNODE*ZNNODE;
	}
	else { //spherical wave
		t.NRBNODE = XNNODE*ZNNODE + //bottom
			XNNODE*(YNNODE - 1) + //back
			(YNNODE - 1)*(ZNNODE - 1) + //left
			(XNNODE - 1)*(YNNODE - 1); //front
		//note the difference between NRBNODE[0] in plane wave and spherical wave case
		t.NRBNODE_loc[0] = XNNODE*ZNNODE; //bottom
		t.NRBNODE_loc[1] = XNNODE*YNNODE; //back
		t.NRBNODE_loc[2] = YNNODE*ZNNODE; //left 
		t.NRBNODE_loc[3] = XNNODE*YNNODE; //front 
	}
	//checkpoint: If this is the correct spherical NRB node number
	
	t.AD = new double*[NNODE];
	t.NRBA = new int[t.NRBNODE]; //the first dimension is initialized to be bigger

	for (i = 0; i < t.NRBNODE; i++) {
		t.NRBA[i] = 0;
	}
	
	for (i = 0; i < NNODE; i++) {
		t.AD[i] = new double[owsfnumber + 1]; //the first dimension is for the total node, the rest of the dimension is for the other surfaces
	}
	for (i = 0; i < NNODE; i++) {
		for (j = 0; j < owsfnumber + 1; j++) {
			t.AD[i][j] = 0.0;
		}
	}
	/*
	t.nrbele = new int[owsfnumber + 1];
	for (i = 0; i < owsfnumber + 1; i++) {
		t.nrbele[i] = 0;
	}
	*/
	t.nrbele = new int[owsfnumber];
	for (i = 0; i < owsfnumber; i++) {
		t.nrbele[i] = 0;
	}
	
	
	for (i = 0; i < owsfnumber; i++) {
		ol[i].DP = new int[NINT*NINT];
	}
	for (i = 0; i < owsfnumber; i++) {
		for (j = 0; j < NINT*NINT; j++) {
			ol[i].DP[j] = 0.0;
		}
	}
	count = 0;
	for (m = 0; m < NINT; m++) {
		for (l = 0; l < NINT; l++) {
			ol[0].DP[count] = LNA[m][0][l]; //bottom
			ol[1].DP[count] = LNA[m][l][0]; //back
			ol[2].DP[count] = LNA[0][m][l]; //left 
			ol[3].DP[count] = LNA[m][l][N]; //front
			count += 1;
		}
	}

	LOBATTOstruct bq;
	bq = LOBATTO(Nq); //one unit higher than the interpolation order
	GLLQUADstruct gq;
	gq = GLLQUAD(bq.Z, bq.WL, Nq, !FEM);
	LOCAL_SHAPEstruct ls;
	ls = LOCAL_SHAPE(LNA, Nq); //one order higher than the interpolation order

	for (i = 0; i < owsfnumber + 1; i++) {
		ol[i].ADMASTER = new double*[NINT*NINT];
		for (j = 0; j < NINT*NINT; j++) {
			ol[i].ADMASTER[j] = new double[NINT*NINT];
		}
	}
	for (i = 0; i < owsfnumber + 1; i++) {
		for (j = 0; j < NINT*NINT; j++) {
			for (k = 0; k < NINT*NINT; k++) {
				ol[i].ADMASTER[j][k] = 0.0;
			}
		}
	}

	if (WAVE == 1) { //if WAVE=1, only bottom fluid domain surface should be NRB 
		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (abs(GCOORD[i][1] + (SY+DY)) < 1e-6) { //error prone!!! could be something like 0.0000001 for boundary nodes (check the equality condition in C++)
				t.NRBA[count] = i + 1; //NRBA node numbering (not sequentially)
				count += 1;
			} //1e-5 is the epsilon
		}
		if (count != t.NRBNODE) {
			std::cout << "total NRB node recognization is wrong!" << std::endl;
		}

		if (AX > 0) {
			t.NRBELE = round(((AX + SX) / XHE)*((2 * BZ + SZ) / ZHE));
		}
		else { //Bleich-Sandler set up 
			t.NRBELE = SXNEL*SZNEL;
		}

		t.NRBELE_ARR = new int*[t.NRBELE];
		for (i = 0; i < t.NRBELE; i++) {
			t.NRBELE_ARR[i] = new int[1];
		}
		for (i = 0; i < t.NRBELE; i++) {
			for (j = 0; j < 1; j++) {
				t.NRBELE_ARR[i][j] = 0;
			}
		}

		count = 0;
		for (i = 0; i < NEL; i++) {
			if (abs(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][1] + (SY + DY)) < 1e-6) {
				t.NRBELE_ARR[count][0] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBELE) {
			std::cout << "NRBELE computation is wrong" << std::endl;
		}
		
		for (i = 0; i < t.NRBELE; i++) {
			for (m = 0; m < NINT*NINT; m++) {
				for (l = 0; l < NINT*NINT; l++) {
					DUNC = 0.0; //accumulator 
					if (m == l) {
						for (j = 0; j < NINT; j++) {
							for (k = 0; k < NINT; k++) {
								DUNC += W[j] * W[k] * (SHL[3][ol[0].DP[m] - 1][j*NINT*NINT+0*NINT+k] / SHOD[0][0][0]) * (SHL[3][ol[0].DP[l] - 1][j*NINT*NINT + 0 * NINT + k] / SHOD[0][0][0]) * (XHE / 2.0)*(ZHE / 2.0);
							}
						}
						t.AD[IEN[ol[0].DP[m] - 1][t.NRBELE_ARR[i][0] - 1] - 1][0] += DUNC;
						if (IEN[ol[0].DP[m] - 1][t.NRBELE_ARR[i][0] - 1]>NNODE) {
							std::cout << "AD out of bound" << std::endl;
						}
					}
				}
			}
		}

		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (j = 0; j < NqINT; j++) {
					for (k = 0; k < NqINT; k++) {
						DUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[0].DP[m] - 1][j*NqINT*NqINT+0*NqINT+k] / SHOD[0][0][0]) * (ls.SHL[3][ol[0].DP[l] - 1][j*NqINT*NqINT + 0 * NqINT + k] / SHOD[0][0][0]) * (XHE / 2.0)*(ZHE / 2.0);
					}
				}
				ol[0].ADMASTER[m][l] += DUNC;
			}
		}
	}

	else { //WAVE==2
		t.nloc = new int*[t.NRBNODE];
		for (i = 0; i < t.NRBNODE;i++) {
			t.nloc[i] = new int[owsfnumber];
		}

		/*
		t.NRBELE_ARR = new int*[NEL];
		for (i = 0; i < NEL; i++) {
			t.NRBELE_ARR[i] = new int[5];
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < 5; j++) {
				t.NRBELE_ARR[i][j] = 0;
			}
		}
		*/
		t.NRBELE_ARR = new int*[NEL];
		for (i = 0; i < NEL; i++) {
			t.NRBELE_ARR[i] = new int[4];
		}
		for (i = 0; i < NEL; i++) {
			for (j = 0; j < 4; j++) {
				t.NRBELE_ARR[i][j] = 0;
			}
		}
		count = 0;
		for (i = 0; i < NNODE; i++) {
			if (abs(GCOORD[i][1] + (SY + DY)) < 1e-6 //bottom
				|| abs(GCOORD[i][2] + SZ / 2 + BZ) < 1e-6     //back
				|| abs(GCOORD[i][2] - SZ / 2 - BZ) < 1e-6 //front
				|| abs(GCOORD[i][0] + (SX + AX)) < 1e-6) {  //left
				t.NRBA[count] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBNODE) {
			std::cout << "total NRBA is wrong!" << std::endl;
		}

		count = 0;
		for (i = 0; i < t.NRBNODE; i++) {
			if (abs(GCOORD[t.NRBA[i] - 1][1] + (SY + DY)) < 1e-6) {  //bottom 
				t.nloc[count][0] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBNODE_loc[0]) {
			std::cout << "bottom NRBA is wrong!" << std::endl;
		}
		count = 0;
		for (i = 0; i < t.NRBNODE; i++) {
			if (abs(GCOORD[t.NRBA[i] - 1][2] + SZ / 2 + BZ) < 1e-6) { //back
				t.nloc[count][1] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBNODE_loc[1]) {
			std::cout << "back NRBA is wrong!" << std::endl;
		}
		count = 0;
		for (i = 0; i < t.NRBNODE; i++) {
			if (abs(GCOORD[t.NRBA[i] - 1][0] + (SX + AX)) < 1e-6) { //left
				t.nloc[count][2] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBNODE_loc[2]) {
			std::cout << "left NRBA is wrong!" << std::endl;
		}
		count = 0;
		for (i = 0; i < t.NRBNODE; i++) {
			if (abs(GCOORD[t.NRBA[i] - 1][2] - SZ / 2 - BZ) < 1e-6) {//front
				t.nloc[count][3] = i + 1;
				count += 1;
			}
		}
		if (count != t.NRBNODE_loc[3]) {
			std::cout << "front NRBA is wrong!" << std::endl;
		}

		count = 0;
		int ele = 0;
		//NRBELE1 = (AXNEL + BXNEL + CXNEL)*AZNEL;
		//Note that the numbering scheme is different from wet surface, in which the back surface is the first one, and then the left, bottom, front, right
		//bottom (without adjacent point with L V R) 
		for (i = 0; i < NEL; i++) {
			if (abs(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][1] + (DY + SY)) < 1e-5) {
				for (m = 0; m < NINT*NINT; m++) {
					for (l = 0; l < NINT*NINT; l++) {
						DUNC = 0.0;
						if (m == l) {
							for (j = 0; j < NINT; j++) {
								for (k = 0; k < NINT; k++) {
									DUNC += W[j] * W[k] * (SHL[3][ol[0].DP[m] - 1][j*NINT*NINT+0*NINT+k] / SHOD[0][0][0]) * (SHL[3][ol[0].DP[l] - 1][j*NINT*NINT + 0 * NINT + k] / SHOD[0][0][0]) * (XHE / 2.0)*(ZHE / 2.0);
								}
							}
							t.AD[IEN[ol[0].DP[m] - 1][i] - 1][1] += DUNC;
							count += 1;
						}
					}
				}
				//t.NRBELE_ARR[ele][1] = i + 1;
				t.NRBELE_ARR[ele][0] = i + 1;
				ele += 1;
			}
		}
		t.nrbele[0] = ele;
		if (ele != XNEL*ZNEL) {
			std::cout << "bottom AD is wrong!" << std::endl;
		}

		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (j = 0; j < NqINT; j++) {
					for (k = 0; k < NqINT; k++) {
						DUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[0].DP[m] - 1][j*NqINT*NqINT+0*NqINT+k] / SHOD[0][0][0]) * (ls.SHL[3][ol[0].DP[l] - 1][j*NqINT*NqINT + 0 * NqINT + k] / SHOD[0][0][0]) * (XHE / 2.0)*(ZHE / 2.0);
					}
				}
				ol[0].ADMASTER[m][l] += DUNC;
			}
		}


		//back (without point shared with L and R NRB)
		count = 0;
		ele = 0;
		for (i = 0; i < NEL; i++) {
			if (abs(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][2] + BZ + SZ / 2) < 1e-5) {
				for (m = 0; m < NINT*NINT; m++) {
					for (l = 0; l < NINT*NINT; l++) {
						DUNC = 0.0;
						if (m == l) {
							for (j = 0; j < NINT; j++) {
								for (k = 0; k < NINT; k++) {
									DUNC += W[j] * W[k] * (SHL[3][ol[1].DP[m] - 1][j*NINT*NINT+k*NINT+0] / SHOD[0][0][0]) * (SHL[3][ol[1].DP[l] - 1][j*NINT*NINT + k*NINT + 0] / SHOD[0][0][0]) * (XHE / 2.0)*(YHE / 2.0);
								}
							}
							t.AD[IEN[ol[1].DP[m] - 1][i] - 1][2] += DUNC;
							count += 1;
						}
					}
				}
				//t.NRBELE_ARR[ele][2] = i + 1;
				t.NRBELE_ARR[ele][1] = i + 1;
				ele += 1;
			}
		}
		t.nrbele[1] = ele;
		if (ele != XNEL*YNEL) {
			std::cout << "back AD is wrong!" << std::endl;
		}

		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (j = 0; j < NqINT; j++) {
					for (k = 0; k < NqINT; k++) {
						DUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[1].DP[m] - 1][j*NqINT*NqINT+k*NqINT+0] / SHOD[0][0][0]) * (ls.SHL[3][ol[1].DP[l] - 1][j*NqINT*NqINT + k*NqINT + 0] / SHOD[0][0][0]) * (XHE / 2.0)*(YHE / 2.0);
					}
				}
				ol[1].ADMASTER[m][l] += DUNC;
			}
		}

		//left (full point)
		count = 0;
		ele = 0;
		for (i = 0; i < NEL; i++) {
			if (abs(GCOORD[IEN[LNA[0][0][0] - 1][i] - 1][0] + (SX + AX)) < 1e-5) {
				for (m = 0; m < NINT*NINT; m++) {
					for (l = 0; l < NINT*NINT; l++) {
						DUNC = 0.0;
						if (m == l) {
							for (j = 0; j < NINT; j++) {
								for (k = 0; k < NINT; k++) {
									DUNC += W[j] * W[k] * (SHL[3][ol[2].DP[m] - 1][0*NINT*NINT+j*NINT+k] / SHOD[0][0][0]) * (SHL[3][ol[2].DP[l] - 1][0 * NINT*NINT + j*NINT + k] / SHOD[0][0][0]) * (YHE / 2.0)*(ZHE / 2.0);
								}
							}
							//if (abs(GCOORD[IEN[ol[2].DP[m] - 1][i] - 1][0] + (SX / 2 + AX)) < 1e-5) {
							t.AD[IEN[ol[2].DP[m] - 1][i] - 1][3] += DUNC;
							count += 1;
							//}
						}
					}
				}
				//t.NRBELE_ARR[ele][3] = i + 1;
				t.NRBELE_ARR[ele][2] = i + 1;
				ele += 1;
			}
		}
		t.nrbele[2] = ele;
		if (ele != YNEL*ZNEL) {
			std::cout << "left AD is wrong!" << std::endl;
		}

		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (j = 0; j < NqINT; j++) {
					for (k = 0; k < NqINT; k++) {
						DUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[2].DP[m] - 1][0 * NqINT*NqINT + j*NqINT + k] / SHOD[0][0][0]) * (ls.SHL[3][ol[2].DP[l] - 1][0 * NqINT*NqINT + j*NqINT + k] / SHOD[0][0][0]) * (YHE / 2.0)*(ZHE / 2.0);
					}
				}
				ol[2].ADMASTER[m][l] += DUNC;
			}
		}

		//front
		count = 0;
		ele = 0;
		for (i = 0; i < NEL; i++) {
			if (abs(GCOORD[IEN[LNA[0][0][N] - 1][i] - 1][2] - SZ / 2 - BZ) < 1e-5) {
				for (m = 0; m < NINT*NINT; m++) {
					for (l = 0; l < NINT*NINT; l++) {
						DUNC = 0.0;
						if (m == l) {
							for (j = 0; j < NINT; j++) {
								for (k = 0; k < NINT; k++) {
									DUNC += W[j] * W[k] * (SHL[3][ol[3].DP[m] - 1][j*NINT*NINT+k*NINT+N] / SHOD[0][N][N]) * (SHL[3][ol[3].DP[l] - 1][j*NINT*NINT + k*NINT + N] / SHOD[0][N][N]) * (XHE / 2.0)*(YHE / 2.0);
								}
							}
							t.AD[IEN[ol[3].DP[m] - 1][i] - 1][4] += DUNC;
							count += 1;
						}
					}
				}
				//t.NRBELE_ARR[ele][4] = i + 1;
				t.NRBELE_ARR[ele][3] = i + 1;
				ele += 1;
			}
		}
		t.nrbele[3] = ele;
		if (ele != YNEL*XNEL) {
			std::cout << "front AD is wrong!" << std::endl;
		}

		for (m = 0; m < NINT*NINT; m++) {
			for (l = 0; l < NINT*NINT; l++) {
				DUNC = 0.0; //accumulator 
				for (j = 0; j < NqINT; j++) {
					for (k = 0; k < NqINT; k++) {
						DUNC += gq.W[j] * gq.W[k] * (ls.SHL[3][ol[3].DP[m] - 1][j*NqINT*NqINT + k*NqINT + Nq] / SHOD[0][N][N]) * (ls.SHL[3][ol[3].DP[l] - 1][j*NqINT*NqINT + k*NqINT + Nq] / SHOD[0][N][N]) * (XHE / 2.0)*(YHE / 2.0);
					}
				}
				ol[3].ADMASTER[m][l] += DUNC;
			}
		}

		for (i = 0; i < NNODE; i++) {
			t.AD[i][0] = t.AD[i][1] + t.AD[i][2] + t.AD[i][3] + t.AD[i][4];
		} //total AD 
	
	}

	t.ADMASTERG = new double[NNODE];
	for (i = 0; i < NNODE; i++) {
		t.ADMASTERG[i] = 0.0;
	}

	for (m = 0; m < owsfnumber; m++) {
		for (i = 0; i < t.nrbele[m]; i++) {
			for (j = 0; j < NINT*NINT; j++) {
				for (k = 0; k < NINT*NINT; k++) {
					t.ADMASTERG[IEN[ol[m].DP[k] - 1][t.NRBELE_ARR[i][m] - 1] - 1] += ol[m].ADMASTER[j][k];
				}
			}
		}
	}

	return t;
}