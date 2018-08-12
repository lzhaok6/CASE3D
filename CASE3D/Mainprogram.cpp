#include "header.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include "data.h"
#include "adapter.h"
#include "mpcci.h"
#include "vector"
#include <string>
#include <sstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <stdlib.h>
#include <omp.h>

/*
This code dedicated for frigate currently only support mapping algorithm 2.  
*/

double time_step_size;
double current_time;
OWETSURF ol[owsfnumber]; //declear the data structure globally
NRBSURF nr[nrbsurfnumber]; 
STRU_WET_SURF ss;
int main()
{	
	double LMAX;
	int h, q, z, e, i, j, k, ii, jj; //error prone: cannot be the same with data structure type (z is the same as Z)
	int TIME = 0;     //CONTROL TIME
	if (FEM == 1 && N != 1) {
		std::cout << "the current code doesn't support high-order FEM" << std::endl;
		system("PAUSE ");
	}
	if (FEM == 1 && Nq != N && (mappingalgo == 1 || mappingalgo == 5)) {
		std::cout << "FEM doesn't need more integration point" << std::endl;
		system("PAUSE ");
	}
	if (FEM == 0 && Nq != N + 1 && (mappingalgo == 1 || mappingalgo == 5)) {
		std::cout << "Are you sure you don't want to use full integration for the boundary forces?" << std::endl;
		system("PAUSE ");
	}
	if (nodeforcemap2 == 0 && mappingalgo != 2) {
		std::cout << "the mapping algorithms other than 2 are not designed to transfer abosolute pressure" << std::endl;
		system("PAUSE ");
	}
	std::cout << "Have you properly configured owsfnumber and nrbsurfnumber?" << std::endl;
	if (tensorfactorization == 1 & element_type == 1) {
		std::cout << "tetrahedral element does not support tensor-product factorization" << std::endl;
		system("PAUSE ");
	}
	if (element_type == 1 && N > 1) {
		std::cout << "The high-order tetrahedral element is not supported" << std::endl;
		system("PAUSE ");
	}
	if (element_type == 1 && FEM != 1) {
		std::cout << "The tet element is only FEM type" << std::endl;
		system("PAUSE ");
	}
	if (element_type == 0 && mappingalgo != 2) {
		std::cout << "The other mapping algorithm is not available for hexahedral element" << std::endl;
		system("PAUSE ");
	}
	if (element_type == 1 && mappingalgo != 4) {
		std::cout << "The other mapping algorithm is not available for tetrahedral element" << std::endl;
		system("PAUSE ");
	}
	if (FEM == 1 && Nq != N) {
		std::cout << "FEM (Gauss-Legendre integration) does not need extra integration nodes" << std::endl;
		system("PAUSE ");
	}
	if (FEM == 0 && Nq == N) {
		std::cout << "SEM (GLL) needs extra integration nodes" << std::endl;
		system("PAUSE ");
	}

	meshgenerationstruct a;
	a = meshgeneration();
	std::cout << "mesh generation done" << std::endl;

	LOBATTOstruct b;
	b = LOBATTO(N);
	std::cout << "LOBATTO() done" << std::endl;

	LOCAL_NODEstruct c;
	c = LOCAL_NODE(N);
	std::cout << "LOCAL_NODE(N) done" << std::endl;

	//Neighborhood_search(a.GCOORD, c.LNA, a.IEN, a.NEL);


	//double** GCOORD, int***LNA, int**IEN_flu, int NEL_flu
	//---------------SHAPE FUNCTION ROUTINE------------------------------//
	//DETERMINE GLL QUADRATURE POINTS AND WEIGHTS
	GLLQUADstruct f;
	f = GLLQUAD(b.Z, b.WL, N, !FEM);
	//std::cout << "GLLQUAD done" << std::endl;
	LOCAL_SHAPEstruct g;
	g = LOCAL_SHAPE(c.LNA, N, N, FEM);
	//std::cout << "LOCAL_SHAPE done" << std::endl;
	//DETERMINE LOCAL GEOMETERY SHAPE FUNCTION AT ELEMENT NODES
	LOCAL_GSHAPEstruct l;
	l = LOCAL_GSHAPE(f.S, c.LNA, NINT);
	//std::cout << "LOCAL_GSHAPE done" << std::endl;
	//DETERMINE JACOBIAN MATRIX 
	JACOBIANstruct m;
	//m = JACOBIAN(a.NEL, a.GCOORD, a.IEN, c.LNA, NINT, f.S);
	m = JACOBIAN(a.NEL, a.GCOORD, a.IEN, c.LNA);
	/*
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 8; j++) {
			for (k = 0; k < NINT; k++) {
				for (h = 0; h < NINT; h++) {
					delete[] l.GSHL[i][j][k][h];
				}
				delete[] l.GSHL[i][j][k];
			}
			delete[] l.GSHL[i][j];
		}
		delete[] l.GSHL[i];
	}
	delete[] l.GSHL;
	*/
	//std::cout << "JACOBIAN done" << std::endl;
	GLOBAL_SHAPEstruct n;
	n = GLOBAL_SHAPE(a.NEL, g.SHL, m.XS, m.JACOB, a.GCOORD, a.IEN, m.JACOB_tet);
	//std::cout << "GLOBAL_SHAPE done" << std::endl;
	MATRIXstruct o;
	o = MATRIX(a.NEL, a.NNODE, g.SHL, f.W, a.IEN, c.LNA, m.XS, n.SHG, m.JACOB, m.JACOB_tet, n.SHG_tet);
	//duration_int = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//std::cout << "time lapse for until this point: " << duration_int << std::endl;

	//std::cout << "MATRIX done" << std::endl;
	//DETERMINE MAXIMUM MESH EIGENVALUE TO FIND CFL TIMESTEP
	LMAX = EIGENMAX(o.QMASTER, o.HMASTER, a.NEL);

	//Derive the integration weight of the wetted surface elements
	FSILINK(c.LNA);

	//read the model file to MpCCI adapter
	char* modelfile = "model.txt";
	readfile(modelfile);

	TIMINTstruct t;
	t = TIMINT(LMAX);

	time_step_size = t.DT;

	double *T; //TIMESTEP ARRAY 
	T = new double[t.NDT + 1];
	for (i = 0; i < t.NDT + 1; i++) {
		T[i] = 0.0;
	}

	//DETERMINE TIMESTEP ARRAY T
	T[0] = 0.0;
	for (i = 1; i < t.NDT; i++) {
		T[i] = T[i - 1] + t.DT;
	}
	//later added
	T[t.NDT] = TTERM;

	if (element_type == 0) { //hexahedral element
		//define mapping functions for interface mapping
		double* basep; //points on base fluid mesh 
		double* origp; //points on original fluid mesh 
		double nomx, nomy, nomz;
		double denomx, denomy, denomz;
		basep = new double[NCINT];
		for (i = 0; i < NCINT; i++) {
			basep[i] = 0.0;
		}
		for (i = 0; i < NCINT; i++) {
			basep[i] = -1.0 + i*(2.0 / (NCINT - 1));
		}

		origp = new double[hpref + 1];
		for (i = 0; i < refine + 1; i++) {
			origp[i*N] = -1 + (2.0 / refine)*i;
		}
		for (i = 0; i < refine; i++) {
			for (j = 1; j < N; j++) {
				origp[i*N + j] = origp[i*N] + (2.0 / refine)*((b.Z[j - 1] + 1) / 2);
			}
		}

		for (z = 0; z < owsfnumber; z++) {
			ol[z].phi_fem = new double**[4];
			for (i = 0; i < 4; i++) { //for all the points in that element
				ol[z].phi_fem[i] = new double*[NINT];
				for (j = 0; j < NINT; j++) {
					ol[z].phi_fem[i][j] = new double[NINT];
				}
			}
			for (i = 0; i < 4; i++) {
				for (j = 0; j < NINT; j++) {
					for (k = 0; k < NINT; k++) {
						ol[z].phi_fem[i][j][k] = 0.0;
					}
				}
			}
			for (h = 0; h < 2; h++) { //stands for every fem point, LNA[u][v][w](shape function is based on those points)
				for (k = 0; k < 2; k++) {
					for (i = 0; i < NINT; i++) {  //i j k are the independent variable in basis function(sem points)
						for (j = 0; j < NINT; j++) {
							nomx = 1.0; nomy = 1.0; //multiplier initialization
							denomx = 1.0; denomy = 1.0; //multiplier initialization
							for (e = 0; e < 2; e++) { //loop through nominator and denominator in basis function expression
								if (e != h) {
									nomx *= (origp[i] - basep[e]);
									denomx *= (basep[h] - basep[e]);
								}
								if (e != k) {
									nomy *= (origp[j] - basep[e]);
									denomy *= (basep[k] - basep[e]);
								}
							}
							ol[z].phi_fem[ol[z].LNA_algo2[h][k] - 1][i][j] = (nomx / denomx)*(nomy / denomy); //tensor product
							//the coordinate definition dof u,v is the same with i,j
						}
					}
				}
			}
		}
	}

	//derive the peak wave pressure and decay rate
	double XC = 0.0; double YC = 0.0; double ZC = 0.0;
	double XO = 0.0; double YO = 0.0; double ZO = 0.0; //stand-off point
	double dist = 0.0;
	double dists = 0.0; //distance from charge center to nearest structural point
	double distf = 0.0; //distance from charge center to nearest freesurface point
	//determine the explosion center
	ZC = stdoff*0.3048 + SZ / 2;
	YC = -depth*0.3048;
	XC = 0.0;

	//determine the stand-off point (nearest structural node or free surface node depending on which one is closer)
	dists = sqrt(pow(XC - 0.0, 2) + pow(YC - (-SY), 2) + pow(ZC - SZ / 2, 2));
	distf = -YC;
	
	if (dists > distf) {
		XO = XC;
		YO = 0;
		ZO = ZC;
	}
	else {
		XO = 0.0;
		YO = -SY;
		ZO = SZ / 2;
	}

	dist = sqrt(pow((XC - XO), 2) + pow((YC - YO), 2) + pow((ZC - ZO), 2));
	if (Bleich == 1 && WAVE == 2) {
		dist = depth;
		ZC = 0;
		YC = -depth;
		XC = -SX / 2;
	}

	double A = 0.0;
	double B = 0.0;
	double K = 0.0;
	if (TNT == 1) { //TNT (1.60 g/cc) Ref: Geers and Hunter (2002)
		A = 0.13;
		B = 0.23;
		K = 52.4e6;
	}
	else { //HBX-1 (1.72 g/cc) Ref: Ref: Geers and Hunter (2002)
		A = 0.15;
		B = 0.29;
		K = 56.7e6;
	}
	double PPEAK = K*pow(dist / pow(W*0.453592, 1.0 / 3.0), -(1 + A)); //pa
	double TAU = pow(W*0.453592, 1.0 / 3.0)*0.084*pow(pow(W*0.453592, 1.0 / 3.0) / dist, -B) / 1000; //sec
	
	double KAPPA = 0.0;

	if (Bleich == 1 && WAVE == 1) {
		PPEAK = 0.712e6;
		TAU = 0.999e-3;
	}

	TIME_INT(a.NNODE, a.GCOORD, c.LNA, a.IEN, a.NEL, TIME, T, t.DT, t.NDT, o.HMASTER, o.Q, KAPPA, PPEAK, TAU, XC, YC, ZC, XO, YO, ZO, g.SHOD, o.gamman, o.gamma_tn, o.Gn, n.SHG,
		o.gamma_t, o.gamma, o.G);

	printf("Cleaning up...\n");
	/* clean up */
	cleanall();
	printf("Bye-bye!\n");
	exitcoupling();
	exit(0);
}
