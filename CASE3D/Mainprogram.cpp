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

/*
**the mesh is for half-symmetric fluid domain (four fluid modules)
this version of code has SEM mesh generation capability. More information could be found in the stand-alone version: mesh_generation_sym.cpp
**This code may not be used to simulate a Bleich-Sandler style 3D barge. The function is redirected to the code "3Dbarge_TFM_BS_Abaqus_EI" in the same folder
**A tabular shock wave is programmed in this code for comparing with Abaqus solution.
Currently, the code only support the standoff point at NRB (xo=AY) for Abaqus wave profile (tabulated).
**This code has fixed the problem of emergence. The free surface in this code is set to be the Y=0 plane,
structure in the Abaqus input file should be moved downward depending on where is the waterline.
For instance, if the waterline in MAESTRO is 0.203, then the structure should be shifted 0.203 downward.
The input file should still be those with _0.203wl token since the wetted surface is changed with waterline location.
**Two mapping algorithms are provided in this code.
The first one uses Farhat's conservative mapping to map the value from SEM mesh
to FEM coupling mesh which is supposedly identical to the structural mesh
The second method breaks the high order SEM mesh into individual 1st order FEM mesh. This method does not generate identical mesh
when order is higher than 2 since the internal SEM nodes are not evenly distributed.
The first and second algorithm should behave the same when N=1.
**The mapping algorithm for Bleich-Sandler case is deleted, please use the code for BS case individually (i.e. 3Dbarge_TFM_BS_Abaqus_EI)
Both TFM and SFM are included in this code as two sub-functions.
Shock properties are automatically calculated in the code based on the stand-off distance and depth of the charge. (Previously, we've been direcly specifying the peak pressure, decay rate, stand-off point and so forth)
11/14/2017 Added the structural displacement predictor for fluid solver. Zhaokuan Lu
(need to fix the displacement predictor in SFM)
11/17/2017 Added the definition for HBX-1 in addition to TNT. Zhaokuan Lu
1/17/2018 Integrate the pure FEM element into the code. Zhaokuan Lu
3/30/2018 Modify the code to take nonuniform element size (Jacobian determinant, Matrix calculation are affected). 
add a comment here
add another comment here
*/

double time_step_size;
double current_time;
OWETSURF ol[4]; //declear the data structure globally
int owsfnumber = 4;
//int output;
int main()
{
	double LMAX;
	int h, i, j, k, q, z, ii, jj; //error prone: cannot be the same with data structure type (z is the same as Z)
	int TIME = 0;     //CONTROL TIME
	//output = 1;
	
	/*
	//ask the user if file output is needed
	std::cout << "Do you want to output the pressure profile? (Y/N): " << std::endl;
	std::string outholder;
	std::getline(std::cin, outholder, '\n');
	if (outholder.compare("Y") != 0) {
		output = 0;
	}
	*/

	if (FEM == 1 && N != 1) {
		std::cout << "the current code doesn't support high-order FEM" << std::endl;
		system("PAUSE ");
	}
	if (FEM == 1 && Nq != N && (mappingalgo == 1 || mappingalgo == 5)) {
		std::cout << "FEM doesn't need more integration point" << std::endl;
		system("PAUSE ");
	}
	if (FEM==0 && Nq!=N+1 && (mappingalgo == 1 || mappingalgo == 5)) {
		std::cout << "Are you sure you don't want to use full integration for the boundary forces?" << std::endl;
		system("PAUSE ");
	}
	if (nodeforcemap2 == 0 && mappingalgo != 2) {
		std::cout << "the mapping algorithms other than 2 are not designed to transfer abosolute pressure" << std::endl;
		system("PAUSE ");
	}
	if (Bleich == 1 && SX == 0.1) {
		std::cout << "Have you divided the WP by SZ and SX to simulate the unit dimension wetted surface?" << std::endl;
		std::cout << "Have you change the peak pressure and the decay rate???" << std::endl;
		system("PAUSE ");
	}
	
	meshgenerationstruct a;
	a = meshgeneration();
	std::cout << "mesh generation done" << std::endl;

	for (i = 0; i < a.NNODE; i++) {
		if (abs(a.GCOORD[i][0]) < 1e-5 && abs(a.GCOORD[i][1] + 2.4384) < 1e-5 && abs(a.GCOORD[i][2] - (2.4384)) < 1e-5) {
			std::cout << " " << std::endl;
		}
	}

	LOBATTOstruct b;
	b = LOBATTO(N);
	std::cout << "LOBATTO() done" << std::endl;
	LOCAL_NODEstruct c;
	c = LOCAL_NODE(N);

	//observation: LNA is identical with the Gmsh user manual 
	TD_LOCAL_NODEstruct ct;
	ct = TD_LOCAL_NODE(NC);

	ELE_GEN(a.NEL, a.GCOORD, a.IEN, c.LNA, b.Z);

	std::cout << "ELE_GENstruct done" << std::endl;
	//---------------SHAPE FUNCTION ROUTINE------------------------------//
	//DETERMINE GLL QUADRATURE POINTS AND WEIGHTS
	GLLQUADstruct f;
	f = GLLQUAD(b.Z, b.WL, N, !FEM);
	//std::cout << "GLLQUAD done" << std::endl;
	LOCAL_SHAPEstruct g;
	g = LOCAL_SHAPE(c.LNA, N); 
	//std::cout << "LOCAL_SHAPE done" << std::endl;
	//DETERMINE LOCAL GEOMETERY SHAPE FUNCTION AT ELEMENT NODES
	LOCAL_GSHAPEstruct l;
	l = LOCAL_GSHAPE(f.S, c.LNA);
	//std::cout << "LOCAL_GSHAPE done" << std::endl;
	//DETERMINE JACOBIAN MATRIX 
	JACOBIANstruct m;
	m = JACOBIAN(a.NEL, l.GSHL, a.GCOORD, a.IEN, c.LNA);

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

	//std::cout << "JACOBIAN done" << std::endl;
	GLOBAL_SHAPEstruct n;
	n = GLOBAL_SHAPE(a.NEL, g.SHL, m.XS, m.JACOB);
	//std::cout << "GLOBAL_SHAPE done" << std::endl;
	MATRIXstruct o;
	o = MATRIX(a.NEL, a.NNODE, g.SHL, f.W, a.IEN, c.LNA, m.XS, n.SHG, m.JACOB);
	//duration_int = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	//std::cout << "time lapse for until this point: " << duration_int << std::endl;

	/*
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			for (k = 0; k < NINT; k++) {
				for (h = 0; h < NINT; h++) {
					delete[] m.XS[i][j][k][h];
				}
				delete[] m.XS[i][j][k];
			}
			delete[] m.XS[i][j];
		}
		delete[] m.XS[i];
	}
	delete[] m.XS;
	for (i = 0; i < NINT; i++) {
		for (j = 0; j < NINT; j++) {
			delete[] m.JACOB[i][j];
		}
		delete[] m.JACOB[i];
	}
	delete[] m.JACOB;
	for (i = 0; i < 4; i++) {
		for (j = 0; j < NINT*NINT*NINT; j++) {
			for (k = 0; k < NINT; k++) {
				for (h = 0; h < NINT; h++) {
					delete[] n.SHG[i][j][k][h];
				}
				delete[] n.SHG[i][j][k];
			}
			delete[] n.SHG[i][j];
		}
		delete[] n.SHG[i];
	}
	delete[] n.SHG;
	*/

	//std::cout << "MATRIX done" << std::endl;
	//DETERMINE MAXIMUM MESH EIGENVALUE TO FIND CFL TIMESTEP
	LMAX = EIGENMAX(o.QMASTER, o.HMASTER, a.NEL);
	FSILINK(f.W, c.LNA, a.IEN, g.SHL, a.GCOORD, a.NNODE, g.SHOD);

	//read the model file to MpCCI adapter
	char* modelfile = "model.txt";
	readfile(modelfile);

	TIMINTstruct t;
	t = TIMINT(LMAX);

	time_step_size = t.DT;

	/*
	double *T; //TIMESTEP ARRAY 
	T = new double[t.NDT];
	for (i = 0; i < t.NDT; i++) {
		T[i] = 0.0;
	}
	*/

	
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

	//define mapping functions for interface mapping
	double ***phi_fem;
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
	phi_fem = new double**[NCINT*NCINT];
	for (i = 0; i < NCINT*NCINT; i++) { //for all the points in that element
		phi_fem[i] = new double*[hpref + 1];
		for (j = 0; j < hpref + 1; j++) {
			phi_fem[i][j] = new double[hpref + 1];
		}
	}
	for (i = 0; i < NCINT*NCINT; i++) {
		for (j = 0; j < hpref + 1; j++) {
			for (k = 0; k < hpref + 1; k++) {
				phi_fem[i][j][k] = 0.0;
			}
		}
	}
	for (h = 0; h < NCINT; h++) { //stands for every fem point, LNA[u][v][w](shape function is based on those points)
		for (k = 0; k < NCINT; k++) {
			for (i = 0; i < hpref + 1; i++) {  //i j k are the independent variable in basis function(sem points)
				for (j = 0; j < hpref + 1; j++) {
					nomx = 1.0; nomy = 1.0; //multiplier initialization
					denomx = 1.0; denomy = 1.0; //multiplier initialization
					for (z = 0; z < NCINT; z++) { //loop through nominator and denominator in basis function expression
						if (z != h) {
							nomx *= (origp[i] - basep[z]);
							denomx *= (basep[h] - basep[z]);
						}
						if (z != k) {
							nomy *= (origp[j] - basep[z]);
							denomy *= (basep[k] - basep[z]);
						}
					}
					phi_fem[ct.LNA[h][k] - 1][i][j] = (nomx / denomx)*(nomy / denomy); //tensor product
					//the coordinate definition dof u,v is the same with i,j
				}
			}
		}
	}
	
	//Define the shape function value at Gauss-Lobatto points on base fluid/structure mesh
	//phi_femg
	double ***phi_femg;
	phi_femg = new double**[NCINT*NCINT];
	for (i = 0; i < NCINT*NCINT; i++) { //for all the points in that element
		phi_femg[i] = new double*[hprefg + 1];
		for (j = 0; j < hprefg + 1; j++) {
			phi_femg[i][j] = new double[hprefg + 1];
		}
	}
	for (i = 0; i < NCINT*NCINT; i++) {
		for (j = 0; j < hprefg + 1; j++) {
			for (k = 0; k < hprefg + 1; k++) {
				phi_femg[i][j][k] = 0.0;
			}
		}
	}
	//use f.S
	if (mappingalgo == 5) {
		LOBATTOstruct bb;
		bb = LOBATTO(hprefg);
		GLLQUADstruct ff;
		ff = GLLQUAD(bb.Z, bb.WL, hprefg, 0); //Gauss-Legendre point
		for (h = 0; h < NCINT; h++) { //stands for every fem point, LNA[u][v][w](shape function is based on those points)
			for (k = 0; k < NCINT; k++) {
				for (i = 0; i < hprefg + 1; i++) {  //i j k are the independent variable in basis function(sem points)
					for (j = 0; j < hprefg + 1; j++) {
						nomx = 1.0; nomy = 1.0; //multiplier initialization
						denomx = 1.0; denomy = 1.0; //multiplier initialization
						for (z = 0; z < NCINT; z++) { //loop through nominator and denominator in basis function expression
							if (z != h) {
								nomx *= (ff.S[i] - basep[z]);
								denomx *= (basep[h] - basep[z]);
							}
							if (z != k) {
								nomy *= (ff.S[j] - basep[z]);
								denomy *= (basep[k] - basep[z]);
							}
						}
						phi_femg[ct.LNA[h][k] - 1][i][j] = (nomx / denomx)*(nomy / denomy); //tensor product
						//the coordinate definition dof u,v is the same with i,j
					}
				}
			}
		}
	}

	//define the linear lagrange shape function defined on Gauss-Legendre points
	double ***phi_fem2;
	phi_fem2 = new double**[NCINT*NCINT];
	for (i = 0; i < NCINT*NCINT; i++) { //for all the points in that element
		phi_fem2[i] = new double*[2];
		for (j = 0; j < 2; j++) {
			phi_fem2[i][j] = new double[2];
		}
	}
	for (i = 0; i < NCINT*NCINT; i++) {
		for (j = 0; j < 2; j++) {
			for (k = 0; k < 2; k++) {
				phi_fem2[i][j][k] = 0.0;
			}
		}
	}
	//use f.S
	LOBATTOstruct bbb;
	bbb = LOBATTO(1);
	GLLQUADstruct fff;
	fff = GLLQUAD(bbb.Z, bbb.WL, 1, 0);
	for (h = 0; h < NCINT; h++) { //stands for every fem point, LNA[u][v][w](shape function is based on those points)
		for (k = 0; k < NCINT; k++) {
			for (i = 0; i < 2; i++) {  //i j k are the independent variable in basis function(sem points)
				for (j = 0; j < 2; j++) {
					nomx = 1.0; nomy = 1.0; //multiplier initialization
					denomx = 1.0; denomy = 1.0; //multiplier initialization
					for (z = 0; z < NCINT; z++) { //loop through nominator and denominator in basis function expression
						if (z != h) {
							nomx *= (fff.S[i] - basep[z]);
							denomx *= (basep[h] - basep[z]);
						}
						if (z != k) {
							nomy *= (fff.S[j] - basep[z]);
							denomy *= (basep[k] - basep[z]);
						}
					}
					phi_fem2[ct.LNA[h][k] - 1][i][j] = (nomx / denomx)*(nomy / denomy); //tensor product
					//the coordinate definition dof u,v is the same with i,j
				}
			}
		}
	}

	//Generation of tabular incident wave
	double x;
	int row;
	int col, col_a, col_t;
	double*timer;
	double*ampt;
	timer = new double[116];
	ampt = new double[116];
	std::string lineA;
	std::ifstream myyfile("shock_amplitude.txt");
	if (!myyfile) {
		std::cout << "can not open the shock amplitude file" << std::endl;
		system("PAUSE ");
	}
	while (myyfile.good()) {
		row = 0;
		while (std::getline(myyfile, lineA)) { //read the file line by line
			std::istringstream streamA(lineA);
			col = 0; col_a = 0; col_t = 0;
			while (streamA >> std::skipws >> x) {
				if (col % 2 == 0) {
					timer[4 * row + col_t] = x;
					col_t += 1;
				}
				else {
					ampt[4 * row + col_a] = x;
					col_a += 1;
				}
				col++;
			}
			row += 1;
		}
	}

	//derive the peak wave pressure and decay rate
	double XC = 0.0; double YC = 0.0; double ZC = 0.0;
	double XO = 0.0; double YO = 0.0; double ZO = 0.0;
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
	if (Bleich == 1) {
		PPEAK = 0.712e6;
		TAU = 0.999e-3;
		//PPEAK = 16.12e6;
		//TAU = 0.423e-3;
	}

	double KAPPA = 0.0;
	//-0.3048
	//-1.2192
	//-0.3048
	//-0.406399995, -0.203199998,           0.
	for (i = 0; i < a.NNODE; i++) {
		if (abs(a.GCOORD[i][0]) < 1e-3 && abs(a.GCOORD[i][1] + SY) < 1e-3 && abs(a.GCOORD[i][2] - SZ / 2) < 1e-3) {
			std::cout << "the node is: " << i << std::endl;
		}
	}
	
	TIME_INT(a.NNODE, a.GCOORD, f.W, ct.LNA, c.LNA, a.IEN, a.NEL, f.S, g.SHL, TIME, T, t.DT, t.NDT, b.Z,
		a.AYIN, o.HMASTER, o.Q, phi_fem, timer, ampt, KAPPA, PPEAK, TAU, XC, YC, ZC, XO, YO, ZO, g.SHOD, o.gamma, o.gamma_t, o.G, o.gamman, o.gamma_tn, o.Gn, phi_femg, phi_fem2);

	printf("Cleaning up...\n");
	/* clean up */
	cleanall();
	printf("Bye-bye!\n");
	exitcoupling();
	exit(0);
}
