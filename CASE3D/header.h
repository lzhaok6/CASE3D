#pragma once
#include <vector>

struct LOBATTOstruct LOBATTO(int n);
struct LOCAL_NODEstruct LOCAL_NODE(int n);
struct TD_LOCAL_NODEstruct TD_LOCAL_NODE(int n);
//local node array is derived from "LNA_info" 
struct meshgenerationstruct meshgeneration();
void Neighborhood_search(double** GCOORD, int***LNA, int*IEN_flu, int NEL_flu);
//read the fem mesh file generated by Gmsh 
//the local node definition by MpCCI kernel 
//generate sem points from fem points by Lagrange interpolation  
//used to generate mesh on wet surface (identical with the model file read by data.c)
struct GLLQUADstruct GLLQUAD(double* Z, double* WL, int n, int SEM);
struct LOCAL_SHAPEstruct LOCAL_SHAPE(int*** LNA, int n, int NQUAD, int fem);
struct LEGENDREstruct LEGENDRE(int Nq, int n);
struct LOCAL_GSHAPEstruct LOCAL_GSHAPE(double* S, int*** LNA, int NINT);
struct JACOBIANstruct JACOBIAN(int NEL, double **GCOORD, int *IEN, int*** LNA);
struct GLOBAL_SHAPEstruct GLOBAL_SHAPE(int NEL, double***SHL, double**GCOORD, int*IEN, double* JACOB_tet, int*** LNA);
struct MATRIXstruct MATRIX(int NEL, int NNODE, double***SHL, double*W, int*IEN, int***LNA, double* JACOB_tet, double*** SHG_tet, double** GCOORD);
double EIGENMAX(double*** QMASTER, double*** HMASTER, int NEL);
struct TIMINTstruct TIMINT(double LMAX);
struct NRBstruct NRB(int NNODE, double **GCOORD, int*** LNA);
double** WAVE_IN(int NNODE, double** GCOORD, double* T, int TIME, double** PIN, double DT, double PPEAK, double TAU, double XC, double YC, double ZC, double XO, double YO, double ZO, std::vector<int> shadow_pts);
void FSILINK(int*** LNA);
struct interface_mappingstruct interface_mapping(int fluid2structure, double ** GCOORD, double* WP, int* IEN, int***LNA);
void TIME_INT(int NNODE, double** GCOORD, int***LNA_3D, int*IEN, int NEL, int TIME, double *T, double DT, int NDT, double* Q, double KAPPA, double PPEAK, double TAU, double XC, double YC, double ZC,
	double XO, double YO, double ZO, double ***SHOD, double gamman[], double gamma_tn[], double***Gn,
	double****gamma_t, double ****gamma, double*****G, double*W, double*** SHL, double*** SHG_tet, double* JACOB_tet, double** HMASTER, std::vector<int> shadow_pts);
//used to map the force value from user defined fluid mesh to MpCCI defined mesh and map the displacement in the opposite way. 
std::vector<int> shadow_pt_identification(double XC, double YC, double ZC, int NNODE, double** GCOORD); 

//extern int** nodesperelem;
//extern int nodesperelem2; 
const int N = 1;    //N is the element order of fluid mesh 
const int NINT = N + 1; //NINT=N+1;
const int hprefg = 1; //The level of Gauss-Legendre integration on the structure mesh (for mapping algorithm 4 and 5)
const int hprefg_flu = N; //The level of Gauss-Legendre integration on the fluid mesh (for mapping algorithm 4 and 5) 
typedef struct owetsurf {
	double *WBS; //wet surface structure force derived from displacement sent back from Nastran 
	double *PSI; //integrated incident pressure 
	double *DI;  //displacement predictor 
	double **DISPI; //incident displacement
					//double **DISP; //total (actual) displacement (the combination of incident displacement and dynamic displacement)
	double *DISP; //total (actual) displacement (the combination of incident displacement and dynamic displacement)
	double **DISP_norm; //normal displacement out of fluid domain for two concecutive time steps
						//double *WP; //wet surface pressure, representing AP, CP, DP previously
	double *WPIN; //wet surface incident pressure, APIN, ...
	double*** FPMASTER;
	double*** FPMASTER_2D;
	int **IEN_2D; //wet surface connectivity matrix (used to write MpCCI model file and information mapping. The wetted surface information from MpCCI is defined this way)
	int **IEN_algo2; //
	int **IEN_gb; //connectivity matrix on 2D surface pointing to global points.
	int **IEN_lc; //
	int FSNEL;
	int *GIDF;   //Coupled fluid element array (numbering)
	int *GIDN_MpCCI; //wet nodes on coupling surfaces
	int GIDNct_MpCCI; //wet nodes number on coupling surfaces
	int *GIDN;
	int GIDNct;
	double** norm; //store the normal direction of linear elements 
	double** Jacob_2D; //the Jacobian value of 2D element on wetted surface
	double** Jacob_test;
	//int ***LNA_2D;
	int *LNA_2D;
	//int LNA_algo2[2][2]; //The local node orientation of the linear element in algorithm2
	//int*** LNA_algo2;
	int* LNA_algo2;
	//int FP[NINT*NINT]; //The 1D version of DP in order to facilitate the calculation of FPMASTER
	int FP_temp[NINT*NINT];
	int **FP;
	//int FP_2D[NINT*NINT];
	//int** FP_2D;
	int* FP_2D;
	//double** JACOB;
	//int Jacob_face[2]; //Used to identify which coordinate dimension (x,y or z/xi,eta,zeta?) to be used for surface 2D Jacobian matrix calculation.  
	int** Jacob_face;
	int* LNA_norm;
	int* LNA_norm_3D; //Used along with LNA_norm to identify GIDF
	int** IEN_py; //the linear element connectivity matrix of 2D physical groups (boundaries) separated from the original 2D high-order elements
				  //double**** xs;
	double****xs_2D; //for the 2D elements on wetted surface
	double ***phi_fem;
	int* LNA_JB2D;
	//LNA_JB2D is the node numbering of the wetted surface 2D element in the linear element numbering scheme. Used to derive the 2D Jacobian matrix. The numbering sequence set to be the same as LNA_norm which is the local numbering of the 2D element used to obtain the coordinate. 
	//LNA_norm is originally devised for mapping algorithm 2. However, it is leveraged in 2D Jacobian determinant derivation (assuming linear geometric mapping) because it numbering is the corner nodes of the 2D wetted surface element.  
	double*****GSHL_2D;
	int dir;
	double *dimension;
	double *dimension_hex;
	double** flu_local;//used to store the local coordinate of the projected fluid point on structural elements
	int* flu_stru;
	//std::vector<int>orphan_flag_flu; 
	int* orphan_flag_flu;
	double** flu_stru_global;
	int** IEN_flu_2D;
	double** GCOORD_flu_gs; //the gauss point on fluid element 
							//double** DISP_gs; //The displacement on gauss nodes  
	double* DISP_gs; //The displacement on gauss nodes  
	double* GCOORD_fs; //The array built for fast access
	int* IEN_flu_3D; //Get the corner point of the fluid elements 
	double** PIN_gs; 
	double** PSI_inc; 
	double* WPIN_gs;
} OWETSURF;

typedef struct stru_wet_surf {
	int*gs_flu; //used to store the corresponding fluid element or fluid point (for orphan node)
	double**gs_flu_global;
	double**gs_flu_local; //used to store the local coordinate of the projected structural gauss point on fluid element
	//std::vector<int>orphan_flag_gs; //the container to store the node numbering of orphan nodes
	int* orphan_flag_gs;
	int ELE_stru; //total number of structural wetted surface elements
	//int** IEN_stru; //Connectivity matrix of the structural wetted surface elements
	int* IEN_stru; //Connectivity matrix of the structural wetted surface elements
	double **GCOORD_stru; //the coordinate of structure nodes
	//int** IEN_stru_MpCCI; //Renumber the node on wetted surface to start from 1
	int* IEN_stru_MpCCI; //Renumber the node on wetted surface to start from 1
	int Node_stru; //total node on structural wetted surface
	int* Node_glob; //Associate the corresponding global node number to the local node in MpCCI model file
	double** Jacob_stru; //the Jacobian determinant of the structural wetted surface element (the value is on quadrature node)
	double*** FPMASTER_stru;
	//int** LNA_stru; 
	int* LNA_stru;
	int** LNA_gs;
	//double** P_gs;
	double* P_gs;
	double** norm_stru;
	double W_stru[hprefg + 1];
	double*** phi_stru;
	double**** GSHL_2D;
	double****xs_2D; //for the 2D elements on wetted surface
	int** IEN_stru_norm;
	double* GCOORD_stru_gs; //the coordinate of structure wetted surface gauss points
	double* dispi_stru; //the incident displacement on the structural mesh points (not gauss points)
	int** IEN_stru_gs;
	int gs_num;
	int** FP_flu; //store the interface node for mapping algorithm 4
	int* elenode; //the number of node of the element (3 means triangular element; 4 means quad element)
	double** GCOORD_stru_fs;
	double* PSI; 
	double* DI; 
	double** DISPI; 
	double** PSI_inc; 
	double** norm_stru_pt; //normal vector on structural wetted surface points
	double **PIN; 
}STRU_WET_SURF;

//Store the properties on NRB surface (currently just one)
//Different from the wetted surface, the element on the NRB surface is high-order
typedef struct nrbsurf {
	double*** ADMASTER; //The boundary force integration weight for each element on NRB surface
						//double* ADMASTERG; //the lumped boundary force integration weight for each node on NRB surface
						//double** JACOB; //3D Jacobian determinant with one additional quadrature point (Dedicated for boundary force integration on boundaries). 
	int** IEN_gb; //connectivity matrix on 2D surface pointing to global points. (From the physical group definition of the mesh file)
	int **IEN_lc;
	int NEL_nrb; //the element number of the 2D element 
				 //int NNODE_nrb; //The total number of node on NRBC
				 //int*** LNA_2D;
	int* LNA_2D;
	//int DP[NINT*NINT]; //The 1D version of DP in order to facilitate the calculation of ADMASTER
	int DP_temp[NINT*NINT];
	int** DP;
	//int DP_2D[NINT*NINT];
	//int **DP_2D;
	int *DP_2D;
	int* NRBA;
	int NRBNODE;
	double**norm;
	int* LNA_norm;
	int* LNA_norm_3D;
	//double**** xs;
	double **XEST_kn;
	double **XEST;
	double **XEST_ukn;
	double **XNRBORG2;
	double ***XNRB_kn;
	double ***XNRB_ukn;
	double ***XNRB;
	double **XCOR_ukn;
	double **XCOR;
	int* NRBELE_ARR;
	int* LNA_JB2D;
	double ***P_dev; //The gradient of the pressure of NRB element points [NEL_nr][NINT*NINT][3]
	double** Jacob_2D;
	double*****GSHL_2D;
	double****xs_2D;
	//int Jacob_face[2];
	int** Jacob_face;
	double **angle_disp1;
	double **angle_disp2;
	double ***disp_mag;
	//double* ADMASTERG; 
	//double *BNRB; 
	double *dimension;
	double **XNRBORG;
	double*** FEEDOT_inc;
} NRBSURF;

struct LOBATTOstruct {
	double* Z;
	double* WL;
};

struct LOCAL_NODEstruct {
	int*** LNA;
	int L[NINT*NINT*NINT][3];
};

struct TD_LOCAL_NODEstruct {
	int** LNA;
};
struct TD_ELE_GENstruct {
	double **GCOORD;
	int **IEN;
};

struct meshgenerationstruct {
	double **GCOORD;
	//int **IEN;
	int *IEN;
	int NNODE;
	int NEL;
	//double** AYIN;
	//int wetnodenum;
	//int*wetnode;
	std::vector<std::vector<std::vector<int>>> BCIEN;
	int** IEN_wt; ////connecvitity matrix of wetted surface (after removing the free surface elements)
};

struct interface_mappingstruct {
	double energy_sent;
	double energy_rec;
};

struct GLLQUADstruct {
	double* S;   //ARRAY OF N+1 GLL QUADRATURE POINTS
	double* W;   //ARRAY OF N+1 GLL QUADRATURE WEGHTS
};

struct LEGENDREstruct {
	double** LN;
};

struct LOCAL_SHAPEstruct {
	double*** SHL;
	double*** SHL_2D;
	double ***SHOD;
};

struct LOCAL_GSHAPEstruct {
	double*** GSHL;
	//double****GSHL_2D; 
	double MCOORD[8][3];
	double MCOORD_2D[4][2];
};

struct JACOBIANstruct {
	//double***XS;
	//double**JACOB; 
	double**XS_tet;
	double*JACOB_tet;
	int dummy;
};

struct GLOBAL_SHAPEstruct {
	//double****SHG;
	double***SHG_tet;
	int dummy;
};

struct MATRIXstruct {
	//double***HMASTER;
	double**HMASTER;
	//double*HMASTER;
	double *Q;  //GLOBAL CAPACITANCE MATRIX
	double ****gamma;
	double gamman[NINT*NINT*NINT * 3];
	double **gamman2;
	double ****gamma_t;
	double gamma_tn[NINT*NINT*NINT * 3];
	double *****G;
	//double ****Gn;
	double ***Gn;
	double LMAX;
};

struct NRBstruct {
	double *ADMASTERG;
	int* NRBA_t;
	int NNODE_nrb;
};

struct TIMINTstruct {
	double DT; //TIME STEP
	int NDT;    //NUMBER OF TIME STEP 
				//int dtscale;
};

//Input values
//const double fs_offset = -6.02674; //Used to draft the free surface to y=0 position.
const double fs_offset = -6.0;
//const double fs_offset = 0.0; //FSP
//const double SX = 0.1;
const double SX = 8.5344 / 2; //14ft (FSP)
//const double SY = 1.2192; //4ft (FSP)
const double SY = -fs_offset; //for DDG case
//const double SY = 3.048; //10ft
//const double SZ = 4.8768; //16ft (FSP)
const double SZ = 0.0; //for DDG case (the stand-off is at the keel)
const int NC = 1;   //NC is the element order on coupling mesh 
const int NCINT = NC + 1; //NCINT=NC+1;
const int Nq = N + 1; //The integration order for boundary nodal force term (exact integration). Should be at least one unit higher than the interpolation order (for algorithm 1, 2 and 5) since the Gauss-Legendre-Lobatto nodes are not accuracy enough. Need not to be used for FEM case since the Gauss-Legendre nodes is accurate enough. 
const int NqINT = Nq + 1;
const int refine = 1; //The refinement rate of fluid mesh against base fluid mesh for h refinement. 
const int hpref = refine*N; //total refinement level of h and p refinement
//const int hprefg = refine*N; //The level of Gauss-Legendre integration on the base mesh (dedicated for mapping algorithm 5) this could integrate the nodal force on the linear base mesh upto the order 2(refine*N)-2
//const int hprefg = 1;
const int mappingalgo = 2; //Mapping algoritm, please refer to the description in the main file (1, 2, 3, 4)
const double RHO = 1025.0; //original
//const double RHO = 989.0; //Bleich-Sandler
const int WAVE = 2; //1 for plane wave; 2 for spherical wave 
const int step = 0;
const int wavdirc[3] = { -1,0,0 }; //the direction of incident plane wave (positive y axis) 
const double C = 1500.0; //original  
//const double C = 1450.0; //Bleich_Sandler	
//const double CFLFRAC = 0.5;  //original 
const double CFLFRAC = 1.0;
const int dtscale = 1;
const double BETA = 0.25;   //original 
const double TTERM = 0.030;    //SIMULATION END TIME 
const int CAV = 1; //1 for cavitation, 0 for non-cavitation 
const double PSAT = 0.0; //saturated pressure 
const double pi = 3.141593;
const double grav = 9.81;
const double PATM = 101.3e3; //pa 
const double stdoff = -3.40 / 0.3048; //ft
//const double stdoff = 0; //xft
//const double depth = 26 / 0.3048; //ft
const double depth = 35.02 / 0.3048; //ft
//const double depth = 70; //ft
//const double x_loc = 74.22;//m for DDG case
const double x_loc = 57.34; //midship (m)
const double y_loc = 0;
const double z_loc = 0;
const double W = 92.53; //charge weight (lb)	

//standoff point in spherical wave case 
//const double XO = 0; //the origin is at the deck center of the barge
//const double YO = -SY;
//const double ZO = SZ / 2;
//the parameter to control whether a tabulated smoothed waveform is used 
const double output_int = 1e-3; //output file time interval (0.5ms)
const int debug = 0; //is the code in debug mode?
const int debug2 = 0;
const int debug3 = 0;
const int debug4 = 0;
const int debug5 = 0;
const int debug_algo5 = 0;
const int debug_PH = 1;
const int debug_IEN = 0;
const int charge_loc = 0; 
const int fsdebug = 0;
const int tfm = 1; //is the total field model used? 
const int tensorfactorization = 1;
const int TNT = 1;
const int output = 0;
const int FEM = 0; //Is this a first order FEM code? 
const int nodeforcemap2 = 1; //If the property to be mapped by MpCCI is nodal force (use 0 if the property is absolute pressure)
const int owsfnumber = 1; //The number of fluid wetted surfaces. 
const int nrbsurfnumber = 1; //The number of fluid NRB surface. 
const int ssnumber = 1; //The number of structural wetted surface (used for algorithm 4 and 5)
const int wt_pys_num[owsfnumber] = { 0 };  //the physical group number that corresponds to the wet surface (physical group 3)
const int nrb_pys_num[nrbsurfnumber] = { 1 };
//const int wt_pys_num[owsfnumber] = { 0,1,2,3 };
//const int nrb_pys_num[nrbsurfnumber] = { 4,5,6,7 };
const double XHE = 0.3048;
//const double XHE = 0.1; //Bleich_Sandler
const double YHE = 0.3048;
//const double YHE = 0.141 / 2; //Bleich_Sandler 
const double DY = 2 * SY;
//const double DY = 3.807; //Bleich_Sandler
//const int SYNEL = 4 * refine;
const int SYNEL = 0;
const double xo = 0.3048;
const int Bleich = 0;
const int improvednrb = 0;
const int element_type = 0; //0 for hexahedral element; 1 for tetrahedral element; 
const int nodeadj = 1; //If the node coordinate needs to be adjusted. 
const int freesurface = 1;
const int hydrostatic = 1;
const int tecplot = 0;
const int incidentdisp_on_fluid = 1; //If the incident pressure displacement is defined on the fluid side rather than the structural side
const int Colewave = 0; 