/****************************************************************************
 *                                                                          *
 *  Very simple program to compute the motion of an elastic foundation      *
 *  ------------------------------------------------------------------      *
 *                                                                          *
 *  This program was written as an example to demonstrate the creation      *
 *  of a code adapter for MpCCI                                             *
 *                                                                          *
 *  Copyright:    Fraunhofer Institute SCAI, http://www.scai.fraunhofer.de  *
 ****************************************************************************/
#if __cplusplus
extern "C" {
#endif

/* maximum possible number of foundations */
#define MAXWETSURF 10 //coupling surfaces A B C D... 

/* structure for data storage */
//only need to define the values on wet surfaces???
//defined by user
typedef struct wetsurf{
  char* name;
  int dim;       //3 for 3D problem. (used for definition of primary var and secondary var)
  int direction; //direction=2 for y direction??? 
  //double stiffness; //seems not useful for the moment
  int nnodes; //number of nodes on coupling surface (e.g. surface A)
  int nelems; //2D element number on coupling surfaces
  double* nodecoord; //node coordinate on coupling surfaces 
  double* nodeforce; //node force/pressure giving to MpCCI service, eventually pass to Nastran./lu
  double* nodepressure;
  //double* nodedisp; //node force get from the MpCCI service computed from Nastran
  //double* elempressure; //may not be useful
  int* elemnodes; //1D connectivity matrix for coupling elements [nelems * nodesperelem]
  double dt;
} WETSURF;

/* list of foundations */
extern WETSURF* wsflist[MAXWETSURF]; //WETSURF=struct wetsurf
//for maximum 10 coupling wet surfaces A B C.../lu
//fndlist is a pointer array to data structure foundation

/* number of wet surface */
extern int wetsurfnumber;

extern double time_step_size;
extern double current_time;

/*  read foundation definitions from a file  */
void readfile(const char* filename);

/*  free all memory  */
void cleanall();

#if __cplusplus
}
#endif