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
//#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "data.h"



int wetsurfnumber;
WETSURF* wsflist[MAXWETSURF];
void readfile(const char* filename){
  FILE* fhandle;
  int ret, nodesperelem;
  int n, i, idx;
  double x,y,z;
  
  fhandle = fopen(filename,"r");
  if(fhandle == 0){ 
    printf("Cannot open >>%s<<.",filename);
    exit(99);
  }
  printf("Reading file >>%s<<...\n",filename);
  
  wetsurfnumber = 0;
  while(! feof(fhandle)){
	//Checks whether the end-of-File indicator associated with stream is set, returning a value different from zero if it is.
	//This indicator is generally set by a previous operation on the stream that attempted to read at or past the end-of-file
	WETSURF* fptr = malloc(sizeof(WETSURF));
    fptr->name=malloc(100);
    //ret = fscanf(fhandle,"EF %s %i %i\n",fptr->name,& fptr->dim,& fptr->direction); //scan from the first line of module definition
	ret = fscanf(fhandle, "EF %s %i %i %i\n", fptr->name, &fptr->dim, &fptr->direction, &fptr->element_type);
	//name is the name of each coupling surface
    if(ret != 4){ 
		printf("Reading error."); exit(-1);    
	}
	printf("foundation >%s<   dim=%i  direction=%i\n",
		fptr->name, fptr->dim, fptr->direction);
    switch(fptr->dim){
      case 2:
        nodesperelem = 2; //linear line element on coupling surface (2D)
        break;
      case 3:
		  if (fptr->element_type == 0) {
			  nodesperelem = 4; //linear quad element on coupling surface (3D)
		  }
		  if (fptr->element_type == 1) {
			  nodesperelem = 3; //linear triangular element for tetrahedral element (3D)
		  }
        break;
      default:
        printf("wrong dimension: %i\n",fptr->dim);
        exit(-1);
    }
    fscanf(fhandle,"\n");
    /*  nodes   */
    printf("  nodes...\n");
    if(fscanf(fhandle,"NODES %i",& fptr->nnodes) != 1){
      printf("input syntax error!\n");
      exit(-1);
    }
	printf("node number: %i\n", fptr->nnodes);
    fptr->nodecoord = malloc(sizeof(double) * fptr->nnodes * fptr->dim);
	//nodecoord is defined in 1D
    for(i=0; i<fptr->nnodes; ++i){
      printf("    %i : ",i);
      if(fptr->dim == 2){
        fscanf(fhandle,"%i %lf %lf\n",& idx, & x, & y);  //proceed to scan all the nodes on that coupling surface
        fptr->nodecoord[fptr->dim*i] = x;
        fptr->nodecoord[fptr->dim*i+1] = y;
        printf(" %g %g\n",x,y);
      }else{ //dim=3
        fscanf(fhandle,"%i %lf %lf %lf\n",& idx, & x, & y, & z);
	    //Since the handle is i, the node numbering in data file should be in sequence. 
        fptr->nodecoord[fptr->dim*i] = x;
        fptr->nodecoord[fptr->dim*i+1] = y;
        fptr->nodecoord[fptr->dim*i+2] = z;
        printf(" %g %g %g\n",x,y,z);
      }
    }
    /*  elements   */
    printf("  elements...\n");
    if(fscanf(fhandle,"ELEMENTS %i",& fptr->nelems) != 1){
      printf("input syntax error!\n");
      exit(-1);
    }
    fptr->elemnodes = malloc(sizeof(int) * fptr->nelems * nodesperelem);
    for(i=0; i<fptr->nelems; ++i){
      printf("    %i :",i);
      fscanf(fhandle,"%i",& idx);
      for(n=0; n<nodesperelem; ++n){
        fscanf(fhandle,"%i", & fptr->elemnodes[nodesperelem * i + n]);
        printf(" %i",fptr->elemnodes[nodesperelem * i + n]);
      }
      printf("\n");
      fscanf(fhandle,"\n");
    }

    //scan all the way to the end of the module. (start from here in the next loop) 
    printf("  allocating memory...\n"); 
	//allocating memory for values being exchanged in this module
	fptr->nodeforce = malloc(sizeof(double)* fptr->nnodes * fptr->dim);
	//std::cout << "surface point number: " << fptr->nnodes << std::endl; 
    //fptr->nodedisp  = malloc(sizeof(double)* fptr->nnodes);
    //fptr->elempressure = malloc(sizeof(double)* fptr->nelems);
	fptr->nodepressure = malloc(sizeof(double)*fptr->nnodes);
	wsflist[wetsurfnumber++] = fptr;
  }
}

/*
 *  free all memory 
 */
void cleanall(){
  int i;
  for (i = 0; i < wetsurfnumber; ++i) {
    free(wsflist[i]->name);
    free(wsflist[i]->nodecoord);
	free(wsflist[i]->nodeforce);
    free(wsflist[i]->elemnodes);
	//free(wsflist[i]->nodepressure);
    free(wsflist[i]);
	//free(wsflist[i]->elempressure);
	//free(wsflist[i]->nodedisp);
  }
}

