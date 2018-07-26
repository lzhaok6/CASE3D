/****************************************************************************
 *   adapter.c                                                              *
 *                                                                          *
 *   MpCCI Code API  -  Code Adapter Source File                            *
 *   Changed for use with the C version of "foundation"                     *
 *                                                                          *
 *   This file contains definitions of all required driver functions.       *
 *   See the MpCCI Programmers Guide for more information.                  *
 *                                                                          *
 *      http://www.scai.fraunhofer.de/mpcci.html                            *
 *                                                                          *
 ****************************************************************************/

/* includes: */

#include <stdio.h>    /* standard C libs */
#include <stdlib.h>
#include <string.h>
#include "adapter.h"   /* code adapter header file */
#include "data.h"     /* CHANGED: added data.h for    data access */

/*#########################################################################
 *############  DATA STRUCTURES                             ###############
 *#########################################################################*/

/* list of driver functions
 *****************************************************/
//the structure is defined in mpcci_types.h/lu
//define the data value in that data structure (communicate with MpCCI server)/lu
static MPCCI_DRIVER MpCCIDriverFunctions = {
  //pointers to those driver functions are passed to MpCCI by ampcci_config/lu
  /******* REQUIRED: information to do a compatibility check  *********/
  (unsigned)sizeof(MPCCI_DRIVER),     /* size of structure MPCCI_DRIVER */ //MPCCI_DRIVER is a data structure/lu 
  MPCCI_CCM_VERSION,                  /* this_version */ //a macro/lu
  0,                                  /* tact_required bitmask */
  /* CHANGED: Name of simulation code and code-specific names */
  {                       /* list of code-specific names (for output) */
     "single point",      /* name for point */
     "Line",              /* name for lines */
     "Wet surface",       /* name for faces */
     "Volume - UNUSED",   /* name for volumes */
     "Global var"         /* name for global quantities */
  },

  /********* REQUIRED: (unsigned)sizeof(float|double): tell the manager
             the size of floating point values that are arguments
             to the set of putXXXValues() functions  *******************/
  (unsigned)sizeof(double),                     /* values_size */


  /********** OPTIONAL: methods called before/after some actions
              Replace by NULL pointers if not needed. **********/
  NULL,                            /* MpCCI_Driver_afterCloseSetup   */
  NULL,                            /* MpCCI_Driver_beforeGetAndSend  */
  NULL,                            /* MpCCI_Driver_afterGetAndSend   */
  NULL,                            /* MpCCI_Driver_beforeRecvAndPut  */
  NULL,                            /* MpCCI_Driver_afterRecvAndPut   */

  NULL,                            /* MpCCI_Driver_partSelect        */
  MpCCI_Driver_partUpdate,         /* MpCCI_Driver_partUpdate        */
  NULL,                            /* MpCCI_Driver_appendPart        */
  NULL,                            /* MpCCI_Driver_getPartRemeshState*/

  /********* REQUIRED: methods to update and check component names and
                       to define the grid (nodes & elements) *********/
  MpCCI_Driver_definePart,         /* MpCCI_Driver_definePart        */

  NULL,                            /* MpCCI_Driver_partInfo          */
  NULL,                            /* MpCCI_Driver_getNodes          */
  NULL,                            /* MpCCI_Driver_getElems          */

  /********* REQUIRED for moving nodes instead of full remeshing:
                      returns the MpCCI error code *******************/
  NULL,                            /* MpCCI_Driver_moveNodes         */

  /********* OPTIONAL: Data exchange functions ***********/
  NULL,                            /* MpCCI_Driver_getPointValues    */
  NULL,                            /* MpCCI_Driver_getLineNodeValues */
  NULL,                            /* MpCCI_Driver_getLineElemValues */
  MpCCI_Driver_getFaceNodeValues,  /* MpCCI_Driver_getFaceNodeValues */
  NULL,                            /* MpCCI_Driver_getFaceElemValues */
  NULL,                            /* MpCCI_Driver_getVoluNodeValues */
  NULL,                            /* MpCCI_Driver_getVoluElemValues */
  //NULL,                            /* MpCCI_Driver_getGlobValues     */
  MpCCI_Driver_getGlobValues,    /* MpCCI_Driver_getGlobValues     */
  NULL,                            /* MpCCI_Driver_putPointValues    */
  NULL,                            /* MpCCI_Driver_putLineNodeValues */
  NULL,                            /* MpCCI_Driver_putLineElemValues */
  MpCCI_Driver_putFaceNodeValues,  /* MpCCI_Driver_putFaceNodeValues */
  NULL,                            /* MpCCI_Driver_putFaceElemValues */
  NULL,                            /* MpCCI_Driver_putVoluNodeValues */
  NULL,                            /* MpCCI_Driver_putVoluElemValues */
  NULL,                            /* MpCCI_Driver_putGlobValues     */
};


static MPCCI_TINFO mpcciTinfo = {0};
//initialize information function.(used to initialize the coupling)/lu
static MPCCI_JOB  *mpcciJob   = NULL;

/*#########################################################################
 *############  ADDITONAL HELPER FUNCTIONS                  ###############
 *#########################################################################*/

/* CHANGED: Implemented some helper functions: */

/* output of adapter messages */
void adapter_output(const char *msg, int len)
{
   if (len) fputs(msg,stdout);
   fflush(stdout);
}

/* print error message */
void error(const char* message, int len){
   if(len) printf("ERROR: %s",message);
   fflush(stdout);
}

/* get ID of a surface if name is given */
int getSurfaceID(const char* name){
   int i, surfaceID;

   surfaceID = MAXWETSURF;

   for(i=0; i<wetsurfnumber; ++i){
      if(strcmp(name, wsflist[i]->name) == 0)
         surfaceID = i;
   }
   if(surfaceID == MAXWETSURF)
      printf("Coupling surface not found!!!");

   return surfaceID;
}

/*#########################################################################
 *############  INTERFACE FUNCTIONS                         ###############
 *#########################################################################*/

/* initialize the coupling process */
//The code has an initcoupling() method available to implement some stuff to initialize the code./lu
//At this point the information about the mesh should be available./lu
void initcoupling(){
   MPCCI_CINFO cinfo; //used to pass information from MpCCI server to user code (MPCCI_CINFO is a data structure type)/lu

   /* already initialized? */
   //if already initialized, then jump through this function
   if (mpcciTinfo.mpcci_state > 0) return;
   //Tinfo is short for transfer information
   printf("Initializing the coupling process!\n");

   /* initialize output routines */
   //It is not necessary to define all functions (give NULL pointers instead), 
   //but you should define either the first set of three functions or the second set of three functions./lu
   umpcci_msg_functs(adapter_output,  /* printFullMessage */
                     adapter_output,  /* printFullWarning */
                     error,           /* printFullFatal   */
                     NULL,            /* printLineMessage */
                     NULL,            /* printLineWarning */
                     NULL,            /* printLineFatal   */
                     NULL);           /* atExitHandler    */


   umpcci_msg_prefix("MpCCI: ");
   //This function is used to set a prefix for the coupled simulation./lu

   /* call ampcci_tinfo_init */
   //This function is used to determine the transfer information selected in the MpCCI GUI./lu
   if (!ampcci_tinfo_init(&mpcciTinfo,NULL))
   {
      mpcciTinfo.mpcci_state = -1; 
      return;
   }
   //The return value of "ampcci_tinfo_init" is 1 if the transfer information was successfully set up./lu
   //if information is not set up correctly, assign mpcci_state to -1 to stop the initialization/lu

  /*
   * initialize the transfer/time info and the code info
   */
   mpcciTinfo.iter = -1;
   //If this iteration information is set to the value -1, this designs the code to be an explicit transient code./lu
   //Otherwise this notifies an imdplicit transient code and enables the code to request data inside the same time step for different iteration steps./lu
   //mpcciTinfo.time = -1.0;
   mpcciTinfo.time = 0.0;
   //time is set to -1.0 for steady state iterative code, the same rule is applied to MPCCI_TINFO/lu
   //During the coupling process, data is exchanged, which must be identified by the simulation time or iteration. This allows the MpCCI server to manage the data exchange and to provide the corresponding data requested by a code.
   //mpcciTinfo.dt   =  0; 
   mpcciTinfo.dt = time_step_size;
   //mpcciTinfo.dt = 1.3826008770376172e-05;
   memset(&cinfo,0,sizeof(cinfo));

   /* CHANGED: Name of simulation code and error handling */
   cinfo.codename = "3dshock";
   cinfo.flags    = MPCCI_CFLAG_TYPE_FEA|MPCCI_CFLAG_GRID_CURR; //The code is a finite element analysis code //The code has the current grid available.
   cinfo.time     = mpcciTinfo.time; //the initial physical time of the code/lu
   cinfo.iter     = mpcciTinfo.iter; //The initial iteration counter of the code./lu
   cinfo.nclients = 1; //The exact number of expected server clients./lu
   cinfo.nprocs   = 1; //The number of parallel processes of the current code./lu

   /* call mpcci_init */
   //The initialization functiodn must be called from your code before the beginning of the iteration or time loop./lu
   //This call should establish a connection to the MpCCI server./lu
   mpcciJob = mpcci_init(NULL, &cinfo);
   //user manual 780
   //function mpcci_init is declared in mpcci.coreapi.h (definition resides in the library)/lu
   //The return value is an MPCCI JOB value and not null if a connection to the server could be established./lu
   //The MPCCI JOB value returned after calling the mpcci init containing the connection information to the MpCCI server./lu
   //mpcciJob=1 if the function "mpcci_init" executed correctly

   if (mpcciJob && ampcci_config(mpcciJob,&MpCCIDriverFunctions) <= 0)
      mpcci_quit(&mpcciJob);
   //&MpCCIDriverFunctions: Pointer to list of driver functions/lu
   //the return value of ampcci_config is 1 if the code successfully defined the mesh and quantities to the server.
   //the details of what ampcci_config does can be found page 781/lu
   //call MpCCI Driver definePart if defined(see Figure 21), otherwise Call MpCCI Driver getNodes and MpCCI...

   if (!mpcciJob)  //mpcciJob should be 1 if correctly configured
   {
      mpcciTinfo.mpcci_state = -1;
      MPCCI_MSG_FATAL0("This is no coupled computation!");
   }
   else
   {
      mpcciTinfo.mpcci_state = 1; //mpcci_state=1 means successfully initialized
      mpcciTinfo.mpcci_used  = 1;
   }
}

/* perform a transfer of quantities */
//The code has an dotransfer() method available to implement a call to the MpCCI data transfer function ampcci transfer/lu
void dotransfer(){
   int ret;
   extern double current_time;
   extern double time_step_size;

   if (!mpcciTinfo.mpcci_state) //mpcciTinfo.mpcci_state=1 if successfully initialized already/lu
	   initcoupling();
   //if the code has not been initialized yet, initialize it./lu  

   if (!mpcciTinfo.mpcci_used)
   {
      MPCCI_MSG_WARNING0("MpCCI is not used.\n");
      return;
   }

   /* Update the time /iteration information for mpcciTinfo */
   /* foundation code is just a QTID code */

   mpcciTinfo.iter = -1;
   mpcciTinfo.time = current_time; //This is the current simulation time 
   //mpcciTinfo.time = -1.0;
   //mpcciTinfo.time = 0.0;
   mpcciTinfo.dt = time_step_size; /* current physical time step size: updated by the code before ampcci_transfer() */
   //mpcciTinfo.dt = 0;
   //may have to link these information with user defined code/lu 
   //First the code must update the MPCCI TINFO data structure (See ? 2.9.6 Transfer Information: MPCCI TINFO ?). The current time, iteration must be set.
   //then the call of ampcci transfer could be performed

   /* do the transfer */ 
   ret = ampcci_transfer(mpcciJob,&mpcciTinfo); //P790,P812
   //return value: -1 error 1 transfer done 0 no transfer occured/lu
   //details of ampcci_transfer is in page 791/lu
   if(ret < 0)
       MPCCI_MSG_WARNING0("Could not complete all requested transfers!");
}

/* exit the coupling process */
//The code has an exitcoupling() method available to implement a call to the MpCCI exit function mpcci quit or mpcci stop./lu
void exitcoupling(){
   mpcci_quit(&mpcciJob);
}
//These functions always return 0/lu


/*###########################################################################
 *############  COUPLING COMPONENT DEFINITIONS                ###############
 *#########################################################################*/

/****************************************************************************
 *  MpCCI_Driver_partUpdate
 *
 *  This function is called during ampcci_config.
 *  You should set the number of nodes and number of elements
 *  of all coupling components here. (number!!!)/lu
 *
 ****************************************************************************/
/****************************************************************************/
int MpCCI_Driver_partUpdate(MPCCI_PART *part, MPCCI_QUANT *quant)
//This function is used to get information of the coupling surface and send to MpCCI server using macro functions/lu
//update PART property (can optionally update those properties in MpCCI_Driver_definePart function)
//This function is called during ampcci_config
/****************************************************************************/
{
   int surfaceID;         /* ID of a coupling surface */

   /* get internal ID of the coupling component */
   surfaceID = getSurfaceID(MPCCI_PART_NAME(part));
   //MPCCI_PART_NAME(part) returns the name of the coupling component
   MPCCI_PART_PARTID(part) = surfaceID;
   /* set number of nodes */
   MPCCI_PART_NNODES(part) = wsflist[surfaceID]->nnodes;
   //wsflist is defined in data.h/lu
   /* set number of elements */
   MPCCI_PART_NELEMS(part) = wsflist[surfaceID]->nelems;
   MPCCI_PART_CSYS(part) = (wsflist[surfaceID]->dim == 2) ? MPCCI_CSYS_C2D : MPCCI_CSYS_C3D;

   /* Define the coordinates system and Model dimension:
      Update this to your current model */

   MPCCI_PART_MDIM(part) = MPCCI_MDIM_FACE; //mesh dimension/lu

   return 0;
}

/***************************************************************************
 * MpCCI_Driver_definePart
 *
 * This function get called during MpCCI_Init.
 *
 * MpCCI needs information about the used mesh for a coupled component.
 * You must define nodes and elements within this driver function.
 *
 ****************************************************************************/
/*****************************************************************************/
int MpCCI_Driver_definePart(MPCCI_SERVER *server, MPCCI_PART *part)
/*****************************************************************************/
{
   /* IMPORTANT: The node definitions may look completely different in your
                 code. This is just an example! */

   double *nodeCoords;      /* node coordinates */
   int surfaceID;           /* surface ID */
   int maxNumNodes;         /* maximum number of nodes */
   int i;                   /* loop counter */
   int dim;                 /* space dimension */
   int * nodeNumbers;       /* array of node IDs */
   //int * elemNodes;         /* array of element nodes */
   //unsigned * elemTypes;         /* array of element types */
   int * elemNodePtr;       /* pointer to element nodes */
   int elemID;              /* element ID */
   int ret;

   /* get internal ID of the coupling component */
   surfaceID = MPCCI_PART_PARTID(part);
   dim = wsflist[surfaceID]->dim;
   switch(dim)
   {
      case 2:
         maxNumNodes = 2;
         break;
      case 3:
         maxNumNodes = 4;
         break;
      default:
        MPCCI_MSG_FATAL0("Wrong space dimension, must be 2 or 3.");
   }
   /***********************
   STEP 1: Define the part
   ***********************/
   MPCCI_MSG_INFO1
   (
      "Coupling grid definition for component \"%s\" ...\n",
      MPCCI_PART_NAME(part)
   );

   //the returned information is updated by driver function "MpCCI_Driver_partUpdate"/lu
   //those values can alternatively defined in this function to eliminate the need of "MpCCI_Driver_partUpdate"/lu
   ret = smpcci_defp
         (
            server, //The MpCCI server pointer passed by the MpCCI Driver definePart
            MPCCI_PART_MESHID(part),
            MPCCI_PART_PARTID(part),
            MPCCI_PART_CSYS  (part),
            MPCCI_PART_NNODES(part),
            MPCCI_PART_NELEMS(part),
            MPCCI_PART_NAME  (part)
         );
   //It returns 0 for a successful operation./lu
   //check if those coupling part properties are already sent to server/lu 

   /***********************
   STEP 2: Collect nodes of the surface in one array
   ***********************/
   if (!ret)  //if step 2 is successful/lu
   {
      /* allocate memory for the node coordinates (3-dimensional)*/
      nodeCoords = (double*) malloc(MPCCI_PART_NNODES(part)*sizeof(double)*dim);

      /* array for node-numbers */
      nodeNumbers = (int*) malloc(MPCCI_PART_NNODES(part)*sizeof(int));

      /* loop over all nodes */
      for(i=0; i<MPCCI_PART_NNODES(part); i++)
      {
          int c; /* coordinate counter */
          /* get node number */
          nodeNumbers[i] = i ;
          /* and the coordinates */
          /* copy node coordinates to buffer */
          for(c=0; c<dim; ++c)
          {
             nodeCoords[dim*i + c] = wsflist[surfaceID]->nodecoord[dim*i + c];
          }
      }

      /* Send the nodes definiton for this coupled component */
	  //smpcci_pnode defines the nodes of a coupling component/lu
      ret = smpcci_pnod
            (
               server,
               MPCCI_PART_MESHID(part),  /* mesh id */
               MPCCI_PART_PARTID(part),  /* part id */
               MPCCI_PART_CSYS  (part),  /* coordinates system */
               MPCCI_PART_NNODES(part),  /* no. of nodes */
               nodeCoords,               /* node coordinates */
               (unsigned)sizeof(double), /* data type of coordinates */
               nodeNumbers,              /* node ids */
               NULL
            );
      free(nodeCoords);
      free(nodeNumbers);
   }
   //It returns 0 for a successful operation/lu 

   /***********************
   STEP 3: Collect elements of the surface in one array
   ***********************/
   if (!ret)  //if step 2 is successful/lu
   {
	  int * elemNodes;         /* array of element nodes */
	  unsigned * elemTypes;         /* array of element types */
      /* allocate memory for the list of nodes for each element */
      elemNodes    = (int*) malloc(MPCCI_PART_NELEMS(part)*maxNumNodes*sizeof(int));
      elemNodePtr = elemNodes;
      elemTypes   = (unsigned*) malloc(MPCCI_PART_NELEMS(part)*sizeof(unsigned));

      /* loop over elements in coupling component */
	  for (i = 0; i < MPCCI_PART_NELEMS(part); i++)
	  {
		  /* get element ID and element information */
		  elemID = i;

		  switch (dim)
		  {
		  case 2:  /* -> line elements */
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 2];
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 2 + 1];
			  elemTypes[i] = MPCCI_ETYP_LINE2;
			  break;
		  case 3:  /* -> quadrilaterals */
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 4];
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 4 + 1];
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 4 + 2];
			  *elemNodePtr++ = wsflist[surfaceID]->elemnodes[i * 4 + 3];
			  //elemTypes[i]   = MPCCI_ETYP_QUAD4;  //change element type by changing macro/lu
			  if (ele_type == 0) {
				  elemTypes[i] = MPCCI_ETYP_QUAD4;
			  }
			  if (ele_type == 1) { //tetrahedral element
				  elemTypes[i] = MPCCI_ETYP_TRIA3;
			  }
			  //linear quad element, need to change to higher order/lu
			  //*elemNodePtr may have to be changed too/lu
			  //element type macros are defined in the file "mpcci elements.h"/lu
			  break;
		  default:
			  MPCCI_MSG_FATAL0("Unsupported element type!");
			  break;
		  }
	  }

      /* Send the element definition for this coupled component */
	  //define the elements of a coupling component/lu
	  //All nodes used in smpcci pels must be already defined with smpcci pnod/lu
      ret = smpcci_pels(
                    server,
                    MPCCI_PART_MESHID(part),    /* mesh id */
                    MPCCI_PART_PARTID(part),    /* part id */
                    MPCCI_PART_NELEMS(part),    /* number of elements */
                    elemTypes[0],               /* first element type */
                    elemTypes,                  /* element types */
                    elemNodes,                  /* nodes of the element */
		            //array of the node numbers (your node numbers) of the element nodes/lu
                    NULL);                      /* element IDs */
	  free(elemTypes);
	  free(elemNodes);
   }
   //It returns 0 for a successful operation/lu
   //The element types can be given in two ways: Either you have a mesh, which consists of elements of the same type. 
   //In this case you can set eltype1 and eltypes[] = NULL, i. e. the list of element types only
   //consists of one type, which is valid for all elements.
   //If you have a mesh which consists of different types you can give the type for each element.
   //Then size of eltypes should have the same value as nelems and the argument eltype1 should be set to 0.
   return ret;
}

/*###########################################################################
 *############  DATA EXCHANGE                                 ###############
 *#########################################################################*/
//driver functions defined below are called by ampcci_transfer function during dotransfer function call (P790)
//The get functions read data from the code and are called before sending data to the partner code via the MpCCI servers. 
//The put functions are called after data was received and must write the data to the code.
//the relay station of that process is the foundation data structure
/*****************************************************************************/
int MpCCI_Driver_getFaceNodeValues( const MPCCI_PART *part,
                                    const MPCCI_QUANT *quant,
                                    void *values)
//get function reads data from user defined code via wsflist/lu 
//details could be found P801
//all information on the coupling component and the requested quantity are given in part and quant/lu
//if your code only supports coupling on faces with nodal values, only MPCCI DRIVER getFaceNodeValues must be defined./lu
/*****************************************************************************/
{
   int surfaceID;                       /* internal ID of coupling component */
   int nnodes;                          /* number of nodes */
   double * valptr = (double*) values; /* we have double values */
   //the values get here is the nodal displacement
   /* CHANGED: Added loop counters and dimension */
   int i, d;                            /* loop counters */
   int dim;                             /* dimension */

   MPCCI_MSG_INFO0("entered get_values...\n");

   /* check whether this is the appropriate method */
   MPCCI_MSG_ASSERT(MPCCI_PART_IS_FACE(part));

   /* get internal ID of the coupling component */
   surfaceID = MPCCI_PART_PARTID(part);

   switch (MPCCI_QUANT_SMETHOD(quant)) //how to store sent/received quantities, one of MPCCI QSM UNDEF invalid send method/lu
   {
      case MPCCI_QSM_DIRECT: /* direct load/store values */ 
		  //MPCCI QSM DIRECT directly written into buffer/lu
		  //that may correspond to the send/receive option in gui.xcf/lu 

         /* number of nodes */
         nnodes = MPCCI_PART_NNODES(part);

         /* distinguish between the different quantities */
         switch (MPCCI_QUANT_QID(quant)) //see page 820 for reference (code API symbol could be found in 943)
         {
			case MPCCI_QID_WALLFORCE:
				dim = wsflist[surfaceID]->dim;
				for (i = 0; i<nnodes; i++)
				{
					for (d = 0; d<3; d++) //force in three directions
					{
						/*
						if (d == wsflist[surfaceID]->direction)
						{
							*valptr++ = wsflist[surfaceID]->nodeforce[i*dim + d];
						} //gives no force value on other direction 
						*/
						*valptr++ = wsflist[surfaceID]->nodeforce[i*dim + d];
						//not sure if nodeforce in MpCCI is defined that way/lu
						//since valptr is a 1D pointer, I think it should be defined just like that used to define node coordinate/lu
					}
				}
				break;
				
			case MPCCI_QID_ABSPRESSURE:
				//dim = wsflist[surfaceID]->dim;
				for (i = 0; i < nnodes; i++)
				{
					*valptr++ = wsflist[surfaceID]->nodepressure[i];
				}
				break;
				
			case MPCCI_QID_RELWALLFORCE:
				dim = wsflist[surfaceID]->dim;
				for (i = 0; i<nnodes; i++)
				{
					for (d = 0; d<3; d++) //force in three directions
					{
						/*
						if (d == wsflist[surfaceID]->direction)
						{
							*valptr++ = wsflist[surfaceID]->nodeforce[i*dim + d];
						} //gives no force value on other direction 
						else
						{
							*valptr++ = 0.0; //2D case: z=0.0 
						}
					    */
						*valptr++ = wsflist[surfaceID]->nodeforce[i*dim + d];
						//not sure if nodeforce in MpCCI is defined that way/lu
						//since valptr is a 1D pointer, I think it should be defined just like that used to define node coordinate/lu
					}
				}
				break;
            default:
               MPCCI_MSG_FATAL0("Quantity not supported\n");
               break;
         }
         break;

      default:
         MPCCI_MSG_FATAL0("Unsupported storage method.");
         break;
   }
   MPCCI_MSG_INFO0("finished get_values...\n");

   return sizeof(double); /* return the size of the value data type */
}

void MpCCI_Driver_putFaceNodeValues(const MPCCI_PART *part,
	const MPCCI_QUANT *quant,
	void *values)
	/*****************************************************************************/
{
	int nnodes;
	int i; int d;
	int surfaceID;                      /* internal ID of coupling component */
	double * valptr = (double*)values;  /* we have double values */
	int dim;

	MPCCI_MSG_INFO0("entered put_values...\n");

	/* check whether this is the appropriate method */

	MPCCI_MSG_ASSERT(MPCCI_PART_IS_FACE(part));

	/* get internal ID of the coupling component */
	surfaceID = MPCCI_PART_PARTID(part);

	/* number of nodes */
	nnodes = MPCCI_PART_NNODES(part);

	switch (MPCCI_QUANT_SMETHOD(quant))
	{
	case MPCCI_QSM_DIRECT: /* direct load/store values */

						   /* distinguish between the different quantities */
		switch (MPCCI_QUANT_QID(quant))
		{
		case MPCCI_QID_NPOSITION:
			dim = wsflist[surfaceID]->dim;
			for (i = 0; i < nnodes; i++) {
				for (d = 0; d < 3; d++)
				{
					wsflist[surfaceID]->nodecoord[i*dim + d] = *valptr++;
					//pointer valptr is the node position, not displacement
				}
			}
				break;
		default:
			MPCCI_MSG_FATAL0("Quantity not supported\n");
			break;
		}
		break;

	default:
		MPCCI_MSG_FATAL0("Unsupported storage method.");
		break;

	}
	MPCCI_MSG_INFO0("finished put_values...\n");
}


//It is important to keep in mind, that sending of data is always possible: The data is stored by MpCCI,
//i. e. it can be received later by the other code. (stored in pointer buffer: valptr)
extern double time_step_size;
static int MpCCI_Driver_getGlobValues(const MPCCI_GLOB *glob, void *values)
{
	static const char module[] = "Get global quantity";
	int qid = MPCCI_QUANT_QID(glob);
	double *valptr = (double *)values;
	MPCCI_MSG_INFO2("%s from variable \"%s\"\n", module, MPCCI_QUANT_NAME(glob));

	switch (qid)
	{
	case MPCCI_QID_TIMESTEP_SIZE:
		//valptr++ = time_step_size; //Provide here the time step size to export 
		*valptr++ = time_step_size;  //Provide here the time step size to export 
		//*valptr++ = 1.3826008770376172e-05;
		break;

	default:
		MPCCI_MSG_FATAL2("%s: variable \"%s\" is not supported.\n", module, MPCCI_QUANT_NAME(glob));
		break;
	}
	return (int) sizeof(double);
}



