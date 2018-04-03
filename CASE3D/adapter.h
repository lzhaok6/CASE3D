/****************************************************************************
 *   adapter.h                                                              *
 *                                                                          *
 *   MpCCI Code API  -  Code Adapter Header File                            *
 *   Changed for use with the C version of "foundation"                     *
 *                                                                          *
 *   This file contains declarations of all driver functions.               *
 *   See the MpCCI Programmers Guide for more information.                  *
 *                                                                          *
 *      http://www.scai.fraunhofer.de/mpcci.html                            *
 *                                                                          *
 ****************************************************************************/
#ifndef __MPCCI_ADAPTER_HEADER_INCLUDED__
#define __MPCCI_ADAPTER_HEADER_INCLUDED__

/* use the mpcci code API */
#include "mpcci.h"

#if __cplusplus
extern "C" {
#endif

/*#########################################################################
 *############  INTERFACE FUNCTIONS                         ###############
 *#########################################################################*/
//interface functions are defined in mpcci.h/lu
/* initialize the coupling procedure */
void initcoupling();

/* transfer quantities */
void dotransfer();

/* exit the coupling procedure */
void exitcoupling();

/*#########################################################################
 *############  DRIVER FUNCTIONS                            ###############
 *#########################################################################*/
/* TODO: Remove those declarations you do not need -- do not modify any! */
//The implementation of those driver functions are defined in adapter.c/lu


/* CHANGED: Removed all methods called before/after some actions */
/* Update MPCCI_PART and /or update MPCCI_QUANT */
int MpCCI_Driver_partUpdate(MPCCI_PART *part, MPCCI_QUANT *quant);

/****************************************************************************
 *       COUPLING COMPONENT DEFINITIONS
 ****************************************************************************/
/* TODO: Remove those declarations you do not need -- do not modify any! */

/* REQUIRED: Method 1: Define the grid */
int MpCCI_Driver_definePart(MPCCI_SERVER *server, MPCCI_PART *part);

/* CHANGED: Removed Method 2: Define the grid */
/* CHANGED: removed MpCCI_Driver_moveNodes */

/****************************************************************************
 *       DATA EXCHANGE
 ****************************************************************************/
/****************************************************************************
 * Required are only those methods which are needed to exchange quantities.
 * TODO: Remove those you do not need.
 ****************************************************************************/
/* CHANGED: Removed all but two exchange functions */

/* get data from nodes of a face - called before send */
int MpCCI_Driver_getFaceNodeValues(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);

/* put data from elements of a face - called after receive */
//void MpCCI_Driver_putFaceElemValues(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);

/* put data from nodes of a face - called after receive */
void MpCCI_Driver_putFaceNodeValues(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);

static int MpCCI_Driver_getGlobValues(const MPCCI_GLOB *glob, void *values);

#if __cplusplus
}
#endif

#endif
