#ifndef __MPCCI4_LIMITS_HEADER_INCLUDED
#define __MPCCI4_LIMITS_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_limits.h
 *
 * Purpose:
 *    Defines some basic limits used to dimension some arrays.
 *
 *    In principle (by design) there are no limits (except memory), neither
 *    for the MpCCI server nor the client library since all objects are managed as chains.
 *    However for performance reasons sometimes the chained objects (or pointers)
 *    temporarily need to be copied/compressd/collected into a simple vector.
 *    In this case the arrays needs to have a dimension which is defined here.
 *
 *    In case of trouble or error messages simply increase the limits and rebuild
 *    the server and client and all adapters.
 *
 * Reviews/changes:
 *    2008/Mar: Carsten Dehning, Initial release
 *    $Id: mpcci_limits.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 *****************************************************************************************
 */

/*
 * Define the max. no. of coupled codes handles be the server.
 * This parameter is only used for dimensioning some temporary arrays used within
 * collect operations.
 */
#define MPCCI_NCODES_MAX       64 /* 64 codes should be enought */


/*
 * Define the max. no. of quantities send or received on a single mesh.
 * This parameter is only used for dimensioning arrays within some structures
 * with fixed length to avoid malloc() and free() and temporary arrays used within
 * the "collect special quantities" operations
 */
#define MPCCI_MESHQ_MAX       32 /* who really exchanges more than 32 quantities? */


/*
 * Define the max. no. of chars used for a mesh/part name.
 * This parameter is only used for dimensioning some char arrays within some structures
 * with fixed length to avoid strdup() and free()
 */
#define MPCCI_MESHNAME_MAX   128


/*
 * Define the max. no. of chars used for quantity names.
 * This parameter is only used for dimensioning some char arrays within some structures
 * with fixed length to avoid strdup() and free()
 */
#define MPCCI_QUANTNAME_MAX   64


/*
 * Define the max. no. of chars used for the auxialiary string in the MPCCI_PINFO.
 * This parameter is only used for dimensioning some char arrays within some structures
 * with fixed length to avoid strdup() and free()
 */
#define MPCCI_AUXSTRING_MAX  256


/*
 * Define the max. no. of faces on a polyhedron cell (PCELL)
 */
#define MPCCI_PCELLFACES_MAX  256


/*
 * Define the max. no. of vertices of a polygonal face (POLY)
 */
#define MPCCI_POLYVERTICES_MAX  256


/*
 * Define the max. no. of user pointers in the MPCCI_QUANT & MPCCI_PART
 */
#define MPCCI_USERPOINTERS_MAX        8

#endif
