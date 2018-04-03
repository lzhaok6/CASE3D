#ifndef __MPCCI4_CCVX_HEADER_INCLUDED
#define __MPCCI4_CCVX_HEADER_INCLUDED
/* mpcci_ccvx.h
 *
 *****************************************************************************************
 *
 * Purpose:
 *    Defines basic structures and CPP macros used with the CCVX trace file format.
 *
 * Author:
 *    Carsten Dehning <carsten.dehning@scai.fhg.de>
 *
 * Reviews/changes:
 *    2010/Jan: Carsten Dehning, Initial release
 *    $Id: mpcci_ccvx.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 *****************************************************************************************
 */

#ifdef __cplusplus
extern "C" {
#endif

#define CCVX_SMAGIC        "CCVX"
#define CCVX_SVERSION_01   "01"
#define CCVX_SVERSION_02   "02"
#define CCVX_SVERSION      CCVX_SVERSION_02

#define CCVX_FLAG_TAGGED   0x01 /* bit is set in case the .ccvx file contains a full implicit step tag */

typedef struct _CCVX_HEADER
{
   char  magic[4];                     /* magic number "CCVX" */
   char  version[2];                   /* '00' .. 'ZZ' */
   char  endianess;                    /* 'I' = Intel or 'M' = Motorola */
   char  flags;                        /* 8 bits (CCVX_FLAG_...) used for some settings */
   int   fileid;                       /* 0=New file, else int id of follower file */
   char  reserved[64-8-sizeof(int)];   /* fill up to 64 bytes */
} CCVX_HEADER;

#ifdef __cplusplus
}
#endif

#endif
