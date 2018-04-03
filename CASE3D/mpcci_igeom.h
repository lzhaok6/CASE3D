#ifndef __MPCCI4_IGEOM_HEADER_INCLUDED__
#define __MPCCI4_IGEOM_HEADER_INCLUDED__
/* mpcci_igeom.h
 *
 *****************************************************************************************
 *
 * Purpose:
 *    Defines the structures and functions for the indexed geometry.
 *    The indexed geometry data structures and functions are used by the morpher and
 *    its file converters.
 *
 * Author:
 *    Carsten Dehning <carsten.dehning@scai.fhg.de>
 *
 * Reviews/changes:
 *    2004/Jan/   Carsten Dehning, Initial release
 *    $Id: mpcci_igeom.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 *****************************************************************************************
 */
#include <stdio.h>
#include "mpcci_elements.h"

#ifdef __cplusplus
extern "C" {
#endif

#define IGEOM_MAJOR_VERSION  1
#define IGEOM_MINOR_VERSION  4


#define IFILE_HEADER_FMT  \
      "<gmd-magic"      \
      " version=%d.%d"  \
      " format=%c"      \
      " partial=%c"     \
      " nodesize=%d"    \
      " facesize=%d"    \
      " cellsize=%d"    \
      " cplesize=%d />\n"

#define IFILE_GROUP_WFMT "<gmd-%s num=%d min=%d max=%d />\n"
#define IFILE_GROUP_RFMT "gmd-%s num=%%d min=%%d max=%%d"
#define IFILE_GROUP_EOF  "\n<gmd-eof />\n"
#define IFILE_GROUP_EFMT "<gmd-%3s />\n"


#define FOREACH_IGEOM(_obj,_list)                  for(_obj=(_list)->head; _obj; _obj=_obj->next)
#define FOREACH_IGEOM_OFTYPE(_obj,_list,_type)     for(_obj=(_type *)((_list)->head); _obj; _obj=_obj->next)
#define FOREACH_IGEOM_NEXT(_obj,_list,_next)       for(_obj=(_list)->head; _obj; _obj=_next)

#define FOREACH_IVERT(_obj,_list)                  FOREACH_IGEOM_OFTYPE(_obj,_list,IVERT)
#define FOREACH_IFACE(_obj,_list)                  FOREACH_IGEOM_OFTYPE(_obj,_list,IFACE)
#define FOREACH_ICELL(_obj,_list)                  FOREACH_IGEOM_OFTYPE(_obj,_list,ICELL)
#define FOREACH_ICPLE(_obj,_list)                  FOREACH_IGEOM_OFTYPE(_obj,_list,ICPLE)


typedef struct _IGEOM         IGEOM;
typedef struct _IVERT         IVERT;
typedef struct _IFACE         IFACE;
typedef struct _ICELL         ICELL;
typedef struct _ICPLE         ICPLE;
typedef struct _IGEOM_LIST    IGEOM_LIST;
typedef void                  IGEOMVIS; /* in fact a socket */

struct _IGEOM /* basic object with first two values */
{
   IGEOM *next;         /* next object in chain */
   int    index;        /* id of the object */
   int    group;        /* a key to identify a group to which the object belongs */
};

#define IVERT_SIZE      (3*sizeof(double)+sizeof(int))
struct _IVERT           /* derived from IGEOM */
{
   IVERT *next;         /* next IVERT in chained list */
   int    index;        /* index>=0 needed if we have gaps in indexing */
   int    group;        /* a key to identify a group to which the object belongs */

   double coord[3];     /* xyz */
};

#define IFACE_SIZE   (4*sizeof(int))
struct _IFACE           /* derived from IGEOM */
{
   IFACE *next;         /* next IFACE in chained list */
   int    index;        /* check only: no gaps allowed, index MUST be 0,1,2,3....nface-1 */
   int    group;        /* a key to identify a group to which the object belongs */

   int   *fnodes;       /* malloc'ed face nodes vector */
   int    ftype;        /* face type */
   int    btype;        /* boundary type IFACE_TYPE_xxxxx */
};

#define ICELL_SIZE    (11*sizeof(int)) /* 8 nodes minimum + 3*int */
struct _ICELL           /* derived from IGEOM */
{
   ICELL *next;         /* next ICELL in chained list */
   int    index;        /* element index */
   int    group;        /* a key to identify a group to which the object belongs */

   int   *cnodes;       /* vertex indices of the nvarr in case of polyhedrons */
   int    ctype;        /* cell type */
};

#define ICPLE_SIZE (2*sizeof(int))
struct _ICPLE           /* derived from IGEOM */
{
   ICPLE *next;         /* next ICPLE in chained list */
   int    index;        /* couple id */
   int    group;        /* a key to identify a group to which the object belongs */

   int    masterFace;   /* index of master face */
   int    childFace;    /* index of child face */
};


struct _IGEOM_LIST
{
   IGEOM   *head;
   IGEOM  **vect;
   void   (*igeom_save)  (const void *g, FILE *fp, int binary);
   int    (*igeom_cmp)   (const void *g1, const void *g2);
   int      numGeom;
   int      minGeomId;
   int      maxGeomId;
};

#define IFACE_TYPE_INTERIOR   'I' /* free moveable interior face */
#define IFACE_TYPE_FIXEDBND   'F' /* fixed boundary face */
#define IFACE_TYPE_FLOATBND   'S' /* sliding boundary (symmetry or cyclic) */

#define IFACE_NNODES(_f)      CAST_INT(MPCCI_ETYP_NNE((_f)->ftype))
#define IFACE_ISPOLY(_f)      MPCCI_ETYP_IS_POLY((_f)->ftype)

#define ICELL_NNODES(_c)   (int)MPCCI_ETYP_NNE((_c)->ctype)
#define ICELL_ISPOLY(_c)        MPCCI_ETYP_IS_PCELL((_c)->ctype)

/*
 * Special treatment of polyhedrons:
 *
 *    INTEGER*4 NVARR(NV) layout in a .cel file is taken from the
 *    pro-STAR V4 documentation
 *
 *    dimension: nvarr[0...nv-1]
 *
 *    nvarr is splitted into two parts
 *       a) nvarr[0 ... nvarr[0]-1]   lookup table (header part of nvarr)
 *       b) nvarr[nvarr[0] .... nv]   vertex index list
 *
 *    a)
 *       nvarr[0]          =>    start index of 1st vertex of face[0] in nvarr
 *       nvarr[1]          =>    start index of 1st vertex of face[1] in nvarr
 *       nvarr[2]          =>    start index of 1st vertex of face[2] in nvarr
 *       .........
 *
 *    b)
 *       nvarr[nvarr[0]  ] =>    vertex id of 1st vertex of face 0
 *       nvarr[nvarr[0]+1] =>    vertex id of 2nd vertex of face 0
 *       nvarr[nvarr[0]+2] =>    vertex id of 3rd vertex of face 0
 *       .........
 *       nvarr[nvarr[1]  ] =>    vertex id of 1st vertex of face 1
 *       nvarr[nvarr[1]+1] =>    vertex id of 2nd vertex of face 1
 *       .........
 *
 *    nvarr[0] - 1      => no. of faces in polyhedron
 */
#define IPOLY_NFACES(_nvarr)     ((_nvarr)[0] - 1 )             /* no. of faces in polyhedron */
#define IPOLY_FNVERT(_nvarr,_i)  ((_nvarr)[_i+1]-(_nvarr)[_i])  /* no. of vertices of face i */
#define IPOLY_ISTART(_nvarr,_i)  (         (_nvarr)[_i] )       /* index   of 1st vertex of face i */
#define IPOLY_VSTART(_nvarr,_i)  ((_nvarr)[(_nvarr)[_i]])       /* value   of 1st vertex of face i */
#define IPOLY_PSTART(_nvarr,_i)  ((_nvarr)+(_nvarr)[_i] )       /* pointer to 1st vertex of face i */



#if !defined(INCLUDE_STATIC) || !INCLUDE_STATIC

   /*
    * All functions are simply declared extern.
    * Therefore we cannot create a MSWin dll which whould require
    *    _decspc(dllexport) during library compilation
    *    _decspc(dllimport) during external usage/compilation
    */

   extern void    IGEOM_Delete         (IGEOM *g);
   extern void    IGEOM_List_init      (IGEOM_LIST *igl,
                                          void (*igeom_save)(const void *g, FILE *fp, int binary),
                                          int  (*igeom_cmp) (const void *g1, const void *g2)
                                       );
   extern void    IGEOM_List_destroy   (IGEOM_LIST *igl);
   extern IGEOM  *IGEOM_List_append    (IGEOM_LIST *igl, IGEOM *g);
   extern void    IGEOM_List_save      (IGEOM_LIST *igl, const char *type, FILE *fp, int binary);
   extern int     IGEOM_List_groups    (IGEOM_LIST *igl, int *grpCount, int size);
   extern void    IGEOM_Vect_create    (IGEOM_LIST *igl, const char *type, int reindex);
   extern void    IGEOM_Vect_reindex   (IGEOM_LIST *igl);
   extern void    IGEOM_List_rebuild   (IGEOM_LIST *igl);
   extern IGEOM  *IGEOM_Find_index     (IGEOM_LIST *igl, int index);
   extern int     IGEOM_compare_index  (const void *c1, const void *c2);


   extern IVERT  *IVERT_New            (int index, double coord[3]);
   extern void    IVERT_Save           (const void *g, FILE *fp, int binary);


   extern IFACE  *IFACE_New            (int index, int ftype, const int *inodes);
   extern IFACE  *IFACE_VNew           (int index, int ftype, ...);
   extern void    IFACE_Save           (const void *g, FILE *fp, int binary);
   extern void    IFACE_Delete_set     (IFACE *facev[], int nfaces);
   extern int     IFACE_List_groups    (IGEOM_LIST *igl, int *grpCount, int size);
   extern int     IFACE_List_types     (IGEOM_LIST *igl, int *vertCount, int maxVert);
   extern int     IFACE_Vect_compress  (IGEOM_LIST *igl);
   extern IFACE  *IFACE_Vect_find      (IGEOM_LIST *igl, int ftype, int *inodes);
   extern int     IFACE_compare        (const void *_f1, const void *_f2);
   extern IFACE  *IFACE_setup          (IFACE *face, const char *caller);


   extern ICELL  *ICELL_New            (int index, int group, int ctype, int *inodes, IGEOM_LIST *vlist);
   extern void    ICELL_Save           (const void *g, FILE *fp, int binary);
   extern void    ICELL_List_types     (IGEOM_LIST *igl, int cellCount[5]);
   extern void    ICELL_List_simplify  (IGEOM_LIST *igl);


   extern ICPLE  *ICPLE_New            (int masterFace, int childFace);
   extern void    ICPLE_Save           (const void *g, FILE *fp, int binary);

   extern size_t  IFILE_Write_bin      (FILE *fp, const void *objBuff, size_t objCount, size_t objSize);
   extern size_t  IFILE_Read_bin       (FILE *fp,       void *objBuff, size_t objCount, size_t objSize, int swap);
   extern void    IFILE_Write_header   (FILE *fp, const char *name, int numVal, int minIdx, int maxIdx);
   extern FILE   *IFILE_Write_open     (const char *modelName, int binary, int partial);
   extern void    IFILE_Write_intv     (FILE *fp, const int *ibuf, int n, const char *trailer);
   extern FILE   *IFILE_Read_open      (const char *modelName, int *valid);
   extern int     IFILE_Check          (const char *modelName, const char *extv[]);

   extern IGEOMVIS  *IGEOMVIS_Open     (int port);
   extern void       IGEOMVIS_Close    (IGEOMVIS *vissock);
   extern void       IGEOMVIS_Step     (IGEOMVIS *vissock, const int step);
   extern void       IGEOMVIS_Boundary (IGEOMVIS *vissock, const char *pname, const char *bname, IGEOM_LIST *vert_list, IGEOM_LIST *face_list, const int group, const int isblk);


   extern int     memindex_int         (const int *buf, size_t n, int seek);
   extern void   *memrot               (      void *buf, size_t n, size_t size, int nrot);
   #if 0 /* these are included locally */
   extern int    *memrev_int           (      int  *buf, size_t n);
   extern void   *memswapb             (      void *buf, size_t n, size_t size);
   extern int     getendianess         (void);
   extern time_t  fxmodtime            (const char *pathname, const char *extension);
   #endif

#endif

#define IVERT_List_init(igl)              IGEOM_List_init(igl,IVERT_Save,IGEOM_compare_index)
#define IFACE_List_init(igl)              IGEOM_List_init(igl,IFACE_Save,IFACE_compare)
#define ICELL_List_init(igl)              IGEOM_List_init(igl,ICELL_Save,IGEOM_compare_index)
#define ICPLE_List_init(igl)              IGEOM_List_init(igl,ICPLE_Save,IGEOM_compare_index)

#define IVERT_List_append(igl,vert)       IGEOM_List_append(igl,(IGEOM *)(vert))
#define IFACE_List_append(igl,face)       IGEOM_List_append(igl,(IGEOM *)(face))
#define ICELL_List_append(igl,cell)       IGEOM_List_append(igl,(IGEOM *)(cell))
#define ICPLE_List_append(igl,cple)       IGEOM_List_append(igl,(IGEOM *)(cple))

#define IVERT_List_save(igl,fp,binary)    IGEOM_List_save(igl,"nodes"  ,fp,binary)
#define IFACE_List_save(igl,fp,binary)    IGEOM_List_save(igl,"faces"  ,fp,binary)
#define ICELL_List_save(igl,fp,binary)    IGEOM_List_save(igl,"cells"  ,fp,binary)
#define ICPLE_List_save(igl,fp,binary)    IGEOM_List_save(igl,"couples",fp,binary)

#ifdef __cplusplus
}
#endif
#endif
