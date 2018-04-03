#ifndef __MPCCI4_SERVER_HEADER_INCLUDED
#define __MPCCI4_SERVER_HEADER_INCLUDED
/* mpcci_server.h
 *
 *****************************************************************************************
 *
 * Purpose:
 *    Defines all classes/structures used within the server and the operators
 *
 * Author:
 *    Carsten Dehning <carsten.dehning@scai.fhg.de>
 *
 * Reviews/changes:
 *    2008/Mar: Carsten Dehning, Initial release
 *    $Id: mpcci_server.h 944 2014-11-05 14:20:01Z pbayrasy $
 *
 *****************************************************************************************
 */
#include "mpcci_types.h"
#include "mpcci_chain.h"
#include "tmat44.h"
#include "ipool.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * QUANT_MAXBUF should be at least 2, but NEVER(!) > 32, since somewhere we use a
 * 32 bit unsigned bitpattern variable which holds a bit per buffer.
 */
#define QUANT_MAXBUF    32


/*
 * these PART flags bits are set whenever a PART is updated due to remeshing or mesh definitions
 */
#define PART_FLAG_DEFINED        0x00000001  /* part was newly defined */
#define PART_FLAG_NODEIDS        0x00000002  /* part got new node ids */
#define PART_FLAG_COORDS         0x00000004  /* part got new coords assigned */
#define PART_FLAG_ELEMIDS        0x00000008  /* part got new element ids */
#define PART_FLAG_ELNODES        0x00000010  /* part got new elements */
#define PART_FLAG_ELTYPES        0x00000020  /* part got new element types */
#define PART_FLAG_ELSIZES        0x00000040  /* part got new element sizes */
#define PART_FLAG_ELNORMS        0x00000080  /* part got new face element normal vectors */
#define PART_FLAG_ELNOFFS        0x00000100  /* part got new element nodes offsets assigned */
#define PART_FLAG_ANGLES         0x00000200  /* part got new anodal point angles assigned */
#define PART_FLAG_MOTION         0x00000400  /* part got a new motion assigned */
#define PART_FLAG_UPDATED      (PART_FLAG_DEFINED|PART_FLAG_NODEIDS|PART_FLAG_COORDS|PART_FLAG_ELTYPES|PART_FLAG_ELNODES|PART_FLAG_ANGLES|PART_FLAG_MOTION)
#define PART_FLAG_REMESHED       0x00000800  /* part was remeshed */

#define PART_FLAG_CHECKED        0x00010000 /* free to use check bit */

/*
 * this bit is set if the part needs to be ...
 */
#define PART_FLAG_TSAVE_COOR     0x00100000 /* save coords in the trace file */
#define PART_FLAG_TSAVE_ELEM     0x00200000 /* save elems in the trace file */
#define PART_FLAG_TSAVE_MESH     0x00300000 /* save coors + elems in the trace file */
#define PART_FLAG_TSAVE_DEFINED  0x00400000 /* true if part was defined to all trace file writers */

#define PART_FLAG_MSEND_COOR     0x01000000 /* send coords to a monitor */
#define PART_FLAG_MSEND_ELEM     0x02000000 /* send elems to a monitor */
#define PART_FLAG_MSEND_MESH     0x03000000 /* send coords + elems to a monitor */
#define PART_FLAG_MSEND_PDEF     0x04000000 /* true if part-plot has to be defined to the monitor */


enum
{
   /* basic states */
   STATE_NULL  = 0,
   STATE_RUN      ,

   /* mesh/part states */
   STATE_DEFINE   ,
   STATE_ASSIGNED ,
   STATE_CONFIGURED ,

   /* various QUANT states after operator->execute() was called */
   STATE_RECEIVED ,  /* IN : all parallel clients finished the PUTQ */
   STATE_PREFILT  ,  /* IN : the STATE_PUTQ quantity was prefiltered: ready for mapping */
   STATE_MAPPED   ,  /* OUT: the quantity was just mapped */
   STATE_FILLED   ,  /* OUT: all node values of the STATE_MAPPED quantity were filled */
   STATE_POSTFILT ,  /* OUT: all node values of the STATE_FILLED quantity were post-filtered */
   STATE_INTEG    ,  /* OUT: all node values of the STATE_POSTFILT quantity were integrated */
   STATE_ADJUST   ,  /* OUT: all node values of the STATE_INTEG quantity were adjusted */
   STATE_READY    ,  /* OUT: limited quantity ready for a GETQ by a client */
   STATE_USED     ,  /* OUT: all QTID type clients(!) finished the GETQ */

   /*
    * all socket blocking sync points or wait states for clients start here!
    */
   STATE_BLOCKING = 100,

   /* sync points: locked until all client have called this function */
   STATE_SYNC     ,  /* sync a code & convergence */
   STATE_SYNT     ,  /* set a time/iter tag and sync */
   STATE_SYNM     ,  /* sync meshes */
   STATE_SYNJ     ,  /* sync the job (all clients) */

   STATE_GETC     ,  /* sync request of partner code information */
   STATE_QTAG     ,  /* sync request of qtag information */
   STATE_LOCK     ,  /* client is temporary locked and must wait */

   STATE_LAST
};

#define CLIENT_IS_BLOCKED(_client)  ((_client)->state > STATE_BLOCKING)


/*
 * cast to avoid compiler complains on data pointers:
 *
 *    (void **)(&(pointer))
 */
#define PPVOID_ADDR(_p)          ( (void *)(&(_p)) )


/* free if pointer is non NULL and assign NULL to the pointer */
#define FREE_NN(_p)        if (_p) { FREE(_p); _p = NULL; }



/*
 * all classes used by the server
 */
typedef struct _JOB           JOB;
typedef struct _MONITOR       MONITOR;
typedef struct _CODE          CODE;
typedef struct _CNAMES        CNAMES;
typedef struct _CLIENT        CLIENT;
typedef struct _GLOB          GLOB;
typedef struct _RMESH         RMESH;
typedef struct _MESH          MESH;
typedef struct _PART          PART;

typedef struct _TAG           TAG;
typedef struct _STAG          STAG;
typedef struct _QTAG          QTAG;
typedef struct _WTAG          WTAG;
typedef struct _SYNC          SYNC;

typedef struct _OPERAT        OPERAT;
typedef struct _QINFO         QINFO;
typedef struct _QUANT         QUANT;
typedef struct _QMON          QMON;
typedef struct _OINFO         OINFO;

/*
 * SUNOS platform special: sys/stream.h #define QNORM 0x00
 * Undef and pray :-)
 */
#ifdef QNORM
   #undef QNORM
#endif
typedef struct _QNORM         QNORM;

typedef struct _RELSET        RELSET;  /* a mesh-mesh relationship */
typedef struct _MAPSET        MAPSET;  /* a mapping set */

typedef struct _WRITER        WRITER;

typedef struct _CODE_QINFO    CODE_QINFO; /* server internal representation of clients MPCCI_QINFO */
typedef struct _CODE_PINFO    CODE_PINFO; /* server internal representation of clients MPCCI_PINFO */
typedef struct _CODE_GINFO    CODE_GINFO; /* server internal representation of clients MPCCI_GLOB */

typedef struct _VNODE         VNODE;      /* helper struct: node coords + values vector */
typedef struct _ECLASS        ECLASS;     /* element class info table */
typedef struct _VELEM         VELEM;      /* element information + values */


typedef struct _CMIRROR       CMIRROR;

/*
 * Structure used to transport part specific periodicity information
 */
typedef struct _PERIODIC  PERIODIC;
struct _PERIODIC
{
   double   	axis[3]; 	/* rotational axis */
   double   	origin[3]; 	/* axis origin */
   int   		n;    		/* number of periodicities*/
   PART			*pOrig;		/*OriginalPart, if NULL, then it is the original Part itself*/
   PART			**pRep;		/*List of children part pointers from Original Part, NULL if Child*/
};


/*
 * Simplified name of a MPCCI_MOTION object + type bitmasks
 * Motion types
 *    nibbles: |0000|0000|0000|0000|0000|0000|0000|0000|
 *             |----|----|----|---Q|---C|-ORD|--AR|--AT|
 *    T(0/1)      : if this bit is set, we have a valid translation (MPCCI_MOTION.velo[3])
 *    R(0/1)      : if this bit is set, we have a valid rotation    (MPCCI_MOTION.omeg/axis[3]/orig[3])
 *    A(0/1)      : if this bit is set, omega/velocity are in fact accelerations
 *    ORD(0/1/2/3): if both R&T are set, these 2 bits (1/2/3) describe the order of R/T
 *    C(0/1)      : if this bit is set, coordinates send to the server are already transformed
 *    Q(0/1)      : if this bit is set, vector quantities send to the server are already transformed
 */
typedef MPCCI_MOTION    MOTION;


/* The bit patterns below MUST MATCH with the MPCCI_MOTION_... bitpatterns */
#define MOTION_MASK_T            0x00000001
#define MOTION_MASK_R            0x00000010
#define MOTION_MASK_ORD          0x00000f00
#define MOTION_MASK_RT           0x00000111  /* rotate, then translate order */
#define MOTION_MASK_TR           0x00000211  /* translate, then rotate order */
#define MOTION_MASK_CO           0x00000311  /* simultanious continuous transformations */
#define MOTION_MASK_TYPE         0x00000fff
#define MOTION_MASK_COORD        0x00001000
#define MOTION_MASK_QUANT        0x00010000

#define MOTION_MASK_ACC_T        MPCCI_MOTION_ACC_T /* velocity is in fact an acceleration */
#define MOTION_MASK_ACC_R        MPCCI_MOTION_ACC_R /* omega is in fact an angular acceleration */
#define MOTION_MASK_ACC          MPCCI_MOTION_ACC

#define MOTION_HAS_MASK(_t,_m)   (  (_t) & _m )
#define MOTION_HAS_T(_t)         (  (_t) & MOTION_MASK_T )
#define MOTION_HAS_R(_t)         (  (_t) & MOTION_MASK_R )
#define MOTION_IS_CO(_t)         ( ((_t) & MOTION_MASK_CO) == MOTION_MASK_CO )


/*
 * Structure which holds the command string and a pointer to a
 * client mirror function inside the server.
 * This mapping between command token and mirror function is(must) static(!) defined
 * inside CLIENT_server().
 */
struct _CMIRROR
{
   /*
    * The client->server 4 bytes command string MPCCI_CMD_XXXX, interpreted as
    * either a 0-terminated char[5] string or a 4 byte unsigned
    */
   union mirrorcmd
   {
      char     s[8]; /* 4 chars + trailing 0 */
      unsigned u;    /* unsigned variant of the first 4 bytes */
   } cmd;

   /* pointer to the client mirror_xxxx() function */
   int (*fnc)(CLIENT *_client, const char *_cmd);

   /* if true, then a halted client must not call the mirror function */
   int  halt_check;
};


/*
 * Helper struct to transport orphan level information
 */
struct _OINFO
{
   char    *olevel;  /* pointer to an array[mesh|part->nnodes|nelems] */
   unsigned flags;   /* currently 0=node base olevel, else element based olevel */
   int      nvals;   /* mesh|part->nnodes|nelems depending on the orphan flags */
};

#define OINFO_FLAG_NODE    0x00000000
#define OINFO_FLAG_ELEM    0x00000001

#define OINFO_FLAG_IS_NODE(_oinfoflags) ((_oinfoflags) == OINFO_FLAG_NODE)
#define OINFO_FLAG_IS_ELEM(_oinfoflags) ((_oinfoflags) == OINFO_FLAG_ELEM)

#define OINFO_IS_NODE(_oinfo) OINFO_FLAG_IS_NODE((_oinfo)->flags)
#define OINFO_IS_ELEM(_oinfo) OINFO_FLAG_IS_ELEM((_oinfo)->flags)

#define OINFO_IS_EMPTY(_oinfo) ((_oinfo)->olevel == NULL)
#define OINFO_IS_VALID(_oinfo) ((_oinfo)->olevel != NULL)


/*
 * Step tag: Just a TAG, but with different namings
 */
struct _STAG
{
   double st_time;   /* time >= 0.0 in case of a transient run */
   int    st_iter;   /* iteration counter, either global for steady or inner iteration count */
   int    st_csid;   /* increasing only step number (=counter) */
};

/*
 * Temprary tag to keep the clients synt tag data
 */
struct _SYNC
{
   double   sync_time;  /* coupling time >= 0.0 in case of a transient explicit/implicit run */
   double   sync_dt;    /* if (dt>0), this is the prospective dt of a code/client */
   int      sync_iter;  /* relative iteration counter (used in case of iterative coupling) */
   int      sync_conv;  /* CLIENT convergence state at this time/iter */
   unsigned sync_flags; /* CLIENT flags */
};


#define TAG_BASE_DATA \
   double time;      /* coupling time >= 0.0 in case of a transient explicit/implicit run */\
   double dt;        /* if (dt>0), this is the prospective dt of a code/client */\
   int    iter;      /* relative iteration counter (used in case of iterative coupling) */\
   int    qtid;      /* increasing only tag id (=counter) */\
   int    csid;      /* only increasing coupling step counter */\
   int    iter_a;    /* absolute iteration counter of a code */

#define TAG_BASE_SIZE            (2*sizeof(double)+4*sizeof(int))
#define TAG_BASE_EQUAL(_t1,_t2)  !memcmp(_t1,_t2,TAG_BASE_SIZE)
#define TAG_BASE_DIFFER(_t1,_t2)  memcmp(_t1,_t2,TAG_BASE_SIZE)
#define TAG_BASE_COPY(_d,_s)      memcpy(_d,_s,TAG_BASE_SIZE)


/* The base TAG type */
struct _TAG
{
   TAG_BASE_DATA
   int  conv;   /* CODE convergence state of this tag */
};

#define TAG_COPY(_d,_s)        memcpy(_d,_s,sizeof(TAG))


/*
 * Quantity tag: TAG + quantity state
 */
struct _QTAG
{
   TAG_BASE_DATA

   int    conv;   /* CODE convergence state of this tag */
   int    state;  /* state of the associated QUANT */

   /*
    * Current type of the quantity (location & integral type)
    * The current type changes due do mapping/filling and integration.
    * The value is QTYPE_...
    */
   int    qtype;
   double inorm;  /* calculated iterative norm */
};

/*
 * Writer tag: TAG + writer state state
 */
struct _WTAG
{
   TAG_BASE_DATA
   int state;  /* state of the associated QUANT */
};


#define TAG_QTID_IS_USED(_qtid)  ( (_qtid) >= 0 ) /* a negative tag id indicates an unused tag */
#define TAG_QTID_IS_FREE(_qtid)  ( (_qtid) <  0 )

#define TAG_TIME_IS_USED(_time)  ( (_time) >= 0 ) /* a negative tag id indicates an unused tag */
#define TAG_TIME_IS_FREE(_time)  ( (_time) <  0 )

#define TAG_ITER_IS_USED(_iter)  ( (_iter) >= 0 ) /* a negative tag id indicates an unused tag */
#define TAG_ITER_IS_FREE(_iter)  ( (_iter) <  0 )

#define TAG_IS_USED(_t)          TAG_QTID_IS_USED((_t)->qtid)
#define TAG_IS_FREE(_t)          TAG_QTID_IS_FREE((_t)->qtid)

#define TAG_HAS_TIME(_t)         TAG_TIME_IS_USED((_t)->time)
#define TAG_NOT_TIME(_t)         TAG_TIME_IS_FREE((_t)->time)

#define TAG_HAS_ITER(_t)         TAG_ITER_IS_USED((_t)->iter)
#define TAG_NOT_ITER(_t)         TAG_ITER_IS_FREE((_t)->iter)

#define TAG_HAS_DATE(_t)         ( TAG_TIME_IS_USED((_t)->time) || TAG_ITER_IS_USED((_t)->iter) )
#define TAG_NOT_DATE(_t)         ( TAG_TIME_IS_FREE((_t)->time) && TAG_ITER_IS_FREE((_t)->iter) )

#define TAG_IS_TIME(_t)          ( TAG_IS_USED(_t) && TAG_HAS_TIME(_t) )  /* is a time tag */
#define TAG_IS_QTID(_t)          ( TAG_IS_USED(_t) && TAG_NOT_DATE(_t) )  /* is a serial tag */
#define TAG_IS_DATE(_t)          ( TAG_IS_USED(_t) && TAG_HAS_DATE(_t) )  /* time || iter */
#define TAG_IS_IMPL(_t)          ( TAG_HAS_TIME(_t) && TAG_HAS_ITER(_t) )

/*
 * Symbolic type of tag depends in the values of tag->time and tag->iter
 */
enum
{
   TAG_TYPE_FREE = 0, /* tag not set */
   TAG_TYPE_QTID = 1, /* simply steps */
   TAG_TYPE_ITER = 2, /* iterative */
   TAG_TYPE_EXPL = 3, /* transient explicit */
   TAG_TYPE_IMPL = 4  /* transient implicit */
};

#define TAG_TYPE_BYVAR(_time,_iter)\
(\
   TAG_TIME_IS_FREE(_time) ? TAG_ITER_IS_FREE(_iter) ? TAG_TYPE_QTID : TAG_TYPE_ITER\
                           : TAG_ITER_IS_FREE(_iter) ? TAG_TYPE_EXPL : TAG_TYPE_IMPL\
)

#define TAG_CTYPE_BYVAR(_time,_iter)\
(\
   TAG_TIME_IS_FREE(_time) ? TAG_ITER_IS_FREE(_iter) ? "QTID" : "ITER"\
                           : TAG_ITER_IS_FREE(_iter) ? "EXPL" : "IMPL"\
)

#define TAG_AXNAME_BYVAR(_time,_iter)\
(\
   TAG_TIME_IS_FREE(_time) ? TAG_ITER_IS_FREE(_iter) ? "Coupling step" : "Iteration number"\
                           : "Time [s]"\
)

#define TAG_TYPE(_t)                ( TAG_IS_FREE(_t) ? TAG_TYPE_FREE : TAG_TYPE_BYVAR ((_t)->time,(_t)->iter) )
#define TAG_CTYPE(_t)               ( TAG_IS_FREE(_t) ? "FREE"        : TAG_CTYPE_BYVAR((_t)->time,(_t)->iter) )
#define TAGTYPE_MISMATCH(_t1,_t2)   ( ((_t1)==TAG_TYPE_FREE || (_t2)==TAG_TYPE_FREE) ? 0 : ((_t1) != (_t2)) )
#define TAGTYPE_MATCH(_t1,_t2)      ( ((_t1)==TAG_TYPE_FREE || (_t2)==TAG_TYPE_FREE || ((_t1) == (_t2)) )

/*
 * Helper struct: nodal point with values + coordinates
 */
struct _VNODE
{
   const real *coords;  /* coordinates[1|2|3] of this node */
   const real *values;  /* pointer to the values[QDIM] at this node */
};


/*
 * Table which describes an element class (name, etyp, function pointers)
 */
struct _ECLASS
{
   /* get the element size, length, area or volume */
   real (*elem_esize)(const VELEM *ve);

   /* get the nodal size, length, area or volume */
   void (*elem_nsize)(const VELEM *ve, real *nsize);

   /* integrate a quantity over the element length, area or volume */
   void (*elem_einte)(const VELEM *ve, double *integral);

   /* get the local element coordinates from global coordinates */
   void (*get_xsi_from_coords)(const real coords[], real xsi[]);

   /* test if local coordinates are EXACTLY inside the element */
   int (*is_xsi_inside_exact)(const real xsi[]);

   /* test if local coordinates are inside a boundary around the element */
   int (*is_xsi_inside_bounds)(const real xsi[], const real esize, const real tol);

   /* force the local coordinates to be inside the element */
   void (*move_xsi_inside)(real xsi[]);

   /* reorder the nodes according to the internal scheme */
   void (*reorder_nodes)(int *nodeid);

   /*
    * to be continued .....
    */

   unsigned etyp;
   unsigned csys;
};


/*
 * All information required to describe an element and possibley nodal values for this element.
 *
 *    This is a VARIABLE ELEMENT structure used to transport element specific information
 *    down into element basic functions like
 *
 *       ECLASS->elem_esize()
 *       ECLASS->elem_nsize()
 *       ECLASS->elem_einte()
 *       etc.
 *
 *    The pre-settings (required information for input) strongly depends on the called function.
 *    The post-settings (output) depends on the called function
 */
struct _VELEM
{
   const int  *elnoffs; /* pointer to the part->elnoffs+? only required for PCELLS! */
   VNODE pv[MPCCI_POLYVERTICES_MAX]; /* help array used internally */
   real        esize;   /* element size (length, aera or volume) calculated only once */
   unsigned    etyp;    /* element type MPCCI_ETYP_*, only requred for PCELL and POLY */
   int         qdim;    /* quantity dimension used by the element integrator */
};

/*
 * Macros mainly to be used in functions which have an argument VELEM or VNODE pv[]
 * VNODE *pv = VELEM->pv;
 */
#define PV_C(_i)      pv[_i].coords
#define PV_X(_i)      pv[_i].coords[0]
#define PV_Y(_i)      pv[_i].coords[1]
#define PV_Z(_i)      pv[_i].coords[2]
#define PV_V(_i,_j)   pv[_i].values[_j]


/*
 * mesh/mesh relationship operator set
 */
struct _RELSET
{
   /*
    * The ->state variable of the CHAIN_DEFINE_CLASSBASE() is in fact
    * the no. of references to this mesh relationship.
    */
   CHAIN_DEFINE_CLASSBASE(RELSET)

   MESH    *smesh;   /* source mesh for this relationship */
   MESH    *tmesh;   /* target mesh for this relationship */
   OPERAT  *o_rel;   /* the relation operator */
   void    *h_rel;   /* returned handle of o_rel->create() */

   /*
    * The orphaned/childless information is created/allocated by the relation operator
    * o_rel->create() and/or o_rel->execute() and can later be used by the
    * filler operator (o_fill->create()) for the quantities.
    * The fille is called just after the mapping.
    *
    * These guys may be NULL if there are no orphans or the relation operator
    * does not deliver this information.
    *
    * The arrays must be free'd() by the o_rel->destroy()
    */
   const char *olevel_n;    /* weak+orphaned nodes on the target mesh: char orphans_n[tmesh->nnodes] */
   const char *childless_n;  /* childless elems on the source mesh    : char childless_n[smesh->nnodes] */
   const char *olevel_e;    /* weak+orphaned elems on the target mesh: char orphans_e[tmesh->nelems] */
   const char *childless_e;  /* childless nodes on the source mesh    : char childless_e[smesh->nelems] */

   /* these gyus are updated in RELSET_validate() */
   int norphans_n;   /* no. of true orphaned nodes */
   int nchildless_n; /* no. of childless nodes */
   int norphans_e;   /* no. of true orphaned elements */
   int nchildless_e; /* no. of childless elements */
};

enum
{
   OLEVEL_PERFECT_MATCH = 0,  /* properly mapped */
   OLEVEL_CLOSE_MATCH,        /* week mapped */
   OLEVEL_FILLED_DIFFUSION,   /* values filled by diffusion filler */
   OLEVEL_FILLED_DISTANCE,    /* values filled by distance filler */
   OLEVEL_FILLED_NEAREST,     /* values filled by nearest filler */
   OLEVEL_FILLED_AVERAGE,     /* values filled by average filler */
   OLEVEL_FILLED_DEFAULT,     /* values filled by default value */
   OLEVEL_FULL_ORPHANED,      /* fully orphaned */
   OLEVEL_MAXDIM              /* always last, used for array dimensions */
};

#define CHILDLESS_LEVEL_U      0x00 /* properly used */
#define CHILDLESS_LEVEL_N      0x01 /* not used */

#define RELSET_FLAG_NODE      0x00000001
#define RELSET_FLAG_ELEM      0x00000002
#define RELSET_FLAG_COORD     0x00000010
#define RELSET_FLAG_FIELD     0x00000020
#define RELSET_FLAG_FLUXI     0x00000040
#define RELSET_FLAG_FLUXD     0x00000080

#define RELSET_FLAG_CONFORM   0x00000100 /* conformal mesh relationship only */

#define RELSET_FLAG_BIDIR     0x00100000 /* bidirectional: can also be used for a tmesh->smesh relation */

#define RELSET_IS_BIDIR(_rs)           ( ((_rs)->flags & RELSET_FLAG_BIDIR) != 0 )
#define RELSET_IS_UNIDIR(_rs)          ( ((_rs)->flags & RELSET_FLAG_BIDIR) == 0 )

#define RELSET_HAS_ROPR(_rs,_ropr)     ( (_rs)->o_rel == _ropr )
#define RELSET_HAS_SMESH(_rs,_mesh)    ( (_rs)->smesh == _mesh || (RELSET_IS_BIDIR(_rs) && (_rs)->tmesh == _mesh) )
#define RELSET_HAS_TMESH(_rs,_mesh)    ( (_rs)->tmesh == _mesh || (RELSET_IS_BIDIR(_rs) && (_rs)->smesh == _mesh) )
#define RELSET_HAS_MESH(_rs,_mesh)     ( (_rs)->smesh == _mesh || (_rs)->tmesh == _mesh )

#define RELSET_HAS_MESHES(_rs,_smesh,_tmesh)\
(\
   (                        (_rs)->smesh == _smesh && (_rs)->tmesh == _tmesh ) ||\
   (RELSET_IS_BIDIR(_rs) && (_rs)->smesh == _tmesh && (_rs)->tmesh == _smesh )\
)

#define RELSET_HAS_SCODE(_rs,_code)  ( (_rs)->smesh->code == _code || (RELSET_IS_BIDIR(_rs) && (_rs)->tmesh->code == _code) )
#define RELSET_HAS_TCODE(_rs,_code)  ( (_rs)->tmesh->code == _code || (RELSET_IS_BIDIR(_rs) && (_rs)->smesh->code == _code) )
#define RELSET_HAS_CODE(_rs,_code)   ( (_rs)->smesh->code == _code || (_rs)->tmesh->code == _code )

#define RELSET_HAS_CODES(_rs,_scode,_tcode)\
(\
   (                        (_rs)->smesh->code == _scode && (_rs)->tmesh->code == _tcode ) ||\
   (RELSET_IS_BIDIR(_rs) && (_rs)->smesh->code == _tcode && (_rs)->tmesh->code == _scode )\
)

#define RELSET_MATCH(_rs,_smesh,_tmesh,_opr) (RELSET_HAS_MESHES(rset,smesh,_tmesh) && RELSET_HAS_ROPR(rset,_opr))


/*
 * mesh/mesh mapping operator set
 */
struct _MAPSET
{
   CHAIN_DEFINE_CLASSBASE(MAPSET)

   OPERAT  *o_map;   /* the mapping operator */
   void    *h_map;   /* returned handle of o_map->create() */

   /*
    * The mapper has an associated MESH-MESH relationship:
    *    source mesh: this->rset->smesh
    *    target mesh: this->rset->tmesh
    */
   RELSET  *rset;
};

#define MAPSET_FLAG_CONFORM   0x00000001


/*
 * monitor object + flags
 */
struct _MONITOR
{
   CHAIN_DEFINE_CLASSBASE(MONITOR)

   JOB            *job;       /* job which is monitored */
   MPCCI_CONNECT  *connect;   /* connection for a monitor process in client mode  */
   OPERAT         *o_monitor; /* up to now just a copy of job->o_monitor */
   const char     *axisname;  /* name of xyplot abszissa : time or iteration */
   STAG            steptag;   /* current step tag */
   STAG            plottag;   /* step tag for xy plots */
};

struct _JOB
{
   CHAIN_DEFINE_CLASSBASE(JOB)

   /* jobid is a unique generated name to identify all clients of this job */
   char     *jobid;

   CODE     *codes;           /* chain of all codes participating in this job */

   RELSET   *rsets;           /* chain of all mesh->mesh relationships sets */
   MAPSET   *msets;           /* chain of all mesh->mesh mapping sets */

   MONITOR  *monitors;        /* chain of connected monitors */
   OPERAT   *o_monitor;       /* the monitor is valid for all meshes & quantities */

   WRITER   *writers;         /* chain of tracefile writers */

   char     *writers_home;    /* root path for all the writers subdirectories */
   WTAG      writers_step;    /* step size between write actions */
   WTAG      writers_next;    /* next write action */

   double    time_tolerance;  /* a value indicating that t_diff/delta_t < time_tolerance means identical time */
   unsigned  monitor_flags;   /* commons flag bits of all monitors */
   int       quant_nbuffer;   /* no. of buffers per quantity */
};


/*
 * job flags
 */
#define JOB_FLAG_MCHECK          0x00000001 /* let the job machine check the meshes */
#define JOB_FLAG_MJOIN           0x00000002 /* let the job machine join the meshes */
#define JOB_FLAG_KEEPPARTS       0x00000010 /* do not delete parts, but keep them as empty parts (used by some writers!) */
#define JOB_FLAG_MESHMOTION      0x00000020 /* if true do the motions for all codes/meshes */
#define JOB_FLAG_BAFFLESHIFT     0x00000040 /* allow the baffle shift for face meshes */
#define JOB_FLAG_TSAVE           0x00000100 /* bit set if we have to write to the trace file */
#define JOB_FLAG_OVERLAPCHECK    0x00000200 /* check the size ratio and overlap of mesh-pairings during join */
#define JOB_FLAG_DOMAINCHECK     0x00000400 /* check the mesh/parts for multi domains or isolated stray elements */
#define JOB_FLAG_MRFCORRECT      0x00000800 /* correct the outgoing vector quantities in case of transient MRF */
#define JOB_FLAG_IZEROORDER      0x00001000 /* zero order interpolation */
#define JOB_FLAG_IFIRSTORDER     0x00002000 /* first order interpolation */
#define JOB_FLAG_IHIGHORDER      0x00004000 /* high order interpolation */
#define JOB_FLAG_EZEROORDER      0x00010000 /* zero order extrapolation */
#define JOB_FLAG_EFIRSTORDER     0x00020000 /* first order extrapolation */
#define JOB_FLAG_EHIGHORDER      0x00040000 /* high order extrapolation */
#define JOB_FLAG_MOTIONUSEINITNBH 0x00080000 /* keep initial relationship for mesh motion */

#define JOB_KEEPS_PARTS(_job)    ( (_job)->flags & JOB_FLAG_KEEPPARTS   )
#define JOB_HAS_MESHMOTION(_job) ( (_job)->flags & JOB_FLAG_MESHMOTION  )

#define JOB_IS_IZEROORDER(_job)  ( (_job)->flags & JOB_FLAG_IZEROORDER  )
#define JOB_IS_IFIRSTORDER(_job) ( (_job)->flags & JOB_FLAG_IFIRSTORDER )
#define JOB_IS_IHIGHORDER(_job)  ( (_job)->flags & JOB_FLAG_IHIGHORDER  )

#define JOB_IS_EZEROORDER(_job)  ( (_job)->flags & JOB_FLAG_EZEROORDER  )
#define JOB_IS_EFIRSTORDER(_job) ( (_job)->flags & JOB_FLAG_EFIRSTORDER )
#define JOB_IS_EHIGHORDER(_job)  ( (_job)->flags & JOB_FLAG_EHIGHORDER  )

/*
 * static table with code specific names
 */
struct _CNAMES
{
   const char *ctyp;
   const char *cobj[4][3];
};


struct _CODE
{
   CHAIN_DEFINE_CLASSBASE(CODE)

   CNAMES      *nametab;   /* pointer to statically defined name table */
   JOB         *job;       /* link to my parent job */
   CLIENT      *clients;   /* chain of all the clients of this code */
   MESH        *meshes;    /* chain of all meshes for this code */
   GLOB        *globs_r;   /* chain of all globals received by the client */
   CODE_QINFO  *cqinfos;   /* chain of code specific quantity properties */
   CODE_PINFO  *cpinfos;   /* chain of code specific coupled parts */
   CODE_GINFO  *cginfos;   /* chain of code specific coupled globals */

   const CODE  *limit_code;/* pointer to the limiting code */

   double       time_curr; /* abolute physical time set by the SYNT cmd */
   double       time_last; /* absolute physical time at the last time step */
   double       time_init; /* absolute physical time when code loggon on */
   int          iter_curr; /* abolute currenten iteration counter */
   int          iter_last; /* absolute iteration counter at last SYNT */
   int          iter_init; /* absolute iteration when code loggend on */

   TAG          tag;       /* current coupling tag of the code */

   /*
    * parameters for the mesh transformation (rotation + translation):
    * relevant if ((flags&MPCCI_CFLAG_TFM44) != 0).
    */
   TMAT44   tmat44_s;   /* send direction: client->server transformation */
   TMAT44   tmat44_r;   /* recv direction: server->client transformation */

   /*
    * parameters for the mesh length scaling:
    * relevant if ((flags&MPCCI_CFLAG_SCALE) != 0).
    */
   real     scale;

   int      nice;
};

#define CODE_print(_code)        CHAIN_print(_code,"Code","ID(%d)",", nice(%d), clients(%d), type(%s).\n",(_code)->nice,CHAIN_getcount((_code)->clients),CODE_C_TYPE(_code))
#define CODE_IS_PARALLEL(_code)  ((_code)->clients->next != NULL)
#define CODE_IS_SERIAL(_code)    ((_code)->clients->next == NULL)
#define CODE_IS_TRANSIENT(_code) ((_code)->time_init >= 0 )
#define CODE_IS_STEADY(_code)    ((_code)->time_init <  0 )
#define CODE_IS_IMPL(_code)      (CODE_IS_TRANSIENT(_code) && (_code)->iter_init >= 0 )
#define CODE_IS_EXPL(_code)      (CODE_IS_TRANSIENT(_code) && (_code)->iter_init <  0 )
#define CODE_IS_ITER(_code)      (CODE_IS_STEADY(_code)    && (_code)->iter_init >= 0 )
#define CODE_IS_QTID(_code)      (CODE_IS_STEADY(_code)    && (_code)->iter_init <  0 )


#define CODE_C_TYPE(_code)             ( (_code)->nametab->ctyp )

#define CODE_C_NODE(_code)             ( (_code)->nametab->cobj[0][2] )
#define CODE_C_NODES(_code)            ( (_code)->nametab->cobj[0][1] )
#define CODE_C_NNODES(_code)           ( (_code)->nametab->cobj[0][0] )

#define CODE_C_ELEM(_code)             ( (_code)->nametab->cobj[3][2] )
#define CODE_C_ELEMS(_code)            ( (_code)->nametab->cobj[3][1] )
#define CODE_C_NELEMS(_code)           ( (_code)->nametab->cobj[3][0] )

#define CODE_C_DELEM(_code,_mdim)      ( (_code)->nametab->cobj[_mdim][2] )
#define CODE_C_DELEMS(_code,_mdim)     ( (_code)->nametab->cobj[_mdim][1] )
#define CODE_C_DNELEMS(_code,_mdim)    ( (_code)->nametab->cobj[_mdim][0] )



#define FOREACH_JOBCLIENT(_client,_code,_job)\
   FOREACH(_code,_job->codes)\
      FOREACH(_client,_code->clients)

struct _CLIENT
{
   CHAIN_DEFINE_CLASSBASE(CLIENT)

   CODE           *code;         /* link to the parent code which has the meshes */
   MPCCI_CONNECT  *connect;      /* link to the client connection */
   GLOB           *globs_s;      /* chain of all globals send by the client */

   /* waiting for a read quantity on a part */
   QUANT          *wait_quant;   /* NULL or a pointer to a QUANT/PART the client is waiting for */
   PART           *wait_part;    /* NULL or a pointer to a QUANT/PART the client is waiting for */
   unsigned        wait_flags;


   /* waiting for a sync'ed tagv request */
   MPCCI_QTAG     *wait_tagv;    /* NULL or a temporary malloc'ed tagv[] the client is waiting for */
   int             wait_tagc;    /* the no. of tags in temporary malloc'ed wait_tagv[] */

   /* halted client */
   const CMIRROR  *wait_mirror;  /* NULL or pointer to the cliens mirror map in case of a halted client */

   SYNC            sync;         /* data to keep clients synch'd data: required in case of a parallel code */

   unsigned        realsize;     /* dynamic defined realsize for the last requested real array */
   int             nprocs;       /* no. of parallel processes assigned to this client (on the client side) */

   /* Current requested server command from client */
   const CMIRROR  *last_cmd;     /* NULL or pointer to the cliens mirror map in case of new request */
};


struct _RMESH
{
   CHAIN_DEFINE_CLASSBASE(RMESH)

   void *mesh;
};

/*
 * Each node has a dynamic chain of neighbor nodes assigned
 */
typedef struct _NBNODE        NBNODE;
typedef struct _NBREL         NBREL;

struct _NBNODE
{
   NBNODE *next;     /* next in chain */
   real    weight;   /* square distance or element size to the neighbor node/element */
   int     index;    /* global index of node/element */
};

struct _NBREL
{
   MESH    *mesh;    /* mesh for the operator */
   NBNODE **nbheads; /* array[mesh->nnodes|nelems] of (NBNODE *) */
   IPOOL    ipool;   /* item pool with NBNODEs */
};


/*
 * Max no. of domains used to store information
 */
#define MPCCI_NDOMAINS_MAX  64

struct _MESH
{
   CHAIN_DEFINE_CLASSBASE(MESH)

   CODE     *code;      /* reference to the parent code */
   PART     *parts;     /* chain of partions */
   QUANT    *quants_s;  /* chain of send(by client)/received(by server) quantities */
   QUANT    *quants_r;  /* chain of received(by client)/send(by server) quantities */

   RMESH    *rmeshes;   /* chain of other representations of this mesh used by specific relationship operators (e.g. maplib relation) */

   NBREL    *rel_nn;    /* optional node-node neighborhood relation used by an an orphan filler */
   NBREL    *rel_ee;    /* optional elem-elem neighborhood relation used by an an orphan filler */
   NBREL    *rel_ne;    /* optional node-elem neighborhood relation used by an an orphan filler */

   ECLASS  **eclasses;  /* pointers to the element class[mesh->nelems] */
   real     *sizes_e;   /* optional: elem size (length/area/volume), may be calculated on demand */
   real     *sizes_n;   /* optional: nodal weight/size, may be calculated on demand */

   /* master-slave node-node part boundary relationships in case of a multipart mesh */
   int      *master_n;  /* global index of master node for each node */

   int      *domains_e;  /* optional: the domain index of each element[mesh->nelems] */
   int       domains_c[MPCCI_NDOMAINS_MAX]; /* no. of elements per domain */
   int       ndomains;   /* no. of element domains */

   /* temporary pointer to the currently accessed values and orphan level */
   OINFO     q_oinfo;   /* temporarily used mesh local copy of QUANT->oinfov[bid] */
   real     *q_values;  /* temporarily used mesh local pointer to QUANT->valuev[bid] */
   real     *l_values;  /* temporarily used mesh local pointer to QUANT->last_values */
   int       q_qtype;   /* temporarily used mesh local pointer to QUANT->qtypev[bid] */

   /*
    * all values below are duplicated or integrated values
    * for all parts of this mesh.
    */
   real      tsize;     /* total volume/area/length of this mesh == sum(part->tsize) */

   real      xmin,xmax; /* bounding box min. and max. coordinates for all parts */
   real      ymin,ymax;
   real      zmin,zmax;
   real      lmin,lmax,lmean; /* min/max/mean distance between nodes within the same element */

   int       nnodes;    /* sum(nnodes)   for all valid parts of this mesh */
   int       nelems;    /* sum(nelems)   for all valid parts of this mesh */
   int       nelnodes;  /* sum(nelnodes) for all valid parts of this mesh */
   int       nelverts;  /* sum(nelverts) for all valid parts of this mesh */
   int       npcells;   /* sum(npcells)  for all valid parts of this mesh */
   int       nelbad;    /* no. of bad/degenerated elements with a size <= 0 */

   unsigned  csys;      /* common parts coordinates system */
   int       mdim;      /* common part dimension: -1=unknown, 0=point, 1=line, 2=surface, 3=volume */
};

#define MESH_IS_NONE(_m)         ( (_m)->mdim == MPCCI_MDIM_NONE  )
#define MESH_IS_PNTS(_m)         ( (_m)->mdim == MPCCI_MDIM_PNTS  )
#define MESH_IS_LINE(_m)         ( (_m)->mdim == MPCCI_MDIM_LINE  )
#define MESH_IS_FACE(_m)         ( (_m)->mdim == MPCCI_MDIM_FACE  )
#define MESH_IS_VOLU(_m)         ( (_m)->mdim == MPCCI_MDIM_VOLU  )
#define MESH_IS_IPNT(_m)         ( (_m)->mdim == MPCCI_MDIM_IPNT  )
#define MESH_IS_PNT1(_m)         ( (_m)->mdim == MPCCI_MDIM_PNTS && (_m)->nnodes == 1 )

/* MESH is any kind of single point, either integration of coordinate point */
#define MESH_IS_POINT1(_m)       (\
                                   ( (_m)->mdim == MPCCI_MDIM_IPNT || (_m)->mdim == MPCCI_MDIM_PNTS )\
                                    && (_m)->nnodes == 1\
                                 )


struct _PART
{
   CHAIN_DEFINE_CLASSBASE(PART)

   CLIENT   *client;       /* the client to which this partition belongs */
   MESH     *mesh;         /* the parent mesh to which this partition belongs */

   int      *nodeids;      /* optional: node id array, if NULL, the ids are 0,1,2,3,... */
   real     *coords;       /* required: 3D=x,y,z,... 2D=x,y,... */
   real     *angles;       /* optional point rotational orientation */

   int      *elemids;      /* optional: element id array, if NULL, the ids are 0,1,2,3,... */
   int      *elnodes;      /* required: element topology = nodeids received */
   int      *elnoffs;      /* derived : mapping of elnodes->local indices, may be calculated in the server */
   unsigned *eltypes;      /* required: contains the type of each element[nelems] */
   real     *normals_f;    /* optional: face normals, may be calculated in the server */

   real     *q_values;     /* temporarily used part local pointer to QUANT->valuev[bid] */
   real     *l_values;     /* temporarily used mesh local pointer to QUANT->last_values */
   OINFO     q_oinfo;      /* temporarily used part local copy of QUANT->oinfov[bid] */

   MOTION    motion;       /* part motion definition */
   PERIODIC  periodic;	   /* part cyclic symmetry definition */

   TAG       ctag;         /* tag when coords were defined */

   real      bthick1;      /* optional shift of a baffle or shell surface */

   real      tsize;        /* total volume/area/length of this part == sum(elsizes) */
   real      xmin,xmax;    /* bounding box min. and max. coordinates */
   real      ymin,ymax;
   real      zmin,zmax;
   real      lmin,lmax,lmean; /* min/max/mean distance between nodes within the same element */

   unsigned  cetype;        /* common element type: !0, then all elements are of the common type */

   /* allocated no. of nodes/elements for this part */
   int       nnodes_a;     /* nodeids[nnodes_a],coords[nnodes_a*cdim],... */
   int       nelems_a;     /* eltypes[nelems_a],elsizes[nelems_a],normals[nelems_a*cdim],... */
   int       nelnodes_a;   /* total number of element nodes: elnodes[nelnodes_a], elnoffs[nelnodes_a] */

   /* helper counters: index offsets for this part compared to all valid parts below in this chain */
   int       nnodes_o;     /* sum(nnodes) for all VALID(!) parts before this part */
   int       nelems_o;     /* sum(nelems) for all VALID(!) parts before this part */
   int       nelnodes_o;   /* sum(nelnodes) for all VALID(!) parts before this part */

   /* currently used no. of nodes/elements */
   int       nnodes;       /* nodeids[0..nnodes],coords[0..nnodes*cdim] */
   int       nelems;       /* eltypes[0..nelems],elsizes[0..nelnodes[0..nelnodes], elnoffs[0..nelnodes] */
   int       nelnodes;     /* total number of element nodes: elems],normals[0..nelems*cdim] */

   /* helper counters */
   int       nelverts;     /* total no. of geometrical vertices for all elements (nelverts <= nelnodes) */
   int       npcells;      /* no. of polyhedrons <= nelems */
   int       nelbad;       /* no. of bad/degenerated elements with a size <= 0 */

   /* helper int IDs for the part */
   int       cpos;         /* local part position within the chain mesh->parts (0,1,2,3,4....) */
   int       guid;         /* globally unique part id (0,1,2,...) set inside PART_new(). never changes */
   int       wuid;         /* writers unique part id changes whenever the part has been remeshed */

   unsigned  csys;         /* coordinates system */
   int       mdim;         /* mesh dimension: -1=unknown, 0=point, 1=line, 2=surface, 3=volume */
};

#define PART_ELEMID(_p,_i)             ( ((_p)->elemids) ? (_p)->elemids[_i] : _i )
#define PART_IS_VALID(_p)              ( (_p)->nnodes >  0 )
#define FOREACH_VALID_PART(_p,_head)   FOREACH_COND(_p,_head,PART_IS_VALID(_p))
#define PART_HAS_MOTION(_p)            ( (_p)->motion.type != 0 )
#define PART_IS_ORIG_PERIODIC(_p)      ( ((_p)->periodic.pOrig == NULL) && ((_p)->periodic.n > 1) )
#define PART_IS_REP_PERIODIC(_p)       ( (_p)->periodic.pOrig != NULL )
#define PART_IS_PERIODIC(_p)           ( (_p)->periodic.n > 1 )

/*
 * global variable
 */
struct _GLOB
{
   CHAIN_DEFINE_CLASSBASE(GLOB)

   /* pointer to the associated QINFO */
   const QINFO *qinfo;

   /*
    * Only valid for globals received by the client: These GLOBs
    * get a qtag assigned to indicate when the glob was last updated.
    */
   QTAG  qtag;

   /*
    * parameters of the transformation:
    * relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   real  scale;
   real  offs;

   /* default fill value if not assigned */
   real  vfill;
   real  values[MPCCI_QDIM_MAX+1];

   int   nrecv; /* countdown: how many clients received this global */
};



/*
 * the operator is the operators config string + pointers to the operator functions
 */
typedef void       *(OPX_CREATE )(const char *config, const void *some_handle, ... );
typedef int         (OPX_EXECUTE)(void *some_object, ...);
typedef int         (OPX_DESTROY)(void *handle);
typedef const char *(OPX_GETINFO)(int ityp);

struct _OPERAT
{
   /*
    * the name of this object is the name of the operator library without extension,
    * the id is irrelevant, can only use CHAIN_getbyname()
    */
   CHAIN_DEFINE_CLASSBASE(OPERAT)

   char        *config;    /* configuration string for the operator */
   OPX_CREATE  *create;
   OPX_EXECUTE *execute;
   OPX_DESTROY *destroy;
};


/*
 * The operator index (OID) is an int containing a bitpattern
 *
 *    nibbles: |0000|0000|0000|0000|0000|0000|0000|0000|
 *             |----|OSIG|ODIM|SDIM|TDIM|-----CSYS-----|
 *             |----|---OTYP--|SDIM|TDIM|-----CSYS-----|
 *
 * The master operator valid for ALL coordinate systems and source and target mesh dimensions
 * has an SDIM=TDIM=CSYS=0000 pattern.
 *
 *    MASTER : |----|OSIG|ODIM|0000|0000|0000|0000|0000|
 */
#define MPCCI_OID_ODIM1                      0x00100000
#define MPCCI_OID_ODIM2                      0x00200000

#define MPCCI_OID_MAKE(_otyp,_csys,_sdim,_tdim)   ( (_otyp) | ((_sdim)<<16) | ((_tdim)<<12) | (_csys) )

#define MPCCI_OID_CSYS(_oid)           (  (_oid)      & 0x00000fff )
#define MPCCI_OID_TDIM(_oid)           ( ((_oid)>>12) & 0x0000000f )
#define MPCCI_OID_SDIM(_oid)           ( ((_oid)>>16) & 0x0000000f )
#define MPCCI_OID_ODIM(_oid)           (  (_oid)      & 0x00f00000 )
#define MPCCI_OID_OTYP(_oid)           (  (_oid)      & 0x0ff00000 ) /* OSIG + ODIM */

#define MPCCI_OID_IS_MASTER(_oid)      ( ((_oid) & 0x000fffff) == 0 )

/*
 * - the bit pattern for OTYP only
 * - the bit pattern for the master operator with (SDIM=TDIM=CSYS=0)
 */
#define MPCCI_OTYP_PRE                 0x01100000  /* single mesh preprocess operator */
#define MPCCI_OTYP_REL                 0x02200000  /* dual mesh relationship operator */
#define MPCCI_OTYP_MAP                 0x03200000  /* dual mesh mapping operator */
#define MPCCI_OTYP_FILL                0x04100000  /* single mesh filler operator */
#define MPCCI_OTYP_POST                0x05100000  /* single mesh postprocess operator */
#define MPCCI_OTYP_INTG                0x06100000  /* single mesh postprocess integrator operator */
#define MPCCI_OTYP_ADJ                 0x07100000  /* single mesh postprocess adjuster operator */
#define MPCCI_OTYP_LIM                 0x08100000  /* single mesh postprocess limiter operator */
#define MPCCI_OTYP_MON                 0x0ff00000  /* monitor operator */

/*
 * call destroy method for an operator of an object
 * if the handle is NOT NULL
 */
#define OPR_DESTROY_NNULL(_obj,_onm)               \
   if ((_obj)->o_##_onm && (_obj)->h_##_onm)       \
   {                                               \
      (_obj)->o_##_onm->destroy((_obj)->h_##_onm); \
      (_obj)->h_##_onm = NULL;                     \
   }


/*
 * The quantity info (QINFO) contains the memory representation of the
 * general quantities configuration setup file.
 */
struct _QINFO
{
   CHAIN_DEFINE_CLASSBASE(QINFO)

   const char *unit; /* optional string with the units of this quantity */
   OPERAT     *ops;  /* chain(!) of configured operators for this quantity */
};


/*
 * The code specific quantity settings (see struct _MPCCI_QINFO).
 */
struct _CODE_QINFO
{
   CHAIN_DEFINE_CLASSBASE(CODE_QINFO)

   CODE *code; /* the owning code */

   /*
    * parameters of the transformation relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   real  scale;
   real  offs;

   /* default orphaned fill value */
   real  vfill;

   /*
    * parameters for relaxation with iterative coupling
    * true-value(n+1) = relax * value(n+1) + (1-relax) * values(n)
    */
   real  relax;      /* < 0 = unused, else 0.00001 < relax < 2 */

   /*
    * parameters for ramping (e.g. forces)
    *
    * Underrelaxation:  true-value(n) = (ramp0 + n*rampd) * value(n)  for (ramp0+n*rampd) <  1
    * Overrelaxation:   true-value(n) = (ramp0 - n*rampd) * value(n)  for (ramp0-n*rampd) >  1
    */
   real  ramp0;      /* < RAMP0_MIN: unused, else small start value RAMP0_MIN...RAMP0_MAX */
   real  rampd;      /* ramp delta within the range [RAMPD_MIN...RAMPD_MAX] */

   real  maxitol;    /* max. tolerance in case of iterative coupling */

   int   mid;        /* id of the coupled mesh for this quantity setting */

   int   method_s;   /* storage locations in index used by the adapters */
   int   index_s;
   int   method_r;
   int   index_r;

   int   maxiter;    /* max. no. of iterations in case of iterative coupling */
};


/*
 * limits on some constants
 */
#define RAMP0_MIN       0.0
#define RAMP0_MAX       2.0
#define RAMPD_MIN       0.001
#define RAMPD_MAX       1.0

#define MAXITOL_MIN     1e-9
#define MAXITOL_MAX     0.999
#define MAXITER_MIN     1
#define MAXITER_MAX     1000

/*
 * The code specific coupled part definitions (see struct _MPCCI_PINFO)
 * The part id pid is the object id.
 */
struct _CODE_PINFO
{
   CHAIN_DEFINE_CLASSBASE(CODE_PINFO)

   char        *paux;                        /* if not NULL, pointer to the part auxiliary string */

   const QINFO *qinfov_s[MPCCI_MESHQ_MAX+1]; /* list of quantities the client should send */
   const QINFO *qinfov_r[MPCCI_MESHQ_MAX+1]; /* list of quantities the client should receive */

   int          mid;                         /* id of the mesh (relevant for the server + the client code)*/
   int          mdim;                        /* mesh dimension MPCCI_MDIM_* */
};


/*
 * The code specific coupled glob definitions (see struct _MPCCI_GLOB)
 * The global id gid is the object id.
 */
struct _CODE_GINFO
{
   CHAIN_DEFINE_CLASSBASE(CODE_GINFO)

   char        *gaux;   /* if not NULL, pointer to the part auxiliary string */
   const QINFO *qinfo;  /* pointer to the associated QINFO */

   /*
    * parameters of the transformation relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   real  scale;
   real  offs;
   real  vfill; /* default fill value */
};

struct _QBUF /* NOT YET USED... LATER */
{
   TAG_BASE_DATA

   OINFO    oinfo;   /* memdup'ed orphan level array after mapping[nalloc_v] */

   /* integral/mean values for this mesh: [QUANT->nbuffer][qdim] */
   double  *mean_s;  /* source integrals */
   double  *mean_t;  /* my target integrals */
   real    *values;  /* malloc'ed buffer[nalloc_v*qdim*sizeof(real)] */
   int     *state_p; /* state of each part portion of the value[nalloc_p] */
   int      state_q; /* state of the associated QUANT */

   /*
    * Current type of the quantity (location & integral type)
    * The current type changes due do mapping/filling and integration.
    * The value is QTYPE_...
    */
   int      qtype;
};


/*
 * quantity is in fact the implementation of QINFO
 */
struct _QUANT
{
   CHAIN_DEFINE_CLASSBASE(QUANT)

   /* pointer to the associated QINFO */
   const QINFO *qinfo;

   /* pointer to the parent mesh */
   MESH  *mesh;

   /*
    * pointers to operators (with the config strings) and handles to the ->o_name->create() results
    * NO chained objects(!): these are just pointers to the operator selected from
    * the QUANT->qinfo->ops and assigned to this quantity based on the mesh and
    * quantity settings - done within QUANT_getoperat().
    */

   OPERAT *o_pre;    /* only valid for a quantity sent by the client */
   void   *h_pre;

   OPERAT *o_fill;   /* only valid for a quantity received by the client */
   void   *h_fill;

   OPERAT *o_post;   /* only valid for a quantity received by the client */
   void   *h_post;

   OPERAT *o_adj;    /* only valid for a quantity received by the client */
   void   *h_adj;

   OPERAT *o_lim;    /* only valid for a quantity received by the client */
   void   *h_lim;

   MAPSET *mset;     /* pointer to the mapping set for this quantity */


   /* integral values for this mesh: [QUANT->nbuffer][qdim] */
   double  *sintegv[QUANT_MAXBUF];  /* source integrals */
   double  *tintegv[QUANT_MAXBUF];  /* my target integrals */

   real    *valuev[QUANT_MAXBUF];   /* malloc'ed buffers: [QUANT->nbuffer][nalloc_v*qdim)] */
   int     *statev[QUANT_MAXBUF];   /* state of each part portion of the valuev [QUANT->nbuffer][nalloc_p] */
   OINFO    oinfov[QUANT_MAXBUF];   /* memdup'ed orphan level array after mapping: [nalloc_v] */

   real    *last_values;            /* optional buffer containing the values of the last iteration */

   /* interpolated values information */
   real    *vint_values;            /* buffer containing the values from the last mirror_getq() for relaxation */
   real    *vint_norelax;           /* buffer containing the values from the last mirror_getq() without relaxation */
   real    *vint_current;           /* buffere containing unrelaxed values (from current getq) */
   OINFO   *vint_oinfo;             /* pointer to the orphans of the last_vint */
   real     vint_relax;             /* relaxation parameter used for the vint_values */
   real     relax_aitken;
   QTAG     vint_qtag;              /* qtag for the last interpolation */
   double  *vint_integ;             /* current integral values for this mesh */

   /*
    * the tagv contains the time tag information per buffer
    * 'nbuffer' is the no. of buffers/tags to be used for this quantity.
    */
   QTAG     qtagv[QUANT_MAXBUF];

   /* ltagv contains the tags of a already mapped and STATE_USED quantity from the source side */
   TAG      ltagv[QUANT_MAXBUF];

   int      nalloc_v;   /* allocated values buffer size: max(this->mesh->nnodes,this->mesh->nelems) */
   int      nalloc_p;   /* allocated no. of parts for the statev[QUANT_MAXBUF][nalloc_p] */

   /*
    * parameters of the transformation:
    * relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   real  scale;
   real  offs;

   /*
    * Default value for orphaned points. Currently relevant only
    * for fields (not fluxes).
    */
   real  vfill;

   /*
    * parameters for relaxation with iterative coupling
    * true-value(n+1) = relax * value(n+1) + (1-relax) * values(n)
    */
   real  relax;   /* < 0 = unused, else 0.00001 < relax < 2 */

   /*
    * parameters for ramping (e.g. forces)
    *
    * Underrelaxation:  true-value(n) = (ramp0 + n*rampd) * value(n)  for (ramp0+n*rampd) <  1
    * Overrelaxation:   true-value(n) = (ramp0 - n*rampd) * value(n)  for (ramp0-n*rampd) >  1
    */
   real  ramp0;   /* < RAMP0_MIN: unused, else small start value RAMP0_MIN...RAMP0_MAX */
   real  rampd;   /* ramp delta within the range [RAMPD_MIN...RAMPD_MAX] */

   real  maxitol; /* max. tolerance in case of iterative coupling */
   int   maxiter; /* max. no. of iterations in case of iterative coupling */

   int   nbuffer; /* no. of history buffers/tags assigned to this quantity */
};


/*
 * Bitpattern for the qtype location
 */
#define QTYPE_LOC_MASK              0x0000000f
#define QTYPE_LOC_NODE              0x00000001
#define QTYPE_LOC_ELEM              0x00000002

/*
 * Bitpattern for the qtype interpolation/mapping/integration
 */
#define QTYPE_INT_MASK              0x000000f0
#define QTYPE_INT_COORD             0x00000010
#define QTYPE_INT_FIELD             0x00000020
#define QTYPE_INT_FLUXI             0x00000040
#define QTYPE_INT_FLUXD             0x00000080

#define QTYPE_UNDEF                 0x00000000

#define QTYPE_N_COORD               ( QTYPE_LOC_NODE | QTYPE_INT_COORD )
#define QTYPE_N_FIELD               ( QTYPE_LOC_NODE | QTYPE_INT_FIELD )
#define QTYPE_N_FLUXI               ( QTYPE_LOC_NODE | QTYPE_INT_FLUXI )
#define QTYPE_N_FLUXD               ( QTYPE_LOC_NODE | QTYPE_INT_FLUXD )

#define QTYPE_E_COORD               ( QTYPE_LOC_ELEM | QTYPE_INT_COORD )
#define QTYPE_E_FIELD               ( QTYPE_LOC_ELEM | QTYPE_INT_FIELD )
#define QTYPE_E_FLUXI               ( QTYPE_LOC_ELEM | QTYPE_INT_FLUXI )
#define QTYPE_E_FLUXD               ( QTYPE_LOC_ELEM | QTYPE_INT_FLUXD )

/*
 * Macros for the plain qtype argument
 */
#define QTYPE_LOC(_qtype)           ( (_qtype) & QTYPE_LOC_MASK )
#define QTYPE_INT(_qtype)           ( (_qtype) & QTYPE_INT_MASK )
#define QTYPE_IS_NODE(_qtype)       ( QTYPE_LOC(_qtype) == QTYPE_LOC_NODE  )
#define QTYPE_IS_ELEM(_qtype)       ( QTYPE_LOC(_qtype) == QTYPE_LOC_ELEM  )
#define QTYPE_IS_COORD(_qtype)      ( QTYPE_INT(_qtype) == QTYPE_INT_COORD )
#define QTYPE_IS_FIELD(_qtype)      ( QTYPE_INT(_qtype) == QTYPE_INT_FIELD )
#define QTYPE_IS_FLUXI(_qtype)      ( QTYPE_INT(_qtype) == QTYPE_INT_FLUXI )
#define QTYPE_IS_FLUXD(_qtype)      ( QTYPE_INT(_qtype) == QTYPE_INT_FLUXD )

#define QTYPE_SET_NODE(_qtype)      { (_qtype) &= ~QTYPE_LOC_MASK; (_qtype) |= QTYPE_LOC_NODE;  }
#define QTYPE_SET_ELEM(_qtype)      { (_qtype) &= ~QTYPE_LOC_MASK; (_qtype) |= QTYPE_LOC_ELEM;  }
#define QTYPE_SET_FLUXD(_qtype)     { (_qtype) &= ~QTYPE_INT_MASK; (_qtype) |= QTYPE_INT_FLUXD; }
#define QTYPE_SET_FLUXI(_qtype)     { (_qtype) &= ~QTYPE_INT_MASK; (_qtype) |= QTYPE_INT_FLUXI; }

/*
 * Macros for the qtype of a quantity buffer _bid
 */
#define QBUF_STATE(_quant,_bid)        (_quant)->qtagv[_bid].state
#define QBUF_QTYPE(_quant,_bid)        (_quant)->qtagv[_bid].qtype
#define QBUF_LOC(_quant,_bid)          ( QBUF_QTYPE(_quant,_bid) & QTYPE_LOC_MASK )
#define QBUF_INT(_quant,_bid)          ( QBUF_QTYPE(_quant,_bid) & QTYPE_INT_MASK )

#define QBUF_IS_NODE(_quant,_bid)      ( QBUF_LOC(_quant,_bid) == QTYPE_LOC_NODE )
#define QBUF_IS_ELEM(_quant,_bid)      ( QBUF_LOC(_quant,_bid) == QTYPE_LOC_ELEM )
#define QBUF_IS_COORD(_quant,_bid)     ( QBUF_INT(_quant,_bid) == QTYPE_INT_COORD )
#define QBUF_IS_FIELD(_quant,_bid)     ( QBUF_INT(_quant,_bid) == QTYPE_INT_FIELD )
#define QBUF_IS_FLUXI(_quant,_bid)     ( QBUF_INT(_quant,_bid) == QTYPE_INT_FLUXI )
#define QBUF_IS_FLUXD(_quant,_bid)     ( QBUF_INT(_quant,_bid) == QTYPE_INT_FLUXD )


#define QBUF_SET_NODE(_quant,_bid)     QTYPE_SET_NODE (QBUF_QTYPE(_quant,_bid))
#define QBUF_SET_ELEM(_quant,_bid)     QTYPE_SET_ELEM (QBUF_QTYPE(_quant,_bid))
#define QBUF_SET_FLUXD(_quant,_bid)    QTYPE_SET_FLUXD(QBUF_QTYPE(_quant,_bid))
#define QBUF_SET_FLUXI(_quant,_bid)    QTYPE_SET_FLUXI(QBUF_QTYPE(_quant,_bid))

/* No. of values for a mesh or individual part */
#define QBUF_NVALS(_quant,_bid,_mp)    ( QBUF_IS_NODE(_quant,_bid)   ? (_mp)->nnodes : (_mp)->nelems )
#define QUANT_NVALS(_quant,_mp)        ( MPCCI_QUANT_IS_NODE(_quant) ? (_mp)->nnodes : (_mp)->nelems )


/* Offset in the values buffer for a intividual part */
#define QBUF_POFFS(_quant,_bid,_p)     ( QBUF_IS_NODE(_quant,_bid)   ? (_p)->nnodes_o : (_p)->nelems_o )
#define QUANT_POFFS(_quant,_p)         ( MPCCI_QUANT_IS_NODE(_quant) ? (_p)->nnodes_o : (_p)->nelems_o )


#define INTEG_FLAG_MESHINTEG        0x00000001  /* also update the mesh integral on the fly */
#define INTEG_FLAG_WRITEBACK        0x00000010  /* write back the values in case of node->element integration */


/*
 * Norm values structure used to exchange the norm values of a quantity or the diff of quantities.
 */
struct _QNORM
{
   double   n_min;
   double   n_max;
   double   n_mean;
   unsigned n_flags;
};

#define QNORM_FLAG_NMIN    0x00000001  /* bit set if the n_min value is valid */
#define QNORM_FLAG_NMAX    0x00000002  /* bit set if the n_max value is valid */
#define QNORM_FLAG_NMEAN   0x00000004  /* bit set if the n_mean value is valid */

/*
 * Quantity monitor object which contains mean/integral values of a quantity on a mesh
 * and can be monitored vs. iteration/time
 */
struct _QMON
{
   CHAIN_DEFINE_CLASSBASE(QMON)

   const QUANT *quant;                 /* my quantity */
   const char  *unit;                  /* my unit */
   QTAG         qtag;
   real         vmin[MPCCI_QDIM_MAX];
   real         vmax[MPCCI_QDIM_MAX];
   real         mean[MPCCI_QDIM_MAX];  /* mean value */
};

#define QMON_FLAG_VMIN     0x00000001  /* bit set if the vmin value is valid */
#define QMON_FLAG_VMAX     0x00000002  /* bit set if the vmax value is valid */
#define QMON_FLAG_MEAN     0x00000004  /* bit set if the mean value is valid */


/*
 * The writer is a special kind of operator which saves data in a trace file
 */
struct _WRITER
{
   CHAIN_DEFINE_CLASSBASE(WRITER)

   char *config;  /* strdup'ed string */
   char *suffix;  /* strdup'ed string */
   char *jobdir;  /* strdup'ed string */

   /*
    * handles job->step->code->mesh->part.
    */
   void *h_job;
   void *h_step;
   void *h_code;
   void *h_mesh;
   void *h_part;

   OPX_GETINFO *getinfo; /* standard getinfo function for all shared libs */

   /* (optional) function pointers */
   void *(*open_job  )(const WRITER *writer, const JOB *job);
   void  (*close_job )(void *jhandle);

   void *(*open_step )(void *jhandle, const TAG *tag);
   void  (*close_step)(void *shandle);

   void *(*open_code )(void *shandle, const CODE *code);
   void  (*close_code)(void *chandle);

   void *(*open_mesh )(void *chandle, const MESH *mesh);
   void  (*close_mesh)(void *mhandle);

   void *(*open_part )(void *mhandle, const PART *part);
   void  (*close_part)(void *phandle);

   int   (*save_quant)(void *phandle, const QUANT *quant, const PART *part);
   int   (*save_glob )(void *chandle, const GLOB *glob);
   int   (*save_qmon )(void *chandle, const QMON *qmon);
};

/*
 * The following flags may be used by a specific writer. The usability
 * depends on the individual file formats which are supported and the
 * and writer specific capabilities.
 */

/* file format/handling flags */
#define WRITER_FLAG_BINARY    0x00000001  /* save binary file */
#define WRITER_FLAG_GZIPPED   0x00000002  /* save a .gz file */
#define WRITER_FLAG_SINGLE    0x00000004  /* save a single file with all time steps */
#define WRITER_FLAG_MIDEDGE   0x00000008  /* save high order quadratic elements with midedge nodes as is */
#define WRITER_FLAG_SFLUSH    0x00000010  /* flush all files within WRITER_close_setp() */

/* additional quantities */
#define WRITER_FLAG_ORPHANS   0x00000100  /* save node orphan level OLEVEL_... */
#define WRITER_FLAG_NODEIDS   0x00000200  /* save node ID's */
#define WRITER_FLAG_ELEMIDS   0x00000400  /* save element ID's */
#define WRITER_FLAG_ELSIZES   0x00000800  /* save element sizes */
#define WRITER_FLAG_GLOBALS   0x00001000  /* save received global variables */
#define WRITER_FLAG_SLAVES    0x00002000  /* save slave node marks (0/1) */
#define WRITER_FLAG_DOMAINS   0x00004000  /* save node domains */
#define WRITER_FLAG_PEAKS     0x00008000  /* save quantity peak & mean values */
#define WRITER_FLAG_INORM     0x00010000  /* save quantity iteration norm peak & mean values */
#define WRITER_FLAG_PERIODIC  0x00020000  /* save periodic part and quantity information */

/* how to distribute the values among various files */
#define WRITER_FLAG_GRMASK    0x0f000000
#define WRITER_FLAG_GROUPA    0x00000000  /* all on one file */
#define WRITER_FLAG_GROUPC    0x01000000  /* one file per code */
#define WRITER_FLAG_GROUPM    0x02000000  /* one file per mesh */
#define WRITER_FLAG_GROUPP    0x03000000  /* one file per part */
#define WRITER_FLAG_GROUPQ    0x04000000  /* one file per quantity */

/*
 * Flags set by the writer only
 */
#define WRITER_FLAG_KEEPPARTS 0x00100000  /* writer always need all parts, even empty parts */

/*
 * Writer status flags
 */
#define WRITER_SFLAG_ANY      0x00000000 /* save anything available: no check */
#define WRITER_SFLAG_FULL     0x00000001 /* full save: only if all data are available */
#define WRITER_SFLAG_MESH     0x00000002 /* force a mesh save for error checks */

#ifndef MPCCI_COMPILE_SERVER

   #if defined(_WINNT) || defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
      #define _DO_EXPORT extern __declspec(dllexport)
   #else
      #define _DO_EXPORT extern
   #endif

   _DO_EXPORT const char *OPERAT_notify(const char *name);
   _DO_EXPORT const char *WRITER_notify(const char *name);

   #undef _DO_EXPORT

#endif


#ifdef __cplusplus
}
#endif

#endif
