#ifndef __MPCCI4_TYPES_HEADER_INCLUDED
#define __MPCCI4_TYPES_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_types.h
 *
 * MpCCI header file with basic typedef'ed structures
 *
 * $Id: mpcci_types.h 1900 2015-09-01 15:56:38Z jkleinert $
 *
 ****************************************************************************************/

#include "mpcci_elements.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BEGIN_FOREACH_JOBPART(_part,_job) \
if (_job)\
{\
   MPCCI_SERVER *_cppmacro_internal_server;\
   FOREACH(_cppmacro_internal_server,(_job)->servers)\
      FOREACH(_part,_cppmacro_internal_server->parts)\

#define END_FOREACH_JOBPART()  }


#define BEGIN_FOREACH_JOBQUANT(_quant,_job) \
{\
   MPCCI_PART *_cppmacro_internal_part;\
   BEGIN_FOREACH_JOBPART(_cppmacro_internal_part,_job)\
      FOREACH(_quant,_cppmacro_internal_part->quants)

#define END_FOREACH_JOBQUANT()   END_FOREACH_JOBPART() }


/*
 * all structures defined below
 */
typedef struct _MPCCI_JOB        MPCCI_JOB;
typedef struct _MPCCI_SERVER     MPCCI_SERVER;
typedef struct _MPCCI_MSERVER    MPCCI_MSERVER; /* USED LATER */
typedef struct _MPCCI_QUANT      MPCCI_QUANT;
typedef struct _MPCCI_PART       MPCCI_PART;
typedef struct _MPCCI_GLOB       MPCCI_GLOB;
typedef struct _MPCCI_CINFO      MPCCI_CINFO;   /* code info */
typedef struct _MPCCI_QINFO      MPCCI_QINFO;   /* quant info */
typedef struct _MPCCI_PINFO      MPCCI_PINFO;   /* part info */
typedef struct _MPCCI_TINFO      MPCCI_TINFO;   /* time/transfer info */
typedef struct _MPCCI_MOTION     MPCCI_MOTION;  /* solid mesh motion transformation */
typedef struct _MPCCI_PERIODIC   MPCCI_PERIODIC;
typedef struct _MPCCI_DRIVER     MPCCI_DRIVER;
typedef struct _MPCCI_QTAG       MPCCI_QTAG;
typedef struct _MPCCI_QSTAT      MPCCI_QSTAT;
typedef struct _MPCCI_EFACES     MPCCI_EFACES;
typedef struct _MPCCI_IDMAP      MPCCI_IDMAP;


/*
 * struct SOCK_t is not defined if we do not compile with the socket stuff.
 * We then need a placeholder.
 */
#ifdef SOCK_T_DEFINED
   typedef SOCK_t    MPCCI_CONNECT;
#else
   typedef void      MPCCI_CONNECT;
#endif



struct _MPCCI_GLOB
{
   MPCCI_GLOB  *next;

   /*
    * parameters of the transformation relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   double   scale;
   double   offs;
   double   vfill;   /* default fill value if not assigned */

   unsigned flags;   /* contains the direction send/receive and unit bits */
   int      gid;     /* id of the global */
   int      qid;     /* quantity id of the global */
   char     name[MPCCI_MESHNAME_MAX];  /* name of the global (may be used by the client code) */
   char     gaux[MPCCI_AUXSTRING_MAX]; /* optional auxiliary information */
};


struct _MPCCI_QUANT
{
   MPCCI_QUANT *next;
   char        *name;      /* name of quantity */
   void        *sbuffer;   /* optional quantities local transfer buffer (used by the adapter) */

   /*
    * The usr_xxx members can be used by the adapter to keep some stuff locally bound to
    * a MPCCI_QUANT, data which need to be exchanged between the driver calls of the
    * adapter. This help to avoids global/static variables inside the adapter code and
    * keep the logic and data structures clear.
    * usr_xxx is NOT used/modified within any of the mpcci_xxx() functions
    */
   void        *qusr_ptr[MPCCI_USERPOINTERS_MAX];
   unsigned     qusr_flags;
   int          qusr_state;

   double   maxitol; /* max. tolerance in case of iterative coupling */

   /*
    * parameters of the transformation relevant if ((flags&MPCCI_QFLAG_UNIT_MASK) != 0).
    * transformation to corresponding SI unit SI = scale*Q + offs
    */
   double   scale;   /* unit scale of the quantity */
   double   offs;    /* offset for unit translation */
   double   vfill;   /* default fill value if not assigned */
   double   relax;   /* < 0 = unused, else  0.00001 < relax < 2.. */
   double   ramp0;   /* < 0 = unused, else small start value 0.00001 ramp0 < 2 */
   double   rampd;   /* ramp delta value < 1 */
   unsigned flags;   /* contains the direction send/receive bits */
   int      qid;     /* quantity id of the quantity */
   int      maxiter; /* max. no. of iterations in case of iterative coupling */
   int      bufnval; /* size of the local transfer buffer (used by the manager) */
   int      state;   /* state of transfer CQ_TSTATE_xxx (used by the manager) */
   int      smethod; /* storage method and storage index are only used by the adapters */
   int      sindex;
};

#define MPCCI_QUANT_SMETHOD(_q)        ( (_q)->smethod     )
#define MPCCI_QUANT_SINDEX(_q)         ( (_q)->sindex      )
#define MPCCI_QUANT_SBUFFER(_q)        ( (_q)->sbuffer     )
#define MPCCI_QUANT_SBFSIZE(_q)        ( (_q)->bufnval     )

#define MPCCI_QUANT_USRFLAGS(_q)       ( (_q)->qusr_flags  )
#define MPCCI_QUANT_USRSTATE(_q)       ( (_q)->qusr_state  )
#define MPCCI_QUANT_USRPTR (_q,_ix)    ( (_q)->qusr_ptr[_ix] )
#define MPCCI_QUANT_USRPTR0(_q)        ( (_q)->qusr_ptr[0] )
#define MPCCI_QUANT_USRPTR1(_q)        ( (_q)->qusr_ptr[1] )
#define MPCCI_QUANT_USRPTR2(_q)        ( (_q)->qusr_ptr[2] )
#define MPCCI_QUANT_USRPTR3(_q)        ( (_q)->qusr_ptr[3] )


struct _MPCCI_PART
{
   MPCCI_PART  *next;      /* next part in a chain of parts */
   char        *name;      /* name of the coupled partition */
   char        *paux;      /* if not NULL, then pointer to any auxiliary information string */
   MPCCI_QUANT *quants;    /* chain of quantities */

   /*
    * The usr_xxx members can be used by the adapter to keep some stuff locally bound to
    * a MPCCI_PART, data which need to be exchanged between the driver calls of the
    * adapter. This help to avoids global/static variables inside the adapter code and
    * keep the logic and data structures clear.
    *    pointer to a fluent thread or Star-CD node index buffer
    *    an array of the (float/double) node coordinates of a part
    *    element/cell information
    *    .... and more
    * usr_xxx is NOT used/modified within any of the mpcci_xxx() functions
    */
   void        *pusr_ptr[MPCCI_USERPOINTERS_MAX];
   unsigned     pusr_flags;
   int          pusr_state;

   int          mid;       /* mesh id */
   int          pid;       /* part id */
   int          mdim;      /* 0=netpoint, 1=line, 2=face, 3=volume */
   unsigned     csys;      /* coordinates system */
   int          nnodes;    /* no. of nodes/points in part */
   int          nelems;    /* no. of faces/elements/cells in part*/
   int          nodmax;    /* max. no. of nodes per element in this part */
   unsigned     flags;     /* temporary helper, for the future */
   int          state;     /* temporary helper, e.g. used with the remesh state */
};

#define MPCCI_PART_NAME(_p)            ( (_p)->name      )
#define MPCCI_PART_PAUX(_p)            ( (_p)->paux      )
#define MPCCI_PART_MDIM(_p)            ( (_p)->mdim      )
#define MPCCI_PART_CSYS(_p)            ( (_p)->csys      )
#define MPCCI_PART_NNODES(_p)          ( (_p)->nnodes    )
#define MPCCI_PART_NELEMS(_p)          ( (_p)->nelems    )
#define MPCCI_PART_NODMAX(_p)          ( (_p)->nodmax    )
#define MPCCI_PART_MESHID(_p)          ( (_p)->mid       )
#define MPCCI_PART_PARTID(_p)          ( (_p)->pid       )
#define MPCCI_PART_QUANTS(_p)          ( (_p)->quants    )

#define MPCCI_PART_USRFLAGS(_p)        ( (_p)->pusr_flags  )
#define MPCCI_PART_USRSTATE(_p)        ( (_p)->pusr_state  )
#define MPCCI_PART_USRPTR(_p,_ix)      ( (_p)->pusr_ptr[_ix] )
#define MPCCI_PART_USRPTR0(_p)         ( (_p)->pusr_ptr[0] )
#define MPCCI_PART_USRPTR1(_p)         ( (_p)->pusr_ptr[1] )
#define MPCCI_PART_USRPTR2(_p)         ( (_p)->pusr_ptr[2] )
#define MPCCI_PART_USRPTR3(_p)         ( (_p)->pusr_ptr[3] )


#define MPCCI_PART_IS_NONE(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_NONE  )
#define MPCCI_PART_IS_PNTS(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_PNTS  )
#define MPCCI_PART_IS_LINE(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_LINE  )
#define MPCCI_PART_IS_FACE(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_FACE  )
#define MPCCI_PART_IS_VOLU(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_VOLU  )
#define MPCCI_PART_IS_IPNT(_p)         ( MPCCI_PART_MDIM(_p) == MPCCI_MDIM_IPNT  )


#define FOREACH_QSEND(_quant,_head)    FOREACH_COND(_quant,_head,MPCCI_QUANT_IS_SEND(_quant))
#define FOREACH_QRECV(_quant,_head)    FOREACH_COND(_quant,_head,MPCCI_QUANT_IS_RECV(_quant))

#define FOREACH_QID_SEND(_quant,_head,_qid)    FOREACH_COND(_quant,_head,MPCCI_QUANT_QID(quant)==_qid && MPCCI_QUANT_IS_SEND(_quant))
#define FOREACH_QID_RECV(_quant,_head,_qid)    FOREACH_COND(_quant,_head,MPCCI_QUANT_QID(quant)==_qid && MPCCI_QUANT_IS_RECV(_quant))

#ifdef MPCCI_COMPILE_SLIMCLIENT

/*
 * single server resp. broker connection definition
 */
struct _MPCCI_SERVER
{
   MPCCI_CONNECT *connect;    /* connection to this server/broker */
   MPCCI_PART    *parts;      /* parts for this server */
   MPCCI_GLOB    *globs;      /* globals for this server */
   unsigned       flags;
   unsigned       realsize;   /* server's const realsize received during smpcci_init() */
   int            nodelay;    /* if true do not delay flush on socket connection */
};

/*
 * multiple servers + broker connection
 */
struct _MPCCI_JOB
{
   MPCCI_SERVER  *server; /* chain of servers, one server per mesh id */
   unsigned       flags;  /* just a bitfield helper for keeping some strange things */
};

#else

/*
 * single server resp. broker connection definition
 */
struct _MPCCI_SERVER
{
   MPCCI_SERVER  *next;       /* pointer to next assigned server */
   MPCCI_CONNECT *connect;    /* connection to this server/broker */
   MPCCI_PART    *parts;      /* parts for this server */
   MPCCI_GLOB    *globs;      /* globals for this server */
   unsigned       flags;
   unsigned       realsize;   /* server's const realsize received during smpcci_init() */
   int            nodelay;    /* if true do not delay flush on socket connection */
};

#if 0 /* USED LATER */
/*
 * single server resp. broker connection definition
 */
struct _MPCCI_MSERVER
{
   MPCCI_MSERVER *next;       /* pointer to next assigned server */
   MPCCI_SERVER  *server;
   MPCCI_CONNECT *connect;    /* connection to this server/broker */
   MPCCI_PART    *parts;      /* parts for this server */
   MPCCI_GLOB    *globs;      /* globals for this server */
   int            state;
   unsigned       flags;
   unsigned       realsize;   /* server's const realsize received during smpcci_init() */
};
#endif

/*
 * multiple servers + broker connection
 */
struct _MPCCI_JOB
{
   MPCCI_SERVER  *servers; /* chain of servers, one server per mesh id */
   MPCCI_SERVER  *broker;  /* connection to the broker daemon */
   unsigned       flags;   /* just a bitfield helper for keeping some strange things */
};

#endif

/*
 * Structure used be smpcci_init() to keep code specific information
 */
struct _MPCCI_CINFO
{
   const char *jobid;
   const char *codename;
   const char *codeid;
   const char *clientid;
   double      time;       /* initial physical time of the code */
   unsigned    flags;      /* feature flags of the code */
   int         iter;       /* initial iteration counter of the code */
   int         nclients;   /* exact no. of expected server clients: <=1 serial,  >1 parallel code */
   int         nprocs;     /* no. of parallel processes, no necessarily the no. of clients */
};


/*
 * Structure used to transport code specific quantity information between
 * client and server. On the server side this is CODE_QINFO
 */
struct _MPCCI_QINFO
{
   double   scale;      /* unit scale of the quantity */
   double   offs;       /* offset for unit translation */
   double   vfill;      /* default fill value if not assigned */
   double   relax;      /* < 0 = unused, else  0.00001 < relax < 0.999999.. */
   double   ramp0;      /* < 0 = unused, else small start value 0.00001 ramp0 < 2 */
   double   rampd;      /* ramp delta value < 1 */
   double   maxitol;    /* max. tolerance in case of iterative coupling */
   unsigned flags;      /* contains the direction send/receive bits */
   int      qid;        /* quantity id of the quantity */
   int      mid;        /* coupled mesh for which these settings are defined */
   int      method_s;   /* send method and storage index */
   int      index_s;
   int      method_r;   /* receive method and storage index */
   int      index_r;
   int      maxiter;    /* max. no. of iterations in case of iterative coupling */
   char     name[MPCCI_QUANTNAME_MAX]; /* name of the quantity (may be used by the client code) */
};


/*
 * Structure used to transport code specific coupled part information between
 * client and server. data member ordering according to best alignment!
 */
struct _MPCCI_PINFO
{
   int    mid;                       /* id of the mesh (relevant for the server + the client code)*/
   int    pid;                       /* id of the part (may be used by the client code)*/
   int    mdim;                      /* assumed(!) coupling dimension the part MPCCI_MDIM_* */
   int    qidv_s[MPCCI_MESHQ_MAX];   /* list of quantity id's the client should send */
   int    qidv_r[MPCCI_MESHQ_MAX];   /* list of quantity id's the client should receive */
   char   name[MPCCI_MESHNAME_MAX];  /* name of the part */
   char   paux[MPCCI_AUXSTRING_MAX]; /* optional auxiliary information for his part */
};


/*
 * Structure used to keep time/iteration information for ampcci_transfer()
 */
struct _MPCCI_TINFO
{
   unsigned this_size;     /* (unsigned)sizeof(MPCCI_TINFO) */
   unsigned this_version;  /* version 3.0.6=306, 3.1.0=310, 4.0.0=400 */

   /*
    * MpCCI state:
    *  -1 : disconnected from the server
    *   0 : nothing done before
    *   1 : initially connected to the server
    *  >1 : init done and at least one data transfer
   */
   int      mpcci_state;
   int      mpcci_used;    /* bool: true, if we use MpCCI */

   unsigned tact;          /* bitpattern: MPCCI_TACT_* */
   unsigned tact_init;     /* bitpattern: same for the initial exchange, see ampcci_tinfo() */
   unsigned tact_last;

   int      mode_t;        /* Coupling scheme: explicit / implicit */
   int      mode_a;        /* Algorithm type:   GS / Jacobi */
   unsigned mode_c;        /* check mode types */
   int      mode_s;        /* various levels of send modes */
   int      mode_r;        /* various levels of receive modes */
   int      mode_q;        /* quantity selection mode */

   unsigned flag_r;        /* bitpattern: remeshing flag */

   /*
    * in case of a time >= 0.0 we are transient and might use
    * the time and iteration for testing
    */
   double   dt;            /* current physical time step size: updated by the code before ampcci_transfer() */
   double   dt_last;       /* physical time step size during the last call */
   double   rdt_min;       /* relative time distance to interpolated values */
   double   rdt_max;       /* relative time distance to interpolated values */

   double   time;          /* current physical time: updated by the code before ampcci_transfer() */

   double   time_cbeg;     /* time when the coupling starts in terms of the code times */
   double   time_cend;     /* time when the coupling ends in terms of the code times */
   double   time_init;     /* code time when the adapter was initialized */
   double   time_last;     /* code time of the last exchange */

   /*
    * implicit coupling parameters
    */
   double   implt_step;     /* time step size between two exchanges: use in implicit coupling */
   double   implt_next;     /* next implicit coupling time information which can be updated by the adapter */
   int      impli_step;     /* no. of iteration between two exchanges */
   int      impli_max;      /* maximum no. of iterations */

   /*
    * in case of a time < 0.0 we are steady/iterative
    */
   int      iter;          /* current iteration counter: updated by the code before ampcci_transfer() */

   int      iter_cbeg;     /* first iteration where we exchange */
   int      iter_cend;     /* last iteration where we exchange */
   int      iter_pend;     /* no. iteration after last iteration where we exchange */
   int      iter_init;     /* code iteration when the adapter was initialized */
   int      iter_last;     /* iteration of the last exchange */

   /*
    * subcycle
    */
   int      sub_nstep;    /* current coupling step counter */
   int      sub_step;     /* no. of steps between two exchanges */
   int      sub_next;     /* counter of next coupling step */

   /*
    * information updated per call of ampcci_transfer()
    */
   int      count_s;       /* no. of ampcci_transfer() send calls (independent of quantities and meshes) */
   int      count_r;       /* no. of ampcci_transfer() receive calls (independent of quantities and meshes) */

   int      nquant_s;      /* no. of send quantities sum( mesh x nquant ) within the last ampcci_transfer() */
   int      nquant_r;      /* no. of received quantities sum( mesh x nquant ) within the last ampcci_transfer() */
   int      nquant_req;    /* no. of requested quantities sum( mesh x nquant ) within the last ampcci_transfer() */

   int      nglob_s;       /* no. of send globals within the last ampcci_transfer() */
   int      nglob_r;       /* no. of received globals within the last ampcci_transfer() */
   int      nglob_req;     /* no. of requested globals within the last ampcci_transfer() */

   /*
    * optional convergence status
    */
   int      conv_code;     /* MPCCI_CONV_STATE_*  input + returned */
   int      conv_job;      /* returned after the transfer */
   int      conv_last;     /* store the last convergence state */
   int      conv_counter;  /* store the number of converged transfers in a row */
};


/* Transfer action bits to execute all steps of a send and receive */
#define MPCCI_TACT_DIR     0x000f
#define MPCCI_TACT_SEND    0x0001   /* send quantities to the server */
#define MPCCI_TACT_RECV    0x0002   /* receive quantities from the server */
#define MPCCI_TACT_XCHG    (MPCCI_TACT_SEND|MPCCI_TACT_RECV)

#define MPCCI_TACT_XTRA    0x00f0
#define MPCCI_TACT_GETQ    0x0010   /* anyway get quantities from the code - call the driver */
#define MPCCI_TACT_PUTQ    0x0020   /* anyway put quantities to the code - call the driver */

/* Coupling modes */
#define MPCCI_TMODE_MASK      0x000f
#define MPCCI_TMODE_IMPLICIT  0x0001
#define MPCCI_TMODE_EXPLICIT  0x0002

/* Solution type */
#define MPCCI_TSOLU_MASK      0x00f0
#define MPCCI_TSOLU_TRANSIENT 0x0010
#define MPCCI_TSOLU_STEADY    0x0020

/* Coupling algorithm modes */
enum
{
   MPCCI_TMODE_JACOBI          = 0,
   MPCCI_TMODE_GAUSSSEIDEL     = 1
};

/* check modes */
#define MPCCI_TMODE_CHECK_NONE        0x0000   /* no check */
#define MPCCI_TMODE_CHECK_SUBCYCLE    0x0001   /* do subcycle */
#define MPCCI_TMODE_CHECK_DURATION    0x0002   /* do control duration */
#define MPCCI_TMODE_CHECK_FORWARD    (MPCCI_TMODE_CHECK_SUBCYCLE|MPCCI_TMODE_CHECK_DURATION)

/* send modes */
enum
{
   MPCCI_TMODE_SEND_ALWAYS       = 0,  /* always send */
   MPCCI_TMODE_SEND_ONEWAIT      = 1,  /* send if at least one partner is waiting */
   MPCCI_TMODE_SEND_ALLWAIT      = 2   /* send if all partners are waiting */
};

/* receive modes */
enum
{
   MPCCI_TMODE_RECV_AVAIL        = 0,  /* do not wait, get what's available */
   MPCCI_TMODE_RECV_COMPLETE     = 1,  /* receive only if all quantities are available, otherwise skip */
   MPCCI_TMODE_RECV_ANY          = 2,  /* receive all if at least one quantity is available, otherwise skip */
   MPCCI_TMODE_RECV_ALL          = 3   /* always wait for all requested quantities */
};

/*
 * Quantity transfer modes
 */
enum
{
   MPCCI_TMODE_QUANT_ALL         = 0,  /* exchange all quantites */
   MPCCI_TMODE_QUANT_COORD       = 1,  /* exchange only coordinate type quantities */
   MPCCI_TMODE_QUANT_OTHER       = 2   /* exchange all non coordinate type quantities */
};


/*
 * Settings based on the adapter code
 */
#define MPCCI_TINFO_IS_TIME(_tinfo)    ( (_tinfo)->time >= 0                        )
#define MPCCI_TINFO_IS_ITER(_tinfo)    ( (_tinfo)->iter >= 0                        )
#define MPCCI_TINFO_IS_QTID(_tinfo)    ( (_tinfo)->time <  0 && (_tinfo)->iter <  0 )
#define MPCCI_TINFO_IS_IMPL(_tinfo)    ( (_tinfo)->time >= 0 && (_tinfo)->iter >= 0 )
#define MPCCI_TINFO_IS_EXPL(_tinfo)    ( (_tinfo)->time >= 0 && (_tinfo)->iter <  0 )

/*
 * Settings based on the ampcci_tinfo_init()
 */
#define MPCCI_TINFO_WANT_IMPL(_tinfo)  ( (_tinfo)->mode_t & MPCCI_TMODE_IMPLICIT )
#define MPCCI_TINFO_WANT_EXPL(_tinfo)  ( (_tinfo)->mode_t & MPCCI_TMODE_EXPLICIT )

#define MPCCI_TSOLU_WANT_TRANS(_tinfo) ( (_tinfo)->mode_t & MPCCI_TSOLU_TRANSIENT )
#define MPCCI_TSOLU_WANT_STEADY(_tinfo)( (_tinfo)->mode_t & MPCCI_TSOLU_STEADY    )

#define MPCCI_TINFO_WANT_JACOBI(_tinfo)      ( (_tinfo)->mode_a == MPCCI_TMODE_JACOBI      )
#define MPCCI_TINFO_WANT_GAUSSSEIDEL(_tinfo) ( (_tinfo)->mode_a == MPCCI_TMODE_GAUSSSEIDEL )

#define MPCCI_TINFO_WANT_PREDICTOR(_tinfo)  (MPCCI_TINFO_WANT_GAUSSSEIDEL(_tinfo) \
    && ( ((_tinfo)->tact_init & MPCCI_TACT_DIR) ==  MPCCI_TACT_XCHG || !(_tinfo)->tact_init ) ) \


/*
 * Structure used to transport mesh specific moving reference frame information
 */
struct _MPCCI_MOTION
{
   double   omeg;    /* omega */
   double   axis[3]; /* rotational axis */
   double   orig[3]; /* axis origin */
   double   velo[3]; /* velocity */
   unsigned type;    /* type/order of rotation/translation: MPCCI_MOTION_XXXX.... */
};
/*
 * Structure used to transport part specific periodicity information
 */
struct _MPCCI_PERIODIC
{
   double   	axis[3]; 	/* rotational axis */
   double   	origin[3]; 	/* axis origin */
   int   		n;    		/* number of periodicities*/
};

/*
 * List of pure code dependent methods to define the grid and
 * to exchange quantities between the code and the coupling manager.
 */
struct _MPCCI_DRIVER
{
   /***************************************************************************************
    * REQUIRED: information for compatibility checks
    ***************************************************************************************/
   unsigned this_size;     /* (unsigned)sizeof(MPCCI_DRIVER) */
   unsigned this_version;  /* version 3.0.6=306, 3.1.0=310, 4.0.0=400 */

   /***************************************************************************************
    * REQUIRED: bitmask for or'ing MPCCCI_TINFO.tact
    ***************************************************************************************/
   unsigned tact_required;   /* bitmask: MPCCI_TACT_* */

   /***************************************************************************************
    * REQUIRED: define symbolic names for the types of parts in code terminology
    *           example:
    *             [0] = "Point(s) with coordinates, single point or point cloud",
    *             [1] = "Line",
    *             [2] = "Surface",
    *             [3] = "Volume domain",
    *             [4] = "Single network integration point",
    *             [5] = "Global variable"
    *             [6] = NULL, not yet used
    *             [7] = NULL, not yet used
    ***************************************************************************************/
   const char *part_description[8];

   /***************************************************************************************
    * REQUIRED: (unsigned)sizeof(real) for all received quantities
    ***************************************************************************************/
   unsigned realsize;

   /***************************************************************************************
    * OPTIONAL: methods called before/after some actions
    ***************************************************************************************/
   void (*afterCloseSetup) (void);  /* method called after def_close */
   void (*beforeGetAndSend)(void);  /* method called before send operation */
   void (*afterGetAndSend) (void);  /* method called after send operation */
   void (*beforeRecvAndPut)(void);  /* method called before receive operation */
   void (*afterRecvAndPut) (void);  /* method called after receive operation */

   void (*partSelect)         (const MPCCI_PART *part);
   int  (*partUpdate)         (      MPCCI_PART *part, MPCCI_QUANT  *quant);
   void (*appendParts)        (const MPCCI_JOB  *job , MPCCI_SERVER *server);
   int  (*getPartRemeshState) (      MPCCI_PART *part, unsigned flags);

   /***************************************************************************************
    * REQUIRED: methods to get information about a part
    ***************************************************************************************/
   /* variante #1: fast & dirty: returns the MpCCI error code */
   int (*definePart)(MPCCI_SERVER *server, MPCCI_PART *part);

   /* variante 2: may be slow but better interface */
   void (*partInfo) (      int  *rDim,
                      unsigned  *csys,
                           int  *nNodes,
                           int  *nElems,
                           int  *maxNodesPerElem,
                           int  *regId,
                     const char *regName,
                           int   len);
   int  (*getNodes) (const MPCCI_PART *part, unsigned *csys, void *coords, int *nodeIds);
   int  (*getElems) (const MPCCI_PART *part, unsigned *elemType1, unsigned *elemTypes,
                     int *elemNodes, int *elemIds);

   /* REQUIRED for moving nodes instead of full remeshing: returns the MpCCI error code */
   int  (*moveNodes)       (MPCCI_SERVER *server, MPCCI_PART *part);

   /* OPTIONAL(depends): methods used to get quantities from the application */
   int  (*getPointValues   )(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getLineNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getLineElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getFaceNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getFaceElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getVoluNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getVoluElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   int  (*getGlobValues    )(                        const MPCCI_GLOB  *glob , void *values);

   /* OPTIONAL(depends): methods used to store quantities into the application */
   void (*putPointValues   )(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putLineNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putLineElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putFaceNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putFaceElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putVoluNodeValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putVoluElemValues)(const MPCCI_PART *part, const MPCCI_QUANT *quant, void *values);
   void (*putGlobValues    )(                        const MPCCI_GLOB  *glob , void *values);
};


struct _MPCCI_QTAG
{
   double   time;    /* time >= 0.0 in case of a transient run */
   int      iter;    /* iteration counter, either global for steady or inner iteration count */
   int      qtid;    /* increasing only tag id (=counter) */
   int      mid;     /* mesh id */
   int      qid;     /* quantity id */
};

struct _MPCCI_QSTAT
{
   MPCCI_QTAG *tagv;    /* NOT YET USED */
   int         ntags;   /* no. of tags from the tag query */
   int         nvalid;  /* no. of valid tags resp. quantity buffers */
   int         nwait;   /* no. of tags resp. quantity buffers we need to wait for */
};

/*
 * utility structure used by ELEM_faces() for splitting an element info faces.
 */
struct _MPCCI_EFACES
{
   int     *fnodes[MPCCI_PCELLFACES_MAX]; /* face nodes[][MPCCI_ETYP_NNE(ftypes[i])] */
   unsigned ftypes[MPCCI_PCELLFACES_MAX]; /* face types */
   int      nfaces;                       /* no. of faces for this element */
};

/*
 * Reverse mapping of a global node or element integer id into a
 * local MpCCI part/id pair.
 *
 * This allows fast access to the MpCCI local buffer of the received quantities
 * in case a code internally uses only global ids.
 */
struct _MPCCI_IDMAP
{
   MPCCI_PART **part;   /* part to which the global id belongs */
   int         *index;  /* local node/element id within this part */
   int          size;   /* allocated size of the mapping arrays */
};




/*
 * Data types used for the external morpher communication
 */
/*
 * type argument to MORPHER_Send_start().
 * must be a negative number to avoid conflics with vertex indices >= 0
 */
#define MORPHER_VTYPE_DISP   -1 /* send displacements only */
#define MORPHER_VTYPE_COOR   -2 /* send displaced absolute coordinates */
#define MORPHER_VTYPE_GRID   -3 /* send local grid definition coordinates */
#define MORPHER_VTYPE_EXEC   -4 /* terminate transfer and execute */
#define MORPHER_VTYPE_EXIT   -5 /* disconnect from morpher */

#define MORPHER_FLAG_SKIPGREAT   0x00000001 /* special flag for StarCD */


typedef struct _MORPHER    MORPHER;

struct _MORPHER
{
   MPCCI_CONNECT *connect;

   /*
    * items coords, nverts and ioffs are only used by the high level functions below
    *
    * coords:
    *    pointer to a buffer of [3*nverts] reals containing the vertex coordinates
    *    in order x1,y1,z1, x2,y2,z2, x3,...
    *
    * nverts:
    *    no. of vertices in the coords buffer
    *
    * ioffs:
    *    Internal vertex indices are C type indices [0..n-1]
    *    External vertex indices are ? type indices [0+ioffs..n+ioffs]
    *    ioffs is e.g. 1 with FORTRAN vertices coord(3,NVERTS)
    */
   void    *vcoords;
   unsigned realsize; /* client   sizeof(float) || sizeof(double) */
   unsigned morpsize; /* morphers sizeof(float) || sizeof(double) */
   int      nverts;
   int      ioffs;

   /*
    * Internal current state of transport:
    *   -1: stop state, client is waiting for results
    *    0: start-magic-sent but no vertices yet
    *   >0: no. of vertices sent to server
    */
   int      nsent;

   /*
    * Total no. of parallel client processes that are connected with the morpher.
    * If there is more than a single process connected with the morpher
    * we need to fulfill the communication protocol to be in sync with all
    * other parallel processes.
    * In serial mode we may reduce the communication overhead.
    */
   int      nclients;


   /*
    * In cases where one of the partner has less precision than the other
    * the precision is reduced to the lower precision
    *   e.g. double->float
    * before sending vertices or displacements
    */
   int      reduce;
};


#ifdef __cplusplus
}
#endif

#endif
