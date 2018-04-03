#ifndef __MPCCI4_DEFS_HEADER_INCLUDED
#define __MPCCI4_DEFS_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_defs.h
 *
 * MpCCI header file with basic flags and bitpattern definitions
 *
 * $Id: mpcci_defs.h 1830 2015-08-18 14:39:37Z scaicae $
 *
 ****************************************************************************************/

#include "mpcci_limits.h"

/* defines the version numbers of the server */
#define MPCCI_SERVER_VERSION_MAJOR    4
#define MPCCI_SERVER_VERSION_MINOR    4
#define MPCCI_SERVER_VERSION_REVIS    2

/* defines the version numbers of the protocol used */
#define MPCCI_PROTOCOL_VERSION_MAJOR  4
#define MPCCI_PROTOCOL_VERSION_MINOR  4

#ifdef __cplusplus
   #define MPCCI_SCAST_INTO(_type,_expr)     static_cast<_type>(_expr)
   #define MPCCI_RCAST_INTO(_type,_expr)     reinterpret_cast<_type>(_expr)
   #define MPCCI_CCAST_INTO(_type,_expr)     const_cast<_type>(_expr)
#else
   #define MPCCI_SCAST_INTO(_type,_expr)     (_type)(_expr)
   #define MPCCI_RCAST_INTO(_type,_expr)     (_type)(_expr)
   #define MPCCI_CCAST_INTO(_type,_expr)     (_type)(_expr)
#endif

#define MPCCI_ICAST(_expr)                MPCCI_SCAST_INTO(int     ,_expr)
#define MPCCI_UCAST(_expr)                MPCCI_SCAST_INTO(unsigned,_expr)
#define MPCCI_CPCAST(_expr)               MPCCI_SCAST_INTO(char *  ,_expr)
#define MPCCI_VPCAST(_expr)               MPCCI_SCAST_INTO(void *  ,_expr)


#define FOREACH(_obj,_head)               for(_obj=(_head); _obj; _obj=_obj->next)
#define FOREACH_NEXT(_obj,_head,_next)    for(_obj=(_head); _obj; _obj=_next)
#define FOREACH_COND(_obj,_head,_cond)    FOREACH(_obj,_head) if (_cond)


/*
 * client-server communication command tokens
 */
#define MPCCI_CMD_PUTQ     "PUTQ" /* put mesh based quantity values */
#define MPCCI_CMD_GETQ     "GETQ" /* get mesh based quantity values */
#define MPCCI_CMD_PUTG     "PUTG" /* put global values */
#define MPCCI_CMD_GETG     "GETG" /* get global values */
#define MPCCI_CMD_QTAG     "QTAG" /* query quantity tags status */

#define MPCCI_CMD_DEFP     "DEFP" /* (re)define a partition */
#define MPCCI_CMD_DEFQ     "DEFQ" /* (re)define a quantity on a mesh */
#define MPCCI_CMD_DEFG     "DEFG" /* (re)define a meshless global quantity */

#define MPCCI_CMD_PMOT     "PMOT" /* define a part motion */
#define MPCCI_CMD_PPER     "PPER" /* define a part periodicity */
#define MPCCI_CMD_PSHF     "PSHF" /* define a surface part (baffle/shell) shift */
#define MPCCI_CMD_PNOD     "PNOD" /* (re)define the partition nodes (id + coords) */
#define MPCCI_CMD_PELS     "PELS" /* (re)define the partition elements (id + type + nodeids) */

#define MPCCI_CMD_SYNC     "SYNC" /* code/job general sync */
#define MPCCI_CMD_SYNT     "SYNT" /* sync time/iteration/convergence */
#define MPCCI_CMD_SYNM     "SYNM" /* sync mesh */

#define MPCCI_CMD_QQIX     "QQIX" /* query quantity by index [0..n-1] of the defined quantities */
#define MPCCI_CMD_QQID     "QQID" /* query quantity by id */
#define MPCCI_CMD_QQNM     "QQNM" /* query quantity by name */
#define MPCCI_CMD_QPIX     "QPIX" /* query part info by index */
#define MPCCI_CMD_QGIX     "QGIX" /* query glob info by index */

#define MPCCI_CMD_GETV     "GETV" /* get a code parameter or environment variable */
#define MPCCI_CMD_GETC     "GETC" /* get partner codes */

#define MPCCI_CMD_DELM     "DELM" /* delete a mesh */
#define MPCCI_CMD_DELP     "DELP" /* delete a part */
#define MPCCI_CMD_DELQ     "DELQ" /* delete a quantity on a mesh (all parts) */
#define MPCCI_CMD_DELG     "DELG" /* delete a global quantity */

#define MPCCI_CMD_SAVJ     "SAVJ" /* save trace file if all job clients called it */
#define MPCCI_CMD_SAVC     "SAVC" /* save trace file if all code clients called it */
#define MPCCI_CMD_NICE     "NICE" /* set the nice value of a code */

#define MPCCI_CMD_QUIT     "QUIT" /* a client quit */
#define MPCCI_CMD_STOP     "STOP" /* a client stop or server stop response */
#define MPCCI_CMD_FERR     "FERR" /* a fatal error message and exit */


/*
 * ETX/^C is the OOB signalling a server shutdown
 */
#define MPCCI_OOB_SIGC     0x03
#define MPCCI_OOB_SIGS     "\x03"


/*
 * Message tokens exchanged during the connection setup between the processes.
 *
 * To avoid any confusion with the connections due to improper definitions of
 * the port@host environment variables MPCCI_SERVER and MPCCI_BROKER each process
 * identifies itsself as CLIENT/SERVER/MONITOR/BROKER, sends this to the
 * expected partner on the remote side, and checks the received answer.
 */
#define MPCCI_HELLO_MAGIC           "MpCCI"
#define MPCCI_HELLO_CLIENT_FMT      "CLIENT:%d.%d"     /* MPCCI_PROTOCOL_VERSION */
#define MPCCI_HELLO_MONITOR_FMT     "MONITOR:%d.%d"    /* MPCCI_PROTOCOL_VERSION */
#define MPCCI_HELLO_SERVER_FMT      "SERVER:%u:%d"     /* sizeof(real): printout level */
#define MPCCI_HELLO_BROKER_FMT      "BROKER:%u:%d"     /* sizeof(real): printout level */


/*
 * Public environment variables used
 */
#define MPCCI_ENV_JOBID             "MPCCI_JOBID"
#define MPCCI_ENV_JOBNAME           "MPCCI_JOBNAME"
#define MPCCI_ENV_CODEID            "MPCCI_CODEID"
#define MPCCI_ENV_DEBUG             "MPCCI_DEBUG"
#define MPCCI_ENV_TMPDIR            "MPCCI_TMPDIR"
#define MPCCI_ENV_RSHTYPE           "MPCCI_RSHTYPE"
#define MPCCI_ENV_BROKER            "MPCCI_BROKER"
#define MPCCI_ENV_SERVER            "MPCCI_SERVER"
#define MPCCI_ENV_MORPHER           "MPCCI_MORPHER"
#define MPCCI_ENV_MONITOR           "MPCCI_MONITOR"
#define MPCCI_ENV_CCVXVIEW          "MPCCI_CCVXVIEW"
#define MPCCI_ENV_NETDEVICE         "MPCCI_NETDEVICE"
#define MPCCI_ENV_TINFO             "MPCCI_TINFO"
#define MPCCI_ENV_LICFILE           "MPCCI_LICENSE_FILE"
#define MPCCI_ENV_TIMEOUT           "MPCCI_TIMEOUT"

/*
 * Default TCP portnumbers used by MpCCI
 */
#define MPCCI_DEFPORT_LICENSE       47000 /* used UDP, but who knows */
#define MPCCI_DEFPORT_BROKER        47001
#define MPCCI_DEFPORT_MONITOR       47002
#define MPCCI_DEFPORT_CCVXVIEW      47003
#define MPCCI_DEFPORT_SERVER        47010 /* in fact the first server */
#define MPCCI_DEFPORT_MORPHER       47020


/*
 * Bitpattern for the flags used with [s]mpcci_sync()
 */
#define MPCCI_SYNC_FLAG_QUANT       0x00000001  /* sync the quantities of a code */
#define MPCCI_SYNC_FLAG_GLOB        0x00000002  /* sync the globals of a code */
#define MPCCI_SYNC_FLAG_CODE        0x00000004  /* sync the code */
#define MPCCI_SYNC_FLAG_MESH        0x00000008  /* sync the meshes of the code */
#define MPCCI_SYNC_FLAG_JOB         0x00000010  /* sync the job (all codes of the job) */
#define MPCCI_SYNC_FLAG_ALL         (\
                                       MPCCI_SYNC_FLAG_CODE\
                                       |MPCCI_SYNC_FLAG_MESH\
                                       |MPCCI_SYNC_FLAG_QUANT\
                                       |MPCCI_SYNC_FLAG_GLOB\
                                    )

/*
 * Defaults for the convergence check parameter
 */
enum
{
   MPCCI_CONV_STATE_NOSEND    = -2, /* server internally only: skip sending the conv state */
   MPCCI_CONV_STATE_ERROR     = -1,
   MPCCI_CONV_STATE_INVALID   =  0,
   MPCCI_CONV_STATE_DIVERGED  =  1,
   MPCCI_CONV_STATE_STOP      =  2,
   MPCCI_CONV_STATE_CONVERGED =  3,
   MPCCI_CONV_STATE_CONTINUE  =  4,

   MPCCI_CONV_STATE_MIN  = 0,
   MPCCI_CONV_STATE_MAX  = 4
};

/*
 * Defaults for the nice parameter:
 * a lower nice value of a code means a higher priority, similar to UNIX process nice()
 */
enum
{
   MPCCI_NICE_INFINITE = -1,
   MPCCI_NICE_MIN      =  0,
   MPCCI_NICE_DEFAULT  =    MPCCI_NCODES_MAX,
   MPCCI_NICE_MAX      = (2*MPCCI_NCODES_MAX)
};

/*
 * Coordinates system (CSYS) definitions: The CSYS consists of the ...
 *       - coordinates system dimension (CDIM)
 *       - angular system dimension (ADIM)
 *       - coordinates system type (CTYP)
 *    and is specified within a single unsigned.
 *
 *    nibbles: |0000|0000|0000|0000|0000|0000|0000|0000|
 *                                      |CTYP|ADIM|CDIM|
 */
#define MPCCI_CSYS(_ctyp,_cdim)        ( (_ctyp) | (_cdim) )

#define MPCCI_CSYS_DIM(_csys)          ((MPCCI_UCAST(_csys)) & 0x0000000f)
#define MPCCI_CSYS_ANG(_csys)          ((MPCCI_UCAST(_csys)) & 0x000000f0)
#define MPCCI_CSYS_TYP(_csys)          ((MPCCI_UCAST(_csys)) & 0x00000f00)

#define MPCCI_CTYP_CART                0x00000100 /* cartesian: x,y,z */
#define MPCCI_CTYP_AXRZ                0x00000200 /* axissymm  : order radial,z (used intenally in the server) */
#define MPCCI_CTYP_AXZR                0x00000300 /* axissymm  : order z,radial (need to swap z<->r) */
#define MPCCI_CTYP_SPHR                0x00000400 /* spherical */
#define MPCCI_CTYP_NONE                0x00000500 /* used by network codes with integration points, no coords */

/* predefines                          ( CTYP            | ADIM       | CDIM) */
#define MPCCI_CSYS_BAD                 0xffffffff
#define MPCCI_CSYS_NET                 ( MPCCI_CTYP_NONE | 0x00000000 | 0 ) /* used by network code with integration points, no coords */
#define MPCCI_CSYS_C1D                 ( MPCCI_CTYP_CART | 0x00000000 | 1 )
#define MPCCI_CSYS_C2D                 ( MPCCI_CTYP_CART | 0x00000010 | 2 )
#define MPCCI_CSYS_C3D                 ( MPCCI_CTYP_CART | 0x00000030 | 3 )
#define MPCCI_CSYS_ARZ                 ( MPCCI_CTYP_AXRZ | 0x00000000 | 2 )
#define MPCCI_CSYS_AZR                 ( MPCCI_CTYP_AXZR | 0x00000000 | 2 )
#define MPCCI_CSYS_CYL                 ( MPCCI_CTYP_AXRZ | 0x00000030 | 3 )
#define MPCCI_CSYS_SPH                 ( MPCCI_CTYP_SPHR | 0x00000030 | 3 )



/*
 * Define the mesh dimensions based on the element types used
 */
enum
{
   MPCCI_MDIM_NONE=-1,
   MPCCI_MDIM_PNTS= 0, /* no elements: mapping on single points or point clouds with coords */
   MPCCI_MDIM_LINE= 1, /* FEA beam elements or 2D faces == edges */
   MPCCI_MDIM_FACE= 2, /* this is a (sur)face mesh */
   MPCCI_MDIM_VOLU= 3, /* this is a volume mesh */
   MPCCI_MDIM_IPNT= 4, /* this is a integration point (for network nodes without coords) */
   MPCCI_MDIM_GLOB= 5  /* symbolic dim of a global variable */
};

/*
 * Defines a special mesh ID which should NOT be used as a true mesh id
 */
#define MPCCI_MESHID_ANY   0x0fffffff

/*
 * Quantity id (QID) definitions: The QID consists of the ...
 *       - dimension                (DIM) [1...64]
 *       - interpolation            (IOK) [0... 1]
 *       - extrapolation            (XOK) [0... 1]
 *       - field interpolation type (INT) [1...15]
 *       - physical meaning         (PHY) [1...15]
 *       - 4 validity bits          (VAL) : 0=global, quantity valid for Volume/Face/Line/Point coupling
 *       - tensor order             (ORD) [0... 7] 0=scalar/1=vector/2=tensor/3=3dr order tensor...
 *       - ramping                  (ROK) [0... 1]
 *       - signature                (SIG) [0..255]
 *    and is specified within a single unsigned.
 *
 *    nibbles: |0000|0000|0|000|0000|0000|0000|00|00|0000|
 *             |---SIG---|R|ORD|-PHY|-VAL|-INT|XI|---DIM-|
 */
#define MPCCI_QID(_qsig,_qval,_qphy,_qtyp,_qdim)  ((_qsig)|(_qval)|(_qphy)|(_qtyp)|(_qdim))

#define MPCCI_QDIM_MAX                                        0x0000003f
#define MPCCI_QID_DIM(_qid)            (            (_qid)  & 0x0000003f)
#define MPCCI_QID_IOK(_qid)            ((MPCCI_UCAST(_qid)) & 0x00000040)
#define MPCCI_QID_XOK(_qid)            ((MPCCI_UCAST(_qid)) & 0x00000080)
#define MPCCI_QID_INT(_qid)            ((MPCCI_UCAST(_qid)) & 0x00000f00)
#define MPCCI_QID_VAL(_qid)            ((MPCCI_UCAST(_qid)) & 0x0000f000)
#define MPCCI_QID_PHY(_qid)            ((MPCCI_UCAST(_qid)) & 0x000f0000)
#define MPCCI_QID_ORD(_qid)            ((MPCCI_UCAST(_qid)) & 0x00700000)
#define MPCCI_QID_ROK(_qid)            ((MPCCI_UCAST(_qid)) & 0x00800000)
#define MPCCI_QID_SIG(_qid)            ((MPCCI_UCAST(_qid)) & 0xff000000)

/* all interpolation types < 0x0a00 are mesh based quantities */
#define MPCCI_QINT_COORD               0x00000100 /* is a grid coordinate/displacement */
#define MPCCI_QINT_FIELD               0x00000200 /* is a field */
#define MPCCI_QINT_FLUXI               0x00000300 /* is an integral flux */
#define MPCCI_QINT_FLUXD               0x00000400 /* is a flux density */

/* all interpolation types >= 0x0a00 are non mesh based global types */
#define MPCCI_GINT_MIN                 0x00000a00 /* minimum value of all */
#define MPCCI_GINT_MAX                 0x00000b00 /* maximum value of all */
#define MPCCI_GINT_SUM                 0x00000c00 /* sum of all values */
#define MPCCI_GINT_PRD                 0x00000d00 /* product of all values */
#define MPCCI_GINT_AVG                 0x00000e00 /* average of all values */
#define MPCCI_GINT_EQU                 0x00000f00 /* all values must be (nearly) identical */

/* validity */
#define MPCCI_QVAL_GLOB                0x00000000 /* global: invalid for all meshes */
#define MPCCI_QVAL_POINT               0x00001000 /* valid for point meshes */
#define MPCCI_QVAL_LINE                0x00002000 /* valid for line meshes */
#define MPCCI_QVAL_FACE                0x00004000 /* valid for face meshes */
#define MPCCI_QVAL_VOLU                0x00008000 /* valid for volume meshes */
#define MPCCI_QVAL_ALL                 0x0000f000 /* valid for all meshes */

/* physical meaning */
#define MPCCI_QPHY_ANY                 0x00000000 /* quantity is not further specified */
#define MPCCI_QPHY_MASS                0x00010000 /* is a mass sink/source */
#define MPCCI_QPHY_MOM                 0x00020000 /* is a momentum sink/source */
#define MPCCI_QPHY_ENTH                0x00030000 /* is a energy sink/source */
#define MPCCI_QPHY_PROP                0x00040000 /* is a property (material) */
#define MPCCI_QPHY_BCVAL               0x00050000 /* is a boundary condition value */
#define MPCCI_QPHY_BCGRAD              0x00060000 /* is a boundary condition wall normal gradient  */
#define MPCCI_QPHY_COORD               0x00070000 /* is a grid coordinate/displacement */
#define MPCCI_QPHY_MASSFR              0x00080000 /* is a chemical component concentration */

/* transformation type */
#define MPCCI_QORD_SCALAR              0x00000000
#define MPCCI_QORD_VECTOR              0x00100000
#define MPCCI_QORD_TENSOR              0x00200000



#define MPCCI_QID_IS_GLOB(_qid)       !MPCCI_QID_VAL(_qid) /* QVAL == 0 is a global quantity */
#define MPCCI_QID_IS_MESH(_qid)        MPCCI_QID_VAL(_qid) /* QVAL != 0 is a mesh quantity */

#define MPCCI_QID_IS_SCALAR(_qid)      ( MPCCI_QID_ORD(_qid) == MPCCI_QORD_SCALAR )
#define MPCCI_QID_IS_VECTOR(_qid)      ( MPCCI_QID_ORD(_qid) == MPCCI_QORD_VECTOR )
#define MPCCI_QID_IS_TENSOR(_qid)      ( MPCCI_QID_ORD(_qid) == MPCCI_QORD_TENSOR )

#define MPCCI_QID_IS_COORD(_qid)       ( MPCCI_QID_INT(_qid) == MPCCI_QINT_COORD )
#define MPCCI_QID_IS_FIELD(_qid)       ( MPCCI_QID_INT(_qid) == MPCCI_QINT_FIELD )
#define MPCCI_QID_IS_FLUXI(_qid)       ( MPCCI_QID_INT(_qid) == MPCCI_QINT_FLUXI )
#define MPCCI_QID_IS_FLUXD(_qid)       ( MPCCI_QID_INT(_qid) == MPCCI_QINT_FLUXD )

#define MPCCI_QID_IS_POINT(_qid)       ( (_qid) & MPCCI_QVAL_POINT )
#define MPCCI_QID_IS_LINE(_qid)        ( (_qid) & MPCCI_QVAL_LINE  )
#define MPCCI_QID_IS_FACE(_qid)        ( (_qid) & MPCCI_QVAL_FACE  )
#define MPCCI_QID_IS_VOLU(_qid)        ( (_qid) & MPCCI_QVAL_VOLU  )


/* translational kind of coordinate */
#define MPCCI_QID_IS_COORD_T(_qid)     ( (_qid) == MPCCI_QID_NPOSITION || (_qid) == MPCCI_QID_POINTPOSITION )

/* rotational coordinate */
#define MPCCI_QID_IS_COORD_R(_qid)     ( (_qid) == MPCCI_QID_ANGULARCOORD )

/*
 * Bitpattern for the flags used with [s]mpcci_getc()
 */
#define MPCCI_GETC_FLAG_SOURCE     0x00000001  /* codes with my source mesh/quantity */
#define MPCCI_GETC_FLAG_TARGET     0x00000002  /* codes with my target mesh/quantity */
#define MPCCI_GETC_FLAG_OTHERS     0x00000004  /* codes that are neither source nor target */
#define MPCCI_GETC_FLAG_ALL        (MPCCI_GETC_FLAG_SOURCE|MPCCI_GETC_FLAG_TARGET|MPCCI_GETC_FLAG_OTHERS)
#define MPCCI_GETC_FLAG_NAMES      0x00000010

/* extra bit: only "waiting" codes with my target mesh/quantity */
#define MPCCI_GETC_FLAG_TWAIT      0x00000100

/* extra bit: this is a SYNC operation and ALL code clients must call smpcci_getc() */
#define MPCCI_GETC_FLAG_SYNC       0x00001000


/*
 * Bitpattern for the flags used with [s]mpcci_qtag()
 */
#define MPCCI_QTAG_FLAG_WAIT        0x00000100
#define MPCCI_QTAG_FLAG_SYNC        0x00001000

/*
 * Special negative tag.qtid return values for the mpcci_qtag() query
 */
enum
{
   MPCCI_QTAG_STATE_EMPTY   = -1, /* this tag is valid, but the quantity buffer is empty */
   MPCCI_QTAG_STATE_BADMID  = -2, /* tag has an invalid mesh ID */
   MPCCI_QTAG_STATE_BADQID  = -3, /* tag has an invalid quantity ID */
   MPCCI_QTAG_STATE_INVALID = -4  /* tag has an invalid qtid */
};


/*
 * Bitpattern for the flags used with [s]mpcci_getv()
 */
#define MPCCI_GETV_FLAG_ENV         0x00000001  /* get environment variable */
#define MPCCI_GETV_FLAG_CODE        0x00000002  /* get code specific variable */


/*
 * Bitpattern for the flags used with [s]mpcci_getq()
 */
#define MPCCI_GETQ_FLAG_WAIT        0x00000001
#define MPCCI_GETQ_FLAG_SYNC        0x00000010

#define MPCCI_GETQ_FLAG_DOT         0x00000100 /* get first time derivative of quantity */
#define MPCCI_GETQ_FLAG_DOTDOT      0x00000200 /* get second time derivative of quantity */

#define MPCCI_GETQ_FLAG_PRED        0x00001000 /* get predictor value for quantity */

/*
 * Bitpattern for the flags used with ampcci_remesh()
 */
#define MPCCI_REMESH_FLAG_CHECK     0x00000000
#define MPCCI_REMESH_FLAG_NODES     0x00000001
#define MPCCI_REMESH_FLAG_ELEMS     0x00000003
#define MPCCI_REMESH_FLAG_ALLPARTS  0x00000010
#define MPCCI_REMESH_FLAG_FULL      (MPCCI_REMESH_FLAG_ELEMS|MPCCI_REMESH_FLAG_ALLPARTS)

/*
 * Various remeshing states
 */
enum
{
   MPCCI_REMESH_STATE_NONE    = 0,
   MPCCI_REMESH_STATE_SMOOTHED,
   MPCCI_REMESH_STATE_REMESHED,
   MPCCI_REMESH_STATE_MIGRATED
};



/*
 * Define the quantity storage method (used by the adapter codes)
 */
enum
{
   MPCCI_QSM_UNDEF=-1,  /* invalid smethod */
   MPCCI_QSM_DIRECT,    /* direct written into codes buffer */
   MPCCI_QSM_BUFFER,    /* quantity is buffered locally and copied later */
   MPCCI_QSM_USRMEM,    /* quantity stored in users indexed memory */
   MPCCI_QSM_SCALAR,    /* quantity stored in users index scalars */
   MPCCI_QSM_SPECIES    /* quantity stored as chemical species */
};

/*
 * Quantities transfer status
 * this status is important, if the quantities are stored in buffers
 * and we may access buffers before they are filled with valid data
 * the transfer status gives us information about the contents of this buffer
 */
enum
{
   MPCCI_QUANT_TSTATE_NONE = 0,  /* no transfer was done before */
   MPCCI_QUANT_TSTATE_SEND =-1,  /* quantity was sent */
   MPCCI_QUANT_TSTATE_RECV = 1   /* quantity was received */
};

/*
 * The flags bitpattern for all kind of objects on the client and server side
 */

/*
 * Client or code flags (CFLAG) resp. attributes
 */
#define MPCCI_CFLAG_GRID_MASK          0x0000000f  /* bit mask to check for grid type */
#define MPCCI_CFLAG_GRID_CURR          0x00000001  /* client has the current grid available */
#define MPCCI_CFLAG_GRID_ORIG          0x00000002  /* client has the original grid available */

#define MPCCI_CFLAG_TYPE_MASK          0x00000ff0  /* bit mask to check for the code type */
#define MPCCI_CFLAG_TYPE_FV            0x00000010  /* client code is FV type (at most CFD) */
#define MPCCI_CFLAG_TYPE_FEA           0x00000020  /* client code is FEA type */
#define MPCCI_CFLAG_TYPE_NET           0x00000030  /* client code is NETWORK type */
#define MPCCI_CFLAG_TYPE_RAD           0x00000040  /* client code is RADIADION type */
#define MPCCI_CFLAG_TYPE_SBM           0x00000050  /* client code is SOLID-BODY-MOTION type */
#define MPCCI_CFLAG_TYPE_MOL           0x00000060  /* client code is MOLECULAR-DYNAMIC type */

#define MPCCI_CFLAG_SIMPLE_PCELL       0x00001000  /* bit set if code request a simplification of polyhedrons */
#define MPCCI_CFLAG_SINGLE_BUF         0x00002000  /* bit set if code request a single buffer only */
#define MPCCI_CFLAG_REPAIR_ELEM        0x00004000  /* bit set if code request a simplification of polyhedrons */

/* flags only for internal use by the server */
#define MPCCI_CFLAG_SCALE              0x00010000  /* bit set if mesh needs to be scaled */
#define MPCCI_CFLAG_TFM44              0x00020000  /* bit set if mesh needs to be transformed */
#define MPCCI_CFLAG_MHALT              0x00100000  /* bit set if client is temporary halted within a mirror function */
#define MPCCI_CFLAG_CHALT              0x00200000  /* bit set if code is temporary halted since not all clients are connected */
#define MPCCI_CFLAG_TSAVE              0x00400000  /* bit set if client called SAVE (trace save) */
#define MPCCI_CFLAG_TUSED              0x00800000  /* bit set if client called PUTQ/GETQ (tag was used) */
#define MPCCI_CFLAG_SYNC               0x01000000  /* bit set if this is a multi client job and needs synched operations */


#define MPCCI_CODE_HAS_CURRG(_c)       ( (_c)->flags & MPCCI_CFLAG_GRID_CURR )
#define MPCCI_CODE_HAS_ORIGG(_c)       ( (_c)->flags & MPCCI_CFLAG_GRID_ORIG )

#define MPCCI_CODE_TYPE(_c)            ( (_c)->flags & MPCCI_CFLAG_TYPE_MASK )
#define MPCCI_CODE_IS_FV(_c)           ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_FV  )
#define MPCCI_CODE_IS_FEA(_c)          ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_FEA )
#define MPCCI_CODE_IS_NET(_c)          ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_NET )
#define MPCCI_CODE_IS_RAD(_c)          ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_RAD )
#define MPCCI_CODE_IS_SBM(_c)          ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_SBM )
#define MPCCI_CODE_IS_MOL(_c)          ( MPCCI_CODE_TYPE(_c) == MPCCI_CFLAG_TYPE_MOL )

#define MPCCI_CODE_IS_MESHBASED(_c)    (MPCCI_CODE_IS_FV (_c)||MPCCI_CODE_IS_FEA(_c)||MPCCI_CODE_IS_RAD(_c))
#define MPCCI_CODE_IS_POINTBASED(_c)   (MPCCI_CODE_IS_NET(_c)||MPCCI_CODE_IS_SBM(_c)||MPCCI_CODE_IS_MOL(_c))

/*
 * Quantity flags (QFLAG) resp. attributes
 */
#define MPCCI_QFLAG_UNIT_MASK          0x00000001
#define MPCCI_QFLAG_UNIT_VAR           0x00000001  /* variable unit */

#define MPCCI_QFLAG_DIR_MASK           0x000000f0
#define MPCCI_QFLAG_DIR_SEND           0x00000010  /* quantity sent by the client (received by the server ) */
#define MPCCI_QFLAG_DIR_RECV           0x00000020  /* quantity received by the client (send by the server) */

#define MPCCI_QFLAG_LOC_MASK           0x00000f00
#define MPCCI_QFLAG_LOC_CODE           0x00000000  /* no bit set, then the code itself decides */
#define MPCCI_QFLAG_LOC_NODE           0x00000100
#define MPCCI_QFLAG_LOC_ELEM           0x00000200
#define MPCCI_QFLAG_LOC_NODE_O         0x00000400  /* special: indicate node orphan level */
#define MPCCI_QFLAG_LOC_ELEM_O         0x00000800  /* special: indicated element orphan level */

#define MPCCI_QFLAG_FILL_AUTO          0x00001000  /* automatic determination of the orphaned fill value */
#define MPCCI_QFLAG_RELAX_AUTO         0x00002000  /* automatic determination of the relaxation value */

/* for network codes */
#define MPCCI_QFLAG_LOC_POINT          MPCCI_QFLAG_LOC_NODE

/* FV code may use these guys */
#define MPCCI_QFLAG_LOC_VERT           MPCCI_QFLAG_LOC_NODE
#define MPCCI_QFLAG_LOC_CELL           MPCCI_QFLAG_LOC_ELEM


#define MPCCI_QFLAG_VISITED            0x01000000  /* adapter helper bit: indicates quant was visited */


/*
 * Lots of wrapper macros used by the coupling manager and the server
 */
#ifdef MPCCI_COMPILE_SERVER
#define MPCCI_QUANT_QID(_q)            ( (_q)->id    )
#else
#define MPCCI_QUANT_QID(_q)            ( (_q)->qid   )
#endif

#define MPCCI_QUANT_NAME(_q)           ( (_q)->name  )
#define MPCCI_QUANT_STATE(_q)          ( (_q)->state )
#define MPCCI_QUANT_FLAGS(_q)          ( (_q)->flags )

#define MPCCI_QUANT_UNIT(_q)           ( (_q)->flags & MPCCI_QFLAG_UNIT_MASK )
#define MPCCI_QUANT_LOC(_q)            ( (_q)->flags & MPCCI_QFLAG_LOC_MASK  )
#define MPCCI_QUANT_DIR(_q)            ( (_q)->flags & MPCCI_QFLAG_DIR_MASK  )
#define MPCCI_QUANT_DIM(_q)            MPCCI_QID_DIM( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IOK(_q)            MPCCI_QID_IOK( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_XOK(_q)            MPCCI_QID_XOK( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_INT(_q)            MPCCI_QID_INT( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_PHY(_q)            MPCCI_QID_PHY( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_ORD(_q)            MPCCI_QID_ORD( MPCCI_QUANT_QID(_q) )

#define MPCCI_QUANT_IS_VAR(_q)         ( MPCCI_QUANT_UNIT(_q) != 0 )
#define MPCCI_QUANT_IS_GLOB(_q)          MPCCI_QID_IS_GLOB  ( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_MESH(_q)          MPCCI_QID_IS_MESH  ( MPCCI_QUANT_QID(_q) )

#define MPCCI_QUANT_IS_SCALAR(_q)        MPCCI_QID_IS_SCALAR( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_VECTOR(_q)        MPCCI_QID_IS_VECTOR( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_TENSOR(_q)        MPCCI_QID_IS_TENSOR( MPCCI_QUANT_QID(_q) )

#define MPCCI_QUANT_IS_COORD(_q)         MPCCI_QID_IS_COORD ( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_FIELD(_q)         MPCCI_QID_IS_FIELD ( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_FLUXI(_q)         MPCCI_QID_IS_FLUXI ( MPCCI_QUANT_QID(_q) )
#define MPCCI_QUANT_IS_FLUXD(_q)         MPCCI_QID_IS_FLUXD ( MPCCI_QUANT_QID(_q) )

#define MPCCI_QUANT_IS_NODE(_q)        ( MPCCI_QUANT_LOC(_q) == MPCCI_QFLAG_LOC_NODE )
#define MPCCI_QUANT_IS_ELEM(_q)        ( MPCCI_QUANT_LOC(_q) == MPCCI_QFLAG_LOC_ELEM )
#define MPCCI_QUANT_IS_SEND(_q)        ( MPCCI_QUANT_DIR(_q) == MPCCI_QFLAG_DIR_SEND ) /* quantity sent by the client (received by the server ) */
#define MPCCI_QUANT_IS_RECV(_q)        ( MPCCI_QUANT_DIR(_q) == MPCCI_QFLAG_DIR_RECV ) /* quantity received by the client (sent by the server) */



/*
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

/* moving reference frame: neither moving grid nor moving quantity */
#define MPCCI_MOTION_FRAME_T           0x00000001  /* translate only */
#define MPCCI_MOTION_FRAME_R           0x00000010  /* rotate only */
#define MPCCI_MOTION_FRAME_RT          0x00000111  /* rotate, then translate order */
#define MPCCI_MOTION_FRAME_TR          0x00000211  /* translate, then rotate order */
#define MPCCI_MOTION_FRAME_CO          0x00000311  /* simultanious continuous transformations */

/* moving grid: grid and quantities are moved, but the grid is not necessarily sent to the server */
#define MPCCI_MOTION_MGRID_T           0x00011001  /* translate only */
#define MPCCI_MOTION_MGRID_R           0x00011010  /* rotate only */
#define MPCCI_MOTION_MGRID_RT          0x00011111  /* rotate, then translate order */
#define MPCCI_MOTION_MGRID_TR          0x00011211  /* translate, then rotate order */
#define MPCCI_MOTION_MGRID_CO          0x00011311  /* simultanious continuous transformations */

#define MPCCI_MOTION_ACC_T             0x00000002 /* velocity is in fact an acceleration */
#define MPCCI_MOTION_ACC_R             0x00000020 /* omega is in fact an angular acceleration */
#define MPCCI_MOTION_ACC               0x00000022

#endif
