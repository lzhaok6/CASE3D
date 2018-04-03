#ifndef __MPCCI4_COREAPI_HEADER_INCLUDED
#define __MPCCI4_COREAPI_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_coreapi.h
 *
 * Include header file with basic function prototypes
 *
 * $Id: mpcci_coreapi.h 322 2014-02-10 09:49:07Z wirth $
 *
 ****************************************************************************************/

#include "mpcci_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * General utility functions which do not require any server/job objects
 */
MPCCI_EXPORT void umpcci_mem_functs
                  (
                     void *(*mallocp )(size_t s),
                     void *(*reallocp)(void *p, size_t s),
                     void  (*freep   )(void *p)
                  );

/* just redefineable function pointers */
MPCCI_EXPORT char      *(*umpcci_getenv)  (const char *name);
MPCCI_EXPORT const char **umpcci_envvars  (void);

MPCCI_EXPORT MPCCI_QINFO *umpcci_quant_to_qinfo(MPCCI_QUANT *quant, MPCCI_QINFO *qinfo);


MPCCI_EXPORT void          umpcci_idmap_alloc   (MPCCI_IDMAP *idmap, int n);
MPCCI_EXPORT void          umpcci_idmap_realloc (MPCCI_IDMAP *idmap, int n);
MPCCI_EXPORT void          umpcci_idmap_free    (MPCCI_IDMAP *idmap);
MPCCI_EXPORT MPCCI_PART   *umpcci_idmap_part    (MPCCI_IDMAP *idmap, int globalid);
MPCCI_EXPORT int           umpcci_idmap_index   (MPCCI_IDMAP *idmap, int globalid);
MPCCI_EXPORT void         *umpcci_idmap_sbuffer (MPCCI_IDMAP *idmap, int globalid, int qid);
MPCCI_EXPORT void         *umpcci_idmap_value   (MPCCI_IDMAP *idmap, int globalid, int qid, size_t realsize);
MPCCI_EXPORT void          umpcci_idmap_assign  (MPCCI_IDMAP *idmap, int globalid, MPCCI_PART *part, int index);

MPCCI_EXPORT void         *umpcci_gmem_ctrl     (const int cmd, size_t rsize);
MPCCI_EXPORT const char   *umpcci_conv_cstate   (const int conv);
MPCCI_EXPORT int           umpcci_conv_common   (const int c1, const int c2);

#define umpcci_gmem_get(_rsize)  umpcci_gmem_ctrl( 1,_rsize)
#define umpcci_gmem_unlock()     umpcci_gmem_ctrl( 0,0)
#define umpcci_gmem_free()       umpcci_gmem_ctrl(-1,0)

/*
 * Functions for a single server communication.
 * This is the lowest level of communication.
 */
MPCCI_EXPORT MPCCI_SERVER *smpcci_init(const char *porthost, const MPCCI_CINFO *cinfo);


MPCCI_EXPORT int  smpcci_quit(MPCCI_SERVER **pserver);
MPCCI_EXPORT int  smpcci_stop(MPCCI_SERVER **pserver, const char *msg);
MPCCI_EXPORT int  smpcci_ferr(MPCCI_SERVER **pserver, const char *msg);
MPCCI_EXPORT int  smpcci_sbye(MPCCI_SERVER **pserver, const char *cmd, const char *msg);
MPCCI_EXPORT int  smpcci_scmd(MPCCI_SERVER   *server, const char *cmd, int getresp);
MPCCI_EXPORT int  smpcci_ckdt(MPCCI_SERVER   *server);
MPCCI_EXPORT int  smpcci_isok(MPCCI_SERVER   *server);


/* (def)ine (p)art, (q)uantity and (g)lobal */
MPCCI_EXPORT int  smpcci_defp(MPCCI_SERVER *server, int mid, int pid, unsigned csys, int nnodes, int nelems, const char *name);
MPCCI_EXPORT int  smpcci_defq(MPCCI_SERVER *server, int mid, const MPCCI_QINFO *qinfo);
MPCCI_EXPORT int  smpcci_defg(MPCCI_SERVER *server, const MPCCI_GLOB *glob);


/* (del)ete (m)mesh, (p)art, (q)uantity and (g) */
MPCCI_EXPORT int  smpcci_delm(MPCCI_SERVER *server, int mid);
MPCCI_EXPORT int  smpcci_delp(MPCCI_SERVER *server, int mid, int pid);
MPCCI_EXPORT int  smpcci_delq(MPCCI_SERVER *server, int mid, int qid);
MPCCI_EXPORT int  smpcci_delg(MPCCI_SERVER *server, int gid, int qid);


/* define (p)art node/normals/sizes/elements/motion */
MPCCI_EXPORT int  smpcci_pnod(MPCCI_SERVER *server, int mid, int pid, unsigned csys, int nnodes, const void *coords , unsigned realsize, const int *nodeids, const void *angles);
MPCCI_EXPORT int  smpcci_pels(MPCCI_SERVER *server, int mid, int pid, int nelems, unsigned eltype1, const unsigned *eltypes, const int *elnodes, const int *elemids);
MPCCI_EXPORT int  smpcci_pmot(MPCCI_SERVER *server, int mid, int pid, const MPCCI_MOTION *motion);
MPCCI_EXPORT int  smpcci_pper(MPCCI_SERVER *server, int mid, int pid, const MPCCI_PERIODIC *periodic);
MPCCI_EXPORT int  smpcci_pshf(MPCCI_SERVER *server, int mid, int pid, double sval);

#if 0 /* used later for server based morphing */
MPCCI_EXPORT int  smpcci_mset(MPCCI_SERVER *server, int mid, int pid, unsigned flags, unsigned realsize);
MPCCI_EXPORT int  smpcci_mdef(MPCCI_SERVER *server, int id, const void *coord, unsigned realsize);
MPCCI_EXPORT int  smpcci_mrun(MPCCI_SERVER *server);
MPCCI_EXPORT int  smpcci_mget(MPCCI_SERVER *server, int mid, int pid, int nnodes, void *coords, unsigned realsize);
#endif


/* (q)uery (tag)s */
MPCCI_EXPORT int  smpcci_qtag(MPCCI_SERVER *server, MPCCI_QTAG *tagv, int ntags, unsigned flags, MPCCI_QSTAT *qstat);


/* put/get: (q)uantity, (g)lobal, (e)nvironment, (p)arameter, (c)odes */
MPCCI_EXPORT int  smpcci_putg(MPCCI_SERVER *server, int gid, int qid, const void *values, unsigned realsize);
MPCCI_EXPORT int  smpcci_getg(MPCCI_SERVER *server, int gid, int qid,       void *values, unsigned realsize);

MPCCI_EXPORT int  smpcci_putq(MPCCI_SERVER *server, int mid, int pid, int qid, int nval, const void *values, unsigned realsize);
MPCCI_EXPORT int  smpcci_getq(MPCCI_SERVER *server, int mid, int pid, int qid, int nval,       void *values, unsigned realsize, unsigned flags, int *conv, double *rdt);

MPCCI_EXPORT int  smpcci_getv(MPCCI_SERVER *server, unsigned flags, const char *var, char *val, size_t size);
MPCCI_EXPORT int  smpcci_getc(MPCCI_SERVER *server, unsigned flags, char ***codev, int *nwait);


/* (q)uery (q)uantity by index/name/id */
MPCCI_EXPORT int  smpcci_qqix(MPCCI_SERVER *server, int mid, int index        , MPCCI_QINFO *qinfo);
MPCCI_EXPORT int  smpcci_qqnm(MPCCI_SERVER *server, int mid, const char *qname, MPCCI_QINFO *qinfo);
MPCCI_EXPORT int  smpcci_qqid(MPCCI_SERVER *server, int mid, int qid          , MPCCI_QINFO *qinfo);

/* (q)uery (p)part and (g)lobal by index */
MPCCI_EXPORT int  smpcci_qpix(MPCCI_SERVER *server, int index, MPCCI_PINFO *pinfo);
MPCCI_EXPORT int  smpcci_qgix(MPCCI_SERVER *server, int index, MPCCI_GLOB *glob);


MPCCI_EXPORT int  smpcci_savj(MPCCI_SERVER *server);
MPCCI_EXPORT int  smpcci_savc(MPCCI_SERVER *server);
MPCCI_EXPORT int  smpcci_nice(MPCCI_SERVER *server, int nice);

/* (syn)chonize (t)ime, (j)ob, (c)ode, (m)esh, (q)uant, (g)lobal */
MPCCI_EXPORT int  smpcci_sync(MPCCI_SERVER *server, unsigned flags); /* general sync */
MPCCI_EXPORT int  smpcci_synt(MPCCI_SERVER *server, double time, double dt, int iter, int conv);
MPCCI_EXPORT int  smpcci_synm(MPCCI_SERVER *server);  /* sync all mesh (all code clients) */

MPCCI_EXPORT void smpcci_list(MPCCI_SERVER *server);

/* find objects in the servers structure */
MPCCI_EXPORT MPCCI_PART  *smpcci_find_partbyid  (MPCCI_SERVER *server, const int mid, const int pid);
MPCCI_EXPORT MPCCI_PART  *smpcci_find_partbyname(MPCCI_SERVER *server, const char *pname);
MPCCI_EXPORT MPCCI_PART  *smpcci_find_partbyusrp(MPCCI_SERVER *server, const void *p, const int uid, MPCCI_PART *(*cmp)(MPCCI_PART *, const int uid, const void *));
MPCCI_EXPORT MPCCI_QUANT *smpcci_find_quantbyid (MPCCI_SERVER *server, const int mid, const int pid, const int qid);
MPCCI_EXPORT void        *smpcci_find_sbuffbyid (MPCCI_SERVER *server, const int mid, const int pid, const int qid);

/* clear or check the vistited bit in the quant->flags */
MPCCI_EXPORT int        smpcci_qvisit_clear  (MPCCI_SERVER *server, unsigned fmask);
MPCCI_EXPORT int        smpcci_qvisit_insync (MPCCI_SERVER *server, unsigned fmask);

/*
 * Functions for multi server communication.
 */
MPCCI_EXPORT MPCCI_JOB *mpcci_init(const char *porthost, const MPCCI_CINFO *cinfo);

MPCCI_EXPORT int  mpcci_quit(MPCCI_JOB **pjob);
MPCCI_EXPORT int  mpcci_stop(MPCCI_JOB **pjob, const char *msg);
MPCCI_EXPORT int  mpcci_ferr(MPCCI_JOB **pjob, const char *msg);
MPCCI_EXPORT int  mpcci_sbye(MPCCI_JOB **pjob, const char *cmd, const char *msg);
MPCCI_EXPORT int  mpcci_scmd(MPCCI_JOB   *job, const char *cmd, int getresp);
MPCCI_EXPORT int  mpcci_ckdt(MPCCI_JOB   *job);
MPCCI_EXPORT int  mpcci_isok(MPCCI_JOB   *job);

MPCCI_EXPORT int  mpcci_pmot(MPCCI_JOB *job, int mid, int pid, const MPCCI_MOTION *motion);
MPCCI_EXPORT int  mpcci_pshf(MPCCI_JOB *job, int mid, int pid, double sval);

MPCCI_EXPORT int  mpcci_defp(MPCCI_JOB *job, int mid, int pid, unsigned csys, int nnodes, int nelems, const char *name);
MPCCI_EXPORT int  mpcci_defq(MPCCI_JOB *job, int mid, const MPCCI_QINFO *qinfo);
MPCCI_EXPORT int  mpcci_defg(MPCCI_JOB *job, const MPCCI_GLOB *glob);

MPCCI_EXPORT int  mpcci_delm(MPCCI_JOB *job, int mid);
MPCCI_EXPORT int  mpcci_delp(MPCCI_JOB *job, int mid, int pid);
MPCCI_EXPORT int  mpcci_delq(MPCCI_JOB *job, int mid, int qid);
MPCCI_EXPORT int  mpcci_delg(MPCCI_JOB *job, int gid, int qid);

MPCCI_EXPORT int  mpcci_putg(MPCCI_JOB *job, int gid, int qid, const void *values, unsigned realsize);
MPCCI_EXPORT int  mpcci_getg(MPCCI_JOB *job, int gid, int qid,       void *values, unsigned realsize);

MPCCI_EXPORT int  mpcci_putq(MPCCI_JOB *job, int mid, int pid, int qid, int nval, const void *values, unsigned realsize);
MPCCI_EXPORT int  mpcci_getq(MPCCI_JOB *job, int mid, int pid, int qid, int nval,       void *values, unsigned realsize, unsigned flags, int *conv, double *rdt);

MPCCI_EXPORT int  mpcci_getv(MPCCI_JOB *job, unsigned flags, const char *var, char *val, size_t size);
MPCCI_EXPORT int  mpcci_getc(MPCCI_JOB *job, unsigned flags, char ***codev, int *nwait);

MPCCI_EXPORT int  mpcci_sync(MPCCI_JOB *job, unsigned flags);
MPCCI_EXPORT int  mpcci_synt(MPCCI_JOB *job, double time, double dt, int iter, int conv);
MPCCI_EXPORT int  mpcci_synm(MPCCI_JOB *job);

MPCCI_EXPORT int  mpcci_savj(MPCCI_JOB *job);
MPCCI_EXPORT int  mpcci_savc(MPCCI_JOB *job);
MPCCI_EXPORT int  mpcci_nice(MPCCI_JOB *job, int nice);

MPCCI_EXPORT int  mpcci_qtag(MPCCI_JOB *job, MPCCI_QTAG *tagv, int ntags, unsigned flags, MPCCI_QSTAT *qstat);

MPCCI_EXPORT void mpcci_list(MPCCI_JOB *job);

/* find objects in the job structure */
MPCCI_EXPORT MPCCI_PART  *mpcci_find_partbyid   (MPCCI_JOB *job, const int mid, const int pid);
MPCCI_EXPORT MPCCI_PART  *mpcci_find_partbyname (MPCCI_JOB *job, const char *pname);
MPCCI_EXPORT MPCCI_PART  *mpcci_find_partbyusrp (MPCCI_JOB *job, const void *p, const int uid, MPCCI_PART *(*cmp)(MPCCI_PART *, const int uid, const void *));
MPCCI_EXPORT MPCCI_QUANT *mpcci_find_quantbyid  (MPCCI_JOB *job, const int mid, const int pid, const int qid);
MPCCI_EXPORT void        *mpcci_find_sbuffbyid  (MPCCI_JOB *job, const int mid, const int pid, const int qid);

/* clear or check the vistited bit in the quant->flags */
MPCCI_EXPORT int          mpcci_qvisit_clear (MPCCI_JOB *job, unsigned fmask);
MPCCI_EXPORT int          mpcci_qvisit_insync(MPCCI_JOB *job, unsigned fmask);

/*
 * Auxiliary adapter level functions.
 *
 * The set of ampcci_xxx() functions handles all the stuff needed in the coupled code.
 *   - intialize all the stuff
 *   - transfer all quantities on all components
 *   - notify the remeshing of mesh based components to MpCCI
 * what else more is needed?
 */
MPCCI_EXPORT int  ampcci_tinfo_init(MPCCI_TINFO *tinfo, const char *config);
MPCCI_EXPORT int  ampcci_config    (MPCCI_JOB *job, const MPCCI_DRIVER *driver);
MPCCI_EXPORT int  ampcci_recvctl   (MPCCI_JOB *job, const MPCCI_TINFO *tinfo);
MPCCI_EXPORT int  ampcci_sendctl   (MPCCI_JOB *job, const MPCCI_TINFO *tinfo);
MPCCI_EXPORT int  ampcci_transfer  (MPCCI_JOB *job, MPCCI_TINFO *tinfo);
MPCCI_EXPORT int  ampcci_remesh    (MPCCI_JOB *job, unsigned flags);


/*
 * Prototypes for the external morpher communication
 *
 * low level core functions
 */
MPCCI_EXPORT MORPHER *MORPHER_connect     (const char *porthost, int nproc, int vtype, unsigned realsize);
MPCCI_EXPORT void     MORPHER_close       (MORPHER **pmorpher);
MPCCI_EXPORT void     MORPHER_exit        (MORPHER **pmorpher);
MPCCI_EXPORT void     MORPHER_send_start  (MORPHER *morpher, int type, ...);
MPCCI_EXPORT void     MORPHER_send_vertex (MORPHER *morpher, int idx, const void *vcoord);
MPCCI_EXPORT void     MORPHER_send_stop   (MORPHER *morpher);
MPCCI_EXPORT int      MORPHER_recv_start  (MORPHER *morpher);
MPCCI_EXPORT int      MORPHER_recv_vertex (MORPHER *morpher, void *vcoord);
MPCCI_EXPORT void     MORPHER_info        (MORPHER *morpher);

/*
 * high level functions: use coords, nverts and ioffs
 */
MPCCI_EXPORT MORPHER *MORPHER_coordv_init (const char *porthost, int nclients, void *vcoords, int nverts, int ioffs, unsigned realsize, unsigned flags);
MPCCI_EXPORT void     MORPHER_coordv_start(MORPHER *morpher, void *vcoords);
MPCCI_EXPORT void     MORPHER_coordv_moved(MORPHER *morpher, int   index  );
MPCCI_EXPORT void     MORPHER_coordv_exec (MORPHER *morpher, void *vcoords);


#ifdef __cplusplus
}
#endif

#endif
