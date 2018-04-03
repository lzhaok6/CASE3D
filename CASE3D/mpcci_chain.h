#ifndef __MPCCI4_CHAIN_HEADER_INCLUDED
#define __MPCCI4_CHAIN_HEADER_INCLUDED
/* mpcci_chain.h
 *
 *****************************************************************************************
 *
 * Purpose:
 *    Defines basic chained objects class and declares the C-methods
 *
 * Author:
 *    Carsten Dehning <carsten.dehning@scai.fhg.de>
 *
 * Reviews/changes:
 *    2008/Mar: Carsten Dehning, Initial release
 *    $Id: mpcci_chain.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 *****************************************************************************************
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _CHAIN      CHAIN;      /* chained object class */
typedef struct _CLASSINFO  CLASSINFO;  /* inherited class info */


/*
 * structure with class specific functions and information
 */
typedef void  (CVF_CLEAN)(      void *obj);
typedef void  (CVF_DUMP )(const void *obj, FILE *fp);
typedef void  (CVF_LOAD )(      void *obj, FILE *fp);

struct _CLASSINFO
{
   const char *class_name;    /* class name */

   /* class virtual function table ... */
   CVF_CLEAN  *class_clean;   /* pointer to cleanup object function */
   CVF_DUMP   *class_dump;    /* pointer to the dump object function */
   CVF_LOAD   *class_load;    /* pointer to the load object function */

   unsigned    class_size;    /* usizeof(struct _CLASS) */
};


#define CHAIN_DEFINE_CLASSINFO(_class_sym,_clean_fnc,_dump_fnc,_load_fnc)\
static const CLASSINFO class_info = \
{\
   #_class_sym,\
   _clean_fnc,\
   _dump_fnc,\
   _load_fnc,\
   (unsigned)sizeof(_class_sym)\
}


/*
 * common (public) data members for all chained objects
 */
#define CHAIN_DEFINE_CLASSBASE(_type) \
   const CLASSINFO *class_info;  /* pointer to class information */ \
   _type           *next;        /* next object in chain */ \
   _type           *prev;        /* previous object in chain */ \
   const char      *name;        /* strdup()'ed name of this object */ \
   unsigned         flags;       /* flags or attribute bits for this object */ \
   int              id;          /* id of this object */ \
   int              state;       /* state of this object */


struct _CHAIN
{
   CHAIN_DEFINE_CLASSBASE(CHAIN)
};


#if !INCLUDE_STATIC

   C_FUNC_PREFIX void *CHAIN_new             (const CLASSINFO *class_info, int id, unsigned flags, const char *fmt, ...);
   C_FUNC_PREFIX void  CHAIN_free            (      void *obj);
   C_FUNC_PREFIX void *CHAIN_dup             (const void *obj);
   C_FUNC_PREFIX void  CHAIN_print           (const void *obj, const char *prefix, const char *idfmt, const char *fmt, ...);
   C_FUNC_PREFIX void  CHAIN_dump            (const void *obj,   FILE *fp);
   C_FUNC_PREFIX void  CHAIN_dumpall         (const void  *head, FILE *fp);

   C_FUNC_PREFIX int   CHAIN_append          (      void **phead, void *obj);
   C_FUNC_PREFIX void  CHAIN_remove          (      void **phead, void *obj, int keep);
   C_FUNC_PREFIX void  CHAIN_destroy         (      void **phead);
   C_FUNC_PREFIX void  CHAIN_movetoend       (      void **phead, void *obj);

   C_FUNC_PREFIX void *CHAIN_getlast         (const void  *head);
   C_FUNC_PREFIX void *CHAIN_getbyname       (const void  *head, const char *name);
   C_FUNC_PREFIX void *CHAIN_getbyid         (const void  *head, const int id);
   C_FUNC_PREFIX void *CHAIN_getbypos        (const void  *head,       int ipos);
   C_FUNC_PREFIX int   CHAIN_getpos          (const void  *head, const void *obj);
   C_FUNC_PREFIX int   CHAIN_getcount        (const void  *head);

   C_FUNC_PREFIX void  CHAIN_setstate        (      void  *head, const int state);
   C_FUNC_PREFIX int   CHAIN_getstate        (const void  *head);
   C_FUNC_PREFIX int   CHAIN_hasstate        (const void  *head, const int state);

   C_FUNC_PREFIX void  CHAIN_flags_set       (      void  *head, unsigned bits);
   C_FUNC_PREFIX void  CHAIN_flags_clear     (      void  *head, unsigned bits);
   C_FUNC_PREFIX int   CHAIN_flags_test      (const void  *head, unsigned bits);
   C_FUNC_PREFIX int   CHAIN_flags_testmask  (const void  *head, unsigned mask, unsigned res);

#endif


#define CHAIN_DESTROY_NULL(_head)   if (_head) CHAIN_destroy(PPVOID_ADDR(_head))

#define CHAIN_APPEND(_head,_apnd)\
{\
   if (_head)\
   {\
      (_apnd)->prev = (_head);\
      if (((_apnd)->next = (_head)->next) != NULL)\
         ((_apnd)->next)->prev = (_apnd);\
      (_head)->next = (_apnd);\
   }\
   else\
   {\
      (_apnd)->next = (_apnd)->prev = NULL;\
      (_head) = (_apnd);\
   }\
}

#ifdef __cplusplus
}
#endif

#endif
