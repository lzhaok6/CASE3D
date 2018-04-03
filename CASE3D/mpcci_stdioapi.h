#ifndef __MPCCI4_STDIOAPI_HEADER_INCLUDED
#define __MPCCI4_STDIOAPI_HEADER_INCLUDED
/******************************************************************************************
 * mpcci_stdioapi.h
 *
 * low level IO-wrapper interface
 *
 * $Id: mpcci_stdioapi.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 *****************************************************************************************/

#include "mpcci_platform.h"
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif


/* just redefineable function pointers */
#if defined(__GNUC__) && !defined(__ICC) /* let gcc/g++ check the format string and arguments */
MPCCI_EXPORT int  (*umpcci_printf)  (int level, const char *fmt, ...) __attribute__ ((format(printf,2,3)));
#else
MPCCI_EXPORT int  (*umpcci_printf)  (int level, const char *fmt, ...);
#endif
MPCCI_EXPORT int  (*umpcci_vprintf) (int level, const char *fmt, va_list ap);


#ifdef __GNUC__
/* let GCC check the format string and arguments */
MPCCI_EXPORT int   umpcci_msg_prefix(const char *fmt, ...) __attribute__ ((format(printf,1,2)));
#else
MPCCI_EXPORT int   umpcci_msg_prefix(const char *fmt, ...);
#endif

MPCCI_EXPORT int   umpcci_msg_level (int level);
MPCCI_EXPORT void  umpcci_msg_functs
                   (
                        void (*printFullMessage)(const char *str, int len),
                        void (*printFullWarning)(const char *str, int len),
                        void (*printFullFatal  )(const char *str, int len),
                        void (*printLineMessage)(const char *str, int len),
                        void (*printLineWarning)(const char *str, int len),
                        void (*printLineFatal  )(const char *str, int len),
                        void (*atExitHandler)   (void)
                   );

#ifdef __cplusplus
}
#endif


#define MPCCI_MSG_LEVEL_FATAL    -1
#define MPCCI_MSG_LEVEL_WARNING   0
#define MPCCI_MSG_LEVEL_INFO      1
#define MPCCI_MSG_LEVEL_ACTION    2
#define MPCCI_MSG_LEVEL_DEBUG     3


/*
 * basic informational messages passed at debug level >= 1
 */
#define MPCCI_MSG_INFO0(f)                               umpcci_printf(MPCCI_MSG_LEVEL_INFO,f)
#define MPCCI_MSG_INFO1(f,a1)                            umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1)
#define MPCCI_MSG_INFO2(f,a1,a2)                         umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2)
#define MPCCI_MSG_INFO3(f,a1,a2,a3)                      umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3)
#define MPCCI_MSG_INFO4(f,a1,a2,a3,a4)                   umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4)
#define MPCCI_MSG_INFO5(f,a1,a2,a3,a4,a5)                umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4,a5)
#define MPCCI_MSG_INFO6(f,a1,a2,a3,a4,a5,a6)             umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4,a5,a6)
#define MPCCI_MSG_INFO7(f,a1,a2,a3,a4,a5,a6,a7)          umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4,a5,a6,a7)
#define MPCCI_MSG_INFO8(f,a1,a2,a3,a4,a5,a6,a7,a8)       umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4,a5,a6,a7,a8)
#define MPCCI_MSG_INFO9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9)    umpcci_printf(MPCCI_MSG_LEVEL_INFO,f,a1,a2,a3,a4,a5,a6,a7,a8,a9)

/*
 * informational activity messages passed at debug level >= 2
 */
#define MPCCI_MSG_ACTION0(f)                             umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f)
#define MPCCI_MSG_ACTION1(f,a1)                          umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1)
#define MPCCI_MSG_ACTION2(f,a1,a2)                       umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2)
#define MPCCI_MSG_ACTION3(f,a1,a2,a3)                    umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3)
#define MPCCI_MSG_ACTION4(f,a1,a2,a3,a4)                 umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4)
#define MPCCI_MSG_ACTION5(f,a1,a2,a3,a4,a5)              umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4,a5)
#define MPCCI_MSG_ACTION6(f,a1,a2,a3,a4,a5,a6)           umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4,a5,a6)
#define MPCCI_MSG_ACTION7(f,a1,a2,a3,a4,a5,a6,a7)        umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4,a5,a6,a7)
#define MPCCI_MSG_ACTION8(f,a1,a2,a3,a4,a5,a6,a7,a8)     umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4,a5,a6,a7,a8)
#define MPCCI_MSG_ACTION9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9)  umpcci_printf(MPCCI_MSG_LEVEL_ACTION,f,a1,a2,a3,a4,a5,a6,a7,a8,a9)

/*
 * detailed debug messages passed at debug level >= 3
 */
#define MPCCI_MSG_DEBUG0(f)                              umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f)
#define MPCCI_MSG_DEBUG1(f,a1)                           umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1)
#define MPCCI_MSG_DEBUG2(f,a1,a2)                        umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2)
#define MPCCI_MSG_DEBUG3(f,a1,a2,a3)                     umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3)
#define MPCCI_MSG_DEBUG4(f,a1,a2,a3,a4)                  umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4)
#define MPCCI_MSG_DEBUG5(f,a1,a2,a3,a4,a5)               umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4,a5)
#define MPCCI_MSG_DEBUG6(f,a1,a2,a3,a4,a5,a6)            umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4,a5,a6)
#define MPCCI_MSG_DEBUG7(f,a1,a2,a3,a4,a5,a6,a7)         umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4,a5,a6,a7)
#define MPCCI_MSG_DEBUG8(f,a1,a2,a3,a4,a5,a6,a7,a8)      umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4,a5,a6,a7,a8)
#define MPCCI_MSG_DEBUG9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9)   umpcci_printf(MPCCI_MSG_LEVEL_DEBUG,f,a1,a2,a3,a4,a5,a6,a7,a8,a9)

/*
 * warning messages always printed
 */
#define MPCCI_MSG_WARNING0(f)                            umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f)
#define MPCCI_MSG_WARNING1(f,a1)                         umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1)
#define MPCCI_MSG_WARNING2(f,a1,a2)                      umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2)
#define MPCCI_MSG_WARNING3(f,a1,a2,a3)                   umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3)
#define MPCCI_MSG_WARNING4(f,a1,a2,a3,a4)                umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4)
#define MPCCI_MSG_WARNING5(f,a1,a2,a3,a4,a5)             umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4,a5)
#define MPCCI_MSG_WARNING6(f,a1,a2,a3,a4,a5,a6)          umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4,a5,a6)
#define MPCCI_MSG_WARNING7(f,a1,a2,a3,a4,a5,a6,a7)       umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4,a5,a6,a7)
#define MPCCI_MSG_WARNING8(f,a1,a2,a3,a4,a5,a6,a7,a8)    umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4,a5,a6,a7,a8)
#define MPCCI_MSG_WARNING9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9) umpcci_printf(MPCCI_MSG_LEVEL_WARNING,f,a1,a2,a3,a4,a5,a6,a7,a8,a9)

/*
 * warning messages and exit
 */
#define MPCCI_MSG_FATAL0(f)                              umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f)
#define MPCCI_MSG_FATAL1(f,a1)                           umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1)
#define MPCCI_MSG_FATAL2(f,a1,a2)                        umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2)
#define MPCCI_MSG_FATAL3(f,a1,a2,a3)                     umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3)
#define MPCCI_MSG_FATAL4(f,a1,a2,a3,a4)                  umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4)
#define MPCCI_MSG_FATAL5(f,a1,a2,a3,a4,a5)               umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4,a5)
#define MPCCI_MSG_FATAL6(f,a1,a2,a3,a4,a5,a6)            umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4,a5,a6)
#define MPCCI_MSG_FATAL7(f,a1,a2,a3,a4,a5,a6,a7)         umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4,a5,a6,a7)
#define MPCCI_MSG_FATAL8(f,a1,a2,a3,a4,a5,a6,a7,a8)      umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4,a5,a6,a7,a8)
#define MPCCI_MSG_FATAL9(f,a1,a2,a3,a4,a5,a6,a7,a8,a9)   umpcci_printf(MPCCI_MSG_LEVEL_FATAL,f,a1,a2,a3,a4,a5,a6,a7,a8,a9)


#define MPCCI_MSG_ASSERT(cond)  if (!(cond)) MPCCI_MSG_FATAL3("Assertion failed (%s): %s(%d)\n",#cond,__FILE__,__LINE__)

#endif
