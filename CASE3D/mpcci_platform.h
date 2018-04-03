#ifndef __MPCCI4_PLATFORM_HEADER_INCLUDED
#define __MPCCI4_PLATFORM_HEADER_INCLUDED
/*****************************************************************************************
 * mpcci_platform.h
 *
 * Include header file to set platform dependent defines
 *
 * $Id: mpcci_platform.h 104 2013-08-16 11:56:59Z pbayrasy $
 *
 ****************************************************************************************/

#if defined(_WINNT) || defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)

   #ifdef MPCCI_COMPILE_LIBMPCCI

      /* we are compiling the client library */
      #define MPCCI_EXPORT extern __declspec(dllexport)

   #else

      /* we are using the client library */
      #ifdef MPCCI_CLIENT_DLL
         #define MPCCI_EXPORT extern __declspec(dllimport)
      #else
         #define MPCCI_EXPORT extern
      #endif

   #endif

#else

   #define MPCCI_EXPORT    extern

#endif

#endif
