#ifndef __MPCCI4_HEADER_INCLUDED
#define __MPCCI4_HEADER_INCLUDED
/*****************************************************
 *  User include file for programs coupled via MpCCI
 *  $Id: mpcci.h 104 2013-08-16 11:56:59Z pbayrasy $
 *****************************************************/
#ifndef MPCCI_CCM_VERSION
   #define MPCCI_CCM_VERSION 430
#endif


#if MPCCI_CCM_VERSION == 430

   #include <stddef.h>
   #include "mpcci_platform.h"
   #include "mpcci_quantities.h"
   #include "mpcci_coreapi.h"
   #include "mpcci_stdioapi.h"

#else

   #error MPCCI_CCM_VERSION must be 430

#endif

#endif

