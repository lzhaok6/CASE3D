#ifndef __MPCCI4_MONITOR_HEADER_INCLUDED
#define __MPCCI4_MONITOR_HEADER_INCLUDED
/* mpcci_monitor.h
 *
 *****************************************************************************************
 *
 * Purpose:
 *    Special defines for the MpCCI Monitor
 *
 * Author:
 *    Carsten Dehning <carsten.dehning@scai.fhg.de>
 *
 * Reviews/changes:
 *    2012/Oct: Carsten Dehning, Initial release
 *    $Id: mpcci_monitor.h 878 2014-10-24 13:07:32Z pbayrasy $
 *
 *****************************************************************************************
 */

/*
 * Bits in monitor->flags
 */
#define MONITOR_FLAG_ORPHANS     0x00000001  /* send node orphan level OLEVEL_... */
#define MONITOR_FLAG_NODEIDS     0x00000002  /* send node ID's */
#define MONITOR_FLAG_ELEMIDS     0x00000004  /* send element ID's */
#define MONITOR_FLAG_ELSIZES     0x00000008  /* send element sizes */
#define MONITOR_FLAG_GLOBALS     0x00000010  /* send received global variables */
#define MONITOR_FLAG_SLAVES      0x00000020  /* send slave node marks (0/1) */
#define MONITOR_FLAG_DOMAINS     0x00000040  /* send node domains */
#define MONITOR_FLAG_PEAKS       0x00000080  /* send quantity peak & mean values */
#define MONITOR_FLAG_INORM       0x00000100  /* send quantity iteration norm peak & mean values */
#define MONITOR_FLAG_NORMALS     0x00000200  /* send face normal vectors */
#define MONITOR_FLAG_ABSDIFF     0x00000400  /* send absolute quantity differences to the monitor */
#define MONITOR_FLAG_RELDIFF     0x00000800  /* send %% quantity differences to the monitor */
#define MONITOR_FLAG_PERIODIC    0x00001000  /* send periodic part and quantities */
#define MONITOR_FLAG_RELAX       0x00002000  /* send Aitken relaxation value */

#define MONITOR_FLAG_UPDATED     0x80000000  /* the flags have been updated by the visualizer GUI */

#endif
