/*
 * -- SuperLU MT routine (version 2.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 * These macros define which machine will be used.
 *
 * Modified:  March 20, 2013  version 2.1
 */

#ifndef __SUPERLU_MACHINES /* allow multiple inclusions */
#define __SUPERLU_MACHINES

#define PTHREAD         0
#define OPENMP		1
#define SUN             6
#define SGI	        2
#define ORIGIN	        3
#define DEC	        4
#define CRAY_PVP	5

#ifdef __PTHREAD
#define MACH PTHREAD
#endif

#ifdef __OPENMP
#define MACH OPENMP
#endif

#ifdef __SOLARIS
#define MACH SUN 
#endif

#ifdef __SGI
#define MACH SGI 
#endif

#ifdef __ORIGIN
#define MACH ORIGIN 
#endif

#ifdef __DEC
#define MACH DEC 
#endif

#ifdef __CRAY
#define MACH CRAY_PVP 
#endif


#endif /* __SUPERLU_MACHINES */
