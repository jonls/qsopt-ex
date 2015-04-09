/* ========================================================================= */
/* EGlib "Efficient General Library" provides some basic structures and
 * algorithms commons in many optimization algorithms.
 *
 * Copyright (C) 2005 Daniel Espinoza and Marcos Goycoolea.
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA 
 * */
/* ========================================================================= */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>
#include <setjmp.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <unistd.h>

#include "eg_macros.h"

#include "qs_config.h"
#include "logging.h"

/** @file
 * @brief implementation of some macros.
 *
 * @version 0.9.2 
 * @par History:
 * 	-2006-09-28
 * 						- First implementation
 *	-2007-12-06
 *						- Add versioning information in header
 * @ingroup EGmacros */
/** @addtogroup EGmacros */
/** @{ */
/* ========================================================================= */
void EGlib_info(void)
{
	int rval;
	struct utsname uts;
	rval = uname(&uts);
	if(rval)
	{
		fprintf(stderr,"Can't get host info\n");
	}
	else
	{
		fprintf(stderr,"Host: %s, %s %s %s\nCurrent process id: %d\n",uts.nodename, uts.sysname, uts.release, uts.machine ,(int)getpid());
	}
}
/** @} */
/* ========================================================================= */
jmp_buf __EGljmp;
/* ========================================================================= */
void EGsighandler(int s)
{
	switch(s)
	{
		case SIGXCPU:
			/* time is over */ 
			fprintf(stderr,"\nTIME_LIMIT_REACHED (ending now)\n");
			longjmp(__EGljmp,s);
			break;
		case SIGINT:
		case SIGTERM:
		case SIGTSTP:
			/* time is over */ 
			fprintf(stderr,"\nUSER_INTERRUPT (ending now)\n");
			longjmp(__EGljmp,s);
			break;
		case SIGSEGV:
			/* Memory is over or segmentation fault */ 
			fprintf(stderr,"\nMEMORY_LIMIT_REACHED (ending now)\n");
			longjmp(__EGljmp,s);
			break;
		case SIGABRT:
			/* something is sending an abort code.... stop execution, report the
			 * signal, do destruction and report and end */
			fprintf(stderr,"\nSIGABRT received (ending now)\n");
			longjmp(__EGljmp,s);
			break;
		default:
			fprintf(stderr,"\nUnkown signal %d\n",s);
			fprintf(stderr,"Ending with status %d\n",s);
			exit(s);
	}
	fflush(stdout);
	fflush(stderr);
}
/* ========================================================================= */
void __EGsigSetSignal(void)
{
	signal(SIGXCPU,EGsighandler);
	signal(SIGINT,EGsighandler);
	signal(SIGSEGV,EGsighandler);
	signal(SIGABRT,EGsighandler);
}
/* ========================================================================= */
void EGsetLimits(double max_rtime, unsigned long memlimit)
{
	struct rlimit mlim;
	WARNIF(getrlimit(RLIMIT_CPU,&mlim));
	MESSAGE(0, "Cur rtime limit %ld, trying to set to %lg", mlim.rlim_cur, max_rtime);
	if(max_rtime > mlim.rlim_max) max_rtime = (double)mlim.rlim_max;
	mlim.rlim_cur = (rlim_t)max_rtime;
	WARNIF(setrlimit(RLIMIT_CPU,&mlim));
	MESSAGE(0, "New rtime limit %ld (%.3lg)", mlim.rlim_cur, max_rtime);
	WARNIF(getrlimit(RLIMIT_DATA,&mlim));
	MESSAGE(0, "Cur data limit %ld,%ld (soft,hard)", mlim.rlim_cur, 
					mlim.rlim_max);
	mlim.rlim_cur = memlimit;				
	WARNIF( setrlimit(RLIMIT_DATA,&mlim));
	WARNIF( getrlimit(RLIMIT_DATA,&mlim));
	MESSAGE(0, "New data limit %ld,%ld (soft,hard)", mlim.rlim_cur, 
					mlim.rlim_max);
	WARNIF( getrlimit(RLIMIT_AS,&mlim));
	MESSAGE(0, "Cur address space limit %ld,%ld (soft,hard)", 
					mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = memlimit;	
	WARNIF( setrlimit(RLIMIT_AS,&mlim));				
	WARNIF( getrlimit(RLIMIT_AS,&mlim));
	MESSAGE(0, "New address space limit %ld,%ld (soft,hard)", 
					mlim.rlim_cur, mlim.rlim_max);
	mlim.rlim_cur = 0;
	WARNIF( setrlimit(RLIMIT_CORE,&mlim));				
	WARNIF( getrlimit(RLIMIT_CORE,&mlim));
	MESSAGE(0, "New core dump space limit %ld,%ld (soft,hard)", 
					mlim.rlim_cur, mlim.rlim_max);
	return;
}
/* ========================================================================= */
/* end of eg_macros.c */
