/****************************************************************************/
/*                                                                          */
/*  This file is part of QSopt_ex.                                          */
/*                                                                          */
/*  (c) Copyright 2006 by David Applegate, William Cook, Sanjeeb Dash,      */
/*  and Daniel Espinoza                                                     */
/*                                                                          */
/*  Sanjeeb Dash ownership of copyright in QSopt_ex is derived from his     */
/*  copyright in QSopt.                                                     */
/*                                                                          */
/*  This code may be used under the terms of the GNU General Public License */
/*  (Version 2.1 or later) as published by the Free Software Foundation.    */
/*                                                                          */
/*  Alternatively, use is granted for research purposes only.               */
/*                                                                          */
/*  It is your choice of which of these two licenses you are operating      */
/*  under.                                                                  */
/*                                                                          */
/*  We make no guarantees about the correctness or usefulness of this code. */
/*                                                                          */
/****************************************************************************/

/* RCSINFO $Id: machdefs.h,v 1.2 2003/11/05 16:47:22 meven Exp $ */
/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __MACHDEFS_H
#define __MACHDEFS_H
#include "qs_config.h"

#ifdef HAVE_STDIO_H
# include <stdio.h>
#endif

#ifdef HAVE_STDARG_H
# include <stdarg.h>
#endif
#ifdef HAVE_STDLIB_H
# include <stdlib.h>
#endif
#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_STRING_H
# include <string.h>
#endif
#ifdef HAVE_STRINGS_H
#  include <strings.h>
#endif
#ifdef HAVE_ERRNO_H
# include <errno.h>
#endif
#ifdef TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  ifdef HAVE_TIME_H
#   include <time.h>
#  endif
# endif
#endif
#ifdef HAVE_STDDEF_H
# include <stddef.h>
#endif
#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif
#ifdef HAVE_MALLOC_H
# include <stdlib.h>
#endif
#ifdef HAVE_SYS_TYPES_H
# include <sys/types.h>
#endif
#ifdef HAVE_SYS_STAT_H
# include <sys/stat.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
# include <sys/resource.h>
#endif
#ifdef HAVE_FCNTL_H
# include <fcntl.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
# include <sys/socket.h>
#endif
#ifdef HAVE_NETDB_H
# include <netdb.h>
#endif
#ifdef HAVE_NETINET_IN_H
# include <netinet/in.h>
#endif

#ifdef HAVE_SOCKET
#define ILL_NETREADY
#endif

#ifdef USEDMALLOC
#include "dmalloc.h"
#endif

#ifdef ILL_ATTRIBUTE
#define ILL_UNUSED __attribute__ ((unused))
#else
#define ILL_UNUSED
#endif

#ifdef ILL_PROTO_PRINTF
/* assume that if you're missing printf, you're missing a bunch */
extern int printf (
	const char *,
	...),
  fprintf (
	FILE *,
	const char *,
	...),
  fflush (
	FILE *),
  scanf (
	const char *,
	...),
  sscanf (
	const char *,
	const char *,
	...),
  fscanf (
	FILE *,
	const char *,
	...),
  fclose (
	FILE *),
  ungetc (
	int,
	FILE *),
  _filbuf (
	FILE *),
  time (
	int *);

#ifdef ILL_NETREADY
extern int socket (
	int,
	int,
	int),
  connect (
	int,
	const struct sockaddr *,
	int),
  accept (
	int,
	struct sockaddr *,
	int *),
  bind (
	int,
	const struct sockaddr *,
	int),
  listen (
	int,
	int);
#endif
extern void *memset (
	void *,
	int,
	size_t),
  perror (
	const char *);
#endif

#ifdef ILL_PROTO_RENAME
extern int rename (
	const char *,
	const char *);
#endif

#ifdef ILL_PROTO_GETHOSTNAME
extern int gethostname (
	char *,
	int);
#endif

#ifdef ILL_PROTO_GETRUSAGE
extern int getrusage (
	int,
	struct rusage *);
#endif

#ifdef ILL_PROTO___VFORK
extern pid_t __vfork (
	void);
#endif

#ifndef NULL
#define NULL (0)
#endif

//#ifndef INT_MAX
//#define INT_MAX ((int) (~(((unsigned) 1) << ((8*sizeof(int))-1))))
//#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#endif /* __MACHDEFS_H */
