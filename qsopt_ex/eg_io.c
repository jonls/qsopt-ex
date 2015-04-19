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
/** @file
 * @ingroup EGio */
/** @addtogroup EGio */
/** @{ */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#ifdef HAVE_LIBBZ2
# include <bzlib.h>
#endif

#ifdef HAVE_LIBZ
# include <zlib.h>
#endif

#include "logging-private.h"

#include "eg_io.h"

/* file-types: */
#define EGIO_PLAIN 0
#define EGIO_ZLIB  1
#define EGIO_BZLIB 2

/* ========================================================================= */
/** @brief Given a string 'input' this function uses EGioParse to separate
 * up to N words in it, we assume that argc is an array of pointers to strings
 * of size N, and note that the input array will be changed. */
void EGioNParse (char *input,
								 int max_argc,
								 const char *delim,
								 const char *comment,
								 int *argc,
								 char **argv)
{
	const size_t clen = strlen(comment);
	const size_t dlen = strlen(delim);
	char __EGiobuff[256] = 
		"20000000000000000000000000000000" /* 000-031 */
		"11111111111111111111111111111111" /* 032-063 */
		"11111111111111111111111111111111" /* 064-095 */
		"11111111111111111111111111111110" /* 096-127 */
		"00000000000000000000000000000000" /* 128-159 */
		"00000000000000000000000000000000" /* 160-191 */
		"00000000000000000000000000000000" /* 192-223 */
		"00000000000000000000000000000000";/* 224-255 */
	int i, EGiostat = 0, cc;
	char *cur;
	/* set convertion table */
	for( i = clen ; i-- ; )
	{
		cc = comment[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '2';
	}
	for( i = dlen ; i-- ; )
	{
		cc = delim[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '0';
	}
	/* now parse */
	*argc = EGiostat = 0;
	cur = input;
	while(cur && (*argc) < max_argc)
	{
		cc = __EGiobuff[(int)(*cur)];
		switch(cc)
		{
			case '1':
				if(!EGiostat) argv[(*argc)++] = cur;
				EGiostat = 1;
				cur++;
				break;
			case '0':
				EGiostat = 0;
				*(cur++) = '\0';
				break;
			case '2':
				*cur = '\0';
				cur = 0;
				break;
			default:
				EXIT(1,"Imposible, cc=%c, cur = %c, pos %zd",cc, *cur, cur-input);
		}
	}
	/* cleanup */
	for( i = clen ; i-- ; )
	{
		cc = comment[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '1';
	}
	for( i = dlen ; i-- ; )
	{
		cc = delim[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '1';
	}
	return;
}

/* given two *pointers 'next' and 'current' and a constant set  
 * of strings (parse delimiters), it store in 'next the next 
 * meaningfull string, and in current the rest of the secuence, 
 * the idea is to iterate over 'next' while it is not true; 
 * you have to store the original pointer to the string stream 
 * elsewere; also, we assume that the original stream is 
 * terminated with '\0'; also, it will discaard any sub-string 
 * that start with #, that is inteded for discard comments.
 * NOTE: this function WILL change the original string!!!!!
 * There is no guarantee on how it will be changed. */
void EGioParse (char **next,
								char **current,
								const char *delim,
								const char *comment)
{
	/* local variables */
	const size_t clen = strlen(comment);
	const size_t dlen = strlen(delim);
	char __EGiobuff[256] = 
		"20000000000000000000000000000000" /* 000-031 */
		"11111111111111111111111111111111" /* 032-063 */
		"11111111111111111111111111111111" /* 064-095 */
		"11111111111111111111111111111110" /* 096-127 */
		"00000000000000000000000000000000" /* 128-159 */
		"00000000000000000000000000000000" /* 160-191 */
		"00000000000000000000000000000000" /* 192-223 */
		"00000000000000000000000000000000";/* 224-255 */
	int i, EGiostat = 0, cc;
	char *cur;
	/* set convertion table */
	for( i = clen ; i-- ; )
	{
		cc = comment[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '2';
	}
	for( i = dlen ; i-- ; )
	{
		cc = delim[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '0';
	}
	/* parse */
	*next = 0;
	cur = *current;
	while(cur)
	{
		cc = __EGiobuff[(int)(*cur)];
		switch(cc)
		{
			case '1':
				if(!EGiostat) *next = cur;
				EGiostat = 1;
				cur++;
				break;
			case '0':
				*(cur++) = '\0';
				if(EGiostat)
				{
					*current = cur;
					cur = 0;
				}
				EGiostat = 0;
				break;
			case '2':
				*cur = '\0';
				*current = cur = 0;
				break;
			default:
				EXIT(1,"Imposible, cc=%c, cur = %c, pos %zd",cc, *cur, cur-*current);
		}
	}
	/* cleanup */
	for( i = clen ; i-- ; )
	{
		cc = comment[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '1';
	}
	for( i = dlen ; i-- ; )
	{
		cc = delim[i];
		if(cc>31 && cc<127) __EGiobuff[cc] = '1';
	}
	return;
}

/* this discard all lines starting with comments and stores the next 
 * leading line in current and the next token in next, we assume that next 
 * does not contain data, and that current store the remainings of the
 * current line */
void EGioDisCom (char **next,
								 char **current,
								 const char *delim,
								 const char *comment,
								 char *store,
								 unsigned int storeSize,
								 FILE * in)
{
	/* local variables */
	int status = 1;

	/* if no current line we read it from the file */
	if (!(*current))
	{
		status = (store == fgets (store, (int) storeSize, in));
		*current = store;
	}
	/* we process the current line, and while the line does
	 * not have a token we look the next line */
	EGioParse (next, current, delim, comment);
	while (!(*next) && status)
	{
		status = (store == fgets (store, (int) storeSize, in));
		*current = store;
		EGioParse (next, current, delim, comment);
	}
	/* ending */
	return;
}

/* ========================================================================= */
void EGdisplayString (void *str,
											FILE * file)
{
	fprintf (file, "%s", (char *) str);
}

/* ========================================================================= */
int EGioReadLine(char*const str,size_t const max_len, FILE*file)
{
	char *rc = fgets(str,max_len,file);
	int rval = !rc;
	FTESTG(rval, CLEANUP,"Nothing to be read");
	CLEANUP:
	return rval;
	/*
	int c=0;
	size_t len = max_len ;
	while(((c=getc(file))!=EOF) && (c!='\n') && --len)
	{
		str[max_len-1-len] = (char)c;
	}
	str[max_len-len] = '\0';
	FTEST((max_len - len)==0 && c!= '\n',"Nothing to be read"); 
	return 0;
	*/
}
/* ========================================================================= */
int EGioReadNamedDblPlusParam(const int argc,
															char**argv,
															const char*const name,
															double*const param,
															int*const has_param,
															int*const rval)
{
	const size_t len = strlen(name)+5;
	int fail=0;
	*rval = 0;
	if(argc<1) return 0;
	if(strncmp(argv[0],name,len)==0)
		{
			FTESTG((fail=(argc!=2)),CLEANUP,"%s has not 2 tokens",name);
			FTESTG((fail=(*has_param)),CLEANUP,"%s keyword repeated",name);
			*param = strtod(argv[1],0);
			FTESTG((fail=((*param)<0)),CLEANUP,"%s should be positive,"
						"is %lf",name,*param);
			*has_param = 1;
			return 1;
		}
	else
		return 0;
	CLEANUP:
	*rval=1;
	return fail;
}
/* ========================================================================= */
int EGioReadNamedIntPlusParam(const int argc,
															char**argv,
															const char*const name,
															int*const param,
															int*const has_param,
															int*const rval)
{
	const size_t len = strlen(name)+5;
	int fail=0;
	*rval = 0;
	if(argc<1) return 0;
	if(strncmp(argv[0],name,len)==0)
		{
			FTESTG((fail=(argc!=2)),CLEANUP,"%s has not 2 tokens",name);
			FTESTG((fail=(*has_param)),CLEANUP,"%s keyword repeated",name);
			*param = atoi(argv[1]);
			FTESTG((fail=((*param)<1)),CLEANUP,"%s should be positive,"
						"is %d",name,*param);
			*has_param = 1;
			return 1;
		}
	else
		return 0;
	CLEANUP:
	*rval=1;
	return fail;
}
/* ========================================================================= */
int EGioReadNamedIntNNParam(const int argc,
														char**argv,
														const char*const name,
														int*const param,
														int*const has_param,
														int*const rval)
{
	const size_t len = strlen(name)+5;
	int fail=0;
	*rval = 0;
	if(argc<1) return 0;
	if(strncmp(argv[0],name,len)==0)
		{
			FTESTG((fail=(argc!=2)),CLEANUP,"%s has not 2 tokens",name);
			FTESTG((fail=(*has_param)),CLEANUP,"%s keyword repeated",name);
			*param = atoi(argv[1]);
			FTESTG((fail=((*param)<0)),CLEANUP,"%s should be positive,"
						"is %d",name,*param);
			*has_param = 1;
			return 1;
		}
	else
		return 0;
	CLEANUP:
	*rval=1;
	return fail;
}
/* ========================================================================= */
int EGioReadNamedStringParam(const int argc,
															char**argv,
															const char*const name,
															char**const param,
															int*const has_param,
															int*const rval)
{
	const size_t len = strlen(name)+5;
	int fail=0;
	*rval = 0;
	if(argc<1) return 0;
	if(strncmp(argv[0],name,len)==0)
		{
			FTESTG((fail=(argc!=2)),CLEANUP,"%s has not 2 tokens",name);
			FTESTG((fail=(*has_param)),CLEANUP,"%s keyword repeated",name);
			if(param) *param = strdup(argv[1]);
			*has_param = 1;
			return 1;
		}
	else
		return 0;
	CLEANUP:
	*rval=1;
	return fail;
}
/* ========================================================================= */
int EGioReadNamedParam( const int argc,
												char**argv,
												const char*const name,
												int*const has_param,
												int*const rval)
{
	const size_t len = strlen(name)+5;
	int fail=0;
	*rval = 0;
	if(argc<1) return 0;
	if(strncmp(argv[0],name,len)==0)
		{
			FTESTG((fail=(argc!=1)),CLEANUP,"%s has not 1 token",name);
			FTESTG((fail=(*has_param)),CLEANUP,"%s keyword repeated",name);
			*has_param = 1;
			return 1;
		}
	else
		return 0;
	CLEANUP:
	*rval=1;
	return fail;
}
/* ========================================================================= */

int EGioReadIntParam( const int argc,
											char**argv,
											int*const param,
											int*const rval)
{
	*rval = 0;
	if(argc<1) return 0;
	FTESTG((argc!=1),CLEANUP,"line has not 1 token");
	*param = atoi(argv[0]);
	return 0;
	CLEANUP:
	*rval=1;
	return 1;
}
/* ========================================================================= */
#define EGio_BUFSIZE 4096
/* ========================================================================= */
struct EGioFile_st {int type; void*file;};
/* ========================================================================= */
int EGioWrite(EGioFile_t*file,const char*const string)
{
	char buf[EGio_BUFSIZE];
	int len;
	buf[EGio_BUFSIZE-1] = 0;
	snprintf(buf,EGio_BUFSIZE,"%s",string);
	len = strlen(buf);
	if(len<=0 || len >= EGio_BUFSIZE || buf[EGio_BUFSIZE-1]!=0) return 0;
	switch(file->type)
	{
		case EGIO_PLAIN:
			return fwrite(buf, (size_t)1, (size_t)len, (FILE*)(file->file));
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			return gzwrite((gzFile)(file->file),buf,(unsigned)len);
#else
			QSlog("no zlib support");
			return 0;
#endif
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			return BZ2_bzwrite((BZFILE*)(file->file),buf,len);
#else
			QSlog("no bzip2 support");
			return 0;
#endif
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			return 0;
	}
}
/* ========================================================================= */
int EGioPrintf(EGioFile_t*file,const char* format, ...)
{
	char buf[EGio_BUFSIZE];
	va_list va;
	buf[EGio_BUFSIZE-1]=0;
	va_start(va,format);
	vsnprintf(buf,EGio_BUFSIZE,format,va);
	va_end(va);
	return EGioWrite(file,buf);
}
/* ========================================================================= */
EGioFile_t* EGioOpenFILE(FILE*ifile)
{
	EGioFile_t* file= (EGioFile_t*)malloc(sizeof(EGioFile_t));
	if(!file) goto CLEANUP;
	file->type = EGIO_PLAIN;
	file->file = ifile;
	if(!file->file)
	{
		free(file);
		file = 0;
	}
	CLEANUP:
	return file;
}
/* ========================================================================= */
EGioFile_t* EGioOpen(const char *path, const char *mode)
{
	char lmode[8];
	int len = strlen(path);
	EGioFile_t* file= (EGioFile_t*)malloc(sizeof(EGioFile_t));
	if(!file) goto CLEANUP;
	if(len>3 && path[len-3] == '.' && path[len-2] == 'g' && path[len-1] == 'z')
	{
		file->type = EGIO_ZLIB;
		if(index(mode,'b')) snprintf(lmode,7,"%s",mode);
		else snprintf(lmode,7,"%s9b",mode);
	}
	else if(len>4 && path[len-4] == '.' && path[len-3] == 'b' && path[len-2] == 'z' && path[len-1] == '2')
	{
		file->type = EGIO_BZLIB;
		if(index(mode,'b')) snprintf(lmode,7,"%s",mode);
		else snprintf(lmode,7,"%sb",mode);
	}
	else 
	{
		file->type = EGIO_PLAIN;
		snprintf(lmode,7,"%s",mode);
	}
	switch(file->type)
	{
		case EGIO_PLAIN:
			file->file = fopen(path,lmode);
			break;
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			file->file = gzopen(path,lmode);
#else
			QSlog("no zlib support");
			file->file = 0;
#endif
			break;
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			file->file = BZ2_bzopen(path,lmode);
#else
			QSlog("no bzip2 support");
			file->file = 0;
#endif
			break;
		default:
			QSlog("UNKNOWN FILE TYPE %d",file->type);
			file->file = 0;
			break;
	}
	if(!file->file)
	{
		const int __EGserrno = errno;
		free(file);
		errno = __EGserrno;
		return NULL;
	}
	CLEANUP:
	return file;
}
/* ========================================================================= */
int EGioClose(EGioFile_t*file)
{
	int rval = 0;
	if(!file) return 0;
	switch(file->type)
	{
		case EGIO_PLAIN:
			if( (((FILE*)(file->file)) != stdin) &&
					(((FILE*)(file->file)) != stdout) &&
					(((FILE*)(file->file)) != stderr))
				rval = fclose((FILE*)(file->file));
			break;
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			rval = gzclose((gzFile)(file->file));
#else
			QSlog("no zlib support");
			rval = EOF;
#endif
			break;
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			BZ2_bzerror((BZFILE*)(file->file),&rval);
			BZ2_bzclose((BZFILE*)(file->file));
#else
			QSlog("no bzip2 support");
			rval =  EOF;
#endif
			break;
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			rval = EOF;
			break;
	}
	free(file);
	return rval;
}
/* ========================================================================= */
int EGioFlush(EGioFile_t*file)
{
	switch(file->type)
	{
		case EGIO_PLAIN:
			return fflush((FILE*)(file->file));
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			/*return gzflush((gzFile)(file->file),Z_FINISH);*/
			return 0;
#else
			QSlog("no zlib support");
			return EOF;
#endif
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			return 0;
#else
			QSlog("no bzip2 support");
			return EOF;
#endif
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			return EOF;
	}
}
/* ========================================================================= */
char* EGioGets(char*buf, int len, EGioFile_t*file)
{
	char*b = buf;
	switch(file->type)
	{
		case EGIO_PLAIN:
			return fgets(buf,len,(FILE*)(file->file));
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			return gzgets((gzFile)(file->file),buf,len);
#else
			QSlog("no zlib support");
			return NULL;
#endif
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			if(buf == 0 || len <=0 ) return NULL;
			while( --len > 0 && BZ2_bzread(((BZFILE*)(file->file)), buf, 1) == 1 && *buf++ != '\n') ;
			*buf = '\0';
			return b == buf && len >0 ? NULL : b ;
#else
			QSlog("no bzip2 support");
			return NULL;
#endif
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			return NULL;
	}
}
/* ========================================================================= */
int EGioEof(const EGioFile_t*const file)
{
	int err;
	switch(file->type)
	{
		case EGIO_PLAIN:
			return feof((FILE*)(file->file));
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			return gzeof((gzFile)(file->file));
#else
			QSlog("no zlib support");
			return 1;
#endif
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			BZ2_bzerror(((BZFILE*)(file->file)),&err);
			return err == BZ_STREAM_END;
#else
			QSlog("no bzip2 support");
			return 1;
#endif
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			return 1;
	}
}
/* ========================================================================= */
int EGioError(const EGioFile_t*const file)
{
	int errnum;
	switch(file->type)
	{
		case EGIO_PLAIN:
			return ferror((FILE*)(file->file));
		case EGIO_ZLIB:
#ifdef HAVE_LIBZ
			gzerror((gzFile)(file->file),&errnum);
			return errnum;
#else
			QSlog("no zlib support");
			return 1;
#endif
		case EGIO_BZLIB:
#ifdef HAVE_LIBBZ2
			BZ2_bzerror(((BZFILE*)(file->file)),&errnum);
			return errnum;
#else
			QSlog("no bzip2 support");
			return 1;
#endif
		default:
			QSlog("UNKNOWN FILE TYPE %d", file->type);
			return 1;
	}
}/* ========================================================================= */
/** @} */
