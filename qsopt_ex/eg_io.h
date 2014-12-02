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
/** @defgroup EGio EGio
 * input/output utilities
 *
 * Version 0.0.2 2003-05-09 (Marcos)
 *
 * Added the function EGioNParse to get a more argc, argv feel.
 * Also: Changed EGioParse so that it ignores multiple sequential delimiters.
 * 
 * Version 0.0.1 2003-04-11
 * - 2004-08-17
 * 					-	Add EGdisplayString function.
 * - 2006-08-16
 * 					- Add EGioReadLine function.
 * - 2007-12-06
 * 					- Add EGioReadxxxParam functions
 * - 2009-10-19
 * 					- Add bzlib/zlib/plain text interface
 * */
/** @file
 * @ingroup EGio */
/** @addtogroup EGio */
/** @{
 * @example eg_ebtree.ex.c
 * This is an example of writing to a zlib-compresed file
 *
 * @example eg_min_cut.ex.c
 * This is an example of reading files, either plain, bz2 or gz
 * */
/* ========================================================================= */

#ifndef __EG_IO_H__
#define __EG_IO_H__

#include "eg_macros.h"

/* ========================================================================= */
/** @brief safe open function, it test that we can open the given file with the
 * given mode, if an error occurs, display it on screen and exit, otherwise,
 * return the open file stream.
 * @param __file file name to open (is a const char*)
 * @param __mode mode to use to open the file (is a const char*)
 * @return FILE* stream */
#define EGsfopen(__file,__mode) ({\
	const char*__EGsfile = (__file);\
	const char*__EGsmode = (__mode);\
	FILE*__EGsFILE = fopen(__EGsfile,__EGsmode);\
	if(!__EGsFILE)\
	{\
		const int __EGserrno = errno;\
		fprintf(stderr,"fopen() failed with error code %d, error\n%s",__EGserrno,strerror(__EGserrno));\
		MESSAGE(0,"Could not open %s with mode %s", __EGsfile, __EGsmode);\
		exit(__EGserrno);\
	}\
	__EGsFILE;})

/* ========================================================================= */
/**@brief type of functions for display, receives a void* to the structure to print,
 * and a *FILE where to output */
typedef void (*EGdisplay_f) (void *, FILE *);
#define EGnullDisplay ((EGdisplay_f)0)

/* ========================================================================= */
/** @brief type of functions for display, receives a void* to the structure to print, a
 * *FILE where to output, and a set of offsets for datas, the length of that
 * array must be know by the user. */
typedef void (*EGdisplayOS_f) (void *,
															 FILE *,
															 size_t *);

/* ========================================================================= */
/** @brief Given a string 'input' this function uses EGioParse to separate
 * up to N words in it, we assume that argc is an array of pointers to strings
 * of size N, and note that the input array will be changed. */
void EGioNParse (char *input,
								 int max_argc,
								 const char *delim,
								 const char *comment,
								 int *argc,
								 char **argv);

/* ========================================================================= */
/** @brief given two *pointers 'next' and 'current' and a constant set 
 * of strings (parse delimiters), it store in 'next the next 
 * meaningfull string, and in current the rest of the secuence, 
 * the idea is to iterate over 'next' while it is not true; 
 * you have to store the original pointer to the string stream 
 * elsewere; also, we assume that the original stream is 
 * terminated with '\0'; also, it will discaard any sub-string 
 * that start with #, that is inteded for discard comments */
void EGioParse (char **next,
								char **current,
								const char *delim,
								const char *coment);

/* ========================================================================= */
/** @brief read a line from an input stream.
 * @param str where to store the line.
 * @param max_len maximum allowed length.
 * @param file stream from where we read the input
 * @return zero on success, non-zero otherwise, errors are cast when we can not
 * read from the given file.
 * */
int EGioReadLine(char*const str,size_t const max_len, FILE*file);

/* ========================================================================= */
/** @brief read a named string parameter, this function checks that we have two
 * parameters, and check if we had previously readed the parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param name named of the parameter
 * @param param where to save the parameter (if null, don't save it)
 * @param has_param if we had the parameter before, should be one, if not,
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadNamedStringParam( const int argc,
															char**argv,
															const char*const name,
															char**const param,
															int*const has_param,
															int*const rval);
/* ========================================================================= */
/** @brief read a named, non-negative, integer parameter, this function checks that we have two
 * parameters, and check if we had previously readed the parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param name named of the parameter
 * @param param where to save the parameter (if null, don't save it)
 * @param has_param if we had the parameter before, should be one, if not,
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadNamedIntNNParam(const int argc,
															char**argv,
															const char*const name,
															int*const param,
															int*const has_param,
															int*const rval);
/* ========================================================================= */
/** @brief read a named, strictly positive, integer parameter, this function checks that we have two
 * parameters, and check if we had previously readed the parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param name named of the parameter
 * @param param where to save the parameter (if null, don't save it)
 * @param has_param if we had the parameter before, should be one, if not,
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadNamedIntPlusParam(const int argc,
															char**argv,
															const char*const name,
															int*const param,
															int*const has_param,
															int*const rval);
/* ========================================================================= */
/** @brief read a named positive double parameter, this function checks that we have two
 * parameters, and check if we had previously readed the parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param name named of the parameter
 * @param param where to save the parameter (if null, don't save it)
 * @param has_param if we had the parameter before, should be one, if not,
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadNamedDblPlusParam(const int argc,
															char**argv,
															const char*const name,
															double*const param,
															int*const has_param,
															int*const rval);
/* ========================================================================= */
/** @brief read a named parameter, this function checks that we have one
 * parameter, and check if we had previously readed the parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param name named of the parameter
 * @param has_param if we had the parameter before, should be one, if not,
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadNamedParam( const int argc,
												char**argv,
												const char*const name,
												int*const has_param,
												int*const rval);
/* ========================================================================= */
/** @brief read an integer parameter, this function checks that we have one
 * parameter, if an error
 * occurs, report it in the given rval 
 * @param argc number of tokens in the current line
 * @param argv array of strings of tokens
 * @param param named of the parameter
 * should be zero, if successfull, then it will be set to one.
 * @param rval return value if an error while reading occurs.
 * @return one if we found the parameter in the line, or if we had an error
 * while reading the named parameter, zero otherwise. */
int EGioReadIntParam( const int argc,
											char**argv,
											int*const param,
											int*const rval);

/* ========================================================================= */
/** @name Zlib, BZlib compability interface:
 * This functions pretend to provide a common interface to work with plain or
 * compresed files. Up to now, we implement zlib and plain text handling, as
 * well as basic input/output handling. */
/* @{ */
/* ========================================================================= */
/* ========================================================================= */
struct EGioFile_st;
typedef struct EGioFile_st EGioFile_t;
/* ========================================================================= */
/* comon functions */
/* ========================================================================= */
/** @brief Converts, formats, and writes the args to the file under control 
 * of the format string, as in fprintf.
 * @return the number of bytes actually written (0 in case of error).
 * @note The number of bytes written is limited to 4095 (inherited from
 * gzprintf). The caller should assure that this limit is not exceeded.
 * If it is exceeded, then EGioprintf() will return return an error (0) with 
 * nothing written. In this case, there may also be a buffer overflow with
 * unpredictable consequences, which is possible only if zlib was compiled
 * with the insecure functions sprintf() or vsprintf() because the secure
 * snprintf() or vsnprintf() functions were not available. For more details,
 * see zlib's manual. */
int EGioPrintf(EGioFile_t*file,const char *format, ...);
int EGioWrite(EGioFile_t*file,const char* buf);
/* ========================================================================= */
/** @brief open a file for read or write as in fopen. if the file name ends
 * with gz, it will (try to) use gzFile mode, if not, it will use regular FILE
 * modes 
 * @return on success, pointer to an EGioFile_t structure, if can not open
 * file, it will print system error on screen and exit execution. */
EGioFile_t* EGioOpen(const char *path, const char *mode);
/* ========================================================================= */
/** @brief open a C stdandard file FILE as an EGioFile_t, it can be used to
 * open stderr, stdout, stdin, it will use plain access to the file. At
 * closing, it will not close the three standard stdout, stderr, stdin, but
 * other files will be closed.
 * @param file pointer to a std-c FILE structure
 * @return pointer to a EGioFile_t structure, linked to the given file.
 * @note we assume that file points to an already opened file.
 * */
EGioFile_t* EGioOpenFILE(FILE*file);
/* ========================================================================= */
/** @brief close file, flush all pending information, free all internally
 * allocated information.
 * @return zero on success; non-zero otherwise
 * */
int EGioClose(EGioFile_t*file);
/* ========================================================================= */
/** @brief flush all un-writed information into the file.
 * @return zero on success, non-zero otherwise. */
int EGioFlush(EGioFile_t*file);
/* ========================================================================= */
/** @brief reads bytes from file until len-1 characters are read, or a newline
 * character is read and transferred to buf, or an end-of-file condition is
 * encountered. The string is then terminated with a null character.
 * @return buf on success, 0 otherwise. 
 * */
char* EGioGets(char*buf, int len, EGioFile_t*file);
/* ========================================================================= */
/** @brief Test if the given file stream pointer point to the end of file.
 * @return non-zero if stream pointer points to end of file, 0 otherwise.
 * */
int EGioEof(const EGioFile_t*const file);
/* ========================================================================= */
/** @brief Test if the given file stream has an internal error flag
 * @return zero if no error detected, otherwise, return the apropiate error
 * code */
int EGioError(const EGioFile_t*const file);
/* ========================================================================= */
/* end interface for zlib-compresed files */
/* @} */
/* ========================================================================= */
/** @name Deprecated Functions
 * All functions here are marked as deprecated, and may be discontinued in
 * sub-sequent releases */
/* @{ */
/* ========================================================================= */
/** @brief this discard all lines starting with comments and stores the next 
 * leading line in current and the next token in next, we assume that next 
 * does not contain data, and that current store the remainings of the
 * current line */
void EGioDisCom (char **next,
								 char **current,
								 const char *delim,
								 const char *coment,
								 char *store,
								 unsigned int storeSize,
								 FILE * in);
/* ========================================================================= */
/** @brief display function for strings
 * @param str pointer to a null terminated string of chars.
 * @param file pointer to a stream where we write the string.
 * @par Description:
 * This function just print the string on the file, it won't add a '\n' at the
 * end. */
void EGdisplayString (void *str,
											FILE * file);

/* @} */
/* ========================================================================= */
/** @} */
#endif
