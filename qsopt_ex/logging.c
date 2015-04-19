/* logging - Functions for logging/error reporting
 *
 * Copyright (C) 2015  Jon Lund Steffensen <jonlst@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "logging.h"
#include "logging-private.h"


static QSlog_func global_log_func = NULL;
static void *global_log_data = NULL;


void QSlog_set_handler(QSlog_func log_func, void *data)
{
	global_log_func = log_func;
	global_log_data = data;
}

void QSlogv(const char *format, va_list args)
{
	va_list args2;
	va_copy(args2, args);

	int n = vsnprintf(NULL, 0, format, args2);
	if (n < 0) {
		perror("vsnprintf");
		abort();
	}

	va_end(args2);

	char *buffer = malloc(n+1);
	if (buffer == NULL) {
		perror("malloc");
		abort();
	}

	n = vsnprintf(buffer, n+1, format, args);
	if (n < 0) {
		perror("vsnprintf");
		free(buffer);
		abort();
	}

	if (global_log_func != NULL) {
		global_log_func(buffer, global_log_data);
	} else {
		fprintf(stderr, "%s\n", buffer);
	}

	free(buffer);
}

void QSlog(const char *format, ...)
{
	va_list args;
	va_start(args, format);
	QSlogv(format, args);
	va_end(args);
}
