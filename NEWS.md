
2.5.10.3
--------
- All output to stdout/stderr now goes through an internal logging function.
	The log messages can optionally be redirected to a user supplied function
	(by default they are written to stderr).
- Some headers that are internal to the library are no longer installed.
- Restored possibility of compressing the output files (Denis Rosset).
- Remove unused support for dmalloc.
- Convert demo_qs into test_qs TAP test suite.

2.5.10.2
--------
* (The patch `_p` version notation used in 2.5.10_p1 is not widely supported so
  this release simply uses the dot notation.)
* Removes Python module (moved to separate repository at
  https://github.com/jonls/python-qsoptex).
* Move library source files into `qsopt_ex` so the header files can be accessed
  in the same way as when installed (`<qsopt_ex/XX.h>`).
* Explicitly include the typename of template types when using a macro,
  variable or function that is template generated. This removes the
  build dependencies on Exuberant Ctags and GNU sed, as well as speeding up
  the template generation significantly.
* Remove util.c/util.h from template generation since the majority of symbols
  are not template type specific.

2.5.10_p1
---------
* Based on original version v2.5.10 published by Daniel Espinoza et al.
* Changed to autotools-based build system.
* Removed external dependency on EGlib.
* Build library using libtool for portability.
* Add Cython-based Python module to interface with libqsopt_ex.
* Fix sprintf calls with missing format string.
* Add missing header declarations, includes.
* Clean up headers to make external use easier.
* Remove some unused functions/macros that caused compiler warnings/errors.
* Workaround: Writing solution to gz-file was broken; temporarily disabled gzip output.
* Add README file with build instructions and code examples.
* Add Travis CI build script.
