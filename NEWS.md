
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
