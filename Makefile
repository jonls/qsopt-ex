# =============================================================================
# This is the Makefile of EGlib
# Revision 2003-4-14
# - 2007-12-27
# 					- Separate Makefile into multiple files
# - 2005-08-17
#						- Update EGlib.so rule
# - 2005-08-01 
#						- Update to generate templates, and compute dependencies on the 
#							fly
# - 2003-11-18 
#						- Update to automatic dependency generation
# =============================================================================

.PHONY: indent doc clean default template library 
DEFAULT := template library
default = $(DEFAULT)

#==============================================================================
# Default rules for each type of file

library: template
	@$(MAKE) -f Makefile.library

template:
	@$(MAKE) -f Makefile.template

doc: template
	@$(MAKE) -f Makefile.library doc

clean:
	@$(MAKE) -f Makefile.library clean
	@$(MAKE) -f Makefile.template clean

indent: template
	@$(MAKE) -f Makefile.library indent

# end of Makefile
# =============================================================================
