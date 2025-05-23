
# -------------------------------------------------------------------------
#    Copyright 2002-2025 National Technology & Engineering Solutions of
#    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
#    NTESS, the U.S. Government retains certain rights in this software.
# 
#    This file is part of Xyce(TM) Parallel Electrical Simulator.
#    
#    Xyce(TM) is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#   
#    Xyce(TM) is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with Xyce(TM) Parallel Electrical Simulator.
#    If not, see <http://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
#
#  Purpose        :
#  This automake Makefile may be used to compile and link a runtime
#  loaded sharable object from a Verilog module.
#
#  It has been generated by admsXml-2.3.5 from the verilog file
#  rlc.va using the xyceMakefile.am.xml template.
#
#  To function properly, the directory containing this Makefile.am and the
#  verilog source code for the model should be placed under the main Xyce
#  source tree.  The "AC_CONFIG_FILES" list in configure.ac must be updated
#  to include this directory's "Makefile", and the "bootstrap" script rerun.
#
#  Once these modifications to the build system are performed, "configure"
#  executed in a build directory will create a subdirectory for this model,
#  in which "make plugin" and "make install-plugin" will work.

AM_CPPFLAGS = @Xyce_INCS@
lib_LTLIBRARIES = libADMSrlc.la
BUILT_SOURCES=N_DEV_ADMSrlc_bootstrap.C
plugin: libADMSrlc.la
install-plugin: install

N_DEV_ADMSrlc.C N_DEV_ADMSrlc.h: $(srcdir)/rlc.va
	admsXml -x -e $(top_srcdir)/utils/ADMS/adms.implicit.xml -e $(top_srcdir)/utils/ADMS/xyceVersion_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceBasicTemplates_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceAnalogFunction_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceHeaderFile_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceImplementationFile_nosac.xml $(srcdir)/rlc.va

N_DEV_ADMSrlc_bootstrap.C: $(srcdir)/rlc.va
	admsXml -x -e $(top_srcdir)/utils/ADMS/adms.implicit.xml -e $(top_srcdir)/utils/ADMS/xyceVersion_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceBasicTemplates_nosac.xml -e $(top_srcdir)/utils/ADMS/xyceBootstrapFile.xml $(srcdir)/rlc.va

N_DEV_ADMSrlc_bootstrap.lo: N_DEV_ADMSrlc.h

CLEANFILES=N_DEV_ADMSrlc.C N_DEV_ADMSrlc.h N_DEV_ADMSrlc_bootstrap.C

libADMSrlc_la_SOURCES = N_DEV_ADMSrlc.C $(BUILT_SOURCES)

