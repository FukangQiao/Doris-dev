#########################################################################
# Makefile for the Doris software.					#
# Created by: configure						#
# 04-Dec-1998                                                           #
#                                                                       #
# Delft University of Technology                                        #
# Delft Institute of Earth Observation and Space Systems                #
# http://doris.tudelft.nl               				#
# See also the user manual, annex installation.				#
#									#
### Usage ### 								#
# The command: "make" creates and optimized version of doris.		#
# The command: "make install" installs it andcleans up.			#
# If this does not work, try a more controlled approach by:		#
#									#
# 0. inspect set up of this Makefile					#
# 1. set CC DEFS LIBS etc., first set these to debug values		#
#    (read NOTE there if problems)					#
# 2. compile software sources to *.o:     "make swobjs"			#
# 3. link object together to executable:  "make doris"			#
# 4. install executable:                  "make install"		#
# 5. remove object files:                 "make clean"			#
#									#
# BTW. "make" does "swobjs doris" by default.				#
#									#
###        								#
# Successfully tested for compilers/platforms:				#
#  GNU g++ version 2.7.2.2 on HP-UX B.10.20 A 9000/785			#
#  GNU g++ version 2.95.2  on HP-UX B.10.20 A 9000/785			#
#  HP  aCC version A.01.07 on HP-UX B.10.20 A 9000/785			#
#  GNU g++ version ?       on Linux X86					#
#  SGI ?								#
#  Sun Studion 9,10,11     on Solaris Sparc and X86/AMD64		#
#  GNU g++ versions v3.4, v4.0, v4.1, v4.2, v4.3 on Linux X86 and AMD64 # 
#  Sun Sudion v11          on Linux AMD64				#
#  Intel Compilers v9, v10, v11 on Linux AMD64				#
#########################################################################
###################################################################
###################################################################
# Please change if required: DEF[1-6], CC, CFLAGS, LIBDIR,      ###
###################################################################
### The shell used by this makefile ###
SHELL   = /bin/sh

### Specify compiler/installation directory ###
INSTALLDIR = /usr/local/bin
CC	   = g++
SCRIPTSDIR = ../bin

### Define statements controlling compilation ###
# Comment out DEF[1-8] statements when appropriate.
DEF1    = -D__DEBUGMAT1	# 		-1- range checking in matrix class
DEF2    = -D__DEBUGMAT2	# 		-2- info matrix class (debug only)
DEF3    = -D__DEBUG #			-3- debug for other files
DEF4    = -Wno-deprecated            #  -4- do not display warnings due to depricated entries 
DEF5    = -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE #  -5- support for large files
#DEF6    = -DNO_FASTTRIG              #  -6-  extra functional switches
#DEF9     = -D__NO_IOS_BINARY__ # -9- comment out if you ios::binary (most will)
#
### Do not change following 2, simply comment DEF1-8 above
DEFSDEBUG 	= $(DEF1) $(DEF2) $(DEF3)
DEFS		= $(DEF4) $(DEF5) $(DEF6)   -D__X86PROCESSOR__  

### NOTE ###
# If you get warnings like the following (I got them with g++, not with aCC)
# then use define __GNUC__ to define these functions
#%// BK 15-Aug-2000
#In file included from matrixbk.h:38,
#		 from processor.cc:23:
#constants.h: In method `real8 cn::dist(cn) const':
#constants.h:106: implicit declaration of function `int sqr(...)'
### END NOTE ###


# LOCAL dir
LOCAL_INSTALLDIR = ../bin

### CFLAGS ###
# Compiler options flag ###
# Specify compiler flags (I selected these for g++)
CFLAGSDEBUG 	= -g -O $(DEFS) $(DEFSDEBUG)
CFLAGSOPT 	= -O $(DEFS)
#CFLAGSOPT 	= -O2 $(DEFS)   # this crashed doris! (?)
#CFLAGSOPT 	= -O3 $(DEFS)   # this crashed doris! (?)

### CCNOTE ###
### Change here ###
### set CFLAGS depending on if you compile debug version or working version ###
#CFLAGS 	= $(CFLAGSDEBUG)
CFLAGS 	= $(CFLAGSOPT)
### END CCNOTE ###

### Library locations flag ###
### -lcl : used for veclib 
### -lm  : used for fftw
LFLAGS  =  


#####################################################
### No need to change anything below here... ####
#####################################################
EXECUTABLE = 	doris

SWSRCS	=	processor.cc \
		utilities.cc \
		ioroutines.cc \
		readinput.cc \
		readdata.cc \
		coregistration.cc \
		filtering.cc \
		referencephase.cc \
		products.cc \
		geocode.cc \
		unwrap.cc \
		matrixspecs.cc \
		exceptions.cc \
		estorbit.cc
SWOBJS	=	$(SWSRCS:.cc=.o)

### scripts of doris.tar in SCRIPTSDIR, used at make install
SCRIPTS	=	helpdoris \
		baseline.doris.sh \
		baseline.doris.csh \
                construct_dem.sh \
		coregpm.doris \
		doris* \
		heightamb \
		phasefilt.doris \
		plotcpm* \
		plotoffsets* \
		run \
		viewanddel \
		cpx2ps \
		lonlathei2ascii \
		ascii2ascii_UTM \
		ascii2ps* \
                tsx* \
                rs2* \
                csk* \
                gammaReadfiles.csh \
                hhmmss2sec.py \
                sec2hhmmss.py




#####################################################
### And no need to change anything below here... ####
#####################################################
### what to do if make w/o arguments ###
# first make all object files, then link them together
default:        $(EXECUTABLE)


### how to compile object code .o from C++ source files .cc (general rule) ###
# space between -o and filename for SUN make (BK 7 july 2000)
.cc.o:		
		$(CC) $(CFLAGS) -c -o $(@) $<

#####################################################
### Make object code from source ###
swobjs:		$(SWOBJS)


### Load objects to executable ###
### The added .cc files contain class definitions
### matrixclass is explicitly included
### in the source code since we had some problems finding out how to do that
### here, since compilers treat the template class differently than normal classes.
### BK 07-Feb-2002
$(EXECUTABLE):	$(SWOBJS) tmp_strptime.cc slcimage.cc productinfo.cc orbitbk.cc
		$(CC) $(CFLAGS) -c -o triangle.o -DTRILIBRARY -DANSI_DECLARATORS triangle.c
		$(CC) $(CFLAGS)\
		tmp_strptime.cc $(SWOBJS) triangle.o\
		slcimage.cc productinfo.cc orbitbk.cc\
		$(LFLAGS) -o $@
		@echo " "
		@echo "*******************************"
		@echo "* ...Compilation finished...  *"
		@echo "*******************************"
		@echo " "


### Install executable in installdir ###
install:	$(EXECUTABLE)
		@echo "* Installing $(EXECUTABLE) in: $(INSTALLDIR)"
		@cp -f $(EXECUTABLE) $(INSTALLDIR)/.
		( cd $(SCRIPTSDIR); cp -f $(SCRIPTS) $(INSTALLDIR)/. )
		$(MAKE) cleaner
		@echo " "
		@echo "*******************************"
		@echo "* ...Installation finished... *"
		@echo "*******************************"
		@echo " "
		@echo "* Check that $(INSTALLDIR) is in your path search: echo \$$PATH ."
		@echo " "


installcb:	$(EXECUTABLE)
		@echo "* Installing $(EXECUTABLE) in: $(LOCAL_INSTALLDIR)"
		@cp -f $(EXECUTABLE) $(LOCAL_INSTALLDIR)/.
		$(MAKE) cleaner
		@echo " "
		@echo "*******************************"
		@echo "* ...Installation finished... *"
		@echo "*******************************"
		@echo " "
		@echo "* Check that $(LOCAL_INSTALLDIR) is in your path search: echo \$$PATH ."
		@echo " "


#####################################################
### Testers ###
test:		testdoris
testdoris:	$(EXECUTABLE)
		@echo " "
		@echo "* Executing command: $(EXECUTABLE) -v"
		$(EXECUTABLE) -v
		@echo " "

### Orbit test program for debugging ###
test-orbit:	ioroutines.o matrixspecs.o utilities.o exceptions.cc slcimage.cc orbitbk.cc matrixbk.cc bk_messages.hh
		$(CC) $(CFLAGS) -D__TESTMAIN__ \
		ioroutines.o matrixspecs.o utilities.o exceptions.cc slcimage.cc orbitbk.cc \
		$(LFLAGS) \
		-o $@
### Matrix test program for debugging ###
### fast_sin defined in utilities.cc, which requires ioroutines, which, etc.
test-matrix:  matrix_test.cc matrixspecs.o utilities.o utilities.o ioroutines.o matrixbk.cc
		$(CC) $(CFLAGS) matrix_test.cc matrixspecs.o \
		utilities.o ioroutines.o exceptions.cc orbitbk.cc \
		$(LFLAGS) \
		-o $@
		@echo " "


#####################################################
### Cleaners ###
clean:		cleanswobjs 

cleaner:	clean cleanprog
cleanest:	clean cleanprog uninstall
cleanprog:
		@rm -f $(EXECUTABLE) $(EXECUTABLE.debug) \
		matrixbk_test testorbit
		@echo "* Removed executables in current dir."
cleanswobjs:
		@rm -f $(SWOBJS) matrixbk_test.o orbitbk.o slcimage.o triangle.o
		@echo "* Removed object files."
uninstall:	
		@rm -f $(INSTALLDIR)/$(EXECUTABLE)
		@echo "* Removed executables in install dir: $(INSTALLDIR)."

### BK 12-Dec-2000 ###

