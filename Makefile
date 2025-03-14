# The rootbeer makefile

# get root libs
ROOTLIBS	= $(shell root-config --libs)
ROOTINCDIR	= $(shell root-config --incdir)

#if ROOTBEER not defined then make it the current working directory
ifndef ROOTBEER
 OBJDIR		= ./obj/Linux
 BINDIR		= ./bin/Linux
 SODIR		= ./slib/Linux
else
 OBJDIR		= $(ROOTBEER_OBJ)
 BINDIR		= $(ROOTBEER_BIN)
 SODIR		= $(ROOTBEER_SLIB)
endif

#define some directories
SRCDIR		= src
INCDIR		= include
DATDIR		= dat
SCRIPTDIR	= scripts
SC		= sample_code
DOCDIR		= doc
PACK		= extra_packages
ELOSS_RB_DIR	= $(PACK)/eloss_rb
ELOSS_RB_OBJ	= $(PACK)/eloss_rb/$(OS_NAME)

#name of file with bankdefs
CLASBANKS     	=  $(INCDIR)/clasbanks.ddl

#only require this if building stuff hacked from clas packages
CLASPACK =$(CLAS_PACK)
#CLASPACK = /home/claslib/builds/DEVELOPMENT/packages
#CLASINC = -I$(CLASPACK)/include -I$(CLASPACK)/inc_derived

CXX		= g++
LD		= g++
CC		= gcc
#G77		= g77
G77		= gfortran

# allow debugging with >make DEBUG=1 
## ifeq ($(DEBUG),1)
ifdef DEBUG
CXXFLAGS      	= -g -Wall -Wno-deprecated -DR__THREAD -fno-exceptions -fno-strict-aliasing \
		-fPIC -I$(ROOTINCDIR) -I./ -I$(INCDIR) $(CLASINC) -DDEBUG
CFLAGS	      	= -g -Wno-deprecated -D_REENTRANT -fno-strict-aliasing -I$(INCDIR) $(CLASINC) -DDEBUG
#FFLAGS	      	= -g -mcpu=pentium -fno-automatic -finit-local-zero -ffixed-line-length-none \
#		-fno-second-underscore -DLinux
FFLAGS 		= -g
LDFLAGS       	= -g

else

CFLAGS		= -O2 -Wno-deprecated -Wno-unused -D_REENTRANT -fno-strict-aliasing -I$(INCDIR) -I$(SC) $(CLASINC)
CXXFLAGS      	= -O2 -Wno-deprecated -Wno-unused -Wall -DR__THREAD \
		-fno-exceptions -fno-strict-aliasing -fPIC -I$(ROOTINCDIR) -I./ -I$(INCDIR) -I$(SC) $(CLASINC)
LDFLAGS       	= -O2 -Wall -Wno-deprecated -Wno-unused
FFLAGS 		= -O2
endif

FFLAGS += -fno-f2c -fno-automatic -funroll-loops -fomit-frame-pointer \
	-fPIC -ffixed-line-length-none -fno-second-underscore

SOFLAGS       	= -shared
LIBS          	= $(ROOTLIBS) -lm -ldl -rdynamic
RBLIB		= -L$(SODIR) -lRootBeer


GAWK 		= gawk

LIBPTHREAD	= -lpthread

#------------------------------ some rules ------------------------------

.SUFFIXES: 	.cxx

#rules for .cxx
$(OBJDIR)/%.o: 	$(SRCDIR)/%.cxx $(INCDIR)/%.h 
		$(CXX) $(CXXFLAGS) -o $@ -c $<

%.o: 		%.cxx
		$(CXX) $(CXXFLAGS) -o $@ -c $<


#rule for .C
$(OBJDIR)/%.o: 	$(SRCDIR)/%.C
		$(CXX) $(CXXFLAGS) -o $@ -c $<

#rule for c
$(OBJDIR)/%.o: 	$(SRCDIR)/%.c
		$(CC) $(CFLAGS) -o $@ -c $<

%.o:		%.c
		$(CC) $(CFLAGS) -o $@ -c $<

%.o:		%.F
		$(G77) $(FFLAGS) -DLinux -o $@ -c $<

$(ELOSS_RB_OBJ)/%.o:	$(CLASPACK)/eloss/%.F
		$(G77) $(FFLAGS) -DLinux -o $@ -c $<

#-------------------------------------------------------------------------
BANKDUMPEXE	= $(BINDIR)/bankdump

all:           	$(ROOTBEERSO) $(BANKDUMPEXE) getbanks ExpTable 

## htmldoc


#---------------- all the basic rootbeer librares and executables----------------------------
# should not need to touch this stuff - long and explicit, no fancy rules here

BANKVARS		= bankvars
BANKVARSO		= $(OBJDIR)/bankvars.o $(OBJDIR)/bankvarsDict.o

$(BANKVARS):		$(BANKVARSO)

$(OBJDIR)/bankvars.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
			$(INCDIR)/TDST.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/bankvars.cxx

$(OBJDIR)/bankvarsDict.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
				$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
				$(INCDIR)/TDST.h

			@echo "Generating dictionary bankvarsDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/bankvarsDict.cxx -c \
			$(INCDIR)/bankvars.h $(INCDIR)/bankvarsLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/bankvarsDict.cxx

## generate headers, banklist etc
$(INCDIR)/bankvars.h:	$(SCRIPTDIR)/bank2struct.gk $(CLASBANKS) $(DATDIR)/singles.dat
			$(GAWK) -v incdir=$(INCDIR) -v srcdir=$(SRCDIR) \
			-f $(SCRIPTDIR)/bank2struct.gk $(CLASBANKS) > $(INCDIR)/bankheader.h




ROOTBEER    		= RootBeer
ROOTBEERO   		=  $(OBJDIR)/TRootBeer.o  $(OBJDIR)/RootBeerDict.o

$(ROOTBEER):		$(ROOTBEERO)

$(OBJDIR)/TRootBeer.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(SRCDIR)/TRootBeer.cxx
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/TRootBeer.cxx


$(OBJDIR)/RootBeerDict.o: $(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/RootBeerLinkDef.h  

			@echo "Generating dictionary RootBeerDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/RootBeerDict.cxx -c \
			$(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/RootBeerDict.cxx




BOS         		= Bos
BOSO        		= $(OBJDIR)/TBos.o $(OBJDIR)/BosDict.o

$(BOS):			$(BOSO)
$(OBJDIR)/TBos.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h $(SRCDIR)/TBos.cxx
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/TBos.cxx

$(OBJDIR)/BosDict.o: 	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h  \
			$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h $(INCDIR)/BosLinkDef.h

			@echo "Generating dictionary BosDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/BosDict.cxx -c $(INCDIR)/TBos.h \
			$(INCDIR)/BosLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/BosDict.cxx




DST         		= Dst
DSTO        		= $(OBJDIR)/TDST.o $(OBJDIR)/DSTDict.o

$(DST):			$(DSTO)
$(OBJDIR)/TDST.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/TDST.h 
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/TDST.cxx

$(OBJDIR)/DSTDict.o: 	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/TDST.h $(INCDIR)/DSTLinkDef.h

			@echo "Generating dictionary BosDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/DSTDict.cxx -c $(INCDIR)/TDST.h \
			$(INCDIR)/DSTLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/DSTDict.cxx




DSTWriter       	= DSTWriter
DSTWriterO      	= $(OBJDIR)/TDSTWriter.o $(OBJDIR)/DSTWriterDict.o

$(DSTWriter):		$(DSTWriterO)
$(OBJDIR)/TDSTWriter.o:	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
			$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
			$(INCDIR)/TDST.h $(SRCDIR)/TDSTWriter.cxx \
			$(INCDIR)/TDSTWriter.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/TDSTWriter.cxx

$(OBJDIR)/DSTWriterDict.o: 	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
				$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
				$(INCDIR)/TDST.h \
				$(INCDIR)/TDSTWriter.h $(INCDIR)/DSTWriterLinkDef.h

			@echo "Generating dictionary BosDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/DSTWriterDict.cxx -c $(INCDIR)/TDSTWriter.h \
			$(INCDIR)/DSTWriterLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/DSTWriterDict.cxx




ROOTBEERUTIL    = 	RootBeerUtil
ROOTBEERUTILO   = 	$(OBJDIR)/RootBeerUtil.o $(OBJDIR)/RootBeerUtilDict.o

$(OBJDIR)/RootBeerUtil.o: 	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
				$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
				$(INCDIR)/TDST.h \
				$(INCDIR)/RootBeerUtil.h $(SRCDIR)/RootBeerUtil.cxx
				$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/RootBeerUtil.cxx

$(OBJDIR)/RootBeerUtilDict.o: 	$(INCDIR)/bankvars.h $(INCDIR)/TRootBeer.h $(INCDIR)/RootBeerTypes.h \
				$(INCDIR)/TBos.h $(INCDIR)/BosTypes.h \
				$(INCDIR)/TDST.h \
				$(INCDIR)/RootBeerUtil.h $(INCDIR)/RootBeerUtilLinkDef.h

				@echo "Generating dictionary BosDict..."
				$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/RootBeerUtilDict.cxx -c \
				$(INCDIR)/RootBeerUtil.h $(INCDIR)/RootBeerUtilLinkDef.h
				$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/RootBeerUtilDict.cxx




ROOTBEERSO  	= $(SODIR)/libRootBeer.so

$(ROOTBEERSO):	$(ROOTBEERO)  $(BOSO)  $(DSTO) $(ROOTBEERUTILO) $(DSTWriterO) $(BANKVARSO)

		$(LD) $(SOFLAGS) $(LDFLAGS) -o $@ $^ $(ROOTLIBS) -lThread

		@echo "$(ROOTBEER) done"




$(BANKDUMPEXE): $(ROOTBEERSO) $(SRCDIR)/bankdump.cxx
		$(CXX) $(CXXFLAGS) -DROOTEXE  -o $(BANKDUMPEXE) $(SRCDIR)/bankdump.cxx $(RBLIB) \
		$(ROOTLIBS) -lThread -lm


EXPTABLESO      = $(SODIR)/libExpTable.so
EXPTABLEO       = $(OBJDIR)/TExpTable.o  $(OBJDIR)/ExpTableDict.o

$(EXPTABLESO):  $(EXPTABLEO)
		$(LD) $(SOFLAGS) $(LDFLAGS) -o $@ $^ $(ROOTLIBS) -lThread

		@echo "$(EXPTABLESO) done"


ExpTable:       $(EXPTABLESO)

$(OBJDIR)/TExpTable.o:  $(INCDIR)/TExpTable.h $(SRCDIR)/TExpTable.cxx
		$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/TExpTable.cxx


$(OBJDIR)/ExpTableDict.o: $(INCDIR)/TExpTable.h $(INCDIR)/ExpTableLinkDef.h
			@echo "Generating dictionary ExpTableDict..."
			$(ROOTSYS)/bin/rootcint -f $(SRCDIR)/ExpTableDict.cxx -c \
			$(INCDIR)/TExpTable.h $(INCDIR)/ExpTableLinkDef.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SRCDIR)/ExpTableDict.cxx


htmldoc:	$(ROOTBEERSO) $(SRCDIR)/TRootBeer.cxx $(SRCDIR)/TDST.cxx $(SRCDIR)/TDSTWriter.cxx\
		$(SRCDIR)/TBos.cxx $(SRCDIR)/TExpTable.cxx
		@echo " Invoking ROOT html automatic documentation"
#		@rm -rf htmldoc
		root -b -n -q sample_code/htmldoc.C


GETBANKS	= $(BINDIR)/getbanks
getbanks:	$(GETBANKS)
$(GETBANKS):	$(SCRIPTDIR)/getbanks
		cp -f $(SCRIPTDIR)/getbanks $(BINDIR)/getbanks

#---------------- end of basic rootbeer librares and executables--------------------------------------------------

#---------------- make tar -----------------------------------------------------------------------------
TD	= $(shell basename `pwd`)
tar:
#rm all obj files in $(PACK) subdirs 
#	rm -f $(PACK)/*/*.o 
	cd ../ ;\
	tar --create --verbose --file $(TD)/$(TD).tar \
	 --exclude='*Dict*' \
	 --exclude='*.so' \
	 --exclude='*.o' \
	 --exclude='*.d' \
	 --exclude='*core*' \
	 --exclude='*~' \
	 --exclude='htmldoc/*' \
	 --exclude='lib/*' \
	 --exclude='slib/*' \
	 --exclude='bin/*' \
	 --exclude='test/*' \
	 --exclude='obj/*' \
	$(TD)

#-------------- end make tar ---------------------------------------------------------------------------

clean:  
	/bin/rm -f $(BOSO) $(DSTO) $(ROOTBEERO) $(ROOTUTILO) $(ROOTBEERSO) $(BANKSERVERO) \
	$(EXPTABLESO):  $(EXPTABLEO) \
	$(BANKSERVEREXE) $(BANKDUMPO) $(DSTWriterO) $(BANKDUMPEXE) $(BANKVARSO) $(ROOTBEERUTILO) $(PERFO)
	/bin/rm -f $(INCDIR)/bankvars.h $(INCDIR)/bankvars.h \
	$(SRCDIR)/bankvars.cxx \
	$(INCDIR)/bankvarsLinkDef.h $(SRCDIR)/*Dict.* $(INCDIR)/*Dict.* core*

#_______________________________________________________________________________________________________________________


# Some more obscure examples on how to turn clas packages into root shared library
# with functions callable in CINT

## incorporation of Eugene's eloss package _____________________________________________________________________________
## Build the library 
ELOSSSO		= $(SODIR)/libeloss.so

eloss: 		$(ELOSSSO)

# build the shared library by using the .F files from the standard eloss distribution, which must exist, or be linked as 
# extra_packages/eloss
# and some extra .F and .cxx and .h in extra_packages eloss_rb
#
$(ELOSSSO):	$(ELOSS_RB_OBJ)/elossDict.o $(INCDIR)/eLoss.h $(ELOSS_RB_DIR)/elossLinkDef.h \
		$(addprefix $(ELOSS_RB_OBJ)/, $(addsuffix .o, $(notdir $(basename $(filter %.F %.cxx, $(wildcard $(ELOSS_RB_DIR)/*)))))) \
		$(addprefix $(ELOSS_RB_OBJ)/, $(addsuffix .o, $(notdir $(basename $(filter %.F, $(wildcard $(CLASPACK)/eloss/*))))))


		$(LD) $(SOFLAGS) $(LDFLAGS) $(ELOSS_RB_OBJ)/*.o  -o $@ -lgfortran

$(INCDIR)/eLoss.h:	$(ELOSS_RB_DIR)/eLoss.h
			cp $(ELOSS_RB_DIR)/eLoss.h $(INCDIR)/

$(ELOSS_RB_OBJ)/elossDict.o: 	$(ELOSS_RB_DIR)/eLoss.h $(ELOSS_RB_DIR)/elossLinkDef.h
				@echo "Generating dictionary $(PACK)/eloss/elossDict.cxx"
				@test -d $(ELOSS_RB_OBJ) || mkdir -p $(ELOSS_RB_OBJ)
				cp $(ELOSS_RB_DIR)/*.cxx $(ELOSS_RB_OBJ)
				cp $(ELOSS_RB_DIR)/*.[hF] $(ELOSS_RB_OBJ)
				$(ROOTSYS)/bin/rootcint -f $(ELOSS_RB_OBJ)/elossDict.cxx \
				-c $(ELOSS_RB_OBJ)/eLoss.h $(ELOSS_RB_OBJ)/elossLinkDef.h
				$(CXX) $(CXXFLAGS) -o $@ -c $(ELOSS_RB_OBJ)/elossDict.cxx

eloss_clean: 
		rm -f $(ELOSS_RB_OBJ)/*.o $(SODIR)/libeloss.so
#_______________________________________________________________________________________________________________________

#shared library PolHandler - an example of a class
#but not very tidy code in the class!

POLHANDLERSO		= $(SODIR)/libPolHandler.so
PolHandler:		$(POLHANDLERSO)
$(POLHANDLERSO):	$(OBJDIR)/PolHandler.o 	$(OBJDIR)/PolHandlerDict.o
			$(LD) $(SOFLAGS) $(LDFLAGS) $^ -o $@ $(ROOTLIBS) $(RBLIB) -L./  -lm

$(OBJDIR)/PolHandler.o: $(ROOTBEERSO) $(SC)/PolHandler.C $(SC)/PolHandler.h
			$(CXX) $(CXXFLAGS) -o $@ -c $(SC)/PolHandler.C

$(OBJDIR)/PolHandlerDict.o: $(ROOTBEERSO) $(SC)/PolHandler.h $(SC)/PolHandlerLinkDef.h
			@echo "Generating dictionary $@"
			$(ROOTSYS)/bin/rootcint -f $(SC)/PolHandlerDict.cxx -c -I$(INCDIR) $(SC)/PolHandler.h \
			$(SC)/PolHandlerLinkDef.h
			$(CXX) $(CXXFLAGS) -o $(OBJDIR)/PolHandlerDict.o -c $(SC)/PolHandlerDict.cxx

PolHandler_clean :
	rm -f $(OBJDIR)/PolHandler*.o $(SC)/PolHandlerDict.cxx

#_______________________________________________________________________________

dstmaker_pass0: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_pass1: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_n2pi: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_KpSm: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_strange: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_Ppim: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_2pos1neg: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_pol: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "

dstmaker_gpid: $(SC)/$< $(BOSO) $(ROOTBEERO) $(DSTO) $(INCDIR)/bankheader.h $(INCDIR)/bankvars.h
	$(CXX) $(CXXFLAGS)  -o $(BINDIR)/$@ -DROOTEXE $(SC)/$@.C $(ROOTLIBS) -L./ -L$(SODIR) -lRootBeer -lm
	@echo "Built $@ "


#_______________________________________________________________________________________________________________________

