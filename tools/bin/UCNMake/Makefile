ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# Linux with egcs
CXX           = g++
#CXX           = icc
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
#LD            = icc
LDFLAGS       = -g -v
SOFLAGS       = -shared


CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lMinuit
GLIBS         = $(ROOTGLIBS) -lMinuit

#------------------------------------------------------------------------------

FILE	    = MakeUCN
TESTO       = $(FILE).$(ObjSuf)
TESTS       = $(FILE).$(SrcSuf)
TEST        = $(FILE)$(ExeSuf)

OBJS          = $(TESTO) 

PROGRAMS      = $(TEST)

all:            $(PROGRAMS)
#$(HWORLD):      $(HWORLDO)
#		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
#		$(MT_EXE)
#		@echo "$@ done"

$(TEST):      	$(TESTO)
		$(LD) $(LDFLAGS) $(TESTO) $(LIBS) $(OutPutOpt) $(TEST)
		$(MT_EXE)
		@echo "$(TEST) done"

clean:
		@rm -f $(OBJS) *Dict.* core

.SUFFIXES: .$(SrcSuf)

###


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
