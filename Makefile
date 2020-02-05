
ROOTSYS = /cvmfs/ams.cern.ch/Offline/root/Linux/root-v5-34-9-icc64.14-slc6
#ROOTSYS = /afs/cern.ch/ams/Offline/root/Linux/root-v5-34-9-icc64-slc6

XRD    = /afs/cern.ch/ams/local2/opt/xrootd-icc64-12

Offline = /cvmfs/ams.cern.ch/Offline
#Offline = /afs/cern.ch/ams/Offline

AMS  = $(Offline)/vdev
AMS  = /afs/cern.ch/work/q/qyan/public/AMS307
AMSL = $(Offline)/vdev/lib/linuxx8664icc5.34
AMSL =  /afs/cern.ch/work/q/qyan/public/AMS307/lib/linuxx8664icc5.34

DEFS = -D_PGTRACK_ -D__ROOTSHAREDLIBRARY__
FLGS = -I$(ROOTSYS)/include -I$(AMS)/include -ICC -Iinclude -fPIC -std=c++11

ifdef DEBUGFLAG
  OPTS = -g
else
  OPTS = -O3
endif

INTEL = /cvmfs/projects.cern.ch/intelsw/psxe/linux/x86_64/2018/compilers_and_libraries_2018.2.199/linux/
#INTELDIR = /afs/cern.ch/ams/local2/opt/intel/
#INTELVER_DEFAULT = Compiler/11.1/073
STAT = -static-intel -static-libgcc -static-libstdc++ -Bstatic
PARA = -D__X8664__ -qopenmp -D__AMSPARALLEL__ 
IOPT = -axsse4.2,ssse3,AVX,CORE-AVX2 -fp-model source -fast-transcendentals
CXX  = $(INTEL)/bin/intel64/icc   $(IOPT) -vec-report0 $(PARA) -m64 
LD   = $(INTEL)/bin/intel64/icpc  $(IOPT) $(PARA)
LDS  = $(INTEL)/bin/intel64/icpc  $(IOPT) $(PARA) $(STAT)
F77  = $(INTEL)/bin/intel64/ifort $(IOPT) $(PARA)

CXXR = g++
LDR  = g++

UNILIB = /afs/cern.ch/work/s/selu/private/lib/unilib
NAGDIR  = /afs/cern.ch/exp/ams/Offline/CERN/NagLib
CERNDIR = /afs/cern.ch/exp/ams/Offline/CERN/2005
#CERNDIR = /cvmfs/sft.cern.ch/lcg/external/cernlib/2005/slc4_ia32_gcc34/lib 
FFDIR = $(INTEL)/lib/intel64

LIBF = -L$(NAGDIR) -lnag64 -lifcore -L$(CERNDIR)/lib -lmathlib -L$(FFDIR) -limf
LIBX = -lMinuit -lMinuit2 -lRFIO -lTMVA -lXMLIO -lMLP -lTreePlayer
LIBR = $(shell $(ROOTSYS)/bin/root-config --libs) -lNetx $(LIBX)

LISX = -llzma -Llib -lshifts -lcrypto -L$(XRDLIB)/lib64 -lXrdClient -lXrdUtils 
#LISR = -L$(ROOTSYS)/lib -lRoot -lfreetype -pthread -lpcre $(LISX) $(LIBF) 
LISR =  $(shell $(ROOTSYS)/bin/root-config --libs) -lRoot -lfreetype -pthread -lpcre $(LISX) $(LIBF) 

LIB4 = $(AMSL)/libntuple_slc4_PG.a
LIB6 = $(AMSL)/libntuple_slc6_PG.a
LIBP = $(LIB6)


all:

obj/%.o: CC/%.C
	@if ! [ -d obj ] ; then mkdir -p obj; fi
	$(CXX) -c -o $@ $(DEFS) $(FLGS) $(OPTS) $<

bin/%: obj/%.o $(LIBP)
	@if ! [ -d bin ] ; then mkdir -p bin; fi
	$(CXX) -o $@ $^ $(LISR)
       #readelf -d $@

F/%.o: F/%.f
	$(F77) -c $< -I$(UNILIB)/src/ -o $@
	


readme:  doc/README.md
	Markdown.pl doc/README.md > doc/README.html
clean:
	rm -f bin/* obj/*
