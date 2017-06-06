
L=-pthread
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}

##########
#
# EDIT the lines below so that they point to your kent source
# and gsl header files and libraries 
#
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a

HG_INC += -I/home/lowec/src/gsl/gsl-1.16_install/include
L += /home/lowec/src/gsl/gsl-1.16_install/lib/libgsl.a /home/lowec/src/gsl/gsl-1.16_install/lib/libgslcblas.a

HG_INC += -I/home/lowec/src/R/R-3.1.2/include
L += /home/lowec/src/R/R-3.1.2/src/nmath/standalone/libRmath.a

#
# END of basic editing
#
##########

##########
#
#
ifeq (${USE_SAMTABIX},1)
    KNETFILE_HOOKS=1
    USE_BAM=1
    USE_TABIX=1
    ifeq (${SAMTABIXINC},)
        SAMTABIXINC = ${SAMTABIXDIR}
    endif
    ifeq (${SAMTABIXLIB},)
        SAMTABIXLIB = ${SAMTABIXDIR}/libsamtabix.a
    endif
    HG_INC += -I${SAMTABIXINC}
    L+=${SAMTABIXLIB} -lz
    HG_DEFS+=-DUSE_SAMTABIX -DUSE_BAM -DUSE_TABIX -DKNETFILE_HOOKS
endif

#
# End of sam/bam and tabix
#
##########

##########
#
# If you compiled your kent source with SSL support
# then you will need to use these lines
#
ifeq (${USE_SSL},1)
	L+=-lssl -lcrypto
	HG_DEFS+=-DUSE_SSL
endif
#
###########


CC=gcc
ifeq (${COPT},)
    COPT=-O -g
endif
ifeq (${CFLAGS},)
    CFLAGS=
endif

SYS = $(shell uname -s)

ifeq (${HG_WARN},)
  ifeq (${SYS},Darwin)
      HG_WARN = -Wall -Wno-unused-variable
      HG_WARN_UNINIT=
  else
    ifeq (${SYS},SunOS)
      HG_WARN = -Wall -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    else
      HG_WARN = -Wall -Werror -Wformat -Wimplicit -Wreturn-type
      HG_WARN_UNINIT=-Wuninitialized
    endif
  endif
  # -Wuninitialized generates a warning without optimization
  ifeq ($(findstring -O,${COPT}),-O)
     HG_WARN += ${HG_WARN_UNINIT}
  endif
endif

CFLAGS += ${HG_WARN}

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

L += -lm -lz

all: copyNumberDiff countKmers countKmersWithTree faToGcStats fastqTrim bamToGcStats

copyNumberDiff: copyNumberDiff.o ${MYLIBS}
	${CC} ${COPT} -o copyNumberDiff${EXE} copyNumberDiff.o ${MYLIBS} $L

countKmers: countKmers.o murmur3.o ${MYLIBS}
	${CC} ${COPT} -o countKmers${EXE} countKmers.o murmur3.o ${MYLIBS} $L

countKmersWithTree: countKmersWithTree.o ${MYLIBS}
	${CC} ${COPT} -o countKmersWithTree${EXE} countKmersWithTree.o ${MYLIBS} $L

faToGcStats: faToGcStats.o ${MYLIBS}
	${CC} ${COPT} -o faToGcStats${EXE} faToGcStats.o ${MYLIBS} $L

fastqTrim: fastqTrim.o ${MYLIBS}
	${CC} ${COPT} -o fastqTrim${EXE} fastqTrim.o ${MYLIBS} $L

bamToGcStats: bamToGcStats.o ${MYLIBS}
	${CC} ${COPT} -o bamToGcStats${EXE} bamToGcStats.o ${MYLIBS} $L

murmur3.o: murmur3.c murmur3.h

clean:
	rm -f copyNumberDiff copyNumberDiff.o countKmers countKmers.o faToGcStats faToGcStats.o bamToGcStats bamToGcStats.o murmur3.o fastqTrim fastqTrim.o countKmersWithTree countKmersWithTree.o

