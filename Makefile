.PHONY: all

BIN=fabe
SRCDIR_C=libfa
SRCDIR_CPP=.

OBJDIR=obj
DEPDIR=dep

CMP_CPP=g++ -std=c++1z
CMP_C=gcc
INC=-I$(shell pwd) -I$(shell pwd)/libfa
WARN=-Wall -Wno-unused-result -Wno-deprecated-declarations -Wno-sign-compare -Wno-maybe-uninitialized -Wno-ignored-attributes -Wno-strict-aliasing 
WARN+=-Wno-misleading-indentation -Wno-format-overflow -Wno-nonnull-compare
OPTIM=-Ofast -march=native -funroll-loops -funsafe-loop-optimizations -falign-functions=16 -falign-loops=16 -fopenmp
NOOPTIM=-O0 -march=native -fopenmp
DBG=-g ${NOOPTIM}
PROF=-g ${OPTIM}
LINK=-lpthread -lprofiler -ltcmalloc

SRC_C=$(shell find "${SRCDIR_C}" -name "*.c")
SRC_CPP=$(shell find "${SRCDIR_CPP}" -name "*.cpp")
OBJF=${SRC_C:.c=.o} ${SRC_CPP:.cpp=.o}
OBJB=$(notdir ${OBJF})
OBJ=$(addprefix ${OBJDIR}/,${OBJB})

BOLD=$$(tput bold)
RED=$$(tput setaf 1)
GREEN=$$(tput setaf 2)
YELLOW=$$(tput setaf 3)
BLUE=$$(tput setaf 4)
MAGENTA=$$(tput setaf 5)
CYAN=$$(tput setaf 6)
NORMAL=$$(tput sgr0)

ECHOCC=>&2 echo "[${BOLD}${YELLOW}  C  ${NORMAL}]"
ECHOCP=>&2 echo "[${BOLD}${YELLOW} C++ ${NORMAL}]"
ECHOLD=>&2 echo "[${BOLD}${CYAN}  L  ${NORMAL}]"

OPT=${PROF} # Put desired optimisation level here

define compilecpp
${ECHOCP} $< ;\
mkdir -p ${DEPDIR} ;\
tmp=`mktemp` ;\
${CMP_CPP} ${DEFS} ${INC} -MM ${OPT} $< >> $$tmp ;\
if [ $$? -eq 0 ] ;\
then echo -n "${OBJDIR}/" > ${DEPDIR}/$(notdir $<).d ;\
cat $$tmp >> ${DEPDIR}/$(notdir $<).d ;\
rm $$tmp ;\
mkdir -p ${OBJDIR} ;\
cd ${OBJDIR} ;\
${CMP_CPP} ${DEFS} -c ${INC} ${OPT} ${WARN} ../$< ;\
else \
ret=$$? ;\
rm $$tmp ;\
exit $$ret ;\
fi
endef

define compilec
${ECHOCC} $< ;\
mkdir -p ${DEPDIR} ;\
tmp=`mktemp` ;\
${CMP_C} ${DEFS} ${INC} -MM ${OPT} $< >> $$tmp ;\
if [ $$? -eq 0 ] ;\
then echo -n "${OBJDIR}/" > ${DEPDIR}/$(notdir $<).d ;\
cat $$tmp >> ${DEPDIR}/$(notdir $<).d ;\
rm $$tmp ;\
mkdir -p ${OBJDIR} ;\
cd ${OBJDIR} ;\
${CMP_C} ${DEFS} -c ${INC} ${OPT} ${WARN} ../$< ;\
else \
ret=$$? ;\
rm $$tmp ;\
exit $$ret ;\
fi
endef

all: ${BIN}
	@true

-include ${DEPDIR}/*.d

${BIN}: ${OBJ}
	@${ECHOLD} ${BIN}
	@${CMP_CPP} ${OPT} ${LDIR} $^ ${LINK} -o ${BIN}

${OBJDIR}/%.o: ${SRCDIR_CPP}/%.cpp
	@$(compilecpp)

${OBJDIR}/%.o: ${SRCDIR_C}/%.c
	@$(compilec)

clean:
	@echo "Removing subdirectories..."
	@rm -rf ${OBJDIR} ${DEPDIR}
