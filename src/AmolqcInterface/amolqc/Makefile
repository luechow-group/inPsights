#
# recursive Makefile for amolqc
#

# set default target to all
.PHONY: all clean debug serialtest
all:
utilsdir := utils
FPPTARGETS :=
FPP90TARGETS :=
FNOTARGETS :=
FOPENMPTARGETS :=

# user config
include rules.mk

ifdef DEBUG
OPTLEVEL = debug
endif

# rules
include vars.mk
include utils/rules.mk
include src/rules.mk

ifdef DEBSYM
FFLAGS += -g
F90FLAGS += -g
FNOFLAGS += -g
LDLAGS += -g
endif

.PHONY: all clean debug tests
all: src/amolqc
	@git log --pretty=format:'%ad %h %d' --abbrev-commit --date=short -1 > version
clean: clean-mod clean-utils clean-src
debug: all
tests: testserial
.PHONY: depend
depend: src/depends.mk utils/depends.mk

.PHONY: clean-mod
clean-mod:
	@echo "Cleaning modules"
	@-rm -f mod/*
	@-rm -f src/*.mod
	@-rm -f utils/*.mod
.PHONY:testserial
testserial:
	cd testsuite/run;chmod +x test.sh;./test.sh  | tee test.log

# Remove implict rules
.SUFFIXES:

# add -fopenmp for files that require it
$(FOPENMPTARGETS): FFLAGS += -fopenmp
$(FOPENMPTARGETS): F90FLAGS += -fopenmp
$(FOPENMPTARGETS): FNOFLAGS += -fopenmp

# use FPP (defined in vars.mk) for targets that require preprocessing
$(FPPTARGETS): FC = $(FPP)
ifeq ($(PARALLEL),yes)
$(FPPTARGETS): FFLAGS += -DPARALLEL
$(FPPTARGETS): F90FLAGS += -DPARALLEL
$(FPPTARGETS): FNOFLAGS += -DPARALLEL
endif

# Pattern rules for compilation
%.o: %.f
	$(FC) -c $(FFLAGS) $(INC) $< -o $@

%.o: %.c
	gcc -std=c99 -O3 -c $< -o $@

%.o: %.f90
	$(FC) -c $(F90FLAGS) $(INC) $< -o $@

# Targets that should not be optimized (f77 only)
$(FNOTARGETS): %.o: %.f
	$(FC) -c $(FNOFLAGS) $(INC) $< -o $@
