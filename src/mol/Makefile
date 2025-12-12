################################################################################
#                                                                              #
#  Libcutils - library of C functions                                          #
#  Version: 3.4                                                                #
#                                                                              #
#  Date:    28/01/2017                                                         #
#                                                                              #
#  Author:  Zdenek Futera                                                      #
#                                                                              #
#  Address: University College London (UCL)                                    #
#           Department of Physics and Astronomy                                #
#           Gower Street, London WC1E 6BT                                      #
#           United Kingdom                                                     #
#                                                                              #
#  E-Mail:  z.futera@ucl.ac.uk                                                 #
#                                                                              #
################################################################################

################################################################################

.PHONY: all
.PHONY: clean

################################################################################

LIBHOME = ../..

LIBDIR = ../../lib
INCDIR = ../../include
OBJDIR = ../../obj/mol
DEPDIR = ../../dep/mol

LIBFILE = ${LIBDIR}/libmol.a

AR      = ar
AFLAGS  = rcs
CC      = gcc
FOPT    = -O2

FLAGS   = -g -c -std=gnu99 -Wall ${FOPT} -I${INCDIR}

HEADERS = ${INCDIR}/mol

################################################################################

SRCS := $(shell find . -name "*.c" | sed 's/^\.\/\(.*\)/\1/')
OBJS := $(addprefix $(OBJDIR)/,$(patsubst %.c,%.o,$(SRCS)))
DEPS := $(addprefix $(DEPDIR)/,$(patsubst %.c,%.d,$(SRCS)))

################################################################################

all: ${LIBFILE}

$(OBJDIR)/%.o: %.c
	@echo "CC  $<"
	@mkdir -p $(@D)
	@$(CC) $(FLAGS) $< -o $@

$(DEPDIR)/%.d: %.c
	@mkdir -p $(@D)
	@$(CC) $(FLAGS) -MM $< | \
	sed "s,\($(*F)\)\.o[ :]*,$(OBJDIR)/$*.o $@ : ," > $@

$(LIBDIR):
	@mkdir -p $(LIBDIR)

${LIBFILE}: $(OBJS) | $(LIBDIR)
	@echo "LD  libmol"
	@${AR} ${AFLAGS} ${LIBFILE} ${OBJS}

################################################################################

clean:
	rm -rf ${OBJDIR} ${DEPDIR} ${LIBFILE}

################################################################################

-include $(DEPS)

################################################################################
