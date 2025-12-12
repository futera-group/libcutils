################################################################################
#                                                                              #
#  Libcutils - library of C functions                                          #
#                                                                              #
#  Author:  Zdenek Futera                                                      #
#                                                                              #
#  Address: University of South Bohemia                                        #
#           Faculty of Science                                                 #
#           Branisovska 1760, 370 05 Ceske Budejovice                          #
#           Czech Republic                                                     #
#                                                                              #
#  E-Mail:  zfutera@prf.jcu.cz                                                 #
#                                                                              #
################################################################################

################################################################################

.PHONY: all
.PHONY: clean
.PHONY: distr

.PHONY: cmn
.PHONY: qmc
.PHONY: mol
.PHONY: prg

################################################################################

LIBNAME=libcutils

DISTDIR=${LIBNAME}

all: cmn qmc mol prg

cmn:
	$(MAKE) -C src/cmn

qmc:
	$(MAKE) -C src/qmc

mol:
	$(MAKE) -C src/mol

prg:
	$(MAKE) -C src/prg

clean:
	rm -rf dep lib obj ${LIBNAME}.tar.bz2
	rm -rf ${DISTDIR}

distr:
	mkdir ${DISTDIR}
	cp Makefile LICENSE.md README.md ${DISTDIR}/
	cp -r include ${DISTDIR}/
	cp -r src ${DISTDIR}/
	tar cjvf ${DISTDIR}.tar.bz2 ${DISTDIR}/
	rm -rf ${DISTDIR}/

################################################################################
