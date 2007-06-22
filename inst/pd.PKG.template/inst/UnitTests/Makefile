TOP=../../..
RUNIT_DIR=$(shell pwd)
INST_PATH=$(shell dirname ${RUNIT_DIR})
PKG_PATH=$(shell dirname ${INST_PATH})
RCMD=R25

all: install test

install:
	cd ${TOP}; \
	${RCMD} CMD INSTALL --clean ${PKG_PATH} && \
	cd ${RUNIT_DIR}

test:
	${RCMD} --slave < runalltests.R

%.R: FORCE
	./runfile.sh $@

FORCE:
