

####
## 
## Makefile for libutils 
##
####

## All user-defineable parameters are set at the top-level; you should
## not need to modify this makefile.

CC?=gcc
RANLIB?=ranlib
AR?=ar
CP?=cp
RM?=rm

ifeq ($(origin ARCH), undefined)
 ARCH=$(subst -,_,$(subst .,_,$(shell uname -s)))
endif


# Uncomment for wx logging support
WXCONFIG=/c/wxWidgets-2.8.9/wx-config

CFLAGS_OPT?= -O -g

INCLUDE_DIRS= -I.
CFLAGS_REQ?= -finline -D$(ARCH)

CFLAGS= $(INCLUDE_DIRS) $(CFLAGS_REQ) $(CFLAGS_OPT)
LDFLAGS= $(LDFLAGS_REQ) $(LDFLAGS_OPT)

LIB_TXTLOG_OBJS=                    \
	logging/txt/die_impl.o      \
	logging/txt/debug_impl.o    \
	logging/txt/info_impl.o     \
	logging/txt/warn_impl.o     \
	logging/txt/pb_impl.o       \
	logging/txt/status_impl.o

LIB_WXLOG_OBJS=                     \
	logging/wx/die_impl.o       \
	logging/wx/debug_impl.o     \
	logging/wx/info_impl.o      \
	logging/wx/warn_impl.o      \
	logging/wx/status_impl.o    \
	logging/wx/pb_impl.o


OBJS=                               \
        logging/set_locus.o         \
	matrix/create_matrix.o      \
	matrix/free_matrix.o        \
	sorting/dcompare.o          \
	sorting/icompare.o          \
	sorting/fcompare.o          \
	text/split.o                \
	text/strip_comments.o       \
	text/trim.o                 \
	text/getline.o              \
        text/parsefloat.o           \
	text/strndup.o		    \
	text/strnlen.o		    \
	text/strchrnul.o	    \
	text/mempcpy.o		    \
	utils_ver.o                 \
	endian.o                    \
	getopt.o		    \
	getopt1.o		    \
	argp-fmtstream.o	    \
	argp-help.o		    \
	argp-parse.o		    \
	argp.o			    
	
		

LIB_TARGETS = libtxtlog.a libutils.a


# If wxWidgets support is enabled, set things up.
ifneq ($(origin WXCONFIG), undefined)
  CXXFLAGS=`$(WXCONFIG) --cxxflags` $(CFLAGS)
  LIB_TARGETS += libwxlog.a
endif

.PHONY: all
all: $(LIB_TARGETS)

endian_config.h: endian
	./endian > endian_config.h

endian: endian.c
	$(CC) -DINTERROGATE_ENDIAN -o $@ $<

typetests: typetests.c
	$(CC) -o $@ $<

.PHONY: run-typetests
run-typetests: typetests
	./typetests 

libutils.a: $(OBJS)
	$(AR) crv $@ $^
	-$(RANLIB) $@

libtxtlog.a: $(LIB_TXTLOG_OBJS)
	$(AR) crv $@ $^
	-$(RANLIB) $@

libwxlog.a: $(LIB_WXLOG_OBJS)
	$(AR) crv $@ $^
	-$(RANLIB) $@

example: utils_example.c libutils.a libtxtlog.a
	$(CC) -o $@ -g $^ 

.PHONY: clean
clean:
	-$(RM) endian_config.h example endian endian.exe .depend $(OBJS) libutils.a libtxtlog.a libwxlog.a

# Dependencies - automatically generated from the object
# files named in OBJS, converted to C equivalents.
.PHONY: 

