
CXXFLAGS=`wx-config --cxxflags` -I../.. 

OBJS=\
debug.o \
die.o \
info.o \
pb.o \
status.o \
warn.o

libwxlog.a: $(OBJS)
	ar r ../../libwxlog.a $(OBJS)
