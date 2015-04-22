# Manual Makefile must be replaced with GBS
LIBS = `pkg-config gtk+-2.0 gthread-2.0 --libs` `pkg-config gsl --libs` \
       `echo $$LD_LIBRARY_PATH | sed -e 's/^/-L/;s/:/ -L/g'` \
       -lpgl -lvmath `pkg-config libxml-2.0 --libs`
GTKFLAGS = `pkg-config gtk+-2.0 --cflags` 
QINC = -I$(HOME)/include $(GTKFLAGS) -ggdb \
       `pkg-config libxml-2.0 --cflags` -Wall
OBJS = 3dvl.o draw.o inter.o mask.o mathp.o plot.o projection.o \
       proj_cyl.o proj_sph.o request.o
probez-gtk: $(OBJS)
	gcc $(LIBS) -o probez-gtk $(OBJS)
%.o:%.c
	gcc $< -c $(QINC)
