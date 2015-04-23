# Manual Makefile must be replaced with GBS
LIBS =  `pkg-config gsl --libs` \
       `echo $$LD_LIBRARY_PATH | sed -e 's/^/-L/;s/:/ -L/g'` \
       -lpgl -lvmath `pkg-config libxml-2.0 --libs`
GTKFLAGS = `pkg-config gtk+-2.0 --cflags` 
GTKLIBS = `pkg-config gtk+-2.0 gthread-2.0 --libs`
QINC = -I$(HOME)/include -ggdb \
       `pkg-config libxml-2.0 --cflags` -Wall
GOBJS = window.o dialogs.o
LOBJS = 3dvl.o mask.o mathp.o plot.o projection.o \
       proj_cyl.o proj_sph.o request.o
probez-gtk: $(GOBJS) libprobez.so
	gcc $(LIBS) -o probez-gtk $(GOBJS) -lprobez -L. 
libprobez.so: $(LOBJS)
	gcc --shared $(LOBJS) -o libprobez.so $(LIBS)
%.o:%.c
	gcc $< -c $(QINC) -fPIC -Wextra
window.o:window.c
	gcc window.c -c $(GTKFLAGS) $(QINC)
dialogs.o:dialogs.c
	gcc dialogs.c -c $(GTKFLAGS) $(QINC)
