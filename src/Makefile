BINDIR = 	../bin
PROD =		$(BINDIR)/compactSchemeTest
LIB =		-lm -lrrnum -lgl -lglu -lglut -lfftw3 #-lgsl -lcblas 
LIBDIR =	-L/usr/local/lib -L/opt/X11/lib -L/sw/lib 
INCDIR =	-I/usr/local/include -I/opt/X11/include -I/sw/include 
DBFLAGS =	#-g
OPTFLAGS = 	-O3
CC =		cc
CFLAGS = 	-Wall $(INCDIR) $(DBFLAGS) $(OPTFLAGS)
LD =		cc
LDFLAGS =	-Wall $(LIBDIR) $(LIB) $(DBFLAGS) $(OPTFLAGS)
OBJ =		main.o \
		rk.o \
		glinterface.o \
		gldisplay.o \
		glkeys.o \
		common.o \
		airprop.o \
		initialconditions.o \
		compact.o \
		domain.o \
		array.o \
		deriv.o
SRCS = 		$(OBJ:.o=.c)

$(PROD): $(OBJ)
	$(LD) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -c $< $(CFLAGS) -o $@

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean:
	rm -f *.o $(PROD)
run:
	./$(PROD)
debug:
	gdb $(PROD)
