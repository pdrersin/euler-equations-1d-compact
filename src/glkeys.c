#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "glinterface.h"
#include "gldisplay.h"
#include "common.h"
#include "compact.h"
#include "domain.h"
#include "glkeys.h"

extern int queueAnother;
extern staggeredGrid* grid;
extern int displayIndex;
extern real min, max;
extern int filterOn;
extern real k;

void keyb(unsigned char key, int x, int y) {
    switch(key) { //respond to keys
	case 'K':
	    grid->t = 0.0;
	    k *= pow(10.0,1.0/5.0);
	    InitializeStaggeredGrid( grid );
	    CalculatePrimitiveVariables( grid );
	    update_plan();
	    update_fft();
	    remax();
	    disp();
	    break;

	case 'k':
	    grid->t = 0.0;
	    k /= pow(10.0,1.0/5.0);
	    InitializeStaggeredGrid( grid );
	    CalculatePrimitiveVariables( grid );
	    update_plan();
	    update_fft();
	    remax();
	    disp();
	    break;
	    
	case 'F':
	    filterOn++;
	    disp();
	    break;
	case 'f':
	    filterOn=MAX(0,filterOn-1);
	    disp();
	    break;
	case 'c':
	    fprintf(stderr,"CFL=%f\n",GetCFLNumber( grid ));
	break;
	case 'm':
	    remax();
	break;
	case 'r':
	    grid->t = 0.0;
	    InitializeStaggeredGrid( grid );
	    CalculatePrimitiveVariables( grid );
	    update_plan();
	    update_fft();
	    remax();
	    disp();
	break;
	case '/':
	    displayIndex = (displayIndex+1)%grid->nVar;
	    update_plan();
	    update_fft();
	    remax();
	    disp();
	break;
	case ' ':
	    iterate = 1-iterate;
	    if(iterate) {
		queueAnother=1;
		redisplay(0);
	    } else {
		queueAnother=0;
	    }
	break;

	case ')': //full screen without window border
	    glutFullScreen();
	    full=1;
	break;

	case '0': //full screen with window border
	    if(full==1) {
		glutReshapeWindow(500,500);
		glutPositionWindow(100,100);
	    } else {
		glutReshapeWindow(1024,768);
		glutPositionWindow(0,-10);
	    }
	    full=!full;
	break;

	case 'q': //quit
	    FreeStaggeredGrid( grid );
	    glutDestroyWindow(win);
	    exit(0);
	break;
    }
}
