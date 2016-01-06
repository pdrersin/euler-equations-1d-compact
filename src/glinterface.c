#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "glinterface.h"
#include "glkeys.h"



int iterate=0;
int win;
int width; //of window
int height; //of window
int full=0; //fullscreen?
int moveMode = VIEW_TRANSLATE;
float mpx=0.0,mpx0=0.0; //variables for past and present window extents
float mpy=0.0,mpy0=0.0;
float xcen0,xcen=1.0;
float ycen0,ycen=-0.15;
float ymin0,ymin,ymax0,ymax;
float xmin0,xmin,xmax0,xmax;
float dx0,dy0,dx,dy=2.7;

void initGlut(int argc, char** argv) {
    //INITIALIZE GLUT WITH SOME OPTIONS
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE );
    glutInitWindowSize(800,800);
    glutInitWindowPosition(0,0);
    win = glutCreateWindow("2D Window");

    //callback functions
    glutDisplayFunc(disp);    //called every time glut thinks a redraw is necessary
    glutReshapeFunc(reshape); //called when window is reshaped
    glutKeyboardFunc(keyb);   //called when keys pressed
    glutMouseFunc(mouseClick);//called on mouse click
    glutMotionFunc(mouseMove);//called when clicked mouse moves

    fprintf(stdout,"\n\nGlut 2d Demo\n\nby Ricky Reusser\n\n");
    fprintf(stdout,"Left click to move view,\nMiddle click to zoom.\n");
}



void reshape(int w, int h) {
    width=w;
    height=h;

    ymin = ycen - dy/2.0; //recalculate the window extents
    ymax = ycen + dy/2.0;
    dx = dy*(float)w/h;
    xmin = xcen-dx/2.0;
    xmax = xcen+dx/2.0;

    glViewport(0, 0, w, h); //tell OpenGL what kind of view to use
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(xmin,xmax,ymin,ymax); //2-D orthographic
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void mouseClick(int button, int state, int x, int y) {
    //function called only when mouse clicked or unclicked
    mpx = xmin+(xmax-xmin)*((float)x/width);
    mpy = ymin+(ymax-ymin)*(1.0-((float)y/height));
    mpx0 = mpx;   mpy0 = mpy;   xmin0 = xmin; ymin0 = ymin;
    xmax0 = xmax; ymax0 = ymax; xcen0 = xcen; ycen0 = ycen;
    dy0 = dy;
    if(state==GLUT_DOWN) { //if button is down
	switch(button) {
	case GLUT_LEFT_BUTTON:
	    moveMode = VIEW_TRANSLATE;
	break;
	case GLUT_MIDDLE_BUTTON:
	    moveMode = VIEW_ZOOM;
	break;
	}
    }
}

void mouseMove(int x, int y) {
    //function called when mouse moved when clicked
    mpx = xmin0+(xmax0-xmin0)*((float)x/width);
    mpy = ymin0+(ymax0-ymin0)*(1.0-((float)y/height));
    switch(moveMode) {
    case VIEW_TRANSLATE:
	xcen = xcen0 - (mpx-mpx0); //move window center by mouse dx, dy
	ycen = ycen0 - (mpy-mpy0);
    break;
    case VIEW_ZOOM:
	dy = dy0*exp(7.0*(mpy-mpy0)/dy0); //multiply window extents by e^dy
	xcen = mpx0 + (xcen0-mpx0)*dy/dy0; //zoom window about inital mouse click point
	ycen = mpy0 + (ycen0-mpy0)*dy/dy0;
    break;
    }
    reshape(width,height); //reshape and redisplay
    disp();
}
