#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>

#include "gldisplay.h"
#include "common.h"
#include "compact.h"
#include "domain.h"
#include "glinterface.h"
#include "deriv.h"

#define c sqrt(1.4*101325.0/1.225)

inline double floatmod(double a, double base) {
    double result = a - (int)(a/base-1)*base;
    return result;
}

fftw_plan plan;

extern real k;
extern staggeredGrid* grid;
double *fftinbuf = NULL;
fftw_complex *fftbuf = NULL;
double *fftabs = NULL;
double fftmin, fftmax;

int queueAnother=0;
int displayIndex = 0;
real min=-10.0, max=10.0;
int frame=0;
int filterOn=0;

void update_plan() {
    if(fftinbuf) free(fftinbuf);
    if(fftbuf) free(fftbuf);
    if(fftabs) free(fftabs);
    fftinbuf = malloc(sizeof(double)*grid->nx);
    fftbuf = malloc(sizeof(fftw_complex)*(grid->nx/2+1));
    fftabs = malloc(sizeof(double)*(grid->nx/2+1));
    plan = fftw_plan_dft_r2c_1d(grid->nx, fftinbuf, fftbuf, FFTW_ESTIMATE);
}

double fft1;
void update_fft() {
    int i;
    memcpy(fftinbuf, grid->W[displayIndex], sizeof(double)*grid->nx);
    double sum=0.0;
    for(i=0; i<grid->nx; i++) {
	sum += fftinbuf[i];
    }
    sum /= grid->nx;
    for(i=0; i<grid->nx; i++) {
	fftinbuf[i] -= sum;
    }
    fftw_execute(plan);
    for(i=0; i<grid->nx/2+1; i++) {
	fftabs[i] = log10(cabs(fftbuf[i]));
    }
}

double calcDB() {
    int i;
    double mean=0.0;
    for(i=0; i<grid->nx; i++) {
	mean += grid->W[2][i];
    }
    mean /= grid->nx;
    double var=0.0;
    for(i=0; i<grid->nx; i++) {
	var += SQR(grid->W[2][i]-mean);
    }
    var /= grid->nx;
    double pref = 2.0e-5; //Pa
    return 10.0*log10(var/SQR(pref));
}

void redisplay(int value) {
    queueAnother=1;
    disp();
}

void remax() {
    int i;
    GetArrayMinMax(grid->W[displayIndex],grid->nx,&min,&max);
    fprintf(stderr,"min = %f, max = %f\n",min,max);
    if(max-min < 1.0e-10) {
	max = min+1.0e-5;
	min = max-2.0e-5;
    }
    fftmin = 1e10;
    fftmax = -1e10;
    for(i=0; i<grid->nx/2+1; i++) {
	fftmin = MIN(fftmin,fftabs[i]);
	fftmax = MAX(fftmax,fftabs[i]);
    }
    fftmin = MAX(fftmin,-2*20);
    fft1 = fftabs[1];
    disp();
}

void disp(void) {
    int i,j;

    if(iterate) {
	for(i=0; i<10; i++) {
	    grid->RK( grid );
	    if(filterOn)
		if((frame++)%filterOn==0)
		    FilterAll( grid );
	}
	CalculatePrimitiveVariables( grid );
	update_fft();
    }

    real mmm, gmmm, vxmin, x, dx, y, dy;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    
    glLineWidth(1.0);
    glBegin( GL_LINES );
    dy = dx = 0.1;
    real gxmin = 0.0;
    real gxmax = 2.0;
    real gymin = 0.0;
    real gymax = 1.0;
    for(x=gxmin;x<gxmax+dx/2;x+=dx) {
	if( ((int)((x+0.5*dx)/dx)) % (int)(1.0/dx) == 0) {
	    glColor3f( 0.8, 0.8, 0.8 );
	} else {
	    glColor3f( 0.2, 0.2, 0.2 );
	}
	glVertex2f( x, gymin );
	glVertex2f( x, gymax );
    }
    for(x=gymin;x<gymax+dx/2;x+=dx) {
	if( ((int)((x+0.5*dx)/dx)) % (int)(1.0/dx) == 0) {
	    glColor3f( 0.8, 0.8, 0.8 );
	} else {
	    glColor3f( 0.2, 0.2, 0.2 );
	}
	glVertex2f( gxmin, x );
	glVertex2f( gxmax, x );
    }
    glEnd();

    glColor3f(0.5,0.5,0.5);
    float xpixel = 1.0/(float)width*(xmax-xmin);
    float ypixel = 1.0/(float)height*(ymax-ymin);
    char buf[128];

    int xpow = MAX((int)log10(ABS(gxmin)), (int)log10(ABS(gxmax)));
    int ypow = MAX((int)log10(ABS(min)), (int)log10(ABS(max)));
    real xsc = pow(10.0, xpow);
    real ysc = pow(10.0, ypow);

    for(x=gxmin;x<gxmax+dx/2;x+=dx*2) {
	sprintf(buf,"%1.2f",x/xsc);
	j=0;
	int len=0;
	while(buf[j++]!='\0') len++;
	glRasterPos2f( x-xpixel*9*len/2, 0.0-ypixel*19 );
	j=0;
	do {
	    glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
	} while(buf[j++]!='\0');
    }
    x-=dx/2;
    glRasterPos2f( x-xpixel*9*2, 0.0-ypixel*10 );
    sprintf(buf,"* (10^%i)",xpow);
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');


    for(y=gymin;y<gymax+dy/2;y+=dy*2) {
	sprintf(buf,"%1.6f",(min + (max-min)*y)/ysc);
	j=0;
	int len=0;
	while(buf[j++]!='\0') len++;
	glRasterPos2f( 0.0-xpixel*(len+1)*9, y - ypixel*7 );
	j=0;
	do {
	    glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
	} while(buf[j++]!='\0');
    }
    y -= dy;
    glRasterPos2f( 0.0-xpixel*5*9, y - ypixel*7 );
    sprintf(buf,"    * (10^%i)",ypow);
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');

    glColor3f(1.0,1.0,1.0);
    switch(displayIndex) {
	case 0: sprintf(buf,"Density, kg/m^3"); break;
	case 1: sprintf(buf,"Velocity, m/s"); break;
	case 2: sprintf(buf,"Pressure, Pa"); break;
    }
    j=0;
    int len=0;
    while(buf[j++]!='\0') len++;
    glRasterPos2f( 1.0-xpixel*len/2*9, 1.0 + ypixel*10 );
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');
    
    glColor3f(0.5,0.5,0.5);
    switch(filterOn) {
	case 0: sprintf(buf,"Filter = OFF"); break;
	default: sprintf(buf,"Filter = %i ",filterOn); break;
    }
    j=0;
    len=0;
    while(buf[j++]!='\0') len++;
    glRasterPos2f( 2.0-xpixel*len*9, 1.0 + ypixel*10 );
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');
    
    sprintf(buf,"%6.3f dB",calcDB());
    j=0;
    len=0;
    while(buf[j++]!='\0') len++;
    glRasterPos2f( 2.0-xpixel*len*9, 1.0 + ypixel*25 );
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');
    
    sprintf(buf,"%6.3f periods",grid->t*c);
    j=0;
    len=0;
    while(buf[j++]!='\0') len++;
    glRasterPos2f( 2.0-xpixel*len*9, 1.0 + ypixel*40 );
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');
    


    //if(1) {
	//GetArrayMinMax(grid->W[displayIndex],grid->nx,&min,&max);
    //} else {
	//min = 1.0;
	//max = 2.0;
    //}
    vxmin = grid->x[0];
    gmmm = 1.0;
    mmm = 1.0/(max-min);

    glLineWidth(1.0);
    glBegin( GL_LINES );
    glColor3f(1.0,0.0,0.0);
    real x1, x2;
    for(i=0;i<grid->nx-1;i++) {
	x1 = floatmod((grid->xm[i  ]-vxmin)*gmmm - grid->t*c,2.0);
	x2 = floatmod((grid->xm[i+1]-vxmin)*gmmm - grid->t*c,2.0);
	if(ABS(x1-x2) < 0.1) {
	    glVertex2f( x1, (grid->W[displayIndex][i  ]-min)*mmm );
	    glVertex2f( x2, (grid->W[displayIndex][i+1]-min)*mmm );
	}
	x1 = floatmod((grid->xm[i  ]-vxmin)*gmmm+1 - grid->t*c,2.0);
	x2 = floatmod((grid->xm[i+1]-vxmin)*gmmm+1 - grid->t*c,2.0);
	if(ABS(x1-x2) < 0.1) {
	    glVertex2f( x1, (grid->W[displayIndex][i  ]-min)*mmm );
	    glVertex2f( x2, (grid->W[displayIndex][i+1]-min)*mmm );
	}
    }
    glEnd();
    glPointSize( 3.0 );
    glBegin( GL_POINTS );
    glColor3f(0.0,0.0,1.0);
    for(i=0;i<grid->nx;i++) {
	glVertex2f( floatmod((grid->xm[i  ]-vxmin)*gmmm - grid->t*c,2.0), (grid->W[displayIndex][i  ]-min)*mmm );
	glVertex2f( floatmod((grid->xm[i  ]-vxmin)*gmmm+1 - grid->t*c,2.0), (grid->W[displayIndex][i  ]-min)*mmm );
    }
    glEnd();

    glLineWidth(1.0);
    glBegin( GL_LINES );
    dy = dx = 0.1;
    gxmin = 0.0;
    gxmax = 2.0;
    gymin = -1.5;
    gymax = -0.5;
    real lsxmin = 1.0;
    real lsxmax = 100.0;
    real llsxmin = log10(lsxmin);
    real llsxmax = log10(lsxmax);
    float lsx;
    int tpmin = (int)(llsxmin + 1.0e-5);
    int tpmax = (int)(llsxmax + 1.0e-5);
    for(i=tpmin+1; i<=tpmax; i++) {
	for(x=0.1; x<0.95; x+=0.1) {
	    lsx = gxmin + (gxmax-gxmin)*(log10(x*pow(10,i)) - llsxmin)/(llsxmax-llsxmin);
	    if( x < 0.15 ) {
		glColor3f( 0.8, 0.8, 0.8 );
	    } else {
		glColor3f( 0.2, 0.2, 0.2 );
	    }
	    glVertex2f( lsx, gymin );
	    glVertex2f( lsx, gymax );
	}
    }
    real lsymin = 0.0001;
    real lsymax = 1.0;
    real llsymin = log10(lsymin);
    real llsymax = log10(lsymax);
    float lsy;
    tpmin = (int)(llsxmin + 1.0e-5);
    tpmax = (int)(llsxmax + 1.0e-5);
    for(i=tpmin-3; i<=tpmax-2; i++) {
	for(x=0.1; x<0.95; x+=0.1) {
	    lsy = gymin + (gymax-gymin)*(log10(x*pow(10,i)) - llsymin)/(llsymax-llsymin);
	    if( x < 0.15 ) {
		glColor3f( 0.8, 0.8, 0.8 );
	    } else {
		glColor3f( 0.2, 0.2, 0.2 );
	    }
	    glVertex2f( gxmin, lsy );
	    glVertex2f( gxmax, lsy );
	}
    }
    glColor3f( 0.8, 0.8, 0.8 );
    glVertex2f( 0.0, -0.5 );
    glVertex2f( 2.0, -0.5 );
    glVertex2f( 2.0, -1.5 );
    glVertex2f( 2.0, -0.5 );
    glEnd();

    glLineWidth(2.0);
    glColor3f(1.0,0.5,0.0);
    glBegin( GL_LINES );
    int nf = grid->nx/2+1;
    mmm = (gymax-gymin)/(fftmax-fftmin);
    double gdx = grid->x[1] - grid->x[0];
    for(i=0; i<nf-1; i++) {
	glVertex2f( gxmin + (gxmax-gxmin)*(log10(1.0*(i  )/(nf-1)*0.5/gdx)-llsxmin)/(llsxmax-llsxmin), gymin + (fftabs[i  ]-fft1-llsymin)/(llsymax-llsymin)*(gymax-gymin));
	glVertex2f( gxmin + (gxmax-gxmin)*(log10(1.0*(i+1)/(nf-1)*0.5/gdx)-llsxmin)/(llsxmax-llsxmin), gymin + (fftabs[i+1]-fft1-llsymin)/(llsymax-llsymin)*(gymax-gymin));
    }
    glEnd();

    xpow = MAX((int)log10(ABS(gxmin)), (int)log10(ABS(gxmax)));
    ypow = MAX((int)log10(ABS(fftmin)), (int)log10(ABS(fftmax)));
    xsc = pow(10.0, xpow);
    ysc = pow(10.0, ypow);

    /*glColor3f(0.5,0.5,0.5);
    for(y=gymin;y<gymax+dy/2;y+=dy*2) {
	sprintf(buf,"%3i",(int)((fftmin + (fftmax-fftmin)*(y-gymin)/(gymax-gymin))));
	j=0;
	int len=0;
	while(buf[j++]!='\0') len++;
	glRasterPos2f( 0.0-xpixel*(len+1)*9, y - ypixel*7 );
	j=0;
	do {
	    glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
	} while(buf[j++]!='\0');
    }*/

    glColor3f(0.5,0.5,0.5);
    sprintf(buf,"Spectrum");
    j=0;
    len=0;
    while(buf[j++]!='\0') len++;
    glRasterPos2f( 1.0-xpixel*len/2*9, -0.5 + ypixel*10 );
    j=0;
    do {
	glutBitmapCharacter(GLUT_BITMAP_9_BY_15,buf[j]);
    } while(buf[j++]!='\0');
    
    glColor3f(0.5,0.5,0.5);


    glutSwapBuffers();

    if(iterate==1 && queueAnother==1) {
	queueAnother=0;
	glutTimerFunc(3, redisplay, 0);
    }

}

