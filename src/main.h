#ifndef __MAIN_H__
#define __MAIN_H__

#ifndef __INTERFACE_EXTERNS__
#define __INTERFACE_EXTERNS__
extern int win;
extern int width; //of window
extern int height; //of window
extern int full; //fullscreen?
extern int moveMode;
extern float mpx,mpx0; //variables for past and present window extents
extern float mpy,mpy0;
extern float xcen0,xcen;
extern float ycen0,ycen;
extern float ymin0,ymin,ymax0,ymax;
extern float xmin0,xmin,xmax0,xmax;
extern float dx0,dy0,dx,dy;
#endif /*__INTERFACE_EXTERNS__*/

int main(int argc, char **argv);

#endif
