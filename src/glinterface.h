#ifndef __INTERFACE_H__
#define __INTERFACE_H__

//moveMode:
#define VIEW_TRANSLATE 0
#define VIEW_ZOOM 1

void reshape(int w, int h); //called each time window is reshaped
void mouseClick(int button, int state, int x, int y);
void mouseMove(int x, int y);

void disp(void); //called each time the screen is displayed
void initGlut(int argc, char** argv); 

#endif /*__INTERFACE_H__*/
