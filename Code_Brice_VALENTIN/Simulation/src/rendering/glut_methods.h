#ifndef GLUT_METHODS_H
#define GLUT_METHODS_H

#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

void reshape(int w, int h);
void display();
void idle();
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyboard(unsigned char key, int x, int y);
void special(int key, int x, int y);
void close();
void menu(int item);
void drawTextFont(float size, bool serif = false);
void drawTextAlign(float h, float v);
void drawTextPos(float x, float y);
float drawText(const char* str);
float drawTextF(const char* fmt, ...);
void init_glut(int* argc, char** argv);
void setup_glut();
void run_glut();

extern int window_width;
extern int window_height;
extern double fps;

#endif
