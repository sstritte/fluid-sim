#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "refRenderer.h"
#include "image.h"
#include "util.h"

#define CELL_DIM 16 // a grid cell is 10x10 pixels

RefRenderer::RefRenderer() {
    image = NULL;
    mousePressedLocation = NULL;
    p = NULL;
    vx = NULL;
    vy = NULL;
}

RefRenderer::~RefRenderer() {

    if (image) delete image;
    if (mousePressedLocation) delete mousePressedLocation;
    if (p) {
        for (int i = 0; i < ny; i++) {
            if (p[i]) delete p[i];
        } 
        delete p;
    }
    if (vx) {
        for (int i = 0; i < ny; i++) {
            if (vx[i]) delete vx[i];
        } 
        delete vx;
    }
    if (vy) {
        for (int i = 0; i < ny + 1; i++) {
            if (vy[i]) delete vy[i];
        } 
        delete vy;
    }

}

const Image*
RefRenderer::getImage() {
    return image;
}

void
RefRenderer::setup() {
   nx = image->width / CELL_DIM;
   ny = image->height / CELL_DIM;
   // p is pressures
   p = new float*[ny];
   for (int i = 0; i < ny; i++) {
       p[i] = new float[nx];
   }
   // vx is left/right velocities
   vx = new float*[ny];
   for (int i = 0; i < ny; i++) {
       vx[i] = new float[nx + 1];
   }   
   // vy is up/down velocities
   vy = new float*[ny + 1];
   for (int i = 0; i < ny + 1; i++) {
       vy[i] = new float[nx];
   }
}

void RefRenderer::setMousePressedLocation(int* mpl) {
    mousePressedLocation = mpl;
}

// allocOutputImage --
//
// Allocate buffer the renderer will render into.  Check status of
// image first to avoid memory leak.
void
RefRenderer::allocOutputImage(int width, int height) {

    if (image)
        delete image;
    image = new Image(width, height);
}

// clearImage --
//
// Clear's the renderer's target image.  
void
RefRenderer::clearImage() {
    image->clear(1.f, 1.f, 1.f, 1.f);
}

void
RefRenderer::render() {

    for (int i = 0; i < 4*image->width*image->height; i+=4) {
        if (mousePressedLocation[i / 4]) {
            image->data[i] = 1.0;
            image->data[i+1] = 1.0;
            image->data[i+2] = 1.0;
            image->data[i+3] = 1.0;
        } else {
            image->data[i] = 0.1798;
            image->data[i+1] = 0.457;
            image->data[i+2] = 0.9063;
            image->data[i+3] = 0.5;
        }
    }
}

