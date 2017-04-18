#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <vector>

#include "refRenderer.h"
#include "image.h"
#include "util.h"

RefRenderer::RefRenderer() {
    image = NULL;
    mousePressedLocation = NULL;
}

RefRenderer::~RefRenderer() {

    if (image) {
        delete image;
    }
}

const Image*
RefRenderer::getImage() {
    return image;
}

void
RefRenderer::setup() {
    // nothing to do here
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

