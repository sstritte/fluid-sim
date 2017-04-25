#include <string>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

#include "cudaRenderer.h"
#include "image.h"
#include "util.h"
#include "cycleTimer.h"

/*
 * THIS IS BASICALLY THE SAME AS refRenderer RIGHT NOW
 *
 */

CudaRenderer::CudaRenderer() {
    image = NULL;
    mousePressedLocation = NULL; 
}

CudaRenderer::~CudaRenderer() {

    if (image) {
        delete image;
    }

}

const Image*
CudaRenderer::getImage() {
    return image;
}

void
CudaRenderer::setup() {
    // Nothing for now because not actually using CUDA yet
}

void 
CudaRenderer::setMousePressedLocation(int* mpl) {
    mousePressedLocation = mpl;
}
 
void 
CudaRenderer::setNewQuantities(double* vxs, double* vys, double* ps) {
    // nothing yet
}
// allocOutputImage --
//
// Allocate buffer the renderer will render into.  Check status of
// image first to avoid memory leak.
void
CudaRenderer::allocOutputImage(int width, int height) {

    if (image)
        delete image;
    image = new Image(width, height);
}

// clearImage --
//
// Clear's the renderer's target image.  The state of the image after
// the clear depends on the scene being rendered.
void
CudaRenderer::clearImage() {
    image->clear(1.f,1.f,1.f,1.f);
}

void
CudaRenderer::render() {
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
