#include <algorithm>
#include <math.h>

#include "circleRenderer.h"
#include "cycleTimer.h"
#include "image.h"
#include "platformgl.h"


void renderPicture();


static struct {
    int width;
    int height;
    bool printStats;
    double lastFrameTime;
    int* mousePressedLocation;
    int prevMouseX;
    int prevMouseY;
    double prevMousePressTime;
    double* newVelocitiesX;
    double* newVelocitiesY;
    double* newPressures;
    CircleRenderer* renderer;

} gDisplay;

// handleReshape --
//
// Event handler, fired when the window is resized
void
handleReshape(int w, int h) {
    gDisplay.width = w;
    gDisplay.height = h;
    glViewport(0, 0, gDisplay.width, gDisplay.height);
    glutPostRedisplay();
}

void
handleDisplay() {

    // simulation and rendering work is done in the renderPicture
    // function below

    renderPicture();

    // the subsequent code uses OpenGL to present the state of the
    // rendered image on the screen.

    const Image* img = gDisplay.renderer->getImage();
    //int* mousePressedLocation = gDisplay.mousePressedLocation;

    int width = std::min(img->width, gDisplay.width);
    int height = std::min(img->height, gDisplay.height);

    glViewport(0, 0, gDisplay.width, gDisplay.height);

    glDisable(GL_DEPTH_TEST);
    glClearColor(0.f, 0.f, 0.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.f, gDisplay.width, 0.f, gDisplay.height, -1.f, 1.f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // copy image data from the renderer to the OpenGL
    // frame-buffer.  This is inefficient solution is the processing
    // to generate the image is done in CUDA.  An improved solution
    // would render to a CUDA surface object (stored in GPU memory),
    // and then bind this surface as a texture enabling it's use in
    // normal openGL rendering
    glRasterPos2i(0, 0);
    glDrawPixels(width, height, GL_RGBA, GL_FLOAT, img->data);//newColors);

    double currentTime = CycleTimer::currentSeconds();

    if (gDisplay.printStats)
        printf("%.2f ms\n", 1000.f * (currentTime - gDisplay.lastFrameTime));

    gDisplay.lastFrameTime = currentTime;

    glutSwapBuffers();
    glutPostRedisplay();
}


// handleKeyPress --
//
// Keyboard event handler
void
handleKeyPress(unsigned char key, int x, int y) {

    switch (key) {
    case 'q':
    case 'Q':
        exit(1);
        break;
    }
}

void handleMouseClick(int button, int state, int x, int y) {
    if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        //printf("MOUSE DOWN\n");
        gDisplay.prevMousePressTime = CycleTimer::currentSeconds();
        gDisplay.prevMouseX = x;
        gDisplay.prevMouseY = y;
    } else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        //printf("MOUSE UP\n");
        //gDisplay.prevMousePressTime = -1; 
    } 

}

// handleMouseMove --
// 
// Mouse event handler
void
handleMouseMove(int x, int y) {
    double mousePressTime = CycleTimer::currentSeconds();
    //double pixelDist = sqrt(pow(x - gDisplay.prevMouseX, 2) + pow(y - gDisplay.prevMouseY, 2)); 
    double distX = x - gDisplay.prevMouseX;
    double distY = y - gDisplay.prevMouseY;
    double d = sqrt(distX * distX + distY * distY);

    int prevIndex = 
        (gDisplay.height - gDisplay.prevMouseY - 1) * gDisplay.width 
        + gDisplay.prevMouseX;
    if (0 <= prevIndex && prevIndex < gDisplay.height * gDisplay.width) {
        gDisplay.mousePressedLocation[prevIndex] = 1;
        gDisplay.newVelocitiesX[prevIndex] = distX / d;
        gDisplay.newVelocitiesY[prevIndex] = -distY / d; //bc screen is upside down
        gDisplay.newPressures[prevIndex] = 1.0;
    }

    /*printf("moved %f pixels, %f x, %f y, from %f ms ago\n", 
            pixelDist, distX, distY, 
            1000 * (mousePressTime - gDisplay.prevMousePressTime));*/
    gDisplay.prevMouseX = x;
    gDisplay.prevMouseY = y;
    gDisplay.prevMousePressTime = mousePressTime;
}

// renderPicture --
//
// At the reall work is done here, not in the display handler
void
renderPicture() {

    double startTime = CycleTimer::currentSeconds();

    // clear screen
    gDisplay.renderer->clearImage();

    double endClearTime = CycleTimer::currentSeconds();

    // SEND RELEVANT INFO BEFORE RENDERING
    //gDisplay.renderer->setMousePressedLocation(gDisplay.mousePressedLocation);
    memset(gDisplay.mousePressedLocation, 0, sizeof(int) * gDisplay.width * gDisplay.height);

    gDisplay.renderer->setNewQuantities(gDisplay.newVelocitiesX, gDisplay.newVelocitiesY, gDisplay.newPressures);
    memset(gDisplay.newVelocitiesX, 0, sizeof(double) * gDisplay.width * gDisplay.height);
    memset(gDisplay.newVelocitiesY, 0, sizeof(double) * gDisplay.width * gDisplay.height);

    // RENDER THE PARTICLES INTO THE IMAGE
    gDisplay.renderer->render();
    double endRenderTime = CycleTimer::currentSeconds();

    if (gDisplay.printStats) {
        printf("Clear:    %.3f ms\n", 1000.f * (endClearTime - startTime));
        printf("Render:   %.3f ms\n", 1000.f * (endRenderTime - endClearTime));
    }
}

void
startRendererWithDisplay(CircleRenderer* renderer) {

    // setup the display
    const Image* img = renderer->getImage();

    gDisplay.renderer = renderer;
    gDisplay.printStats = true;
    gDisplay.lastFrameTime = CycleTimer::currentSeconds();
    gDisplay.width = img->width;
    gDisplay.height = img->height;
    gDisplay.mousePressedLocation = new int[img->width * img->height];
    memset(gDisplay.mousePressedLocation, 0, sizeof(int) * img->width * img->height);
    printf("before\n");
    gDisplay.newVelocitiesX = new double[img->width * img->height];
    memset(gDisplay.newVelocitiesX, 0, sizeof(double) * img->width * img->height);
    gDisplay.newVelocitiesY = new double[img->width * img->height];
    memset(gDisplay.newVelocitiesY, 0, sizeof(double) * img->width * img->height);
    gDisplay.newPressures = new double[img->width * img->height];
    memset(gDisplay.newPressures, 0, sizeof(double) * img->width * img->height);
    printf("after\n");

    // configure GLUT
    glutInitWindowSize(gDisplay.width, gDisplay.height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow("CMU 15-418 Final Project - Fluid Simulator");
    glutDisplayFunc(handleDisplay);
    glutKeyboardFunc(handleKeyPress);
    glutMouseFunc(handleMouseClick);
    glutMotionFunc(handleMouseMove);
    glutMainLoop();
}
