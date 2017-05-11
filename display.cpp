#include <algorithm>
#include <math.h>
#include <utility>
#include <vector>

#include "circleRenderer.h"
#include "cycleTimer.h"
#include "image.h"
#include "platformgl.h"

int number;
void renderPicture();


static struct {
    int width;
    int height;
    bool printStats;
    double lastFrameTime;
    std::vector<std::pair<int,int> > mousePressedLocations;

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
    if (number == 100 || number == 150) {
        printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    } else {
        printf("NUMBER %d ",number);
    }
    number++;
    // simulation and rendering work is done in the renderPicture
    // function below

    renderPicture();

    // the subsequent code uses OpenGL to present the state of the
    // rendered image on the screen.

    const Image* img = gDisplay.renderer->getImage();

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
        //gDisplay.isMouseDown = true;
        std::pair<int, int> coords = std::make_pair(x,gDisplay.height - y);
        gDisplay.mousePressedLocations.push_back(coords);
    } else if(button == GLUT_LEFT_BUTTON && state == GLUT_UP) {
        //gDisplay.isMouseDown = false;
    } 
}

// handleMouseMove --
// 
// Mouse event handler
void
handleMouseMove(int x, int y) {
    //printf("handleMouseMove %d, %d\n", x,y);
    int index = (gDisplay.height - y - 1) * gDisplay.width + x;
    if (0 <= index && index < gDisplay.height * gDisplay.width) {
        std::pair<int, int> coords = std::make_pair(x,gDisplay.height - y);
        gDisplay.mousePressedLocations.push_back(coords);
    }
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
    gDisplay.renderer->setNewQuantities(gDisplay.mousePressedLocations);

    gDisplay.mousePressedLocations.clear();

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

    // configure GLUT
    glutInitWindowSize(gDisplay.width, gDisplay.height);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow("CMU 15-418 Final Project - Fluid Simulator");
    glutDisplayFunc(handleDisplay);
    glutKeyboardFunc(handleKeyPress);
    glutMouseFunc(handleMouseClick);
    glutMotionFunc(handleMouseMove);
    number = 0;
    glutMainLoop();
}
