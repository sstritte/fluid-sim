#ifndef __REF_RENDERER_H__
#define __REF_RENDERER_H__

#include "circleRenderer.h"
#include <list>
#include <utility> // std::pair

class RefRenderer : public CircleRenderer {

private:

    bool isMouseDown;

    Image* image;
    int* mousePressedLocation;
    std::vector<std::pair<int,int> > mousePressedLocations; 

    float** velocitiesX;
    float** velocitiesY;
    int cells_per_side;
    float*** color;
    float*** colorCopy;
    float** pressures;
    float** advectionCopyX;
    float** advectionCopyY;
    float** advectionCopy;
    float** divergence;
    float** vorticity;

    double distanceToSegment(double ax, double ay, double bx, double by, 
        double px, double py, double* fp);
    double distanceToNearestMousePoint(double px, double py, double *fp);
    int isBoundary(int i, int j);
    void advectColor();
    void advectColorForward();
    void advectColorBackward();
    void advectQuantity(float** q);
    void advectForward(float** q);
    void advectBackward(float** q);
    void applyPressure();
    void advectVelocityForward();
    void advectVelocityBackward();
    void applyDivergence();
    void pressureSolve();
    void pressureGradient();
    void applyVorticity();
    void applyVorticityForce();
public:

    RefRenderer();
    virtual ~RefRenderer();

    const Image* getImage();

    void setup();
    
    void setNewQuantities(double* vxs, double* vys, int* mpl, bool mouseDown);

    void allocOutputImage(int width, int height);

    void clearImage();

    void render();
};


#endif
