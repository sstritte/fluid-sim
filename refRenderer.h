#ifndef __REF_RENDERER_H__
#define __REF_RENDERER_H__

#include "circleRenderer.h"
#include <list>
#include <utility> // std::pair
#include <vector>

#include "platformgl.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
class RefRenderer : public CircleRenderer {

private:

    Image* image;
    std::vector<std::pair<int,int> > mousePressedLocations; 

    float** velocitiesX;
    float** velocitiesY;
    int cells_per_side;
    float*** color;
    float*** colorCopy;
    float** pressures;
    float** advectionCopyX;
    float** advectionCopyY;
    float** divergence;
    float** vorticity;

    double distanceToSegment(double ax, double ay, double bx, double by, 
        double px, double py, double* fp);
    double distanceToNearestMouseSegment(double px, double py, double *fp, 
            std::pair<double,double>* mouseSegmentVelocity);
    int isBoundary(int i, int j);
    void advectColor();
    void advectColorForward();
    void advectColorBackward();
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
    
    void setNewQuantities(std::vector<std::pair<int, int> > mpls);

    void allocOutputImage(int width, int height);

    void clearImage(cudaSurfaceObject_t s);

    void render(cudaSurfaceObject_t s);
};


#endif
