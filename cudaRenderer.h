#ifndef __CUDA_RENDERER_H__
#define __CUDA_RENDERER_H__

#ifndef uint
#define uint unsigned int
#endif

#include "circleRenderer.h"


class CudaRenderer : public CircleRenderer {

private:
    Image* image;

    int cells_per_side;
    int mplsSize;

    float* VX;
    float* VY;
    float* color;
    float* colorCopy;
    float* pressures;
    float* pressuresCopy;
    float* VXCopy;
    float* VYCopy;
    float* divergence;
    float* vorticity;

    int* mpls;

    float* cdVX;
    float* cdVY;
    float* cdColor;
    float* cdColorCopy;
    float* cdPressures;
    float* cdPressuresCopy;
    float* cdVXCopy;
    float* cdVYCopy;
    float* cdDivergence;
    float* cdVorticity;
    float* cdImageData;

    int* cdMpls;

    double distanceToSegment(double ax, double ay, double bx, double by, 
        double px, double py, double* fp);
    double distanceToNearestMouseSegment(double px, double py, double *fp, 
            std::pair<double,double>* mouseSegmentVelocity);
    /*int isBoundary(int i, int j);
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
    void applyVorticityForce();*/


   
public:

    CudaRenderer();
    virtual ~CudaRenderer();

    const Image* getImage();

    void setup();
    
    void setNewQuantities(std::vector<std::pair<int, int> > mpls);

    void allocOutputImage(int width, int height);

    void clearImage();

    void render();
};


#endif
