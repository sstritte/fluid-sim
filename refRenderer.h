#ifndef __REF_RENDERER_H__
#define __REF_RENDERER_H__

#include "circleRenderer.h"


class RefRenderer : public CircleRenderer {

private:

    Image* image;
    int* mousePressedLocation; 

    float** velocitiesX;
    float** velocitiesY;
    int cells_per_side;
    float*** color;
    float*** colorCopy;
    float** pressures;
    float** advectionCopy;

    void advectColor();
    void advectColorForward();
    void advectColorBackward();
    void advectQuantity(float** q);
    void advectForward(float** q);
    void advectBackward(float** q);
    void applyPressure();
public:

    RefRenderer();
    virtual ~RefRenderer();

    const Image* getImage();

    void setup();
    
    void setMousePressedLocation(int* mpl);

    void setNewQuantities(double* vxs, double* vys, double* ps);

    void allocOutputImage(int width, int height);

    void clearImage();

    void render();
};


#endif
