#ifndef __CIRCLE_RENDERER_H__
#define __CIRCLE_RENDERER_H__

struct Image;

class CircleRenderer {

public:

    virtual ~CircleRenderer() { };

    virtual const Image* getImage() = 0;

    virtual void setup() = 0;
    
    virtual void setNewQuantities(double* vxs, double* vys, int* mpl, bool mouseDown) = 0;

    virtual void allocOutputImage(int width, int height) = 0;

    virtual void clearImage() = 0;

    virtual void render() = 0;

};


#endif
