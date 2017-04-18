#ifndef __REF_RENDERER_H__
#define __REF_RENDERER_H__

#include "circleRenderer.h"


class RefRenderer : public CircleRenderer {

private:

    Image* image;
    int* mousePressedLocation; 

public:

    RefRenderer();
    virtual ~RefRenderer();

    const Image* getImage();

    void setup();
    
    void setMousePressedLocation(int* mpl);

    void allocOutputImage(int width, int height);

    void clearImage();

    void render();
};


#endif
