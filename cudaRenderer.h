#ifndef __CUDA_RENDERER_H__
#define __CUDA_RENDERER_H__

#ifndef uint
#define uint unsigned int
#endif

#include "circleRenderer.h"


class CudaRenderer : public CircleRenderer {

private:

    Image* image;
    int* mousePressedLocation;
   
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
