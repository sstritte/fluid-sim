#include <algorithm>
#include <math.h>
#include <float.h>
#include <utility>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <unistd.h>
#include "cudaRenderer.h"
#include "image.h"
#include "util.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

#define CELL_DIM 1
#define TIME_STEP 1


///////////////////////////CUDA CODE BELOW////////////////////////////////
struct GlobalConstants {
    int cells_per_side;
    int width;
    int height;

    float* VX;
    float * VY;
    float* pressures;
    float* VXCopy;
    float* VYCopy;
    float* divergence;
    float* vorticity;
    float* color;
    float* colorCopy;
    float* imageData;
};

__constant__ GlobalConstants cuParams;

// kernelClearImage
__global__ void kernelClearImage(float r, float g, float b, float a) {
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;

    int width = cuParams.width;
    int height = cuParams.height;

    if (imageX >= width || imageY >= height) return;

    int offset = 4 * (imageY * width + imageX);
    float4 value = make_float4(r,g,b,a);
    
    // Write to global memory.
    *(float4*)(&cuParams.imageData[offset]) = value;

}


///////////////////////////HOST CODE BELOW////////////////////////////////

CudaRenderer::CudaRenderer() {
    image = NULL;

    VX = NULL;
    VY = NULL;
    color = NULL;
    colorCopy = NULL;
    pressures = NULL;
    VXCopy = NULL;
    VYCopy = NULL;
    divergence = NULL;
    vorticity = NULL;

    cdVX = NULL;
    cdVY = NULL;
    cdColor = NULL;
    cdColorCopy = NULL;
    cdPressures = NULL;
    cdVXCopy = NULL;
    cdVYCopy = NULL;
    cdDivergence = NULL;
    cdVorticity = NULL;
    cdImageData = NULL;
}

CudaRenderer::~CudaRenderer() {

    if (image) delete image;

    if (VX) {
        delete VX;
        delete VY;
        delete pressures;
        delete VXCopy;
        delete VYCopy;
        delete divergence;
        delete vorticity;
        delete color;
        delete colorCopy;
    }

    if (cdVX) {
        cudaFree(cdVX);
        cudaFree(cdVY);
        cudaFree(cdPressures);
        cudaFree(cdVXCopy);
        cudaFree(cdVYCopy);
        cudaFree(cdDivergence);
        cudaFree(cdVorticity);
        cudaFree(cdColor);
        cudaFree(cdColorCopy);
        cudaFree(cdImageData);
    }
}

const Image*
CudaRenderer::getImage() {
    printf("Copying image data from device\n");

    cudaMemcpy(image->data, cdImageData, 
            4 * sizeof(float) * image->width * image->height,
            cudaMemcpyDeviceToHost);

    return image;
}

/*int
CudaRenderer::isBoundary(int i, int j) {
    if (j == 0) return 1; // left 
    if (i == 0) return 2; // top
    if (j == cells_per_side) return 3; // right
    if (i == cells_per_side) return 4; // bottom
    return 0;
}*/

void
CudaRenderer::setup() {
   cells_per_side = image->width / CELL_DIM;

   cudaMalloc(&cdVX, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdVY, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdPressures, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdVXCopy, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdVYCopy, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdDivergence, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdVorticity, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdColor, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdColorCopy, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMalloc(&cdImageData, 4 * sizeof(float) * image->width * image->height);

   cudaMemcpy(cdVX, VX, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdVY, VY, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdPressures, pressures, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdVXCopy, VXCopy, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdVYCopy, VYCopy, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdDivergence, divergence, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdVorticity, vorticity, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdColor, color, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);
   cudaMemcpy(cdColorCopy, colorCopy, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1), cudaMemcpyHostToDevice);

    GlobalConstants params;
    params.cells_per_side = cells_per_side;
    params.width = image->width;
    params.height = image->height;
    params.VX = cdVX;
    params.VY = cdVY;
    params.pressures = cdPressures;
    params.VXCopy = cdVXCopy;
    params.VYCopy = cdVYCopy;
    params.divergence = cdDivergence;
    params.vorticity = cdVorticity;
    params.color = cdColor;
    params.colorCopy = cdColorCopy;
    params.imageData = cdImageData;

    cudaMemcpyToSymbol(cuParams, &params, sizeof(GlobalConstants));
}

// a is prev mouse point, b is cur mouse point, p is point to consider,
// fp is fraction projection to be populated as output
/*double 
CudaRenderer::distanceToSegment(double ax, double ay, double bx, double by, 
        double px, double py, double* fp) {
    double dx = px - ax; //vec2 d = p - a;
    double dy = py - ay;
    double xx = bx - ax; //vec2 x = b - a;
    double xy = by - ay;
    *fp = 0.0; // fractional projection, 0 - 1 in the length of b-a
    double lx = sqrt(xx*xx + xy*xy); //length(x)
    double ld = sqrt(dx*dx + dy*dy); //length(d)
    if (lx <= 0.0001) return ld;
    double projection = dx*(xx/lx) + dy*(xy/lx); //dot(d, x/lx)
    *fp = projection / lx;
    if (projection < 0.0) return ld;
    else if (projection > lx) return sqrt((px-bx) * (px-bx) +
            (py-by) * (py-by)); //length(p - b)
    return sqrt(abs(dx*dx + dy*dy - projection * projection));
}

double 
CudaRenderer::distanceToNearestMouseSegment(double px, double py, double *fp,
        std::pair<double,double> *mouseSegmentVelocity) {
    double minLen = DBL_MAX;
    double fpResult = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    for (std::vector<std::pair<int,int> >::iterator 
            it = mousePressedLocations.begin() 
            ; it + 1 != mousePressedLocations.end(); ++it) {
        std::pair<int,int> coords1 = *it;
        int grid_col1 = coords1.first;
        int grid_row1 = coords1.second;
        std::pair<int,int> coords2 = *(it + 1);
        int grid_col2 = coords2.first;
        int grid_row2 = coords2.second; 
        double len = distanceToSegment(grid_col1, grid_row1, grid_col2, grid_row2, px, py, fp);
        if (len < minLen) {
            minLen = len;
            fpResult = *fp;
            vx = grid_col2 - grid_col1;
            vy = grid_row2 - grid_row1;
        }         
    }
    *fp = fpResult;
    std::pair<double,double> msvResult = std::make_pair(vx,vy);
    *mouseSegmentVelocity = msvResult;
    return minLen;
}
*/
void CudaRenderer::setNewQuantities(std::vector<std::pair<int, int> > mpls) {
    /*mousePressedLocations.clear();

    for (std::vector<std::pair<int,int> >::iterator it = mpls.begin() 
            ; it != mpls.end(); ++it) {
        std::pair<int,int> c = *it;
        mousePressedLocations.push_back(std::make_pair(c.first, c.second));
        //printf("putting in %d,%d\n", c.first, c.second);
    }
    
    for (int i = 0; i < image->height * image-> width; i++) {
        int grid_row = (i / image->width) / (CELL_DIM);
        int grid_col = (i % image->width) / (CELL_DIM);

        cdVX[grid_row][grid_col] *= 0.999;
        cdVY[grid_row][grid_col] *= 0.999;
        if (mpls.size() > 0) { //isMouseDown) {
            //double d = sqrt(vxs[i] * vxs[i] + vys[i] * vys[i]);
            double projection;
            std::pair<double,double> mouseSegmentVelocity;
            double l = distanceToNearestMouseSegment(grid_col, grid_row, 
                    &projection, &mouseSegmentVelocity);
            //printf("velocity %f,%f\n", mouseSegmentVelocity.first, mouseSegmentVelocity.second);
            double vx = mouseSegmentVelocity.first;
            double vy = mouseSegmentVelocity.second;
            double taperFactor = 0.6;
            double projectedFraction = 1.0 - std::min(1.0, std::max(projection, 0.0)) * taperFactor;
            double R = 10;
            double m = exp(-l/R); //drag coefficient
            m *= projectedFraction * projectedFraction;
            double targetVelocityX = vx * 1 * 1.4; 
            double targetVelocityY = vy * 1 * 1.4; 

            //if (vxs[i] != 0.0 || vys[i] != 0.0) {
                cdVX[grid_row][grid_col] += 
                    (targetVelocityX - cdVX[grid_row][grid_col]) * m;
                cdVY[grid_row][grid_col] += 
                    (targetVelocityY - cdVY[grid_row][grid_col]) * m;
                //printf("setting velocity of row %d col %d to [%f,%f]\n", grid_row, grid_col, cdVX[grid_row][grid_col], cdVY[grid_row][grid_col]);
            //}
        }

    }*/

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
// Clear's the renderer's target image.  
void
CudaRenderer::clearImage() {
    dim3 blockDim(16,16,1);
    dim3 gridDim(
            (image->width + blockDim.x - 1) / blockDim.x,
            (image->height + blockDim.y - 1) / blockDim.y);
    kernelClearImage<<<gridDim, blockDim>>>(1.f,1.f,1.f,1.f);
    cudaDeviceSynchronize();
}

/*void CudaRenderer::advectColorBackward() {
     //Advecting the values in cdColor
     for (int row = 0; row < cells_per_side; row++) {
         for (int col = 0; col < cells_per_side; col++) {
            int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
            int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
            int prevPixelRow = round(pixelRow - TIME_STEP * cdVY[row][col] * CELL_DIM);
            int prevPixelCol = round(pixelCol - TIME_STEP * cdVX[row][col] * CELL_DIM);
            int prevCellCol = prevPixelCol / CELL_DIM;
            int prevCellRow = prevPixelRow / CELL_DIM;
 
            if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
                    && prevCellCol >= 0 && prevCellRow >= 0) {
                 cdColor[row][col][0] = cdColorCopy[prevCellRow][prevCellCol][0];
                 cdColor[row][col][1] = cdColorCopy[prevCellRow][prevCellCol][1];
                 cdColor[row][col][2] = cdColorCopy[prevCellRow][prevCellCol][2];
                 cdColor[row][col][3] = cdColorCopy[prevCellRow][prevCellCol][3];
            } 
         }
     }
 }
 
 void CudaRenderer::advectColorForward() {
     //Advecting the values in cdColor
     for (int row = 0; row < cells_per_side; row++) {
         for (int col = 0; col < cells_per_side; col++) {
            int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
            int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
            int nextPixelRow = round(pixelRow + TIME_STEP * cdVY[row][col] * CELL_DIM);
            int nextPixelCol = round(pixelCol + TIME_STEP * cdVX[row][col] * CELL_DIM);
            int nextCellCol = nextPixelCol / CELL_DIM;
            int nextCellRow = nextPixelRow / CELL_DIM;
 
            if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
                    && nextCellCol >= 0 && nextCellRow >= 0) {
                 cdColor[nextCellRow][nextCellCol][0] = cdColorCopy[row][col][0];
                 cdColor[nextCellRow][nextCellCol][1] = cdColorCopy[row][col][1];
                 cdColor[nextCellRow][nextCellCol][2] = cdColorCopy[row][col][2];
                 cdColor[nextCellRow][nextCellCol][3] = cdColorCopy[row][col][3];
            } 
         }
     }
 }
 

void
CudaRenderer::advectColor() {
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            for (int k = 0; k < 4; k++) {
                cdColorCopy[i][j][k] = cdColor[i][j][k];
            }
        }
    }
    advectColorForward();
    advectColorBackward();
}


void
CudaRenderer::applyPressure() {
    for (int i = 0; i < cells_per_side; i ++) {
        for (int j = 0; j < cells_per_side; j++) {
            float force_x = (cdPressures[i][j] - cdPressures[i][j+1]);
            float force_y = (cdPressures[i][j] - cdPressures[i+1][j]);
            cdVX[i][j] += force_x;
            cdVX[i][j+1] += force_x;
            cdVY[i][j] += force_y;
            cdVY[i+1][j] += force_y;
            //if (force_x != 0) printf("adding force_x %f to (%d,%d)\n", force_x, i, j);
        }
    }
}



void CudaRenderer::advectVelocityForward() {
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            cdVXCopy[i][j] = cdVX[i][j];
            cdVYCopy[i][j] = cdVY[i][j];
        }
    }

    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int nextPixelRow = round(pixelRow + TIME_STEP * cdVY[row][col] * CELL_DIM);
           int nextPixelCol = round(pixelCol + TIME_STEP * cdVX[row][col] * CELL_DIM);
           int nextCellCol = nextPixelCol / CELL_DIM;
           int nextCellRow = nextPixelRow / CELL_DIM;

           if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
                   && nextCellCol >= 0 && nextCellRow >= 0) {
                cdVX[nextCellRow][nextCellCol] = cdVXCopy[row][col];
                cdVY[nextCellRow][nextCellCol] = cdVYCopy[row][col];
           } 
        }
    }

}

void
CudaRenderer::advectVelocityBackward() {
    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int prevPixelRow = round(pixelRow - TIME_STEP * cdVY[row][col] * CELL_DIM);
           int prevPixelCol = round(pixelCol - TIME_STEP * cdVX[row][col] * CELL_DIM);
           int prevCellCol = prevPixelCol / CELL_DIM;
           int prevCellRow = prevPixelRow / CELL_DIM;

           if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
                   && prevCellCol >= 0 && prevCellRow >= 0) {
                cdVX[row][col] = cdVXCopy[prevCellRow][prevCellCol];
                cdVY[row][col] = cdVYCopy[prevCellRow][prevCellCol];
           } 
           if (prevCellCol == col && prevCellRow == row) {
                // you don't move so just disappear
                cdVX[row][col] = 0;
                cdVY[row][col] = 0;
           }
        }
    }
}

// Divergence of velocity: This computes how divergent the velocity field is
// (how much in/out flow there is at every point).  Used as input to the 
// pressure solve below.
void
CudaRenderer::applyDivergence() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            
            if (i > 0) T = cdVY[i-1][j];
            if (i < cells_per_side) B = cdVY[i+1][j];
            if (j < cells_per_side) R = cdVX[i][j+1];
            if (j > 0) L = cdVX[i][j-1];
            cdDivergence[i][j] = 0.5*((R-L) + (T-B));   
        }
    }
}

void
CudaRenderer::pressureSolve() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    float** tempPressure = new float*[cells_per_side + 1];
    for(int i = 0; i < cells_per_side+1; i++) {
        tempPressure[i] = new float[cells_per_side+1];
    }

    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) T = cdPressures[i-1][j];
            if (i < cells_per_side) B = cdPressures[i+1][j];
            if (j < cells_per_side) R = cdPressures[i][j+1];
            if (j > 0) L = cdPressures[i][j-1];
            tempPressure[i][j] = (L + R + B + T + -1 * cdDivergence[i][j]) * .25;
                //if (L+R+B+T > 0.0) printf("L+R+B+T is %f, -1*cdDivergence is %f\n", L+R+B+T, -1*cdDivergence[i][j]);
        }
    }
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) { 
            if (isBoundary(i,j)) continue;
            cdPressures[i][j] = tempPressure[i][j];
        }
    }

}

void
CudaRenderer::pressureGradient() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) T = cdPressures[i-1][j];
            if (i < cells_per_side) B = cdPressures[i+1][j];
            if (j < cells_per_side) R = cdPressures[i][j+1];
            if (j > 0) L = cdPressures[i][j-1];

            //if (cdVY[i][j] > 10) printf("doing cdVY = %f - 0.5*(%f)\n", cdVY[i][j], T-B);
            cdVX[i][j] = cdVX[i][j] - 0.5*(R - L);
            cdVY[i][j] = cdVY[i][j] - 0.5*(T - B);
            //if (cdVY[i][j] != 0.0) printf("cdVY is %f\n", cdVY[i][j]);
        }
    }
}

void
CudaRenderer::applyVorticity() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) T = cdVX[i-1][j];
            if (i < cells_per_side) B = cdVX[i+1][j];
            if (j < cells_per_side) R = cdVY[i][j+1];
            if (j > 0) L = cdVY[i][j-1];
            cdVorticity[i][j] = 0.5 * ((R - L) - (T - B));
        }
    }
}

void
CudaRenderer::applyVorticityForce() {
    float vortConfinementFloat = 0.035f;
    float vortL = 0.0;
    float vortR = 0.0;
    float vortB = 0.0;
    float vortT = 0.0;
    float vortC = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) vortT = cdVorticity[i-1][j];
            if (i < cells_per_side) vortB = cdVorticity[i+1][j];
            if (j < cells_per_side) vortR = cdVorticity[i][j+1];
            if (j > 0) vortL = cdVorticity[i][j-1];
            vortC = cdVorticity[i][j];
            float forceX = 0.5 * (abs(vortT) - abs(vortB));
            float forceY = 0.5 * (abs(vortR) - abs(vortL));
            float EPSILON = pow(2,-12);
            float magSqr = std::max(EPSILON, forceX * forceX + forceY * forceY);
            forceX = forceX * (1/sqrt(magSqr));
            forceY = forceY * (1/sqrt(magSqr));
            forceX *= vortConfinementFloat * vortC * 1;
            forceY *= vortConfinementFloat * vortC * -1;
            cdVX[i][j] += forceX;
            cdVY[i][j] += forceY;
        }
    }
}
*/
void
CudaRenderer::render() {
/*    advectVelocityForward();
    advectVelocityBackward();
    applyVorticity();
    applyVorticityForce();
    applyDivergence();
    pressureSolve();
    pressureGradient();

    // Draw
    for (int i = 0; i < 4*image->width*image->height; i+=4) {
            int grid_row = ((i/4) / image->width) / (CELL_DIM);
            int grid_col = ((i/4) % image->width) / (CELL_DIM);
            double vx = cdVX[grid_row][grid_col];
            double vy = cdVY[grid_row][grid_col];
            double v = sqrt(vx * vx + vy * vy);

            if (abs(v) < 0.00001) {
                // make the cdColor go away faster
                cdColor[grid_row][grid_col][0] *= 0.9;
                cdColor[grid_row][grid_col][1] *= 0.9;
                cdColor[grid_row][grid_col][2] *= 0.9;
                cdColor[grid_row][grid_col][3] = 1.0;
            } 
            cdColor[grid_row][grid_col][0] *= 0.9494; 
            cdColor[grid_row][grid_col][1] *= 0.9494; 
            cdColor[grid_row][grid_col][2] *= 0.9696; 
        
            if (mousePressedLocations.size() > 0) {
                //double d = sqrt(vx * vx + vy * vy);
                double projection;
                std::pair<double,double> mouseSegmentVelocity;
                double l = distanceToNearestMouseSegment(grid_col, grid_row, 
                        &projection, &mouseSegmentVelocity);

                //if (l < 1.0) printf("l is %f\n", l);
                float taperFactor = 0.6;
                double projectedFraction = 1.0 - std::min(1.0, 
                        std::max(projection, 0.0)) * taperFactor;
                double R = 12; //0.025; // the bigger, the more stuff gets cdColored
                double m = exp(-l/R); //drag coefficient
                //double speed = d;
                double vx = cdVX[grid_row][grid_col];
                double vy = cdVY[grid_row][grid_col];
                double speed = sqrt(vx * vx + vy * vy);

                //printf("l is %f, m is %f, projection is %f\n", l, m, projection);

                double x = std::min(1.0, std::max(fabs((speed * speed * 0.02 - 
                            projection * 5.0) * projectedFraction), 0.0));

                double r = (2.4 / 60.0) * x + (0.2 /30.0) * (1-x) + (1.0 * pow(x, 9.0));
                double g = (0.0 / 60.0) * x + (51.8 / 30.0) * (1-x) + (1.0 * pow(x, 9.0));
                double b = (5.9 / 60.0) * x + (100.0 / 30.0) * (1-x) + (1.0 * pow(x, 9.0));

                cdColor[grid_row][grid_col][0] += m * r;
                cdColor[grid_row][grid_col][1] += m * g;
                cdColor[grid_row][grid_col][2] += m * b;
                cdColor[grid_row][grid_col][3] = 1.0;
            }

            image->data[i] = cdColor[grid_row][grid_col][0];
            image->data[i+1] = cdColor[grid_row][grid_col][1];
            image->data[i+2] = cdColor[grid_row][grid_col][2];
            image->data[i+3] = cdColor[grid_row][grid_col][3];
    }
    advectColor();*/
}

