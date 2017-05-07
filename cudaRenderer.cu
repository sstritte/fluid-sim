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
    float* VY;
    float* pressures;
    float* VXCopy;
    float* VYCopy;
    float* divergence;
    float* vorticity;
    float* color;
    float* colorCopy;
    float* imageData;

    int* mpls;
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


__device__ __inline__ int
isBoundary(int i, int j) {
    int cells_per_side = cuParams.cells_per_side;
    if (j == 0) return 1; // left 
    if (i == 0) return 2; // top
    if (j == cells_per_side) return 3; // right
    if (i == cells_per_side) return 4; // bottom
    return 0;
}

// a is prev mouse point, b is cur mouse point, p is point to consider,
// fp is fraction projection to be populated as output
__device__ __inline__ double 
distanceToSegment(double ax, double ay, double bx, double by, 
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

__device__ __inline__ double 
distanceToNearestMouseSegment(double px, double py, double *fp,
        double* vx, double *vy) {
    double minLen = DBL_MAX;
    double fpResult = 0.0;
    double vxResult = 0.0;
    double vyResult = 0.0;
    for (int i = 0; i < 400 - 2; i += 2) {

        int grid_col1 = cuParams.mpls[i];
        int grid_row1 = cuParams.mpls[i + 1];
        int grid_col2 = cuParams.mpls[i + 2];
        int grid_row2 = cuParams.mpls[i + 3];
        if (grid_col2 == 0 & grid_row2 == 0) break;
        double len = distanceToSegment(grid_col1, grid_row1, grid_col2, grid_row2, px, py, fp);
        if (len < minLen) {
            minLen = len;
            fpResult = *fp;
            vxResult = grid_col2 - grid_col1;
            vyResult = grid_row2 - grid_row1;
        }        

    }
    *fp = fpResult;
    *vx = vxResult;
    *vy = vyResult;
    return minLen;
}

//kernelFadeVelocities
__global__ void kernelFadeVelocities() {
    int grid_col = blockIdx.x * blockDim.x + threadIdx.x;
    int grid_row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;
    int height = cuParams.height;

    if (grid_col >= width || grid_row >= height) return;
    if (grid_row * width + grid_col >= width * height) return; 
    
    cuParams.VX[grid_row * width + grid_col] *= 0.999;
    cuParams.VY[grid_row * width + grid_col] *= 0.999;
}

//kernelSetNewVelocities
__global__ void kernelSetNewVelocities() {
    int grid_col = blockIdx.x * blockDim.x + threadIdx.x;
    int grid_row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;
    int height = cuParams.height;

    if (grid_col >= width || grid_row >= height) return;
    if (grid_row * width + grid_col >= width * height) return; 
    
    int imageX = grid_col;
    int imageY = grid_row;
    int offset = 4 * (imageY * width + imageX);
    float4 value = make_float4(1.f,0.f,1.f,1.f);

    // Write to global memory.
    *(float4*)(&cuParams.imageData[offset]) = value;
   
    cuParams.VX[grid_row * width + grid_col] *= 0.999;
    cuParams.VY[grid_row * width + grid_col] *= 0.999;
    double projection;
    double vx;
    double vy;
    double l = distanceToNearestMouseSegment(grid_col, grid_row, 
            &projection, &vx, &vy);
    //printf("velocity %f,%f\n", mouseSegmentVelocity.first, mouseSegmentVelocity.second);
    double taperFactor = 0.6;
    double projectedFraction = 1.0 - fminf(1.0, fmaxf(projection, 0.0)) * taperFactor;
    double R = 10;
    double m = exp(-l/R); //drag coefficient
    m *= projectedFraction * projectedFraction;
    double targetVelocityX = vx * 1 * 1.4; 
    double targetVelocityY = vy * 1 * 1.4; 

    cuParams.VX[grid_row * width + grid_col] += 
        (targetVelocityX - cuParams.VX[grid_row * width + grid_col]) * m;
    cuParams.VY[grid_row * width + grid_col] += 
        (targetVelocityY - cuParams.VY[grid_row * width + grid_col]) * m;

}

//kernelAdvectVelocityForward
__global__ void kernelAdvectVelocityForward() {
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    cuParams.VXCopy[row * width + col] = cuParams.VX[row * width + col];
    cuParams.VYCopy[row * width + col] = cuParams.VY[row * width + col];

   int pixelRow = row * CELL_DIM;
   int pixelCol = col * CELL_DIM;
   int nextPixelRow = round(pixelRow + TIME_STEP * cuParams.VY[row * width + col] * CELL_DIM);
   int nextPixelCol = round(pixelCol + TIME_STEP * cuParams.VX[row * width + col] * CELL_DIM);
   int nextCellCol = nextPixelCol / CELL_DIM;
   int nextCellRow = nextPixelRow / CELL_DIM;

   if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
           && nextCellCol >= 0 && nextCellRow >= 0) {
        cuParams.VX[nextCellRow * width + nextCellCol] = cuParams.VXCopy[row * width + col];
        cuParams.VY[nextCellRow * width + nextCellCol] = cuParams.VYCopy[row* width + col];
   } 

}

//kernelAdvectVelocityBackward
__global__ void kernelAdvectVelocityBackward() {
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

   int pixelRow = row * CELL_DIM;
   int pixelCol = col * CELL_DIM;
   int prevPixelRow = round(pixelRow - TIME_STEP * cuParams.VY[row * width + col] * CELL_DIM);
   int prevPixelCol = round(pixelCol - TIME_STEP * cuParams.VX[row * width + col] * CELL_DIM);
   int prevCellCol = prevPixelCol / CELL_DIM;
   int prevCellRow = prevPixelRow / CELL_DIM;

   if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
           && prevCellCol >= 0 && prevCellRow >= 0) {
        cuParams.VX[row * width + col] = cuParams.VXCopy[prevCellRow * width + prevCellCol];
        cuParams.VY[row * width + col] = cuParams.VYCopy[prevCellRow * width + prevCellCol];
   } 
   if (prevCellCol == col && prevCellRow == row) {
        // you don't move so just disappear
        cuParams.VX[row * width + col] = 0;
        cuParams.VY[row * width + col] = 0;
   }
}

//kernelApplyVorticity
__global__ void kernelApplyVorticity(){
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    if (isBoundary(row,col)) return;

    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;

    if (row > 0) T = cuParams.VX[(row-1) * width + col];
    if (row < cells_per_side) B = cuParams.VX[(row+1) * width + col];
    if (col < cells_per_side) R = cuParams.VY[row * width + (col+1)];
    if (col > 0) L = cuParams.VY[row * width + (col-1)];
    cuParams.vorticity[row * width + col] = 0.5 * ((R - L) - (T - B));
}

//kernelApplyVorticityForce
__global__ void kernelApplyVorticityForce(){
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    float vortConfinementFloat = 0.035f;
    float vortL = 0.0;
    float vortR = 0.0;
    float vortB = 0.0;
    float vortT = 0.0;
    float vortC = 0.0;
            
    if (isBoundary(row,col)) return;

    if (row > 0) vortT = cuParams.vorticity[(row-1) * width + col];
    if (row < cells_per_side) vortB = cuParams.vorticity[(row+1) * width + col];
    if (col < cells_per_side) vortR = cuParams.vorticity[row * width + (col+1)];
    if (row > 0) vortL = cuParams.vorticity[row * width + (col-1)];
    vortC = cuParams.vorticity[row * width + col];
    
    float forceX = 0.5 * (fabsf(vortT) - fabsf(vortB));
    float forceY = 0.5 * (fabsf(vortR) - fabsf(vortL));
    float EPSILON = powf(2,-12);
    float magSqr = fmaxf(EPSILON, forceX * forceX + forceY * forceY);
    forceX = forceX * (1/sqrtf(magSqr));
    forceY = forceY * (1/sqrtf(magSqr));
    forceX *= vortConfinementFloat * vortC * 1;
    forceY *= vortConfinementFloat * vortC * -1;
    cuParams.VX[row * width + col] += forceX;
    cuParams.VY[row * width + col] += forceY;
}

//kernelApplyDivergence
__global__ void kernelApplyDivergence(){
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    if (isBoundary(row,col)) return;

    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;

    if (row > 0) T = cuParams.VY[(row-1) * width + col];
    if (row < cells_per_side) B = cuParams.VY[(row+1) * width + col];
    if (col < cells_per_side) R = cuParams.VX[row * width + (col+1)];
    if (col > 0) L = cuParams.VX[row * width + (col-1)];
    cuParams.divergence[row * width + col] = 0.5*((R-L) + (T-B));
}

//kernelPressureSolve
__global__ void kernelPressureSolve(){
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    if (isBoundary(row,col)) return;

    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;

    if (row > 0) T = cuParams.pressures[(row-1) * width + col];
    if (row < cells_per_side) B = cuParams.pressures[(row+1) * width + col];
    if (col < cells_per_side) R = cuParams.pressures[row * width + (col+1)];
    if (col > 0) L = cuParams.pressures[row * width + (col-1)];
    // I FEEL LIKE MAYBE WE NEED A SYNCTHREADS() HERE!!!!!!!
    // BECAUSE IN THE REF WE HAVE A TEMPPRESSURES ARRAY AND COPY EVERYTHING AT THE END
    // SO THAT EVERYONE IS ALWAYS READING FROM THE OLD PRESSURE ARRAY AND THINGS AREN'T BEING
    // UPDATED AS THEY GO
    // ALTERNATIVELY WE COULD HAVE THIS KERNEL MAKE A cuParams.tempPressures IN THIS WAY
    // AND THEN HAVE  SEPARATE KERNEL (TO CALL AFTER THIS ONE) THAT COPIES cuParams.tempPressures
    // OVER TO the real cuParams.pressures
    cuParams.pressures[row * width + col] = (L + R + B + T + -1 * cuParams.divergence[row * width + col]) * .25;
}

//kernelPressureGradient
__global__ void kernelPressureGradient(){
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    if (isBoundary(row,col)) return;

    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;

    if (row > 0) T = cuParams.pressures[(row-1) * width + col];
    if (row < cells_per_side) B = cuParams.pressures[(row+1) * width + col];
    if (col < cells_per_side) R = cuParams.pressures[row * width + (col+1)];
    if (col > 0) L = cuParams.pressures[row * width + (col-1)];
    cuParams.VX[row * width + col] = cuParams.VX[row * width + col] - 0.5*(R - L);
    cuParams.VY[row * width + col] = cuParams.VY[row * width + col] - 0.5*(T - B);
}

//kernelAdvectColorForward
__global__ void kernelAdvectColorForward() {
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    // I HAVE A SIMILAR CONCERN HERE AS THE COMMENT I WROTE ABOVE...
    // DO WE NEED ALL THE KERNELS TO FINISH POPULATING cuParams.colorCopy
    // BEFORE THEY START THE REST OF THE LOGIC??? SINCE EACH ACCESSES OTHER PLACES IN
    // cuParams.colorCopy BELOW. (THIS SAME QUESTION WILL HAVE TO APPLY TO 
    // kernelAdvectVelocityForward()...)
    cuParams.colorCopy[row * width + col] = cuParams.color[row * width + col];

    int pixelRow = row * CELL_DIM;
    int pixelCol = col * CELL_DIM;
    // THERE ARE A ZILLION VERSION OF ROUNDING FUNCTIONS IN THE CUDA MATH API
    // http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html#group__CUDA__MATH__SINGLE
    // I THINK WE WANT TO CHANGE round(x) TO rint(x).... BUT NOT SURE
    // SAME CHANGE WOULD NEED TO BE APPLIED IN kernelAdvectVelocityForward() and kernelAdvectVelocityBackward()
    // and kernelAdvectColorBackward()
    int nextPixelRow = round(pixelRow + TIME_STEP * cuParams.VY[row * width + col] * CELL_DIM);
    int nextPixelCol = round(pixelCol + TIME_STEP * cuParams.VX[row * width + col] * CELL_DIM);
    int nextCellCol = nextPixelCol / CELL_DIM;
    int nextCellRow = nextPixelRow / CELL_DIM;

   if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
           && nextCellCol >= 0 && nextCellRow >= 0) {
        // I SPENT LIKE 10 MINUTES TRYING TO CONVINCE MYSELF IF THE col * 4 THING IS RIGHT 
        // BUT I'M STILL NOT CONVINCED SO WE SHOULD DOUBLE/TRIPLE CHECK THIS
        cuParams.color[nextCellRow * width + nextCellCol * 4 + 0] = cuParams.colorCopy[row * width + col * 4 + 0];
        cuParams.color[nextCellRow * width + nextCellCol * 4 + 1] = cuParams.colorCopy[row * width + col * 4 + 1];
        cuParams.color[nextCellRow * width + nextCellCol * 4 + 2] = cuParams.colorCopy[row * width + col * 4 + 2];
        cuParams.color[nextCellRow * width + nextCellCol * 4 + 3] = cuParams.colorCopy[row * width + col * 4 + 3];

   } 
}

//kernelAdvectColorBackward
__global__ void kernelAdvectColorBackward() {
    int cells_per_side = cuParams.cells_per_side;
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y; 
    int width = cuParams.width;

    int pixelRow = row * CELL_DIM;
    int pixelCol = col * CELL_DIM;
    int prevPixelRow = round(pixelRow - TIME_STEP * cuParams.VY[row * width + col] * CELL_DIM);
    int prevPixelCol = round(pixelCol - TIME_STEP * cuParams.VX[row * width + col] * CELL_DIM);
    int prevCellCol = prevPixelCol / CELL_DIM;
    int prevCellRow = prevPixelRow / CELL_DIM;

    if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
            && prevCellCol >= 0 && prevCellRow >= 0) {
         cuParams.color[row * width + col * 4 + 0] = cuParams.colorCopy[prevCellRow * width + prevCellCol * 4 + 0];
         cuParams.color[row * width + col * 4 + 1] = cuParams.colorCopy[prevCellRow * width + prevCellCol * 4 + 1];
         cuParams.color[row * width + col * 4 + 2] = cuParams.colorCopy[prevCellRow * width + prevCellCol * 4 + 2];
         cuParams.color[row * width + col * 4 + 3] = cuParams.colorCopy[prevCellRow * width + prevCellCol * 4 + 3];
   } 
}


//////////////////////////////////////////////////////////////////////////
///////////////////////////HOST CODE BELOW////////////////////////////////
//////////////////////////////////////////////////////////////////////////

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

    mpls = NULL;

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

    cdMpls = NULL;
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
        delete mpls;
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
        cudaFree(cdMpls);
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


void
CudaRenderer::setup() {
   cells_per_side = image->width / CELL_DIM - 1;

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
   cudaMalloc(&cdMpls, 400 * sizeof(float) * image->width * image->height);

   cudaMemset(cdVX, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdVY, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdPressures, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdVXCopy, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdVYCopy, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdDivergence, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdVorticity, 0, sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdColor, 0, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));
   cudaMemset(cdColorCopy, 0, 4 * sizeof(float) * (cells_per_side + 1) * (cells_per_side + 1));

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
    params.mpls = cdMpls;

    cudaMemcpyToSymbol(cuParams, &params, sizeof(GlobalConstants));
}

// Called after clear, before render
void CudaRenderer::setNewQuantities(std::vector<std::pair<int, int> > mpls) {

    int mplsSize = mpls.size();
    if (mplsSize == 0) {
        // if mpls.size is 0, then call kernel that decreases VX,VY by 0.999
        dim3 blockDim(16,16,1);
        dim3 gridDim(
                (image->width + blockDim.x - 1) / blockDim.x,
                (image->height + blockDim.y - 1) / blockDim.y);
        kernelFadeVelocities<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();

    } else {
        int* mplsArray = new int[mplsSize * 2 + 1];
        int count = 0;
        for (std::vector<std::pair<int,int> >::iterator it = mpls.begin() 
                ; it != mpls.end(); ++it) {
            std::pair<int,int> c = *it;
            mplsArray[count] = c.first;
            mplsArray[count + 1] = c.second;
            count += 2;
        }
        cudaMemset(cdMpls, 0, 400 * sizeof(int));
        cudaMemcpy(cdMpls, mplsArray, (mplsSize * 2 + 1) * sizeof(int), 
                cudaMemcpyHostToDevice);

        dim3 blockDim(16,16,1);
        dim3 gridDim(
                (image->width + blockDim.x - 1) / blockDim.x,
                (image->height + blockDim.y - 1) / blockDim.y);
        kernelSetNewVelocities<<<gridDim, blockDim>>>();
        cudaDeviceSynchronize();
    }
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
    kernelClearImage<<<gridDim, blockDim>>>(1.f,1.f,0.f,1.f);
    cudaDeviceSynchronize();
}

void
CudaRenderer::render() {
//    usleep(1000000);
    dim3 blockDim(16,16,1);
    dim3 gridDim(
            (image->width + blockDim.x - 1) / blockDim.x,
            (image->height + blockDim.y - 1) / blockDim.y);
    
    kernelAdvectVelocityForward<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelAdvectVelocityBackward<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelApplyVorticity<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelApplyVorticityForce<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelApplyDivergence<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelPressureSolve<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelPressureGradient<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
   
    //TO DO: DRAW STUFF
  
    kernelAdvectColorForward<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    kernelAdvectColorBackward<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();


    //advectVelocityForward();
    //advectVelocityBackward();
    //applyVorticity();
    //applyVorticityForce();
    //applyDivergence();
    //pressureSolve();
    //pressureGradient();

    /*// Draw
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

