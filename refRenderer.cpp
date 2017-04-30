#include <algorithm>
#include <math.h>
#include <float.h>
#include <utility>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <unistd.h>
#include "refRenderer.h"
#include "image.h"
#include "util.h"

#define CELL_DIM 1
#define TIME_STEP 1

RefRenderer::RefRenderer() {
    image = NULL;
    mousePressedLocation = NULL;
    //mousePressedLocations = NULL; gives an error
    velocitiesX = NULL;
    velocitiesY = NULL;
    color = NULL;
    colorCopy = NULL;
    pressures = NULL;
    advectionCopyX = NULL;
    advectionCopyY = NULL;
    advectionCopy = NULL;
    divergence = NULL;
    vorticity = NULL;
}

RefRenderer::~RefRenderer() {

    if (image) delete image;
    if (mousePressedLocation) delete mousePressedLocation;
    //if (mousePressedLocations) delete mousePressedLocations; gives error

    if (velocitiesX) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (velocitiesX[i]) delete velocitiesX[i];

        }
        delete velocitiesX;
    }
    if (velocitiesY) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (velocitiesY[i]) delete velocitiesY[i];
        }
        delete velocitiesY;
    }

    if (pressures) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (pressures[i]) delete pressures[i];
        }
        delete pressures;
    }
    if (advectionCopyX) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (advectionCopyX[i]) delete advectionCopyX[i];
        }
        delete advectionCopyX;
    }
    if (advectionCopyY) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (advectionCopyY[i]) delete advectionCopyY[i];
        }
        delete advectionCopyY;
    }
    if (advectionCopy) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (advectionCopy[i]) delete advectionCopy[i];
        }
        delete advectionCopy;
    }
    if (divergence) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (divergence[i]) delete divergence[i];
        }
        delete divergence;
    }
    if (vorticity) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (vorticity[i]) delete vorticity[i];
        }
        delete vorticity;
    }

    if (color) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (color[i]) {
                for (int j = 0; j < cells_per_side; j++) {
                    if (color[i][j]) delete color[i][j];
                }
            }
            delete color[i];
        }
        delete color;
    }
    if (colorCopy) {
        for (int i = 0; i < cells_per_side + 1; i++) {
            if (colorCopy[i]) {
                for (int j = 0; j < cells_per_side; j++) {
                    if (colorCopy[i][j]) delete colorCopy[i][j];
                }
            }
            delete colorCopy[i];
        }
        delete colorCopy;
    }

}

const Image*
RefRenderer::getImage() {
    return image;
}

int
RefRenderer::isBoundary(int i, int j) {
    if (j == 0) return 1; // left 
    if (i == 0) return 2; // top
    if (j == cells_per_side) return 3; // right
    if (i == cells_per_side) return 4; // bottom
    return 0;
}

void
RefRenderer::setup() {
   cells_per_side = image->width / CELL_DIM;

   mousePressedLocation = new int[image->width * image->height];

   velocitiesX = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    velocitiesX[i] = new float[cells_per_side + 1];
    memset(velocitiesX[i], 0, sizeof(float) * (cells_per_side + 1));
   }

   velocitiesY = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    velocitiesY[i] = new float[cells_per_side + 1];
    memset(velocitiesY[i], 0, sizeof(float) * (cells_per_side + 1));
   }

   pressures = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    pressures[i] = new float[cells_per_side + 1];
    memset(pressures[i], 0, sizeof(float) * (cells_per_side + 1));
   }
   advectionCopy = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    advectionCopy[i] = new float[cells_per_side + 1];
   }
   advectionCopyX = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    advectionCopyX[i] = new float[cells_per_side + 1];
   }
   advectionCopyY = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    advectionCopyY[i] = new float[cells_per_side + 1];
   }
   divergence = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    divergence[i] = new float[cells_per_side + 1];
   }
   vorticity = new float*[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
    vorticity[i] = new float[cells_per_side + 1];
   }

   color = new float**[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
       color[i] = new float*[(cells_per_side + 1)];
       for (int j = 0; j < cells_per_side + 1; j++) {
           color[i][j] = new float[4];
           color[i][j][0] = 0;
           color[i][j][1] = 0.0392;
           color[i][j][2] = 0.1098;
           color[i][j][3] = 1.0;
       }
   }
   colorCopy = new float**[cells_per_side + 1];
   for (int i = 0; i < cells_per_side + 1; i++) {
       colorCopy[i] = new float*[(cells_per_side + 1)];
       for (int j = 0; j < cells_per_side + 1; j++) {
           colorCopy[i][j] = new float[4];
       }
    }
}

// a is prev mouse point, b is cur mouse point, p is point to consider,
// fp is fraction projection to be populated as output
double 
RefRenderer::distanceToSegment(double ax, double ay, double bx, double by, 
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
    return sqrt(abs(dx*dx + dy*dy) - projection * projection);
}

// THIS FUNCTION IS SOOOOO SLOW LOL cause when we draw the pixels (and when
// the mouse is pressed) we iterate through the whole grid for each cell
double 
RefRenderer::distanceToNearestMousePoint(double px, double py, double *fp) {
    double minLen = DBL_MAX;
    double fpResult = 0.0;
    for (std::vector<std::pair<int,int> >::iterator 
            it = mousePressedLocations.begin() 
            ; it != mousePressedLocations.end(); ++it) {
        std::pair<int,int> coords = *it;
        int grid_row = coords.first;
        int grid_col = coords.second; 
        double len = distanceToSegment(grid_col, grid_row, grid_col, grid_row, px, py, fp);
        if (len < minLen) {
            minLen = len;
            fpResult = *fp;
        }            
    }
    *fp = fpResult;
    return minLen;
}

void RefRenderer::setNewQuantities(double* vxs, double* vys, int* mpl, bool mouseDown) {
    isMouseDown = mouseDown;
    if (mouseDown) printf("MOUSE DOWN\n");
    mousePressedLocations.clear();
    for (int i = 0; i < image->height * image-> width; i++) {
        int grid_row = (i / image->width) / (CELL_DIM);
        int grid_col = (i % image->width) / (CELL_DIM);

        mousePressedLocation[i] = mpl[i]; // unncessary now pretty sure
        
        if (mpl[i]) {
            std::pair<int,int> coords = std::make_pair(grid_row, grid_col);
            mousePressedLocations.push_back(coords);
        }
        velocitiesX[grid_row][grid_col] *= 0.99;
        velocitiesY[grid_row][grid_col] *= 0.99;
        if (isMouseDown) {
            /*double d = sqrt(vxs[i] * vxs[i] + vys[i] * vys[i]);
            double projection;
            double l = distanceToNearestMousePoint(grid_col, grid_row, &projection);
            double taperFactor = 0.6;
            double projectedFraction = 1.0 - std::min(1.0, std::max(projection, 0.0)) * taperFactor;
            double R = 50;
            double m = exp(-l/R); //drag coefficient
            m *= projectedFraction * projectedFraction;*/
            double targetVelocityX = vxs[i] * 1 * 1.4; 
            double targetVelocityY = vys[i] * 1 * 1.4; 

            if (vxs[i] != 0.0 || vys[i] != 0.0) {
                velocitiesX[grid_row][grid_col] += 
                    (targetVelocityX - velocitiesX[grid_row][grid_col]);// * m;
                velocitiesY[grid_row][grid_col] += 
                    (targetVelocityY - velocitiesY[grid_row][grid_col]);// * m;
                //printf("setting velocity of row %d col %d to [%f,%f]\n", grid_row, grid_col, velocitiesX[grid_row][grid_col], velocitiesY[grid_row][grid_col]);
            }
        }

    }

}

// allocOutputImage --
//
// Allocate buffer the renderer will render into.  Check status of
// image first to avoid memory leak.
void
RefRenderer::allocOutputImage(int width, int height) {

    if (image)
        delete image;
    image = new Image(width, height);
}

// clearImage --
//
// Clear's the renderer's target image.  
void
RefRenderer::clearImage() {
    image->clear(1.f, 1.f, 1.f, 1.f);
}

void RefRenderer::advectBackward(float** quantity) {
    //Advecting the values
    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int prevPixelRow = round(pixelRow - TIME_STEP * velocitiesY[row][col] * 1/CELL_DIM);
           int prevPixelCol = round(pixelCol - TIME_STEP * velocitiesX[row][col] * 1/CELL_DIM);
           int prevCellCol = prevPixelCol / CELL_DIM;
           int prevCellRow = prevPixelRow / CELL_DIM;

           if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
                   && prevCellCol >= 0 && prevCellRow >= 0) {
                quantity[row][col] = advectionCopy[prevCellRow][prevCellCol];
           } else {
                //quantity[row][col] = 0.0;
           }
           if (prevCellCol == col && prevCellRow == row) {
                // you don't move so just disappear
                quantity[row][col] = 0;
           }
        }
    }
}

void RefRenderer::advectForward(float** quantity) {
    //Advecting the values
    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int nextPixelRow = round(pixelRow + TIME_STEP * velocitiesY[row][col] * 1/CELL_DIM);
           int nextPixelCol = round(pixelCol + TIME_STEP * velocitiesX[row][col] * 1/CELL_DIM);
           int nextCellCol = nextPixelCol / CELL_DIM;
           int nextCellRow = nextPixelRow / CELL_DIM;

           if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
                   && nextCellCol >= 0 && nextCellRow >= 0) {
                quantity[nextCellRow][nextCellCol] = advectionCopy[row][col];
           } else {
                //quantity[nextCellRow][nextCellCol] = 0.0;
           }
        }
    }
}

void
RefRenderer::advectQuantity(float** quantity) {
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            advectionCopy[i][j] = quantity[i][j];
        }
    }
    advectForward(quantity);
    advectBackward(quantity);
}

void RefRenderer::advectColorBackward() {
     //Advecting the values in color
     for (int row = 0; row < cells_per_side; row++) {
         for (int col = 0; col < cells_per_side; col++) {
            int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
            int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
            int prevPixelRow = round(pixelRow - TIME_STEP * velocitiesY[row][col] * CELL_DIM);
            int prevPixelCol = round(pixelCol - TIME_STEP * velocitiesX[row][col] * CELL_DIM);
            int prevCellCol = prevPixelCol / CELL_DIM;
            int prevCellRow = prevPixelRow / CELL_DIM;
 
            if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
                    && prevCellCol >= 0 && prevCellRow >= 0) {
                 color[row][col][0] = colorCopy[prevCellRow][prevCellCol][0];
                 color[row][col][1] = colorCopy[prevCellRow][prevCellCol][1];
                 color[row][col][2] = colorCopy[prevCellRow][prevCellCol][2];
                 color[row][col][3] = colorCopy[prevCellRow][prevCellCol][3];
            } 
         }
     }
 }
 
 void RefRenderer::advectColorForward() {
     //Advecting the values in color
     for (int row = 0; row < cells_per_side; row++) {
         for (int col = 0; col < cells_per_side; col++) {
            int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
            int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
            int nextPixelRow = round(pixelRow + TIME_STEP * velocitiesY[row][col] * CELL_DIM);
            int nextPixelCol = round(pixelCol + TIME_STEP * velocitiesX[row][col] * CELL_DIM);
            int nextCellCol = nextPixelCol / CELL_DIM;
            int nextCellRow = nextPixelRow / CELL_DIM;
 
            int prevPixelRow = round(pixelRow - TIME_STEP * velocitiesY[row][col] * CELL_DIM);
            int prevPixelCol = round(pixelCol - TIME_STEP * velocitiesX[row][col] * CELL_DIM);
            int prevCellCol = prevPixelCol / CELL_DIM;
            int prevCellRow = prevPixelRow / CELL_DIM;
 
            if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
                    && nextCellCol >= 0 && nextCellRow >= 0) {
                 color[nextCellRow][nextCellCol][0] = colorCopy[row][col][0];
                 color[nextCellRow][nextCellCol][1] = colorCopy[row][col][1];
                 color[nextCellRow][nextCellCol][2] = colorCopy[row][col][2];
                 color[nextCellRow][nextCellCol][3] = colorCopy[row][col][3];
            } 
         }
     }
 }
 

void
RefRenderer::advectColor() {
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            for (int k = 0; k < 4; k++) {
                colorCopy[i][j][k] = color[i][j][k];
            }
        }
    }
    advectColorForward();
    advectColorBackward();
}


void
RefRenderer::applyPressure() {
    for (int i = 0; i < cells_per_side; i ++) {
        for (int j = 0; j < cells_per_side; j++) {
            float force_x = (pressures[i][j] - pressures[i][j+1]);
            float force_y = (pressures[i][j] - pressures[i+1][j]);
            velocitiesX[i][j] += force_x;
            velocitiesX[i][j+1] += force_x;
            velocitiesY[i][j] += force_y;
            velocitiesY[i+1][j] += force_y;
            //if (force_x != 0) printf("adding force_x %f to (%d,%d)\n", force_x, i, j);
        }
    }
}



void RefRenderer::advectVelocityForward() {
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            advectionCopyX[i][j] = velocitiesX[i][j];
            advectionCopyY[i][j] = velocitiesY[i][j];
        }
    }

    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int nextPixelRow = round(pixelRow + TIME_STEP * velocitiesY[row][col] * CELL_DIM);
           int nextPixelCol = round(pixelCol + TIME_STEP * velocitiesX[row][col] * CELL_DIM);
           int nextCellCol = nextPixelCol / CELL_DIM;
           int nextCellRow = nextPixelRow / CELL_DIM;

           if (nextCellCol < cells_per_side && nextCellRow < cells_per_side 
                   && nextCellCol >= 0 && nextCellRow >= 0) {
                velocitiesX[nextCellRow][nextCellCol] = advectionCopyX[row][col];
                velocitiesY[nextCellRow][nextCellCol] = advectionCopyY[row][col];
           } 
        }
    }

}

void
RefRenderer::advectVelocityBackward() {
    for (int row = 0; row < cells_per_side; row++) {
        for (int col = 0; col < cells_per_side; col++) {
           int pixelRow = row * CELL_DIM;// + CELL_DIM/2;
           int pixelCol = col * CELL_DIM;// + CELL_DIM/2;
           int prevPixelRow = round(pixelRow - TIME_STEP * velocitiesY[row][col] * CELL_DIM);
           int prevPixelCol = round(pixelCol - TIME_STEP * velocitiesX[row][col] * CELL_DIM);
           int prevCellCol = prevPixelCol / CELL_DIM;
           int prevCellRow = prevPixelRow / CELL_DIM;

           if (prevCellCol < cells_per_side && prevCellRow < cells_per_side 
                   && prevCellCol >= 0 && prevCellRow >= 0) {
                velocitiesX[row][col] = advectionCopyX[prevCellRow][prevCellCol];
                velocitiesY[row][col] = advectionCopyY[prevCellRow][prevCellCol];
           } 
           if (prevCellCol == col && prevCellRow == row) {
                // you don't move so just disappear
                velocitiesX[row][col] = 0;
                velocitiesY[row][col] = 0;
           }
        }
    }
}

// Divergence of velocity: This computes how divergent the velocity field is
// (how much in/out flow there is at every point).  Used as input to the 
// pressure solve below.
void
RefRenderer::applyDivergence() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            
            if (i > 0) T = velocitiesY[i-1][j];
            if (i < cells_per_side) B = velocitiesY[i+1][j];
            if (j < cells_per_side) R = velocitiesX[i][j+1];
            if (j > 0) L = velocitiesX[i][j-1];
            divergence[i][j] = 0.5*((R-L) + (T-B));   
        }
    }
}

void
RefRenderer::pressureSolve() {
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
            if (i > 0) T = pressures[i-1][j];
            if (i < cells_per_side) B = pressures[i+1][j];
            if (j < cells_per_side) R = pressures[i][j+1];
            if (j > 0) L = pressures[i][j-1];
            tempPressure[i][j] = (L + R + B + T + -1 * divergence[i][j]) * .25;
                //if (L+R+B+T > 0.0) printf("L+R+B+T is %f, -1*divergence is %f\n", L+R+B+T, -1*divergence[i][j]);
        }
    }
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) { 
            if (isBoundary(i,j)) continue;
            pressures[i][j] = tempPressure[i][j];
        }
    }

}

void
RefRenderer::pressureGradient() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) T = pressures[i-1][j];
            if (i < cells_per_side) B = pressures[i+1][j];
            if (j < cells_per_side) R = pressures[i][j+1];
            if (j > 0) L = pressures[i][j-1];

            //if (velocitiesY[i][j] > 10) printf("doing velocitiesY = %f - 0.5*(%f)\n", velocitiesY[i][j], T-B);
            velocitiesX[i][j] = velocitiesX[i][j] - 0.5*(R - L);
            velocitiesY[i][j] = velocitiesY[i][j] - 0.5*(T - B);
            //if (velocitiesY[i][j] != 0.0) printf("velocitiesY is %f\n", velocitiesY[i][j]);
        }
    }
}

void
RefRenderer::applyVorticity() {
    float L = 0.0;
    float R = 0.0;
    float B = 0.0;
    float T = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) T = velocitiesX[i-1][j];
            if (i < cells_per_side) B = velocitiesX[i+1][j];
            if (j < cells_per_side) R = velocitiesY[i][j+1];
            if (j > 0) L = velocitiesY[i][j-1];
            vorticity[i][j] = 0.5 * ((R - L) - (T - B));
        }
    }
}

void
RefRenderer::applyVorticityForce() {
    float vortConfinementFloat = 0.035f;
    float vortL = 0.0;
    float vortR = 0.0;
    float vortB = 0.0;
    float vortT = 0.0;
    float vortC = 0.0;
    for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            if (isBoundary(i,j)) continue;
            if (i > 0) vortT = vorticity[i-1][j];
            if (i < cells_per_side) vortB = vorticity[i+1][j];
            if (j < cells_per_side) vortR = vorticity[i][j+1];
            if (j > 0) vortL = vorticity[i][j-1];
            vortC = vorticity[i][j];
            float forceX = 0.5 * (abs(vortT) - abs(vortB));
            float forceY = 0.5 * (abs(vortR) - abs(vortL));
            float EPSILON = pow(2,-12);
            float magSqr = std::max(EPSILON, forceX * forceX + forceY * forceY);
            forceX = forceX * (1/sqrt(magSqr));
            forceY = forceY * (1/sqrt(magSqr));
            forceX *= vortConfinementFloat * vortC * 1;
            forceY *= vortConfinementFloat * vortC * -1;
            velocitiesX[i][j] += forceX;
            velocitiesY[i][j] += forceY;
        }
    }
}

void
RefRenderer::render() {
    advectVelocityForward();
    advectVelocityBackward();
    // set boundary conditions
    /*for (int i = 0; i < cells_per_side + 1; i++) {
        for (int j = 0; j < cells_per_side + 1; j++) {
            int side = isBoundary(i,j);
            if (side) {
                divergence[i][j] = 0.0;
                pressures[i][j] = 0.0;
                if (side == 1) {// left
                    velocitiesX[i][j] = -1*velocitiesX[i][j+1];
                }
                if (side == 2) { //top
                    velocitiesY[i][j] = -1*velocitiesY[i+1][j];
                }
                if (side == 3) { //right
                    velocitiesX[i][j] = -1*velocitiesX[i][j-1];
                }
                if (side == 4) { //bottom
                    velocitiesY[i][j] = -1*velocitiesY[i-1][j];
                }
            }
        }
    }*/
    applyVorticity();
    applyVorticityForce();
    applyDivergence();
    pressureSolve();
    pressureGradient();

    // Draw
    for (int i = 0; i < 4*image->width*image->height; i+=4) {
            int grid_row = ((i/4) / image->width) / (CELL_DIM);
            int grid_col = ((i/4) % image->width) / (CELL_DIM);
            double vx = velocitiesX[grid_row][grid_col];
            double vy = velocitiesY[grid_row][grid_col];
            double v = sqrt(vx * vx + vy * vy);
            double s = 2 * (1.0 / (1.0 + exp(-1.0 * v))) - 1;


            if (abs(v) < 0.00001) {
                color[grid_row][grid_col][0] = 0.0;
                color[grid_row][grid_col][1] = 0.0;//0.0392;
                color[grid_row][grid_col][2] = 0.0;//0.1098;
                color[grid_row][grid_col][3] = 1.0;
                image->data[i] = color[grid_row][grid_col][0];
                image->data[i+1] = color[grid_row][grid_col][1];
                image->data[i+2] = color[grid_row][grid_col][2];
                image->data[i+3] = color[grid_row][grid_col][3];
                continue;
            } 
            color[grid_row][grid_col][0] *= 0.9494; 
            color[grid_row][grid_col][1] *= 0.9494; 
            color[grid_row][grid_col][2] *= 0.9696; 
        
            if (isMouseDown){
                double d = sqrt(vx * vx + vy * vy);
                double projection;
                double l = distanceToNearestMousePoint(grid_col, grid_row, &projection);
                float taperFactor = 0.6;
                double projectedFraction = 1.0 - std::min(1.0, 
                        std::max(projection, 0.0)) * taperFactor;
                double R = 50; //0.025; // the bigger, the more stuff gets colored? not really, m is 
                                        // already maxed out to 1
                double m = exp(-l/R); //drag coefficient
                double speed = v;

                //printf("l is %f, m is %f, projection is %f\n", l, m, projection);

                double x = std::min(1.0, std::max(fabs((speed * speed * 0.02 - 
                            projection * 5.0) * projectedFraction), 0.0));

                double r = (2.4 / 60.0) * x + (0.2 /30.0) * (1-x) + (1.0 * pow(x, 5.0));
                double g = (0.0 / 60.0) * x + (51.8 / 30.0) * (1-x) + (1.0 * pow(x, 5.0));
                double b = (5.9 / 60.0) * x + (100.0 / 30.0) * (1-x) + (1.0 * pow(x, 5.0));

                color[grid_row][grid_col][0] += m * r;
                color[grid_row][grid_col][1] += m * g;
                color[grid_row][grid_col][2] += m * b;
                color[grid_row][grid_col][3] = 1.0;
            }

            image->data[i] = color[grid_row][grid_col][0];
            image->data[i+1] = color[grid_row][grid_col][1];
            image->data[i+2] = color[grid_row][grid_col][2];
            image->data[i+3] = color[grid_row][grid_col][3];

            //if (image->data[i] < 1.0) printf("RGB: %f %f %f\n", image->data[i], image->data[i+1], image->data[i+2]);
            //printf("m rgb: %f %f %f %f\n", m, r,g,b);

            /*if ((vx != 0.0 || vy != 0.0) && s >= 0.1) {
                image->data[i] = 0.0;
                image->data[i+1] = 0.0;
                image->data[i+2] = s;
                image->data[i+3] = 1.0;
                //if (velocitiesY[grid_row][grid_col] < 0.0) image->data[i] = 1;
            } else {
                image->data[i] = 0.0;
                image->data[i+1] = 0.0392;
                image->data[i+2] = 0.1098;
                image->data[i+3] = 1.0;

            }*/
    }
    advectColor();

}

