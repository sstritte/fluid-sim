#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>

#include "cudaRenderer.h"
#include "refRenderer.h"
#include "platformgl.h"

void usage(const char* progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Program Options:\n");
    printf("  -b  --baseline  Run the sequential baseline simulation\n");
    printf("  -h  --help                 This message\n");
}

void startRendererWithDisplay(CircleRenderer* renderer); // From display.cpp

int main(int argc, char** argv)
{
    int imageSize = 512;
    bool useBaseline = false;

    // parse commandline options ////////////////////////////////////////////
    int opt;
    static struct option long_options[] = {
        {"help",     0, 0,  'h'},
        {"baseline", 0, 0,  'b'},
        {0 ,0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "hb", long_options, NULL)) != EOF) {
        switch (opt) {
        case 'b':
            useBaseline = true;
            printf("Running baseline simulation...\n");
            break;
        case 'h':
        default:
            usage(argv[0]);
            return 1;
        }
    }
    // end parsing of commandline options //////////////////////////////////////

    printf("Rendering %dx%d simulation\n", imageSize, imageSize);

    CircleRenderer* renderer;
   
    if (useBaseline) 
        renderer = new RefRenderer();
    else
        renderer = new CudaRenderer();

    renderer->allocOutputImage(imageSize, imageSize);
    renderer->setup();

    glutInit(&argc, argv);
    startRendererWithDisplay(renderer);

    return 0;
}
