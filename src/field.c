#include "SDL2/SDL.h"
#include <time.h>
#include <stdio.h>

const unsigned WIDTH    = 512;
const unsigned HEIGHT   = 512;

const unsigned N        = 1000;
const unsigned iter     = 10000;

void init_surface(Uint32 *pixels)
{
    //init RNG
    time_t t;
    srand((unsigned) time(&t));
    
    for (unsigned i = 0; i < N; i++)
    {
        float f = ((float) rand()) / ((float) RAND_MAX);
        unsigned p = (unsigned) (f * WIDTH * HEIGHT);
        //printf("f: %f and p: %i at iteration %i\n", f, p, i);
        pixels[p] = 255*256*256;
    }
}

void update(Uint32 *source, Uint32 *target, const unsigned WIDTH, const unsigned HEIGHT)
{
    for (unsigned x = 0; x < WIDTH; x++)
    {
        for (unsigned y = 0; y < HEIGHT; y++)
        {
            target[y*WIDTH + x] = source[y * WIDTH + x];
        }
    }
}

void copy(Uint32 *source, Uint32 *target, const unsigned LENGTH)
{
    for (unsigned i = 0; i < LENGTH; i++)
    {
        target[i] = source[i];
    }
    
}

int main(int argc, char** argv)
{
    // --SET UP SDL -- //
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window *window          = SDL_CreateWindow("field", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, 0);
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex        = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Surface *surface1       = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    SDL_Surface *surface2       = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);

    float *field1               =(float *)malloc(WIDTH * HEIGHT * sizeof(float));

    init_surface(surface1->pixels); 

    for (unsigned i = 0; i < iter; i++)
    {
        update(surface1->pixels, surface2->pixels, WIDTH, HEIGHT);

        SDL_UpdateTexture(gpu_tex, NULL, surface2->pixels, WIDTH * sizeof(Uint32));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, gpu_tex, NULL, NULL);
        SDL_RenderPresent(renderer);

        copy(surface2->pixels, surface1->pixels, WIDTH * HEIGHT);
    }
    
    free(field1);
    SDL_FreeSurface(surface1);
    SDL_FreeSurface(surface2);
    SDL_DestroyTexture(gpu_tex);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}