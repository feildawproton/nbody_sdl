#include "SDL2/SDL.h"
#include <time.h>
#include <stdio.h>

#include <math.h>

const unsigned WIDTH    = 512;
const unsigned HEIGHT   = 512;

const unsigned N        = 1000;
const unsigned iter     = 10000;

void init_field(float *field)
{
    //init RNG
    time_t t;
    srand((unsigned) time(&t));
    
    for (unsigned i = 0; i < N; i++)
    {
        float f = ((float) rand()) / ((float) RAND_MAX);
        unsigned p = (unsigned) (f * WIDTH * HEIGHT);
        //printf("f: %f and p: %i at iteration %i\n", f, p, i);
        field[p] = -1.0f;
    }
}

int imax(const int a, const int b)
{
    if (a > b)
        return a;
    else
        return b;
}

int imin(const int a, const int b)
{
    if (a < b)
        return a;
    else   
        return b;
}

void update_field(const float *source, float *target, const unsigned WIDTH, const unsigned HEIGHT)
{
    for (int x = 0; x < WIDTH; x++)
    {
        for (int y = 0; y < HEIGHT; y++)
        {
            //need to check bounds
            //the way I check here creates hard boundaries
            int x_l = imax(0, x - 1);
            int x_r = imin(WIDTH - 1, x + 1);
            int y_u = imax(0, y - 1);
            int y_d = imin(HEIGHT - 1, y + 1);
            
            float lu = source[y_u * WIDTH + x_l];
            float up = source[y_u * WIDTH + x];
            float ru = source[y_u * WIDTH + x_r];

            float l = source[y * WIDTH + x_l];
            float c = source[y * WIDTH + x];
            float r = source[y * WIDTH + x_r];

            float ld = source[y_d * WIDTH + x_l];
            float dn = source[y_d * WIDTH + x];
            float rd = source[y_d * WIDTH + x_r];

            float sum = lu + up + ru + l + c + r + ld + dn + rd;
            float val = sum / 9.0f;

            target[y*WIDTH + x] = val;
        }
    }
}

void copy_field(const float *source, float *target, const unsigned LENGTH)
{
    for (unsigned i = 0; i < LENGTH; i++)
    {
        target[i] = source[i];
    }
    
}

void draw_surface(const float *field, Uint32 *pixels, const unsigned LENGTH)
{
    for (unsigned i = 0; i < LENGTH; i++)
    {
        pixels[i] = (field[i] * -1.0 * 255 * 256 * 256);
    }
    
}

int main(int argc, char** argv)
{
    // --SET UP SDL -- //
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window *window          = SDL_CreateWindow("field", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, 0);
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex        = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Surface *surface        = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    memset(surface->pixels, 0, WIDTH * HEIGHT * sizeof(Uint32));

    float *field1               =(float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(field1, 0.0f, WIDTH * HEIGHT * sizeof(float));
    float *field2               =(float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(field2, 0.0f, WIDTH * HEIGHT * sizeof(float));

    init_field(field1);

    for (unsigned i = 0; i < iter; i++)
    {
        update_field(field1, field2, WIDTH, HEIGHT);

        draw_surface(field2, surface->pixels, WIDTH*HEIGHT);

        SDL_UpdateTexture(gpu_tex, NULL, surface->pixels, WIDTH * sizeof(Uint32));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, gpu_tex, NULL, NULL);
        SDL_RenderPresent(renderer);

        copy_field(field2, field1, WIDTH * HEIGHT);
    }
    
    free(field2);
    free(field1);
    SDL_FreeSurface(surface);
    SDL_DestroyTexture(gpu_tex);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}