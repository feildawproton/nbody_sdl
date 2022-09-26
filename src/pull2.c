#include "SDL2/SDL.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>

#include <omp.h>

const unsigned WIDTH    = 512;
const unsigned HEIGHT   = 512;

const unsigned N        = 1000;
const unsigned ITER     = 10000;

const SDL_Rect dstrect1 = {0, 0, WIDTH, HEIGHT};
const SDL_Rect dstrect2 = {WIDTH, 0, WIDTH, HEIGHT};

const double InvSqrt2 = 0.70710678118;

struct Particle{
    float x;
    float dx;
    float y;
    float dy;
};
typedef struct Particle Particle;

struct Vec2{
    float x;
    float y;
};
typedef struct Vec2 Vec2;

Particle *create_particles()
{
    //init RNG
    time_t t;
    srand((unsigned) time(&t));

    Particle *P = (Particle *)malloc(N * sizeof(Particle));
    for(unsigned i = 0; i < N; i++)
    {
        P[i].x  = ((float) rand()) / ((float) RAND_MAX);
        P[i].y  = ((float) rand()) / ((float) RAND_MAX);
        //put them more centered
        //P[i].x  = (P[i].x / 4) + (3.0f / 8.0f);
        //P[i].y  = (P[i].y / 4) + (3.0f / 8.0f);
        P[i].x  = (P[i].x / 2) + (1.0f / 4.0f);
        P[i].y  = (P[i].y / 2) + (1.0f / 4.0f);
        P[i].dx = 0.0;
        P[i].dy = 0.0;
    }
    return P;
}
void destroy_particles(Particle *P)
{
    free(P);
}

//set field doesn't check bounds
double particle2field(const Particle *P, const unsigned N, float *F, const unsigned WIDTH, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);

    memset(F, 0.0, WIDTH * HEIGHT * sizeof(float));

    // -- LOOK OUT FOR RACE CONDITION -- //
    for (unsigned i = 0; i < N; i++)
    {
        unsigned x = (unsigned)(P[i].x * WIDTH);
        unsigned y = (unsigned)(P[i].y * HEIGHT);
        if (x >= 0 && x < WIDTH)
        {
            if (y >= 0 && y < HEIGHT)
            {
                F[y * WIDTH + x] += 1.0;
            }
            
        }       
    }

    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

double draw_field(const float *field, Uint32 *pixels, const unsigned LENGTH)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    #pragma omp parallel for
    for (unsigned i = 0; i < LENGTH; i++)
    {
        pixels[i] = (field[i] * 255 * 256 * 256);
    }
    
    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX

    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

float sample_hardboundary(const float *source, const int x, const int y, const unsigned WIDTH, const unsigned HEIGHT)
{
    int x_this = x;
    int y_this = y;

    if (x < 0)
        x_this = 0;
    else if (x >= WIDTH)
        x_this = WIDTH - 1;
    
    if (y < 0)
        y_this = 0;
    else if (y >= HEIGHT)
        y_this = HEIGHT - 1;

    float val = source[y_this * WIDTH + x_this];

    return val;
}
float sample_zeropad(const float *source, const int x, const int y, const unsigned WIDTH, const unsigned HEIGHT)
{
    float val = 0.0f;
    if ( x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT)
        val = source[y * WIDTH + x];
    return val;
}

/*
In my mental model the entries of the field are connected by edges.
The edges are what convey information (such as pull) between elements.
This can be summarized as a vecter at each element
That is the sum of the effects on the edges
These calculations don't consider the value at the entry, because that value does not pull on the entry
*/
double set_pull(const float *field, Vec2 *pullfield, const unsigned WIDHT, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    memset(pullfield, 0.0f, WIDTH * HEIGHT * sizeof(Vec2));

    #pragma omp parallel for
    for (int  x = 0; x < WIDTH; x++)
    {
        for (int y = 0; y < HEIGHT; y++)
        {
            float lu = sample_zeropad(field, x - 1, y - 1, WIDTH, HEIGHT);
            float up = sample_zeropad(field, x    , y - 1, WIDTH, HEIGHT);
            float ru = sample_zeropad(field, x + 1, y - 1, WIDTH, HEIGHT);

            float l  = sample_zeropad(field, x - 1, y    , WIDTH, HEIGHT);
            float c  = sample_zeropad(field, x    , y    , WIDTH, HEIGHT);
            float r  = sample_zeropad(field, x + 1, y    , WIDTH, HEIGHT);

            float ld = sample_zeropad(field, x - 1, y + 1, WIDTH, HEIGHT);
            float dn = sample_zeropad(field, x    , y + 1, WIDTH, HEIGHT);
            float rd = sample_zeropad(field, x + 1, y + 1, WIDTH, HEIGHT);

            pullfield[y * WIDTH + x].x = - (lu * InvSqrt2) - l - (ld * InvSqrt2) + (ru * InvSqrt2) + r + (rd * InvSqrt2);
            pullfield[y * WIDTH + x].y = - (lu * InvSqrt2) - up - (ru * InvSqrt2) + (ld * InvSqrt2) + dn + (rd * InvSqrt2);
        }
        
    }
    
    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

double update_field(const float *source, float *target, const unsigned WIDTH, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX
    
    #pragma omp parallel for
    for (int x = 0; x < WIDTH; x++)
    {
        for (int y = 0; y < HEIGHT; y++)
        {
            //need to check bounds
            /*
            //the way I check here creates hard boundaries

            float lu = sample_hardboundary(source, x - 1, y - 1, WIDTH, HEIGHT);
            float up = sample_hardboundary(source, x    , y - 1, WIDTH, HEIGHT);
            float ru = sample_hardboundary(source, x + 1, y - 1, WIDTH, HEIGHT);

            float l  = sample_hardboundary(source, x - 1, y    , WIDTH, HEIGHT);
            float c  = sample_hardboundary(source, x    , y    , WIDTH, HEIGHT);
            float r  = sample_hardboundary(source, x + 1, y    , WIDTH, HEIGHT);

            float ld = sample_hardboundary(source, x - 1, y + 1, WIDTH, HEIGHT);
            float dn = sample_hardboundary(source, x    , y + 1, WIDTH, HEIGHT);
            float rd = sample_hardboundary(source, x + 1, y + 1, WIDTH, HEIGHT);
            */

            //float lu = sample_zeropad(source, x - 1, y - 1, WIDTH, HEIGHT);
            float up = sample_zeropad(source, x    , y - 1, WIDTH, HEIGHT);
            //float ru = sample_zeropad(source, x + 1, y - 1, WIDTH, HEIGHT);

            float l  = sample_zeropad(source, x - 1, y    , WIDTH, HEIGHT);
            float c  = sample_zeropad(source, x    , y    , WIDTH, HEIGHT);
            float r  = sample_zeropad(source, x + 1, y    , WIDTH, HEIGHT);

            //float ld = sample_zeropad(source, x - 1, y + 1, WIDTH, HEIGHT);
            float dn = sample_zeropad(source, x    , y + 1, WIDTH, HEIGHT);
            //float rd = sample_zeropad(source, x + 1, y + 1, WIDTH, HEIGHT);

            //float sum = lu + up + ru + l + c + r + ld + dn + rd;
            float sum = up + l + c + r + dn;

            //float val = sum / 9.0f;

            float val = sum / 5.0f;

            target[y*WIDTH + x] = val;
        }
    }

    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX

    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

double copy_field(const float *source, float *target, const unsigned LENGTH)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    #pragma omp parallel for
    for (unsigned i = 0; i < LENGTH; i++)
    {
        target[i] = source[i];
    }
    
    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX

    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}



int main(int argc, char** argv)
{
    // -- TEST OMP -- //
    int n_threads;
    int tid;
    #pragma omp parallel private(n_threads, tid)
    {
        tid = omp_get_thread_num();
        printf("hello from thread %i\n", tid);
        if (tid == 0)
        {
            n_threads = omp_get_num_threads();
            printf("total of %i threads\n", n_threads);
        }
        
    }
    // --SET UP SDL -- //
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window *window          = SDL_CreateWindow("field", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH * 2, HEIGHT, 0);
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex1       = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Texture *gpu_tex2       = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Surface *surface1       = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    memset(surface1->pixels, 0, WIDTH * HEIGHT * sizeof(Uint32));
    SDL_Surface *surface2       = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    memset(surface2->pixels, 0, WIDTH * HEIGHT * sizeof(Uint32));

    // -- Generate particles -- //
    Particle *particles = create_particles();
    float *particlefield        = (float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(particlefield, 0.0f, WIDTH * HEIGHT * sizeof(float));

    Vec2 *pullfield             = (Vec2 *)malloc(WIDTH * HEIGHT * sizeof(Vec2));
    memset(pullfield, 0.0f,WIDTH * HEIGHT * sizeof(Vec2));

    double part2field_time      = 0.0f;
    double particledraw_time    = 0.0f;
    double setpull_time         = 0.0f;
    double sdldraw_time         = 0.0f;
    for (unsigned i = 0; i < ITER; i++)
    {
        part2field_time     += particle2field(particles, N, particlefield, WIDTH, HEIGHT);
        particledraw_time   += draw_field(particlefield, surface1->pixels, WIDTH * HEIGHT);
        setpull_time        += set_pull(particlefield, pullfield, WIDTH, HEIGHT);

        struct timespec wall_start;              
        clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

        SDL_UpdateTexture(gpu_tex1, NULL, surface1->pixels, WIDTH * sizeof(Uint32));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, gpu_tex1, NULL, &dstrect1);
        SDL_RenderPresent(renderer);

        struct timespec wall_stop;
        clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
        sdldraw_time        += (wall_stop.tv_sec - wall_start.tv_sec);
        sdldraw_time        += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    }

    printf("total time spent converting particles to a field: %f over %i iterations\n", part2field_time, ITER);
    printf("total time spent drawing particle field to a cpu surface: %f over %i iterations\n", particledraw_time, ITER);
    printf("total time spent setting pull field: %f over %i iterations\n", setpull_time, ITER);
    printf("total time spent drawing with sdl (GPU): %f over %i iterations\n", sdldraw_time, ITER);
    
    free(particlefield);
    destroy_particles(particles);
    SDL_FreeSurface(surface2);
    SDL_FreeSurface(surface1);
    SDL_DestroyTexture(gpu_tex2);
    SDL_DestroyTexture(gpu_tex1);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}