#include "SDL2/SDL.h"
#include <time.h>
#include <stdio.h>

#include <math.h>
#include <time.h>

#include <omp.h>

const unsigned WIDTH    = 512;
const unsigned HEIGHT   = 512;

const unsigned N        = 1000;
const unsigned ITER     = 10000;

struct Particle{
    float x;
    float dx;
    float y;
    float dy;
};
typedef struct Particle Particle;

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
double set_field(const Particle *P, const unsigned N, float *F, const unsigned WIDTH, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    //memset(F, 0.0, WIDTH * HEIGHT * sizeof(float));

    //firsty clear the locations where the particles will be.
    //this is not right
    #pragma omp parallel for
    for (unsigned i = 0; i < N; i++)
    {
        unsigned x = (unsigned)(P[i].x * WIDTH);
        unsigned y = (unsigned)(P[i].y * HEIGHT);
        if (x >= 0 && x < WIDTH)
        {
            if (y >= 0 && y < HEIGHT)
            {
                F[y * WIDTH + x] = 0.0;
            }
            
        }       
    }

    // -- LOOK OUT FOR RACE CONDITION -- //
    //need to reduce sum
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
            //float c  = sample_zeropad(source, x    , y    , WIDTH, HEIGHT);
            float r  = sample_zeropad(source, x + 1, y    , WIDTH, HEIGHT);

            //float ld = sample_zeropad(source, x - 1, y + 1, WIDTH, HEIGHT);
            float dn = sample_zeropad(source, x    , y + 1, WIDTH, HEIGHT);
            //float rd = sample_zeropad(source, x + 1, y + 1, WIDTH, HEIGHT);

            //float sum = lu + up + ru + l + c + r + ld + dn + rd;
            float sum = up + l + r + dn;

            //float val = sum / 9.0f;

            float val = sum / 4.0f;

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

double draw_surface(const float *field, Uint32 *pixels, const unsigned LENGTH)
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

    SDL_Window *window          = SDL_CreateWindow("field", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH * 2, HEIGHT * 2, 0);
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex        = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Surface *surface        = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    memset(surface->pixels, 0, WIDTH * HEIGHT * sizeof(Uint32));

    float *field1               =(float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(field1, 0.0f, WIDTH * HEIGHT * sizeof(float));
    float *field2               =(float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(field2, 0.0f, WIDTH * HEIGHT * sizeof(float));

    Particle *particles = create_particles();

    double setfield_time    =  0.0f;
    double update_time      = 0.0f;
    double cpudraw_time     = 0.0f;
    double fieldcopy_time   = 0.0f;
    for (unsigned i = 0; i < ITER; i++)
    {
        setfield_time   += set_field(particles, N, field1, WIDTH, HEIGHT); 
        update_time     += update_field(field1, field2, WIDTH, HEIGHT);

        cpudraw_time    += draw_surface(field2, surface->pixels, WIDTH*HEIGHT);

        SDL_UpdateTexture(gpu_tex, NULL, surface->pixels, WIDTH * sizeof(Uint32));
        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, gpu_tex, NULL, NULL);
        SDL_RenderPresent(renderer);

        fieldcopy_time  += copy_field(field2, field1, WIDTH * HEIGHT);
    }
    printf("Time spent setting field: %f for %i iterations\n", update_time, ITER);
    printf("Time spent on updating the field: %f\n", update_time);
    printf("Time spent on cpu draw: %f\n", cpudraw_time);
    printf("Time spent on field copy: %f\n", fieldcopy_time);
    
    destroy_particles(particles);
    free(field2);
    free(field1);
    SDL_FreeSurface(surface);
    SDL_DestroyTexture(gpu_tex);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}