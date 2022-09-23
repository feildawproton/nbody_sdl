
#include "SDL2/SDL.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <omp.h>

const unsigned WIDTH = 512;
const unsigned HEIGHT = 512;

const unsigned N = 1000;
const float dt  = 0.0001f;

const unsigned iter = 10000;

const SDL_Rect dstrect1 = {0, 0, WIDTH, HEIGHT};
const SDL_Rect dstrect2 = {WIDTH, 0, WIDTH, HEIGHT};

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
        P[i].x  = (P[i].x / 4) + (3.0f / 8.0f);
        P[i].y  = (P[i].y / 4) + (3.0f / 8.0f);
        P[i].dx = 0.0;
        P[i].dy = 0.0;
    }
    return P;
}
void destroy_particles(Particle *P)
{
    free(P);
}

void attract_vel(Particle *P, const float dt, const unsigned N)
{
    //unsigned i = 0;
    #pragma omp parallel for
    for (unsigned i = 0; i < N; i++)
    {
        float x_force = 0.0f;
        float y_force = 0.0f;
        for (unsigned j = 0; j < N; j++)
        {
            float x_dist    = P[j].x - P[i].x;
            float y_dist    = P[j].y - P[i].y;
            float sqrdDist  = x_dist * x_dist + y_dist * y_dist;           
            /*
            float invDist      = 1 / sqrtf(sqrdDist);           //cuda reference used rsqrtf() instead of 1 / sqrtf()
            float invDist2  = invDist * invDist;
            x_force += x_dist * invDist2;
            y_force += y_dist * invDist2;
            //for 3d the reference is
            float invDist = rsqrt(sqrdDist);                    //means 1 / sqrtf()
            float invDist3 = invDist * invDist * invDist;
            x_force = x_dist * invDist3;
            y_force = y_dist * invDist3;
            z_force = z_dist * invDist3;
            */
           //need to check if we are evaluating ourselvses
           //or if the particles are on top of eachother (unlickely but possible i suppose)
            if (sqrdDist == 0.0f)
            {
                x_force += 0.0f;
                y_force += 0.0f;
            }
            else
            {
                x_force += x_dist / sqrdDist;
                y_force += y_dist / sqrdDist;
            }           
        }
        P[i].dx += x_force * dt;
        P[i].dy += y_force * dt;
    }
}

//should check bounds here
//this enforces a taurus geometry (it's also square)
void integrate(Particle *P, const float dt, const unsigned N)
{
    #pragma omp parallel for
    for (unsigned i = 0; i < N; i++)
    {
        P[i].x += P[i].dx * dt;
        P[i].y += P[i].dy * dt;

        //need to do loops in case it's way out there
        /*
        while (P[i].x > 1.0f)
        {
            P[i].x = P[i].x - 1.0f;
        }
        while (P[i].x < 0.0f)
        {
            P[i].x = P[i].x + 1.0f;
        }
        while (P[i].y > 1.0f)
        {
            P[i].y = P[i].y - 1.0f;
        }
        while (P[i].y < 0.0f)
        {
            P[i].y = P[i].y + 1.0f;
        }
        */
    }
}

//set field doesn't check bounds
void set_fields(const Particle *P, const unsigned N, float *F_1, float *F_2, const unsigned WIDTH, const unsigned HEIGHT)
{
    memset(F_1, 0.0, WIDTH * HEIGHT * sizeof(float));
    //don't reset F_2
    //care should be when parallelizing this to avoid race conditions
    //this is like reduce-sum but with a field instead of a single variable
    for (unsigned i = 0; i < N; i++)
    {
        unsigned x = (unsigned)(P[i].x * WIDTH);
        unsigned y = (unsigned)(P[i].y * HEIGHT);
        if (x >= 0 && x < WIDTH)
        {
            if (y >= 0 && y < HEIGHT)
            {
                F_1[y * WIDTH + x] += 1.0;
                F_2[y * WIDTH + x] += 1.0;
            }
            
        }       
    }
}

void draw_fields(const float *F_1, const float *F_2, Uint32 *cpu_pix1, Uint32 *cpu_pix2, const unsigned WIDTH, const unsigned HEIGHT)
{
    #pragma omp parallel for
    for (unsigned i = 0; i < WIDTH * HEIGHT; i++)
    {
        cpu_pix1[i]     = (Uint32) (F_1[i] * 256 * 255); 

        Uint32 red      = (Uint32) (F_2[i]);
        Uint32 green    = (Uint32) (F_2[i] / 256);
        Uint32 blue     = (Uint32) (F_2[i] / (256 * 256));
        cpu_pix2[i]     = red + green + blue; 
        //cpu_pixels[i] = 0xFF010101;
    }
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

    SDL_Window *window          = SDL_CreateWindow("base nbody", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH * 2, HEIGHT, 0);   //the window
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex1       = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    SDL_Texture *gpu_tex2       = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    Uint32 *cpu_pix1            = (Uint32 *)malloc(WIDTH * HEIGHT * sizeof(Uint32));
    memset(cpu_pix1, 0, WIDTH * HEIGHT * sizeof(Uint32)); //make it black
    Uint32 *cpu_pix2            = (Uint32 *)malloc(WIDTH * HEIGHT * sizeof(Uint32));
    memset(cpu_pix2, 0, WIDTH * HEIGHT * sizeof(Uint32)); //make it black

    // --SET UP SIM -- //
    float *F_1                  = (float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(F_1, 0.0, WIDTH * HEIGHT * sizeof(float));
    float *F_2                  = (float *)malloc(WIDTH * HEIGHT * sizeof(float));
    memset(F_2, 0.0, WIDTH * HEIGHT * sizeof(float));
    
    SDL_PixelFormat *fmt;                    
    fmt->format = SDL_PIXELFORMAT_ABGR8888;

    Particle *P                 = create_particles();

    double sim_time     = 0;
    double cpudraw_time = 0;
    double sdldraw_time = 0;
    for (unsigned i = 0; i < iter; i++)
    {
        // -- SIM -- //
        struct timespec sim_start;              
        clock_gettime(CLOCK_MONOTONIC, &sim_start);     //CLOCK_MONOTONIC requires POSIX

            attract_vel(P, dt, N);
            integrate(P, dt, N);

        struct timespec sim_stop;
        clock_gettime(CLOCK_MONOTONIC, &sim_stop);      //requires POSIX
        sim_time += (sim_stop.tv_sec - sim_start.tv_sec);
        sim_time += (sim_stop.tv_nsec - sim_start.tv_nsec) / 1000000000.0;

        // -- CPU DRAW -- //
        struct timespec cpudraw_start;
        clock_gettime(CLOCK_MONOTONIC, &cpudraw_start); 

            set_fields(P, N, F_1, F_2, WIDTH, HEIGHT);
            draw_fields(F_1, F_2, cpu_pix1, cpu_pix2, WIDTH, HEIGHT);

        struct timespec cpudraw_stop;
        clock_gettime(CLOCK_MONOTONIC, &cpudraw_stop);
        cpudraw_time += (cpudraw_stop.tv_sec - cpudraw_start.tv_sec);
        cpudraw_time += (cpudraw_stop.tv_nsec - cpudraw_start.tv_nsec) / 1000000000.0;

        // -- SDL DRAW -- //
        struct timespec sdldraw_start;
        clock_gettime(CLOCK_MONOTONIC, &sdldraw_start);

            SDL_UpdateTexture(gpu_tex1, NULL, cpu_pix1, WIDTH * sizeof(Uint32));
            SDL_UpdateTexture(gpu_tex2, NULL, cpu_pix2, WIDTH * sizeof(Uint32));
            SDL_RenderClear(renderer);
            SDL_RenderCopy(renderer, gpu_tex1, NULL, &dstrect1);
            SDL_RenderCopy(renderer, gpu_tex2, NULL, &dstrect2);
            SDL_RenderPresent(renderer);

        struct timespec sdldraw_stop;
        clock_gettime(CLOCK_MONOTONIC, &sdldraw_stop);
        sdldraw_time += (sdldraw_stop.tv_sec - sdldraw_start.tv_sec);
        sdldraw_time += (sdldraw_stop.tv_nsec - sdldraw_start.tv_nsec) / 1000000000.0;
    }
    printf("Total time in sim: %f for %u iterations\n", sim_time, iter);
    printf("Total time in cpu draw: %f for %u iterations\n", cpudraw_time, iter);
    printf("Total time in sdl draw: %f for %u iterations\n", sdldraw_time, iter);

    destroy_particles(P);
    free(F_2);
    free(F_1);
    free(cpu_pix2);
    free(cpu_pix1);
    SDL_DestroyTexture(gpu_tex2);
    SDL_DestroyTexture(gpu_tex1);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}