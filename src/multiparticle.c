
#include <SDL2/SDL.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <omp.h>

const unsigned WIDTH    = 500;
const unsigned HEIGHT   = 500;

const unsigned N_red    = 200;
const unsigned N_green  = 100;
const unsigned N_blue   = 200;

const float dt          = 0.0001f;

const unsigned iter     = 10000;

/*
const SDL_Rect dstrect1 = {0, 0, WIDTH * 2, HEIGHT * 2};
const SDL_Rect dstrect2 = {WIDTH * 2, 0, WIDTH * 2, HEIGHT * 2};
*/

enum Type{red, green, blue};
typedef enum Type Type;

struct Colorf{
    float r;
    float g;
    float b;
};
typedef struct Colorf Colorf;

struct Particle{
    float   x;
    float   dx;
    float   y;
    float   dy;
    Type    PType;
};
typedef struct Particle Particle;

Particle *create_particles(const unsigned N_red, const unsigned N_green, const unsigned N_blue)
{
    //init RNG
    time_t t;
    srand((unsigned) time(&t));

    unsigned N = N_red + N_green + N_blue;
    Particle *P = (Particle *)malloc(N * sizeof(Particle));

    for(unsigned i = 0; i < N_red; i++)
    {
        P[i].x      = ((float) rand()) / ((float) RAND_MAX);
        P[i].y      = ((float) rand()) / ((float) RAND_MAX);
        //put them more centered
        P[i].x      = (P[i].x / 4) + (3.0f / 8.0f);
        P[i].y      = (P[i].y / 4) + (3.0f / 8.0f);
        P[i].dx     = 0.0;
        P[i].dy     = 0.0;

        P[i].PType  = red;
    }
    for(unsigned i = N_red; i < N_red + N_green; i++)
    {
        P[i].x      = ((float) rand()) / ((float) RAND_MAX);
        P[i].y      = ((float) rand()) / ((float) RAND_MAX);
        //put them more centered
        P[i].x      = (P[i].x / 4) + (3.0f / 8.0f);
        P[i].y      = (P[i].y / 4) + (3.0f / 8.0f);
        P[i].dx     = 0.0;
        P[i].dy     = 0.0;

        P[i].PType  = green;
    }
    for(unsigned i = N_red + N_green; i < N; i++)
    {
        P[i].x      = ((float) rand()) / ((float) RAND_MAX);
        P[i].y      = ((float) rand()) / ((float) RAND_MAX);
        //put them more centered
        P[i].x      = (P[i].x / 4) + (3.0f / 8.0f);
        P[i].y      = (P[i].y / 4) + (3.0f / 8.0f);
        P[i].dx     = 0.0;
        P[i].dy     = 0.0;

        P[i].PType  = blue;
    }
    return P;
}
void destroy_particles(Particle *P)
{
    free(P);
}

float rockpaperscissors(Type iType, Type jType)
{
    //itype
    if (iType == red)
    {
        if (jType == red)
        {
            // red on red
            return 0.0f;
        }
        else if (jType == green)
        {
            // red chases green
            return 1.0f;
        }
        else if (jType == blue)
        {
            //red runs aways from blue
            return -1.0f;
        }        
    }
    else if (iType == green)
    {
        if (jType == red)
        {
            // green runs away from red
            return -1.0f;
        }
        else if (jType == green)
        {
            // greeen on green
            return 0.0f;
        }
        else if (jType == blue)
        {
            // green chases blue
            return 1.0f;
        }   
    }
    else if (iType == blue)
    {
        if (jType == red)
        {
            //blue chases red
            return 1.0f;
        }
        else if (jType == green)
        {
            // blue runs away from green
            return -1.0f;
        }
        else if (jType == blue)
        {
            //blue on blue
            return 0.0f;
        }   
    }
    else
        return 0.0f;   
    
}

float cycle(Type iType, Type jType)
{
    //itype
    if (iType == red)
    {
        if (jType == red)
        {
            return -1.0f;
        }
        else if (jType == green)
        {
            return 1.0f;
        }
        else if (jType == blue)
        {
            return 1.0f;
        }        
    }
    else if (iType == green)
    {
        if (jType == red)
        {
            return 0.0f;
        }
        else if (jType == green)
        {
            return 1.0f;
        }
        else if (jType == blue)
        {
            return 0.0f;
        }   
    }
    else if (iType == blue)
    {
        if (jType == red)
        {
            return 1.0f;
        }
        else if (jType == green)
        {
            return 1.0f;
        }
        else if (jType == blue)
        {
            return -1.0f;
        }   
    }
    else
        return 0.0f;   
    
}

double attract_vel(Particle *P, const float dt, const unsigned N)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

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
                Type iType = P[i].PType;
                Type jType = P[j].PType;
                float rps = cycle(iType, jType);
                x_force += rps * x_dist / sqrdDist;
                y_force += rps * y_dist / sqrdDist;
            }      
            //printf("we have types %i and %i with rps: %f\n", iType, jType, rps);    
        }
        P[i].dx += x_force * dt;
        P[i].dy += y_force * dt;
    }

    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

//should check bounds here
//this enforces a taurus geometry (it's also square)
double integrate(Particle *P, const float dt, const unsigned N)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    #pragma omp parallel for
    for (unsigned i = 0; i < N; i++)
    {
        P[i].x += P[i].dx * dt;
        P[i].y += P[i].dy * dt;
    }

    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

//set field doesn't check bounds
double set_field(const Particle *P, const unsigned N, Colorf *F, const unsigned WIDTH, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    memset(F, 0.0, WIDTH * HEIGHT * sizeof(Colorf));

    // -- LOOK OUT FOR RACE CONDITION -- //
    for (unsigned i = 0; i < N; i++)
    {
        unsigned x = (unsigned)(P[i].x * WIDTH);
        unsigned y = (unsigned)(P[i].y * HEIGHT);
        if (x >= 0 && x < WIDTH)
        {
            if (y >= 0 && y < HEIGHT)
            {
                if (P[i].PType == red)
                    F[y * WIDTH + x].r += 1.0;   
                else if (P[i].PType == green)
                    F[y * WIDTH + x].g += 1.0;
                else
                    F[y * WIDTH + x].b += 1.0;
            }
        }       
    }

    struct timespec wall_stop;
    clock_gettime(CLOCK_MONOTONIC, &wall_stop);      //requires POSIX
    double wall_time    = (wall_stop.tv_sec - wall_start.tv_sec);
    wall_time           += (wall_stop.tv_nsec - wall_start.tv_nsec) / 1000000000.0;
    return wall_time;
}

double draw_field(const Colorf *F, Uint32 *cpu_pix1, const unsigned WIDTH, const unsigned HEIGHT)
{
    struct timespec wall_start;              
    clock_gettime(CLOCK_MONOTONIC, &wall_start);     //CLOCK_MONOTONIC requires POSIX

    #pragma omp parallel for
    for (unsigned i = 0; i < WIDTH * HEIGHT; i++)
    {
        Uint32 redpart      = (Uint32) (F[i].r * 255);
        Uint32 greenpart    = (Uint32) (F[i].g * 256 * 255); 
        Uint32 bluepart     = (Uint32) (F[i].b * 256 * 256 * 255);
        cpu_pix1[i]         = redpart + greenpart + bluepart;
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
    #pragma omp parallel private(tid) shared(n_threads)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            n_threads = omp_get_num_threads();
            printf("total of %i threads\n", n_threads);
        }
    }

    omp_set_dynamic(0);
    omp_set_num_threads((n_threads * 7)/ 8);

    #pragma omp parallel private(tid) shared(n_threads)
    {
        tid = omp_get_thread_num();
        if (tid == 0)
        {
            n_threads = omp_get_num_threads();
            printf("total of %i threads\n", n_threads);
        }
    }

    // --SET UP SDL -- //
    SDL_Init(SDL_INIT_VIDEO);

    SDL_Window *window          = SDL_CreateWindow("base nbody", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH * 2, HEIGHT * 2, 0);   //the window
    SDL_Renderer *renderer      = SDL_CreateRenderer(window, -1, 0);    //pick the driver, -1 means init the first one supported with the flags 0
    SDL_Texture *gpu_tex1       = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ABGR8888, SDL_TEXTUREACCESS_STATIC, WIDTH, HEIGHT);
    Uint32 *cpu_pix1            = (Uint32 *)malloc(WIDTH * HEIGHT * sizeof(Uint32));
    memset(cpu_pix1, 0, WIDTH * HEIGHT * sizeof(Uint32)); //make it black

    // --SET UP SIM -- //
    Colorf *F                   = (Colorf *)malloc(WIDTH * HEIGHT * sizeof(Colorf));
    memset(F, 0.0, WIDTH * HEIGHT * sizeof(Colorf));

    Particle *P                 = create_particles(N_red, N_green, N_blue);

    const unsigned N            = N_red + N_green + N_blue;

    double sim_time     = 0;
    double cpudraw_time = 0;
    double sdldraw_time = 0;
    for (unsigned i = 0; i < iter; i++)
    {
        // -- SIM -- //
        sim_time += attract_vel(P, dt, N);
        sim_time += integrate(P, dt, N);

        // -- CPU DRAW -- //
        cpudraw_time += set_field(P, N, F, WIDTH, HEIGHT);
        cpudraw_time += draw_field(F, cpu_pix1, WIDTH, HEIGHT);

        // -- SDL DRAW -- //
        struct timespec sdldraw_start;
        clock_gettime(CLOCK_MONOTONIC, &sdldraw_start);

            SDL_UpdateTexture(gpu_tex1, NULL, cpu_pix1, WIDTH * sizeof(Uint32));
            SDL_RenderClear(renderer);
            SDL_RenderCopy(renderer, gpu_tex1, NULL, NULL);
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
    free(F);
    free(cpu_pix1);
    SDL_DestroyTexture(gpu_tex1);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    
    return 0;
}