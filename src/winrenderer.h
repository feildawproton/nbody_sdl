#ifndef WINRENDERER_H
#define WINRENDERER_H

#include "SDL2/SDL.h"

struct WinRenderer
{
    SDL_Window *window;
    SDL_Renderer *renderer;
    SDL_Texture *gpu_tex1;
};
typedef struct WinRenderer WinRenderer;

__inline__ WinRenderer create_WinRenderer()
{
    
}

#endif