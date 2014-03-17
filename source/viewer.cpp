#include "viewer.h"
#include "xychain.h"

#ifdef VISUALS
 #ifdef LINUX
  #include <SDL/SDL.h>
  #include <SDL/SDL_video.h>
  #include <GL/gl.h>
  #include <GL/glu.h>
 #endif  // LINUX

 #ifdef OSX
  #include <SDL.h>
  #include <SDL_video.h>
  #include <OpenGL/gl.h>
  #include <OpenGL/glu.h>
 #endif  // OSX

#endif  // VISUALS

#ifdef PNG_VIS
 #include <png++/png.hpp>
#endif  // PNG_VIS

#include <memory>
#include <string>
using std::string;

XYChain * Viewer::fgLattice = 0;
float Viewer::fgBoxSize = 0;
int Viewer::fgWidth = 0;
int Viewer::fgHeight = 0;
std::unique_ptr<unsigned char[]> Viewer::fgBuffer = 0;
int Viewer::fgPngFrameNum = 0;

void Viewer::Init(XYChain * l, int width, int height, float bsize)
{
   fgLattice = l;
   fgBoxSize = bsize;
   fgWidth = width;
   fgHeight = height;

   fgBuffer.reset(new unsigned char[3*width*height]);

   // New window
   if (SDL_Init(SDL_INIT_VIDEO) < 0) return;
   const SDL_VideoInfo * info = SDL_GetVideoInfo();
   if (!info) return;
   int bpp = info->vfmt->BitsPerPixel;

   SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 5 );
   SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 5 );
   SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 5 );
   SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 16 );
   SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
   if ( SDL_SetVideoMode( fgWidth, fgHeight, bpp, SDL_OPENGL ) == 0 ) {
      // Error
      return;
   }

   // By now GL should be available for use
   glClearColor(0,0,0,0); // Black
   glViewport(0,0,fgWidth, fgHeight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, fgWidth, fgHeight, 0, 0, 1);
   glMatrixMode(GL_MODELVIEW);
   glDisable(GL_DEPTH_TEST);
}

void Viewer::Render()
{
   if (!fgLattice) return;

   glClear(GL_COLOR_BUFFER_BIT);
   fgLattice->GLRender(fgBoxSize);
   SDL_GL_SwapBuffers();

   SDL_Event e;
   bool pause = false;
   while(SDL_PollEvent(&e)) {
      if (e.key.keysym.sym == SDLK_ESCAPE || e.type == SDL_QUIT) Destroy();
      if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_SPACE) {
         pause = true;
      }
      while (pause == true) {
         while (SDL_PollEvent(&e)) {
            if (e.type == SDL_KEYUP && e.key.keysym.sym == SDLK_SPACE)
               pause = false;
         }
      }
   }
}

void Viewer::RenderToPng(const string& prefix) {
   char fname[256];
   sprintf(fname, "%s_%i.png", prefix.c_str(), fgPngFrameNum++);
   
   // Read the GL window RGB pixels into fgBuffer
   glReadPixels(0, 0, fgWidth, fgHeight, GL_RGB, GL_UNSIGNED_BYTE, fgBuffer.get());

   png::image< png::rgb_pixel > img(fgWidth, fgHeight);
   for (int i = 0; i < fgWidth; i++) {
      for (int j = 0; j < fgHeight; j++) {
         int offset = 3*(fgWidth*j + i);
         img[fgHeight-j-1][i] = png::rgb_pixel(fgBuffer[offset+0],  // R
                                               fgBuffer[offset+1],  // G
                                               fgBuffer[offset+2]); // B
      }
   }

   img.write(fname);
}

void Viewer::Resize(int w, int h)
{
   fgWidth = w;
   fgHeight = h;
   glViewport(0,0,fgWidth, fgHeight);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(0, fgWidth, fgHeight, 0, 0, 1);
   glMatrixMode(GL_MODELVIEW);
   glDisable(GL_DEPTH_TEST);

   fgBuffer.reset(new unsigned char[3*w*h]);
}

void Viewer::Destroy()
{
   fgLattice = 0;
   SDL_Quit();
}
