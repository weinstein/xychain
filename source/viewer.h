#ifndef VIEWER_H
#define VIEWER_H

#include <memory>

class XYChain;
//class std::string;

class Viewer {
   public:
      /**
       * Initialize SDL window and OpenGL context for rendering a given lattice
       * @param l Pointer to a TLattice to draw
       * @param w Width of the window
       * @param h Height of the window
       */
      static void Init(XYChain * l, int w, int h);
      
      /**
       * Initialize SDL window and OpenGL context for rendering a given lattice
       * @param l Pointer to a TLattice to draw
       * @param w Width of the window
       * @param h Height of the window
       * @param scale Drawing scale to use when rendering (default 1).
       * Larger scales mean more "zoomed in" rendering.
       */
      static void Init(XYChain * l, int w, int h, float scale);
      
      /**
       * Set the current lattice to be rendered. This allows a new lattice to be rendered
       * without recreating the window.
       * @param l Pointer to the new TLattice to render
       */
      static void SetLattice(XYChain * l) { fgLattice = l; }
      
      /**
       * Return a pointer to the currently selected lattice for rendering/
       */
      static XYChain * GetLattice() { return fgLattice; }
      
      /**
       * Render the currently selected lattice to the current window.
       * @param field If true, an arrow representing the external field value is also drawn.
       */
      static void Render();

      static void RenderToPng(const std::string& prefix);
      
      /**
       * Resize the current window
       * @param w New width
       * @param h New height
       */
      static void Resize(int w, int h);

      
      /**
       * Destroy the window and clean up
       */
      static void Destroy();

      static float fgBoxSize;

   private:
      static XYChain * fgLattice;
      static int fgWidth;
      static int fgHeight;
      static std::unique_ptr<unsigned char[]> fgBuffer;
      static int fgPngFrameNum;
};

#endif
