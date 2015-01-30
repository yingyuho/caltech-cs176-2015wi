#ifndef DDG_SUBDIVISION_H
#define DDG_SUBDIVISION_H

namespace DDG {
   class Mesh;
   class Vector;
   class Face;
   class Edge;
   class Vertex;

   class Subdivision
   {
   public:
      static void build( const Mesh& src, Mesh& dest );
   protected:
      static Vector computeControlPoint ( const Face& face );
      static Vector computeControlPoint ( const Edge& edge );
      static Vector computeControlPoint ( const Vertex& vertex );
   };
}

#endif
