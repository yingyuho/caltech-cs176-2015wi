#include <cmath>
#include <vector>
#include <set>

#include "Subdivision.h"
#include "Mesh.h"
#include "MeshIO.h"
#include "Vector.h"

using namespace std;

namespace DDG
{
   static const double PI = 3.141592653589793238463;

   static double theta_k ( int k )
   {
      switch (k)
      {
         case 2:
            return 0.0;
         case 3:
            return 0.5;
         default:
            return cos( PI / k );
      }
   }

   void Subdivision :: build( const Mesh& src, Mesh& dest )
   {
      const vector<Face>& meshF = src.faces;
      const vector<Edge>& meshE = src.edges;
      const vector<Vertex>& meshV = src.vertices;
      const FaceCIter meshF0 = meshF.begin();
      const EdgeCIter meshE0 = meshE.begin();
      const VertexCIter meshV0 = meshV.begin();
      const HalfEdgeCIter meshH0 = src.halfedges.begin();
      const HalfEdgeCIter meshH1 = src.halfedges.end();
      const unsigned int 
         sizeF = meshF.size(), 
         sizeE = meshE.size(),
         sizeV = meshV.size();

      MeshData data;
      vector<Vector>& ctlpts = data.positions;

      ctlpts.resize(sizeF + sizeE + sizeV);

      vector<int> degreesF(sizeF, 0);
      // vector<int> degreesV2(sizeV, 0);
      vector<int> degreesV(sizeV, 0);      

      Vector
         *ctlptsF = &ctlpts[0],
         *ctlptsE = &ctlpts[sizeF],
         *ctlptsV = &ctlpts[sizeF + sizeE];

      // count degrees and compute unaveraged face points
      for( HalfEdgeCIter he  = meshH0;
                         he != meshH1; 
                         he++ )
      {
         if( !he->onBoundary )
         {
            unsigned int idxF = he->face - meshF0;
            unsigned int idxV = he->vertex - meshV0;
            ctlptsF[idxF] += he->vertex->position;
            ++degreesF[idxF];
            ++degreesV[idxV];
         }
      }

      // compute averaged face points
      for( unsigned int idxF = 0; idxF < sizeF; idxF++ )
      {
         ctlptsF[idxF] /= degreesF[idxF];
      }

      for( HalfEdgeCIter he  = meshH0;
                         he != meshH1; 
                         he++ )
      {
         unsigned int idxF = he->face - meshF0;
         unsigned int idxE = he->edge - meshE0;
         unsigned int idxV = he->vertex - meshV0;

         ctlptsE[idxE] += he->vertex->position;

         if( !he->onBoundary )
         {
            HalfEdgeCIter he2 = he->next;
            if( !he->flip->onBoundary ) 
            {
               ctlptsE[idxE] += ctlptsF[idxF];
               ctlptsV[idxV] += ctlptsF[idxF] + he2->vertex->position;
               // degreesV2[idxV] += 2;
            }
            else
            {
               // corrections for edges near boundaries
               unsigned int idxE = he2->edge - meshE0;
               unsigned int idxV = he2->vertex - meshV0;
               ctlptsE[idxE] += 
                  theta_k( degreesV[idxV] ) * 
                  (he2->vertex->position - he2->next->vertex->position);
            }
         }
      }

      for( HalfEdgeCIter he  = meshH0;
                         he != meshH1; 
                         he++ )
      {
         if( !he->onBoundary )
         {
            // compute averaged edge points
            unsigned int idxE = he->edge - meshE.begin();
            ctlptsE[idxE] /= 2;
         }
         else
         {
            VertexCIter v0 = he->vertex;
            VertexCIter v1 = he->next->vertex;
            VertexCIter v2 = he->next->next->vertex;
            unsigned int idxV = v1 - meshV0;
            ctlptsV[idxV] = 
               (v0->position + 6 * v1->position + v2->position) / 8;
            degreesV[idxV] *= -1;
         }
      }


      for( unsigned int idxV = 0; idxV < sizeV; idxV++ )
      {
         const int k = degreesV[idxV];
         if( k > 0 )
         {
            ctlptsV[idxV] /= k * k;
            ctlptsV[idxV] += 1.0 * (k - 2) / k * meshV[idxV].position;
         }
      }

      // cout << "V" << endl;

      // for( unsigned int i = 0; i < sizeV; i++ )
      // {
      //    cout << ctlptsV[i] << " " << degreesV[i] << endl;
      // }

      // cout << "E" << endl;

      // for( unsigned int i = 0; i < sizeE; i++ )
      // {
      //    cout << ctlptsE[i] << endl;
      // }

      // cout << "F" << endl;

      // for( unsigned int i = 0; i < sizeF; i++ )
      // {
      //    cout << ctlptsF[i] << " " << degreesF[i] << endl;
      // }

      for( unsigned int i = 0; i < sizeF; i++ )
      {
         HalfEdgeIter he = meshF[i].he, he1 = he;
         do
         {
            data.indices.push_back(vector<Index>());
            data.indices.back().push_back(
               Index(i, -1, -1));
            data.indices.back().push_back(
               Index(sizeF + he->edge - meshE0, -1, -1));
            he = he->next;
            data.indices.back().push_back(
               Index(sizeF + sizeE + he->vertex - meshV0, -1, -1));
            data.indices.back().push_back(
               Index(sizeF + he->edge - meshE0, -1, -1));
         }
         while( he != he1 );
         // for ( unsigned int j = 0; j < data.indices[i].size(); j ++ )
         //    cout << data.indices[i][j].position << " ";
         // cout << endl;
      }

      MeshIO::buildMesh(data, dest);
   }
}
