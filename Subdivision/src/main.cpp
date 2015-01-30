#include <iostream>
#include <unistd.h>
using namespace std;

#include "Viewer.h"
#include "DenseMatrix.h"
using namespace DDG;

int main( int argc, char** argv )
{
   if( argc != 2 )
   {
      cerr << "usage: " << argv[0] << " in.obj" << endl;
      return 1;
   }

   Viewer viewer;
   viewer.mesh.read( argv[1] );

   string path(argv[0]);

   size_t found = path.find_last_of('/');

   if (found != string::npos) {
      path = path.substr(0, found);
      chdir(path.c_str());
   }

   viewer.init();

   return 0;
}

