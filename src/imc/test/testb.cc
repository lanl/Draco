//----------------------------------*-C++-*----------------------------------//
// testb.cc
// Thomas M. Evans
// Wed Feb 18 10:45:23 1998
//---------------------------------------------------------------------------//
// @> test driver for OS_Builder
//---------------------------------------------------------------------------//

#include "OS_Builder.hh"
#include "OS_Mesh.hh"
#include "SP.hh"
#include <iostream>
#include <string>

main()
{
    using IMC::OS_Builder;
    using IMC::OS_Mesh;
    using namespace std;

    SP<OS_Mesh> mesh;

  // scoping blocks
    {
	string infile;
	cout << "Name the input file" << endl;
	cin >> infile;
	
      // initialize the builder
	OS_Builder build(infile);
	mesh = build.Build_Mesh();
    }

    cout << "Coordinate System: " << mesh->Coord().Get_coord() << endl;
    int cell;
    cout << "Give a cell" << endl;
    cin >> cell;
    mesh->Print(cell);
    cout << endl;
    cout << "Mesh Size: " << mesh->Num_cells() << endl;
    cout << "Position of Cells:" << endl;
    for (int i = 1; i <= mesh->Num_cells(); i++)
	cout << i << " " << "(" << mesh->Pos(1,i) << "," << mesh->Pos(2,i)
	     << ")" << endl;
	    
}

//---------------------------------------------------------------------------//
//                              end of testb.cc
//---------------------------------------------------------------------------//
