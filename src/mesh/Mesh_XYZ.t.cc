//----------------------------------*-C++-*----------------------------------//
// Mesh_XYZ.t.cc
// Geoffrey M. Furnish
// Wed May 27 11:02:08 1998
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#include "c4/global.hh"
#include "c4/C4_Req.hh"

template<class T>
void dump( const Mesh_XYZ::cctf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::cctf: " << name << endl;
    {
    //	HTSyncSpinLock h;
	char buf[80];
	for( int i=0; i < data.size(); i++ ) {
        //	    sprintf( buf, "node %d, cell %d, value=%lf \n",
        //		     C4::node(), i, data(i) );
	    sprintf( buf, "cell %d, value=%lf \n",
		     i, data(i) );
	    cout << buf;
	}
    }
}

template<class T>
void dump( const Mesh_XYZ::fcdtf<T>& data, char *name )
{
    cout << "dumping a Mesh_XYZ::fcdtf: " << name << endl;
    {
    //	HTSyncSpinLock h;
	char buf[80];
	for( int i=0; i < data.size()/6; i++ ) {
            for( int j=0; j < 6; j++ ) {
	        sprintf( buf, "cell %d, face %d, value=%lf \n",
		         i, j, data(i,j) );
                cout << buf;
            }
	}
    }

}

template<class T>
Mesh_XYZ::gcctf<T>& 
Mesh_XYZ::gcctf<T>::operator=( const Mesh_XYZ::cctf<T>& c )
{
    for( int k=zoff; k < zoff + nczp; k++ )
        for( int j=0; j < ncy; j++ )
            for( int i=0; i < ncx; i++ )
                data(i,j,k) = c(i,j,k);

    update_guard_cells();
    
    return *this;
}

template<class T>
void Mesh_XYZ::gcctf<T>::update_guard_cells()
{
    using namespace C4;
    C4_Req lrcv, rrcv;

//     Mat2<T> lrbuf( &data(

//     if (node > 0)
//         AsyncRecv( 
}

template <class Op>
void Mesh_XYZ::scatter( Mesh_XYZ::fcdsf& to, const Mesh_XYZ::ccsf& from )
{
    Op::thisIsAnError();
}


//---------------------------------------------------------------------------//
//                              end of Mesh_XYZ.t.cc
//---------------------------------------------------------------------------//
