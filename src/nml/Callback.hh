//----------------------------------*-C++-*----------------------------------//
// Callback.hh
// Geoffrey Furnish
// 15 January 1994
//---------------------------------------------------------------------------//
// @> A class for managing callbacks.
//---------------------------------------------------------------------------//

#ifndef __nml_Callback_hh__
#define __nml_Callback_hh__

#include "ds++/String.hh"

class NML_Callback {

    dsxx::String name;
    void *data;
    void (*f) ( const dsxx::String& arg, void *data );

  public:
    NML_Callback() : data(0), f(0) {}

    NML_Callback( void *_d, void (*_f) ( const dsxx::String&, void * ) )
	: data(_d), f(_f)
    {}

    void invoke( const dsxx::String& current_val )
    {
	if (!f) throw "Cannot invoke NULL Callback.";

	(*f) ( current_val, data );
    }

    bool operator==( const NML_Callback& cb ) const
    {
	return (data == cb.data) && (f == cb.f);
    }
};

#endif				// __nml_Callback_hh__

//---------------------------------------------------------------------------//
//                              end of nml/Callback.hh
//---------------------------------------------------------------------------//
