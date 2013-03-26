//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstSP.cc
 * \author Thomas M. Evans
 * \date   Wed Feb  5 17:29:59 2003
 * \brief  SP test.
 * \note   Copyright (C) 2003-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../SP.hh"
#include "../Assert.hh"

#ifdef DRACO_ENABLE_CXX11
#include <memory> // shared_ptr
#endif
#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>
#include <sstream>

using namespace std;
using rtt_dsxx::SP;

#define PASSMSG(a) ut.passes(a)
#define ITFAILS    ut.failure(__LINE__)
#define FAILURE    ut.failure(__LINE__, __FILE__)
#define FAILMSG(a) ut.failure(a)

//---------------------------------------------------------------------------//
// TEST HELPERS
//---------------------------------------------------------------------------//

int nfoos = 0;
int nbars = 0;
int nbazs = 0;
int nbats = 0;

#define CHECK_0_OBJECTS if (nfoos != 0) ITFAILS; if (nbars != 0) ITFAILS; if (nbazs != 0) ITFAILS; if (nbats != 0) ITFAILS;

#define CHECK_N_OBJECTS(nf, nb, nbz, nbt) if (nfoos != nf) ITFAILS; if (nbars != nb) ITFAILS; if (nbazs != nbz) ITFAILS; if (nbats != nbt) ITFAILS;

//---------------------------------------------------------------------------//

class List
{
  public:
    List() : next() {/*empty*/};
    SP<List> next;
};

// forward declaration
class ListD;

class ListWithDerived
{
  public:
    ListWithDerived(void);
    virtual ~ListWithDerived(void);
    SP<ListD> next;
};

class ListD : public ListWithDerived
{
  public:
    ListD(void);
    ~ListD(void);
};

ListWithDerived::ListWithDerived( void) : next(NULL) {/*empty*/}
ListWithDerived::~ListWithDerived(void)              {/*empty*/}

ListD::ListD( void) : ListWithDerived() {/*empty*/}
ListD::~ListD(void)                     {/*empty*/}
      

class Foo
{
  private:
    int v;

  public:
    Foo(void)           : v(0)   { nfoos++; }
    explicit Foo(int i) : v(i)   { nfoos++; }
    Foo(const Foo &f)   : v(f.v) { nfoos++; }
    virtual ~Foo(void) {  nfoos--; }
    virtual int vf() { return v; }
    int f(void) { return v+1; }
};


//---------------------------------------------------------------------------//

class Bar : public Foo
{
  private:
    Bar(const Bar &);

  public:
    explicit Bar(int i) : Foo(i) { nbars++; }
    virtual ~Bar(void)  { nbars--; }
    virtual int vf()    { return Foo::f() + 1; }
    int f(void)         { return Foo::f() + 2; }
};

//---------------------------------------------------------------------------//

class Baz : public Bar
{
  private:
    Baz(const Baz &);

  public:
    explicit Baz(int i) : Bar(i) { nbazs++; }
    virtual ~Baz() { nbazs--; }
    virtual int vf(void) { return Bar::f() + 1; }
    int f(void) { return Bar::f() + 2; }
};

//---------------------------------------------------------------------------//

class Wombat 
{
  private:
    Wombat(const Wombat &);

  public:
    Wombat() { nbats++; }
    virtual ~Wombat() { nbats--; }
};

//---------------------------------------------------------------------------//
 
SP<Foo> get_foo()
{
    SP<Foo> f(new Foo(10));
    return f;
}

//---------------------------------------------------------------------------//

SP<Bar> get_bar()
{
    SP<Bar> b(new Bar(20));
    return b;
}

//---------------------------------------------------------------------------//

void test_foobar( rtt_dsxx::UnitTest & ut, SP<Foo> f, int v)
{
    if (f->vf() != v) ITFAILS;
    return;
}

//---------------------------------------------------------------------------//

void kill_SPBar(SP<Bar> &b)
{
    b = SP<Bar>();
    return;
}

//---------------------------------------------------------------------------//

void temp_change_SP( rtt_dsxx::UnitTest & ut, SP<Foo> f)
{
    CHECK_N_OBJECTS(1, 1, 0, 0);

    // this is a temporary change
    f = new Foo(100);

    CHECK_N_OBJECTS(2, 1, 0, 0);

    if (f->vf() != 100) ITFAILS;

    if( ut.numFails == 0 )
        PASSMSG("SP<Bar> successfully (temporarily) reassigned to SP<Foo>.");

    return;
}

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
// here we test the following SP members:
//
//    SP();
//    SP(T *);
//    SP(const SP<T> &);
//    SP<T>& operator=(T *);
//    SP<T>& operator=(const SP<T> &);
//    T* operator->() const;
//    bool operator==(const T *) const;
//    bool operator!=(const T *) const;
//    bool operator==(const SP<T> &) const;
//    bool operator!=(const SP<T> &) const;
// 
// plus we test 
//
//    bool operator==(const T *, const SP<T> &);
//    bool operator!=(const T *, const SP<T> &);
//
void type_T_test( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    // test explicit constructor for type T *
    {
        // make a Foo, Bar, and Baz
        SP<Foo> spfoo(new Foo(1));
        SP<Bar> spbar(new Bar(2));
        SP<Baz> spbaz(new Baz(3));
        
        // there should be 3 Foos, 2 Bars and 1 Baz
        CHECK_N_OBJECTS(3, 2, 1, 0);
    }

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Explicit constructor for type T * ok.");

    // test copy constructor for type T *
    {
        SP<Foo> rspfoo;
        SP<Bar> rspbar;
        SP<Baz> rspbaz;
        {
            // no objects yet
            CHECK_0_OBJECTS;

            Foo *f  = new Foo(1);
            Bar *b  = new Bar(2);
            Baz *bz = new Baz(3);
            
            SP<Foo> spfoo(f);
            SP<Bar> spbar(b);
            SP<Baz> spbaz(bz);
            
            // there should be 3 Foos, 2 Bars and 1 Baz
            CHECK_N_OBJECTS(3, 2, 1, 0);

            // now assign
            rspfoo = spfoo;
            rspbar = spbar;
            rspbaz = spbaz;
            
            // there are no additional objects made because the SP will make
            // additional references
            CHECK_N_OBJECTS(3, 2, 1, 0);

            if (ut.numFails == 0)
                PASSMSG("Assignment of SP<T> ok.");

            // now copy construct
            SP<Foo> ispfoo = rspfoo;
            SP<Bar> ispbar = rspbar;

            // still no new foos created
            CHECK_N_OBJECTS(3, 2, 1, 0);
            if (ispfoo->f() != 2) ITFAILS;
            if (spfoo->f()  != 2) ITFAILS;
            if (rspfoo->f() != 2) ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Copy construct of SP<T> ok.");

            // now make a foo pointer and assign
            Foo *ff = new Foo(10);
            ispfoo  = ff;

            // still no new foos created
            CHECK_N_OBJECTS(4, 2, 1, 0);
            if (ispfoo->f() != 11) ITFAILS;
            if (spfoo->f()  != 2)  ITFAILS;
            if (rspfoo->f() != 2)  ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment of T* ok.");

            // now we can check equality 
            if (rspfoo == spfoo)
            {
                PASSMSG("Equality operation ok.");
            }
            else
            {
                FAILMSG("Equality operation failed.");
            }
            
            // now check inequality
            if (rspfoo != spfoo)
            {
                FAILMSG("Equality operation failed.");
            }
            else
            {
                PASSMSG("Equality operation ok.");
            }

            if (rspbar != spbar) ITFAILS;
            if (rspbaz != spbaz) ITFAILS;

            if (spfoo != f)      ITFAILS;
            if (spbar != b)      ITFAILS;
            if (spbaz != bz)     ITFAILS;

            if (spfoo == b)      ITFAILS; // this is ok because a Bar * can
                                          // be passed to Foo *

            if (spbar == dynamic_cast<Bar *>(f)) ITFAILS;

            if (f  != spfoo)     ITFAILS;
            if (b  != spbar)     ITFAILS;
            if (bz != spbaz)     ITFAILS;

            if (f == spfoo)
            {
                PASSMSG("Overloaded equality operators ok.");
            }
            else
            {
                FAILMSG("Overloaded equality operators failed.");
            }

            if (ut.numFails == 0)
                PASSMSG("Equality/Inequality operations ok.");
        }
        
        // we should still have objects left even because we still have
        // viable SPs in scope
        CHECK_N_OBJECTS(3, 2, 1, 0);
    }

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Operations on type T ok");

    return;
}

//---------------------------------------------------------------------------//

#ifdef HAS_CXX11_SHARED_PTR

void type_T_test_shared_ptr( rtt_dsxx::UnitTest & ut )
{
    using std::shared_ptr;
    
    CHECK_0_OBJECTS;

    // test explicit constructor for type T *
    {
        // make a Foo, Bar, and Baz
        shared_ptr<Foo> spfoo(new Foo(1));
        shared_ptr<Bar> spbar(new Bar(2));
        shared_ptr<Baz> spbaz(new Baz(3));
        
        // there should be 3 Foos, 2 Bars and 1 Baz
        CHECK_N_OBJECTS(3, 2, 1, 0);
    }

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Explicit constructor for type T * ok. (std::shared_ptr)");

    // test copy constructor for type T *
    {
        shared_ptr<Foo> rspfoo;
        shared_ptr<Bar> rspbar;
        shared_ptr<Baz> rspbaz;
        {
            // no objects yet
            CHECK_0_OBJECTS;

            Foo *f  = new Foo(1);
            Bar *b  = new Bar(2);
            Baz *bz = new Baz(3);
            
            shared_ptr<Foo> spfoo(f);
            shared_ptr<Bar> spbar(b);
            shared_ptr<Baz> spbaz(bz);
            
            // there should be 3 Foos, 2 Bars and 1 Baz
            CHECK_N_OBJECTS(3, 2, 1, 0);

            // now assign
            rspfoo = spfoo;
            rspbar = spbar;
            rspbaz = spbaz;
            
            // there are no additional objects made because the shared_ptr will make
            // additional references
            CHECK_N_OBJECTS(3, 2, 1, 0);

            if (ut.numFails == 0)
                PASSMSG("Assignment of shared_ptr<T> ok. (std::shared_ptr)");

            // now copy construct
            shared_ptr<Foo> ispfoo = rspfoo;
            shared_ptr<Bar> ispbar = rspbar;

            // still no new foos created
            CHECK_N_OBJECTS(3, 2, 1, 0);
            if (ispfoo->f() != 2) ITFAILS;
            if (spfoo->f()  != 2) ITFAILS;
            if (rspfoo->f() != 2) ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Copy construct of shared_ptr<T> ok. (std::shared_ptr)");

            // now make a foo pointer and assign
            Foo *ff = new Foo(10);
            ispfoo.reset( ff );  // for rtt_dsxx::SP we do 'ispfoo=ff;'

            // still no new foos created
            CHECK_N_OBJECTS(4, 2, 1, 0);
            if (ispfoo->f() != 11) ITFAILS;
            if (spfoo->f()  != 2)  ITFAILS;
            if (rspfoo->f() != 2)  ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment of T* ok. (std::shared_ptr)");

            // now we can check equality 
            if (rspfoo == spfoo)
                PASSMSG("Equality operation ok. (std::shared_ptr)");
            else
                FAILMSG("Equality operation failed. (std::shared_ptr)");
            
            // now check inequality
            if (rspfoo != spfoo)
                FAILMSG("Equality operation failed. (std::shared_ptr)");
            else
                PASSMSG("Equality operation ok. (std::shared_ptr)");

            if (rspbar != spbar)   ITFAILS;
            if (rspbaz != spbaz)   ITFAILS;

            if (spfoo.get() != f)  ITFAILS;
            if (spbar.get() != b)  ITFAILS;
            if (spbaz.get() != bz) ITFAILS;

            if (spfoo.get() == b)  ITFAILS; // this is ok because a Bar * can
                                            // be passed to Foo *

            if (spbar.get() == dynamic_cast<Bar *>(f)) ITFAILS;

            if (f  != spfoo.get()) ITFAILS;
            if (b  != spbar.get()) ITFAILS;
            if (bz != spbaz.get()) ITFAILS;

            if (f == spfoo.get())
                PASSMSG("Overloaded equality operators ok. (std::shared_ptr)");
            else
                FAILMSG("Overloaded equality operators failed. (std::shared_ptr)");

            if (ut.numFails == 0)
                PASSMSG("Equality/Inequality operations ok. (std::shared_ptr)");
        }
    
        // we should still have objects left even because we still have
        // viable shared_ptrs in scope
        CHECK_N_OBJECTS(3, 2, 1, 0);
    }

    // now all should be destroyed
    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Operations on type T ok (std::shared_ptr)");

    return;
}
#endif

//---------------------------------------------------------------------------//
// here we test the following SP members:
//
//    SP();
//    SP(X *);
//    SP<const SP<X> &);
//    SP<T>& operator=(T *);
//    SP<T>& operator=(X *);
//    SP<T>& operator=(const SP<X> &);
//    T* operator->() const;
//    T& operator*() const;
//    T* bp() const;
//    operator bool() const;
//    bool operator!() const;
// 
void type_X_test( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    // check explicit constructor
    {
        // make a foo pointer
        SP<Foo> spfoo(new Bar(10));
        CHECK_N_OBJECTS(1, 1, 0, 0);

        if (spfoo->vf() != 12) ITFAILS;
        if (spfoo->f() != 11)  ITFAILS;

        Foo &f = *spfoo;
        if (f.f() != 11)  ITFAILS;
        if (f.vf() != 12) ITFAILS;

        Foo ff = *spfoo;
        if (ff.vf() != 10) ITFAILS;

        Bar *b = dynamic_cast<Bar *>(spfoo.bp());
        if (b->vf() != 12) ITFAILS;
        if (b->f() != 13)  ITFAILS;

        if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
        if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;

        CHECK_N_OBJECTS(2, 1, 0, 0);
    }

    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Explicit constructor for type X * ok.");

    // check SP<X> constructor and assignment
    {
        // make some objects
        SP<Foo> spfoo;
        SP<Bar> spbar;
        SP<Foo> spfoo2;

        if (spfoo)  ITFAILS;
        if (spbar)  ITFAILS;
        if (spfoo2) ITFAILS;
        {
            spbar = new Bar(50);
            CHECK_N_OBJECTS(1, 1, 0, 0);

            if (spbar->f() != 53)  ITFAILS;
            if (spbar->vf() != 52) ITFAILS;

            // now assign to base class SP
            spfoo = spbar;
            CHECK_N_OBJECTS(1, 1, 0, 0);

            // check reassignment
            spfoo = spbar;
            CHECK_N_OBJECTS(1, 1, 0, 0);

            if (spfoo->f() != 51)  ITFAILS;
            if (spfoo->vf() != 52) ITFAILS;

            if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
            if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;
            if (typeid(spbar.bp()) != typeid(Bar *)) ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment with SP<X> ok.");

            // now do copy construction
            SP<Foo> rspfoo(spbar);
            CHECK_N_OBJECTS(1, 1, 0, 0);

            if (rspfoo->f() != 51)  ITFAILS;
            if (rspfoo->vf() != 52) ITFAILS;

            if (typeid(rspfoo.bp()) != typeid(Foo *)) ITFAILS;
            if (typeid(*rspfoo.bp()) != typeid(Bar))  ITFAILS;
            
            if (ut.numFails == 0)
                PASSMSG("Copy constructor with SP<X> ok.");

            // now check assignment with X *
            rspfoo = new Bar(12);
            CHECK_N_OBJECTS(2, 2, 0, 0);

            if (rspfoo->f() != 13)  ITFAILS;
            if (rspfoo->vf() != 14) ITFAILS;

            if (typeid(rspfoo.bp()) != typeid(Foo *)) ITFAILS;
            if (typeid(*rspfoo.bp()) != typeid(Bar))  ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment with X * ok.");

            // assign SPfoo2 to a bar
            spfoo2 = new Bar(20);
            CHECK_N_OBJECTS(3, 3, 0, 0);

            // assign SPfoo2 to itself
            spfoo2 = spfoo2;
            CHECK_N_OBJECTS(3, 3, 0, 0);
            
            // assign SPfoo2 to itself (underlying object)
            spfoo2 = spfoo2.bp();
            CHECK_N_OBJECTS(3, 3, 0, 0);
            
        }
        // still have 2 object
        CHECK_N_OBJECTS(2, 2, 0, 0);

        // assign spfoo to a baz
        spfoo2 = new Baz(45);
        CHECK_N_OBJECTS(2, 2, 1, 0);

        if (spfoo2->f() != 46)  ITFAILS;
        if (spfoo2->vf() != 49) ITFAILS;

        if (typeid(*spfoo2.bp()) != typeid(Baz)) ITFAILS;

        // assign spbar to NULL
        spbar = SP<Bar>();
        CHECK_N_OBJECTS(2, 2, 1, 0);

        // spfoo should still point to the same bar
        if (spfoo->f() != 51)  ITFAILS;
        if (spfoo->vf() != 52) ITFAILS;
        
        if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;
        if (typeid(*spfoo.bp()) != typeid(Bar))  ITFAILS;
        if (typeid(spbar.bp()) != typeid(Bar *)) ITFAILS;
        
        if (ut.numFails == 0)
            PASSMSG("Set to SP<>() releases pointer.");

        // assign spfoo to NULL
        spfoo = SP<Foo>();
        CHECK_N_OBJECTS(1, 1, 1, 0);

        if (spfoo) ITFAILS;
        if (spbar) ITFAILS;

        if (ut.numFails == 0)
            PASSMSG("Overloaded bool ok.");

        if (!spfoo2) ITFAILS;

        if (ut.numFails == 0)
            PASSMSG("Overloaded ! (not) ok."); 
    }

    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Operations on type X ok");

    return;
}

//---------------------------------------------------------------------------//
#ifdef HAS_CXX11_SHARED_PTR
void type_X_test_shared_ptr( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    // check explicit constructor
    {
        // make a foo pointer
        shared_ptr<Foo> spfoo(new Bar(10));
        CHECK_N_OBJECTS(1, 1, 0, 0);

        if (spfoo->vf() != 12) ITFAILS;
        if (spfoo->f()  != 11) ITFAILS;

        Foo &f = *spfoo;
        if (f.f()  != 11) ITFAILS;
        if (f.vf() != 12) ITFAILS;

        Foo ff = *spfoo;
        if (ff.vf() != 10) ITFAILS;

        Bar *b = dynamic_cast<Bar *>(spfoo.get()); // base_class ptr?
        if (b->vf() != 12) ITFAILS;
        if (b->f()  != 13) ITFAILS;

        if (typeid(spfoo.get()) != typeid(Foo *)) ITFAILS;
        if (typeid(*spfoo.get()) != typeid(Bar))  ITFAILS;

        CHECK_N_OBJECTS(2, 1, 0, 0);
    }

    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Explicit constructor for type X * ok (std::shared_ptr<X>).");

    // check shared_ptr<X> constructor and assignment
    {
        // make some objects
        shared_ptr<Foo> spfoo;
        shared_ptr<Bar> spbar;
        shared_ptr<Foo> spfoo2;

        if (spfoo)  ITFAILS;
        if (spbar)  ITFAILS;
        if (spfoo2) ITFAILS;
        {
            spbar.reset( new Bar(50) ); // spbar = new Bar(50);
            CHECK_N_OBJECTS(1, 1, 0, 0);
            
            if (spbar->f()  != 53) ITFAILS;
            if (spbar->vf() != 52) ITFAILS;

            // now assign to base class shared_ptr
            spfoo = spbar;
            CHECK_N_OBJECTS(1, 1, 0, 0);

            // check reassignment
            spfoo = spbar;
            CHECK_N_OBJECTS(1, 1, 0, 0);

            if (spfoo->f()  != 51) ITFAILS;
            if (spfoo->vf() != 52) ITFAILS;

            if (typeid(spfoo.get())  != typeid(Foo *)) ITFAILS;
            if (typeid(*spfoo.get()) != typeid(Bar))   ITFAILS;
            if (typeid(spbar.get())  != typeid(Bar *)) ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment with shared_ptr<X> ok.");

            // now do copy construction
            shared_ptr<Foo> rspfoo(spbar);
            CHECK_N_OBJECTS(1, 1, 0, 0);

            if (rspfoo->f()  != 51) ITFAILS;
            if (rspfoo->vf() != 52) ITFAILS;

            if (typeid(rspfoo.get())  != typeid(Foo *)) ITFAILS;
            if (typeid(*rspfoo.get()) != typeid(Bar))   ITFAILS;
            
            if (ut.numFails == 0)
                PASSMSG("Copy constructor with shared_ptr<X> ok.");

            // now check assignment with X *
            rspfoo.reset( new Bar(12) ); // rspfoo = new Bar(12);
            CHECK_N_OBJECTS(2, 2, 0, 0);

            if (rspfoo->f()  != 13) ITFAILS;
            if (rspfoo->vf() != 14) ITFAILS;

            if (typeid(rspfoo.get())  != typeid(Foo *)) ITFAILS;
            if (typeid(*rspfoo.get()) != typeid(Bar))   ITFAILS;

            if (ut.numFails == 0)
                PASSMSG("Assignment with X * ok (std::shared_ptr<X>).");

            // assign shared_ptrfoo2 to a bar
            spfoo2.reset( new Bar(20) ); // spfoo2 = new Bar(20);
            CHECK_N_OBJECTS(3, 3, 0, 0);

            // assign spfoo2 to itself
            spfoo2 = spfoo2;
            CHECK_N_OBJECTS(3, 3, 0, 0);
            
            // assign spfoo2 to itself (underlying object)
            // !!! THIS BEHAVIOR IS DIFFERENT FROM rtt_dsxx:SP<X> !!! //
            // spfoo2.reset( spfoo2.get() ); // spfoo2 = spfoo2.bp();
            // CHECK_N_OBJECTS(2, 2, 0, 0);  // CHECK_N_OBJECTS(3, 3, 0, 0);
        }
        // still have 2 object
        CHECK_N_OBJECTS(2, 2, 0, 0);

        // assign spfoo to a baz
        spfoo2.reset( new Baz(45) ); // spfoo2 = new Baz(45);
        CHECK_N_OBJECTS(2, 2, 1, 0);

        if (spfoo2->f()  != 46) ITFAILS;
        if (spfoo2->vf() != 49) ITFAILS;

        if (typeid(*spfoo2.get()) != typeid(Baz)) ITFAILS;

        // assign spbar to NULL
        spbar = shared_ptr<Bar>();
        CHECK_N_OBJECTS(2, 2, 1, 0);

        // spfoo should still point to the same bar
        if (spfoo->f()  != 51) ITFAILS;
        if (spfoo->vf() != 52) ITFAILS;
        
        if (typeid(spfoo.get())  != typeid(Foo *)) ITFAILS;
        if (typeid(*spfoo.get()) != typeid(Bar))   ITFAILS;
        if (typeid(spbar.get() ) != typeid(Bar *)) ITFAILS;
        
        if (ut.numFails == 0)
            PASSMSG("Set to shared_ptr<>() releases pointer.");

        // assign spfoo to NULL
        spfoo = shared_ptr<Foo>();
        CHECK_N_OBJECTS(1, 1, 1, 0);

        if (spfoo) ITFAILS;
        if (spbar) ITFAILS;

        if (ut.numFails == 0)
            PASSMSG("Overloaded bool ok.");

        if (!spfoo2) ITFAILS;

        if (ut.numFails == 0)
            PASSMSG("Overloaded ! (not) ok (std::shared_ptr<X>)."); 
    }
    
    CHECK_0_OBJECTS;

    if (ut.numFails == 0)
        PASSMSG("Operations on type X ok (std::shared_ptr<X>).");

    return;
}
#endif

//---------------------------------------------------------------------------//

void fail_modes_test( rtt_dsxx::UnitTest & ut )
{
    // make an object and try to reference it
    SP<Foo>    spfoo;
    SP<Bar>    spbar;
    SP<Baz>    spbaz;
    SP<Wombat> spbat;

    if (spfoo) ITFAILS;
    if (spfoo) ITFAILS;
    if (spfoo) ITFAILS;

    // try to reference a function
    bool caught = false;
    if( ut.dbcOn() && ! ut.dbcNothrow() )
    {
        try
        {
            spfoo->f();
        }
        catch (rtt_dsxx::assertion & /* assertion */ )
        {
            caught = true;
            ostringstream m;
            m << "Good, caught null access on member function " << endl;
            PASSMSG(m.str());
        }
        if (!caught)
            FAILMSG("Failed to catch illegal access exception.");

        // try assigning a derived NULL to a base; the spfoo base pointer is
        // still a foo in the case 
        spfoo = spbar;
        if (typeid(spfoo.bp()) != typeid(Foo *)) ITFAILS;

        CHECK_0_OBJECTS;

        // now try assigning to a non-derived class of Foo that is NULL,
        // unfortunately, even though this shouldn't be allowed we get away with
        // it because wombat has some virtual functions and is NULL; however,
        // this isn't really dangerous because spfoo still doesn't point to
        // anything 
        spfoo  = spbat;
        caught = false;
        try
        {
            spfoo->f();
        }
        catch (rtt_dsxx::assertion & /* assertion */ )
        {
            caught = true;
            ostringstream m;
            m << "Good, caught null access on member function " << endl;
            PASSMSG(m.str());
        }
        if (!caught)
            FAILMSG("Failed to catch illegal access exception.");
    }
    // now make a wombat and try
    spbat  = new Wombat;
    CHECK_N_OBJECTS(0, 0, 0, 1);

    caught = false;
    try
    {
        spfoo = spbat;
    }
    catch (rtt_dsxx::assertion &assert)
    {
        caught = true;
        ostringstream m;
        m << "Good, caught the following exception, " << endl
          << assert.what();
        PASSMSG(m.str());
    }
    if (!caught)
        FAILMSG("Failed to catch illegal assignment.");

    // now try copy construction
    caught = false;
    try
    {
        SP<Foo> s(spbat);
    }
    catch (rtt_dsxx::assertion &assert)
    {
        caught = true;
        ostringstream m;
        m << "Good, caught the following exception, " << endl
          << assert.what();
        PASSMSG(m.str());
    }
    if (!caught)
        FAILMSG("Failed to catch illegal assignment.");

    // now try copy and assignment on X *
    Wombat *bat = new Wombat();
    CHECK_N_OBJECTS(0, 0, 0, 2);

    caught = false;
    try
    {
        spfoo = bat;
    }
    catch (rtt_dsxx::assertion &assert)
    {
        caught = true;
        ostringstream m;
        m << "Good, caught the following exception, " << endl
          << assert.what();
        PASSMSG(m.str());
    }
    if (!caught)
        FAILMSG("Failed to catch illegal assignment.");

    // now try copy construction
    caught = false;
    try
    {
        SP<Foo> s(bat);
    }
    catch (rtt_dsxx::assertion &assert)
    {
        caught = true;
        ostringstream m;
        m << "Good, caught the following exception, " << endl
          << assert.what();
        PASSMSG(m.str());
    }
    if (!caught)
        FAILMSG("Failed to catch illegal assignment.");
    
    // assign wombat to a pointer to clean it up
    spbat = bat;
    CHECK_N_OBJECTS(0, 0, 0, 1);
    
    if (ut.numFails == 0)
        PASSMSG("Failure modes work ok.");

    return;
}

//---------------------------------------------------------------------------//

void equality_test( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    // try some more advanced stuff
    SP<Foo> f1;
    SP<Foo> f2;

    Foo *f  = new Foo(5);
    Foo *ff = new Foo(5);

    f1 = f;
    f2 = f1;

    if (f2 != f1) ITFAILS;

    // now f and ff are equivalent, but the smart pointers won't be because
    // they don't point to the same instance of Foo *
    f2 = ff;

    if (f2 == f1) ITFAILS;

    CHECK_N_OBJECTS(2, 0, 0, 0);

    if (ut.numFails == 0)
        PASSMSG("Equality tests work ok.");

    return;
}

//---------------------------------------------------------------------------//

void get_test( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    // get a foo and bar
    {
        
        SP<Foo> f  = get_foo();
        SP<Foo> fb = get_bar();
        SP<Bar> b  = get_bar();
        
        CHECK_N_OBJECTS(3, 2, 0, 0);
        
        if (fb == b) ITFAILS;
        
        if (f->f() != 11)   ITFAILS;
        if (fb->vf() != 22) ITFAILS;
        if (b->vf() != 22)  ITFAILS;
        
        if (fb->f() != 21)  ITFAILS;
        if (b->f() != 23)   ITFAILS;
    }

    if (ut.numFails == 0)
        PASSMSG("Get/factory tests work ok.");

    return;
}

//---------------------------------------------------------------------------//

void access_test( rtt_dsxx::UnitTest & ut )
{
    CHECK_0_OBJECTS;

    SP<Bar> b(new Bar(10));
    CHECK_N_OBJECTS(1, 1, 0, 0);

    test_foobar(ut, b, 12);
    CHECK_N_OBJECTS(1, 1, 0, 0);

    kill_SPBar(b);
    CHECK_0_OBJECTS;

    b = new Bar(12);
    temp_change_SP(ut, b); // this temporarily changes to a Foo

    if (b->vf() != 14)                   ITFAILS;
    if (typeid(b.bp()) != typeid(Bar *)) ITFAILS;

    CHECK_N_OBJECTS(1, 1, 0, 0);

    if (ut.numFails == 0)
        PASSMSG("Accessor/set-style tests work ok.");

    return;
}

//---------------------------------------------------------------------------//

void list_test( rtt_dsxx::UnitTest & ut )
{
    {
        // This test was borrowed from Boost's shared_ptr_test.cpp
        
        SP<List> p(new List);
        p->next = SP<List>(new List);
        p = p->next;
        if ( p->next ) ITFAILS;
    }

    {
        // Test of a derived class.
        
        SP<ListWithDerived> p(new ListWithDerived);
        p->next = SP<ListD>(new ListD);
        p = p->next;
        if ( p->next ) ITFAILS;
    }

    if (ut.numFails == 0)
        PASSMSG("Linked-list test works ok.");

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
    try
    {
        // >>> UNIT TESTS

        CHECK_0_OBJECTS;

        type_T_test(ut);     cout << endl;
        type_X_test(ut);     cout << endl;
        fail_modes_test(ut); cout << endl;
        equality_test(ut);   cout << endl;
        get_test(ut);        cout << endl;
        access_test(ut);     cout << endl;
        list_test(ut);       cout << endl;
        CHECK_0_OBJECTS;

#ifdef HAS_CXX11_SHARED_PTR
        type_T_test_shared_ptr(ut);
        type_X_test_shared_ptr(ut);
        CHECK_0_OBJECTS;        
#endif
        
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstSP.cc
//---------------------------------------------------------------------------//
