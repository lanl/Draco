//----------------------------------*-C++-*----------------------------------//
// LSquad.cc
// William D. Hawkins
// Mon Jun 14 10:07:47 1999
// $Id$
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

// Level-Symmetric 3-D Angular Quadrature Class for Isotropic Scattering
// Returns quadrature points and weights for S2 - S24
// LSquad.hh: Declaration of the LSquad class
// LSquad.cc: Member function definitions

#include <assert.h>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <iterator>

#include "LSquad.hh"

using std::cout;
using std::endl;
using std::setprecision;
using std::iterator;
using std::ios;

namespace rtt_quad
{

// Constructor Function
// Calls setLevels, setNorm, and computeQuadSet member functions
// Default number of levels is 2, default normalization is to 4.0*Pi
LSquad::LSquad(int inp_levels, double inp_norm)
{
  computeQuadSet(inp_levels,inp_norm);
}

LSquad::LSquad(int levels_, int dirs_, int odirs_, double norm_,
	       const vector<double> &mu_, const vector<double> &eta_,
	       const vector<double> &xi_, const vector<double> &wt_)
    : levels(levels_), dirs(dirs_), odirs(odirs_), norm(norm_), mu(mu_),
      eta(eta_), xi(xi_), wt(wt_)
{}

// Destructor Function
LSquad::~LSquad()
{
  mu.clear();
  eta.clear();
  xi.clear();
  wt.clear();
}

// Checks and sets the number of quadrature levels
void LSquad::setLevels(int inp_levels)
{
  if(inp_levels >= 2 && inp_levels <= 24 && inp_levels % 2 == 0)
    levels = inp_levels;
  else
  {
    cout << inp_levels << " is an invalid level value.\n"
         << "Valid values are: 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24.\n";
    exit(1);
  }
}

// Checks and sets the normalization constant
void LSquad::setNorm(double inp_norm)
{
  if(inp_norm > 0.0)
    norm = inp_norm;
  else
  {
    cout << inp_norm << " is an invalid normalization value.\n"
         << "A normalization value must be nonzero and positive.\n";
    exit(1);
  }
}

// Returns the value of mu for a specific quadrature direction
double LSquad::getMu(int dir) const
{
  if(dir < 0 || dir > dirs-1)
  {
    cout << dir << " is an invalid direction.\n"
         << "Valid direction values are 0 - " << dirs-1 << ".";
    exit(1);
  }
  return(mu[dir]);
}

// Returns the value of eta for a specific quadrature direction
double LSquad::getEta(int dir) const
{
  if(dir < 0 || dir > dirs-1)
  {
    cout << dir << " is an invalid direction.\n"
         << "Valid direction values are 0 - " << dirs-1 << ".";
    exit(1);
  }
  return(eta[dir]);
}

// Returns the value of xi for a specific quadrature direction
double LSquad::getXi(int dir) const
{
  if(dir < 0 || dir > dirs-1)
  {
    cout << dir << " is an invalid direction.\n"
         << "Valid direction values are 0 - " << dirs-1 << ".";
    exit(1);
  }
  return(xi[dir]);
}

// Returns the weight for a specific quadrature direction
double LSquad::getWt(int dir) const
{
 if(dir < 0 || dir > dirs-1)
 {
    cout << dir << " is an invalid direction.\n"
         << "Valid direction values are 0 - " << dirs-1 << ".";
    exit(1);
 }
 return(wt[dir]);
}

// Returns the omega vector for a specific quadrature direction
vector<double> LSquad::getOmega(int dir) const
{
 vector<double> omega;

 omega.resize(3);
 if(dir < 0 || dir > dirs-1)
 {
    cout << dir << " is an invalid direction.\n"
         << "Valid direction values are 0 - " << dirs-1 << ".";
    exit(1);
 }
 omega[0] = mu[dir];
 omega[1] = eta[dir];
 omega[2] = xi[dir];
 return(omega);
}

// Writes a quadrature table to standard out
void LSquad::printQuad() const
{
  cout << "\n";
  cout << "LQn S-" << levels << " Quadrature Set\n";
  cout << "-----------------------\n";
  cout << "\n";
  cout << "Number of quadrature levels: " << levels << "\n";
  cout << "Number of quadrature directions: " << dirs << "\n";
  cout << "Number of quadrature directions per octant: " << odirs << "\n";
  cout << "Value of normalization constant: " << setprecision(10) << norm
       << "\n\n";
  cout << "      Mu               Eta                Xi              Weight\n"
       << " ------------      ------------      ------------      ------------\n";
  for(int idx = 0;idx < dirs;++idx)
  {
    cout.setf(ios::fixed,ios::floatfield);
    cout.width(13);
    cout << mu[idx] << "     ";
    cout.width(13);
    cout << eta[idx] << "     ";
    cout.width(13);
    cout << xi[idx] << "     ";
    cout.width(13);
    cout << wt[idx] << "     " << "\n";
  }
  cout << "\n";
}

// Computes the quadrature set (mu, eta, xi, and weights)
void LSquad::computeQuadSet(int inp_levels, double inp_norm)
{
  int i,j,m,n;

  int wp4[3]   = {0,0,0};
  int wp6[6]   = {0,1,0,1,1,0};
  int wp8[10]  = {0,1,1,0,1,2,1,1,1,0};
  int wp10[15] = {0,1,2,1,0,1,3,3,1,2,3,2,1,1,0};
  int wp12[21] = {0,1,2,2,1,0,1,3,4,3,1,2,4,4,2,2,3,2,1,1,0};
  int wp14[28] = {0,1,2,3,2,1,0,1,4,5,5,4,1,2,5,6,5,2,3,5,5,3,2,4,2,1,1,0};
  int wp16[36] = {0,1,2,3,3,2,1,0,1,4,5,6,5,4,1,2,5,7,7,5,2,3,6,7,6,3,3,5,5,3,2,
                  4,2,1,1,0};
  int wp18[45] = {0,1,2,3,4,3,2,1,0,1,5,6,7,7,6,5,1,2,6,8,9,8,6,2,3,7,9,9,7,3,4,
                  7,8,7,4,3,6,6,3,2,5,2,1,1,0};
  int wp20[55] = {0,1,2,3,4,4,3,2,1,0,1,5,6,7,8,7,6,5,1,2,6,9,10,10,9,6,2,3,7,
                  10,11,10,7,3,4,8,10,10,8,4,4,7,9,7,4,3,6,6,3,2,5,2,1,1,0};
  int wp22[66] = {0,1,2,3,4,5,4,3,2,1,0,1,6,7,8,9,9,8,7,6,1,2,7,10,11,12,11,10,
                  7,2,3,8,11,13,13,11,8,3,4,9,12,13,12,9,4,5,9,11,11,9,5,4,8,10,
		  8,4,3,7,7,3,2,6,2,1,1,0};
  int wp24[78] = {0,1,2,3,4,5,5,4,3,2,1,0,1,6,7,8,9,10,9,8,7,6,1,2,7,11,12,13,
                  13,12,11,7,2,3,8,12,14,15,14,12,8,3,4,9,13,15,15,13,9,4,5,10,
		  13,14,13,10,5,5,9,12,12,9,5,4,8,11,8,4,3,7,7,3,2,6,2,1,1,0};

  setLevels(inp_levels);
  setNorm(inp_norm);

  dirs = (levels*(levels+2));
  odirs = dirs/8;

  mu.resize(dirs);
  eta.resize(dirs);
  xi.resize(dirs);
  wt.resize(dirs);

  vector<double> att;
  vector<double> wtt;

  switch (levels)
  {
    case 2: // LQn S-2 Quadrature Set
      att.resize(1);
      att[0] = 0.577350269189625764509149;
      wt[0]  = 1.000000000000000000000000;
    break;

    case 4: // LQn S-4 Quadrature Set
      att.resize(2);
      wtt.resize(1);
      att[0] = 0.350021174581540677777041;
      att[1] = 0.868890300722201205229788;
      wtt[0] = 0.333333333333333333333333;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp4[m]];
    break;

    case 6: // LQn S-6 Quadrature Set
      att.resize(3);
      wtt.resize(2);
      att[0] = 0.266635401516704720331535;
      att[1] = 0.681507726536546927403750;
      att[2] = 0.926180935517489107558380;
      wtt[0] = 0.176126130863383433783565;
      wtt[1] = 0.157207202469949899549768;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp6[m]];
    break;

    case 8: // LQn S-8 Quadrature Set
      att.resize(4);
      wtt.resize(3);
      att[0] = 0.218217890235992381266097;
      att[1] = 0.577350269189625764509149;
      att[2] = 0.786795792469443145800830;
      att[3] = 0.951189731211341853132399;
      wtt[0] = 0.120987654320987654320988;
      wtt[1] = 0.0907407407407407407407407;
      wtt[2] = 0.0925925925925925925925926;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp8[m]];
    break;

    case 10:  // LQn S-10 Quadrature Set
      att.resize(5);
      wtt.resize(4);
      att[0] = 0.189321326478010476671494;
      att[1] = 0.508881755582618974382711;
      att[2] = 0.694318887594384317279217;
      att[3] = 0.839759962236684758403029;
      att[4] = 0.963490981110468484701598;
      wtt[0] = 0.0893031479843567214704325;
      wtt[1] = 0.0725291517123655242296233;
      wtt[2] = 0.0450437674364086390490892;
      wtt[3] = 0.0539281144878369243545650;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp10[m]];
    break;

    case 12: // LQn S-12 Quadrature Set */
      att.resize(6);
      wtt.resize(5);
      att[0] = 0.167212652822713264084504;
      att[1] = 0.459547634642594690016761;
      att[2] = 0.628019096642130901034766;
      att[3] = 0.760021014833664062877138;
      att[4] = 0.872270543025721502340662;
      att[5] = 0.971637719251358378302376;
      wtt[0] = 0.0707625899700910439766549;
      wtt[1] = 0.0558811015648888075828962;
      wtt[2] = 0.0373376737588285824652402;
      wtt[3] = 0.0502819010600571181385765;
      wtt[4] = 0.0258512916557503911218290;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp12[m]];
    break;

    case 14: // LQn S-14 Quadrature Set
      att.resize(7);
      wtt.resize(7);
      att[0] = 0.151985861461031912404799;
      att[1] = 0.422156982304796966896263;
      att[2] = 0.577350269189625764509149;
      att[3] = 0.698892086775901338963210;
      att[4] = 0.802226255231412057244328;
      att[5] = 0.893691098874356784901111;
      att[6] = 0.976627152925770351762946;
      wtt[0] = 0.0579970408969969964063611;
      wtt[1] = 0.0489007976368104874582568;
      wtt[2] = 0.0227935342411872473257345;
      wtt[3] = 0.0394132005950078294492985;
      wtt[4] = 0.0380990861440121712365891;
      wtt[5] = 0.0258394076418900119611012;
      wtt[6] = 0.00826957997262252825269908;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp14[m]];
    break;

    case 16: // LQn S-16 Quadrature Set
      att.resize(8);
      wtt.resize(8);
      att[0] = 0.138956875067780344591732;
      att[1] = 0.392289261444811712294197;
      att[2] = 0.537096561300879079878296;
      att[3] = 0.650426450628771770509703;
      att[4] = 0.746750573614681064580018;
      att[5] = 0.831996556910044145168291;
      att[6] = 0.909285500943725291652116;
      att[7] = 0.980500879011739882135849;
      wtt[0] = 0.0489872391580385335008367;
      wtt[1] = 0.0413295978698440232405505;
      wtt[2] = 0.0203032007393652080748070;
      wtt[3] = 0.0265500757813498446015484;
      wtt[4] = 0.0379074407956004002099321;
      wtt[5] = 0.0135295047786756344371600;
      wtt[6] = 0.0326369372026850701318409;
      wtt[7] = 0.0103769578385399087825920;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp16[m]];
    break;

    case 18: // LQn S-18 Quadrature Set
      att.resize(9);
      wtt.resize(10);
      att[0]  = 0.129344504545924818514086;
      att[1]  = 0.368043816053393605686086;
      att[2]  = 0.504165151725164054411848;
      att[3]  = 0.610662549934881101060239;
      att[4]  = 0.701166884252161909657019;
      att[5]  = 0.781256199495913171286914;
      att[6]  = 0.853866206691488372341858;
      att[7]  = 0.920768021061018932899055;
      att[8]  = 0.983127661236087115272518;
      wtt[0]  = 0.0422646448843821748535825;
      wtt[1]  = 0.0376127473827281471532380;
      wtt[2]  = 0.0122691351637405931037187;
      wtt[3]  = 0.0324188352558815048715646;
      wtt[4]  = 0.00664438614619073823264082;
      wtt[5]  = 0.0312093838436551370068864;
      wtt[6]  = 0.0160127252691940275641645;
      wtt[7]  = 0.0200484595308572875885066;
      wtt[8]  = 0.000111409402059638628382279;
      wtt[9]  = 0.0163797038522425240494567;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp18[m]];
    break;

    case 20: // LQn S-20 Quadrature Set
      att.resize(10);
      wtt.resize(12);
      att[0]  =  0.120603343036693597409418;
      att[1]  =  0.347574292315847257336779;
      att[2]  =  0.476519266143665680817278;
      att[3]  =  0.577350269189625764509149;
      att[4]  =  0.663020403653288019308783;
      att[5]  =  0.738822561910371432904974;
      att[6]  =  0.807540401661143067193530;
      att[7]  =  0.870852583760463975580977;
      att[8]  =  0.929863938955324566667817;
      att[9]  =  0.985347485558646574628509;
      wtt[0]  =  0.0370210490604481342320295;
      wtt[1]  =  0.0332842165376314841003910;
      wtt[2]  =  0.0111738965965092519614021;
      wtt[3]  =  0.0245177476959359285418987;
      wtt[4]  =  0.0135924329650041789567081;
      wtt[5]  =  0.0318029065936585971501960;
      wtt[6]  =  0.00685492401402507781062634;
      wtt[7]  =  0.0308105481755299327227893;
      wtt[8]  = -0.000139484716502602877593527;
      wtt[9]  =  0.00544675187330776223879437;
      wtt[10] =  0.00474564692642379971238396;
      wtt[11] =  0.0277298541009064049325246;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp20[m]];
    break;

    case 22: // LQn S-22 Quadrature Set
      att.resize(11);
      wtt.resize(14);
      att[0]  =  0.113888641383070838173488;
      att[1]  =  0.330271760593086736334651;
      att[2]  =  0.452977095507524183904005;
      att[3]  =  0.548905330875560154226714;
      att[4]  =  0.630401360620980621392149;
      att[5]  =  0.702506006153654989703184;
      att[6]  =  0.767869456282208576047898;
      att[7]  =  0.828089557415325768804621;
      att[8]  =  0.884217805921983001958912;
      att[9]  =  0.936989829997455780115072;
      att[10] =  0.986944149751056870330152;
      wtt[0]  =  0.0329277718552552308051381;
      wtt[1]  =  0.0309569328165031538543025;
      wtt[2]  =  0.00577105953220643022391829;
      wtt[3]  =  0.0316834548379952775919418;
      wtt[4]  = -0.00669350304140992494103696;
      wtt[5]  =  0.0368381622687682466526634;
      wtt[6]  =  0.0273139698006629537455404;
      wtt[7]  =  0.0100962716435030437817055;
      wtt[8]  =  0.0195181067555849392224199;
      wtt[9]  =  0.0117224275470949786864925;
      wtt[10] = -0.00442773155233893239996431;
      wtt[11] =  0.0156214785078803432781324;
      wtt[12] = -0.0101774221315738297143270;
      wtt[13] =  0.0135061258938431808485310;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp22[m]];
    break;

    case 24: // LQn S-24 Quadrature Set
      att.resize(12);
      wtt.resize(16);
      att[0]  =  0.107544208775147285552086;
      att[1]  =  0.315151630853896874875332;
      att[2]  =  0.432522073446742487657060;
      att[3]  =  0.524242441631224399254880;
      att[4]  =  0.602150256328323868809286;
      att[5]  =  0.671073561381361944701265;
      att[6]  =  0.733549261041044861004094;
      att[7]  =  0.791106384731321324814121;
      att[8]  =  0.844750913317919895113069;
      att[9]  =  0.895186516397704814461305;
      att[10] =  0.942928254285052510917188;
      att[11] =  0.988366574868785749937406;
      wtt[0]  =  0.0295284942030736546025272;
      wtt[1]  =  0.0281530651743695026834932;
      wtt[2]  =  0.00519730128072174996473824;
      wtt[3]  =  0.0259897467786242920448933;
      wtt[4]  =  0.00146378160153344429844948;
      wtt[5]  =  0.0166609651269037212368055;
      wtt[6]  =  0.0281343344093849194875108;
      wtt[7]  =  0.00214364311909247909952968;
      wtt[8]  =  0.0331943417648083019611294;
      wtt[9]  = -0.0142483904822400753741381;
      wtt[10] =  0.0416812529998231580614934;
      wtt[11] =  0.00323522898964475022578598;
      wtt[12] =  0.000813552611571786631179287;
      wtt[13] =  0.00228403610697848813660369;
      wtt[14] =  0.0338971925236628645848112;
      wtt[15] = -0.00644725595698339499416262;
      for(m=0;m<=odirs-1;++m)
        wt[m] = wtt[wp24[m]];
    break;
  }

  // Evaluate mu and eta for octant 1
  m = 0;
  for(i=0;i<=levels/2-1;++i)
    for(j=0;j<=(levels/2)-i-1;++j)
    {
      mu[m]  = att[i];
      eta[m] = att[j];
      ++m;
    }

  // Evaluate mu and eta for octants 2-4
  for(int octant=2;octant<=4;++octant)
    for(n=0;n<=odirs-1;++n)
    {
      m = (octant-1)*odirs+n;
      switch (octant)
      {
        case 2:
	  mu[m]  = -mu[n];
          eta[m] =  eta[n];
          wt[m]  =  wt[n];
	break;

	case 3:
	  mu[m]  = -mu[n];
          eta[m] = -eta[n];
          wt[m]  =  wt[n];
	break;

	case 4:
	  mu[m]  =  mu[n];
          eta[m] = -eta[n];
          wt[m]  =  wt[n];
	break;
      }
    }

  // Evaluate mu and eta for octants 5-8
  for(n=0;n<=dirs/2-1;++n)
  {
    mu[n+dirs/2]  = mu[n];
    eta[n+dirs/2] = eta[n];
    wt[n+dirs/2]  = wt[n];
  }

  // Evaluate xi for all octants
  for(n=0;n<=dirs/2-1;++n)
    xi[n] = sqrt(1.0-mu[n]*mu[n]-eta[n]*eta[n]);

  for(n=0;n<=dirs/2-1;++n)
    xi[n+dirs/2] = -xi[n];

  // Normalize the quadrature set
  double wsum = 0.0;
  for(n=0;n<=dirs-1;++n)
   wsum = wsum + wt[n];
  for(n=0;n<=dirs-1;++n)
   wt[n] = wt[n]*(norm/wsum);

  att.clear();
  wtt.clear();
}

}

//---------------------------------------------------------------------------//
//                              end of LSquad.cc
//---------------------------------------------------------------------------//
