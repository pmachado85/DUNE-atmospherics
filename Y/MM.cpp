#include <iostream>              // Initūs canonicī c++
#include <string>                // Interpositio lineae
#include <fstream>               // Initūs exitūsque
#include <cstdlib>               // Bibliothecae canonicae ūtilitātum generalium
#include <stdio.h>               // Indecēs
#include <math.h>                // Mathēmatica prior
#include <vector>                // Vectōrēs 
#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <chrono>
#include <ctime>
//
//    GNU Bibliothecae scientificae
//
#include <gsl/gsl_sf_bessel.h>   // Fūnctiōnēs Bessēlis
#include <gsl/gsl_errno.h>       // Audītiō erratōrum
#include <gsl/gsl_spline.h>      // Interpolātiō
#include <gsl/gsl_math.h>        // Mathēmatica prior
#include <gsl/gsl_integration.h> // Integrātiōnēs
#include <gsl/gsl_sf_log.h>      // Logarithmī
#include <gsl/gsl_vector.h>      // Vectōrēs
#include <gsl/gsl_matrix.h>      // Mātrīcēs
#include <gsl/gsl_rng.h>         // Numerī āleātōriī
#include <gsl/gsl_randist.h>     // Distribūtiōnēs numerōrum āleātōriōrum
#include <gsl/gsl_sf_pow_int.h>  // Potentiae
#include <gsl/gsl_sf_gamma.h>    // Fūnctio gamma
#include <gsl/gsl_monte.h>       // Montecarlo Integrātiō
#include <gsl/gsl_monte_plain.h> // Montecarlo - Plain Integrātiō
#include <gsl/gsl_monte_miser.h> // Montecarlo - Miser Integrātiō
#include <gsl/gsl_monte_vegas.h> // Montecarlo - Vegas Integrātiō
#include <gsl/gsl_poly.h>        // Polynomia
#include <gsl/gsl_sf_erf.h>      // Fūnctio erratōrum
#include <gsl/gsl_multimin.h>    // Minima invenire
#include <gsl/gsl_histogram2d.h> // 2D histogrammae
//
#include <armadillo>             // Bibliothecae Armadillo
//
//  Dēfīnītiōnēs
//
#define NEtru 109   // Numerus lacuum energiae verae
#define NErec 60 //100    // Numerus lacuum energiae reficientis
#define Ntrec 40 //200   // Numerus lacuum anguli reficientis 

using namespace std;
using namespace arma;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//                                                                                                                                  //
//                                         Migration Matrix from different Generators                                               //
//                                           Mātrix Migrātionis ē Variīs Generatrīs                                                 //
//                                                                                                                                  //
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//                       Angle smearing of trimomentum                     //
//                        Trimomenti Anguli Tendendum                      //
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

arma::mat R(3, 3, fill::eye), Rphi(3, 3, fill::eye), Rth(3, 3, fill::eye);

std::vector<double> angsmp(std::vector<double> p, double sTh, gsl_rng * r){

  arma::vec pt = {p[0], p[1], p[2]};

  double the = acos(pt(2)/sqrt(pt(0)*pt(0) + pt(1)*pt(1) + pt(2)*pt(2)));
  double phi = atan(pt(1)/pt(0));
  
  if(pt(0) < 0.) phi += M_PI;

  R = {{sin(phi)*sin(phi) + cos(the)*cos(phi)*cos(phi), sin(phi)*cos(phi)*(cos(the) - 1.), -cos(phi)*sin(the)},
       {sin(phi)*cos(phi)*(cos(the) - 1.), cos(phi)*cos(phi) + cos(the)*sin(phi)*sin(phi), -sin(phi)*sin(the)},
       {cos(phi)*sin(the), sin(phi)*sin(the), cos(the)}}; // Mātrix rotātiōnis in axem z trimomentum ponere
  
  double thr = gsl_ran_gaussian(r, sTh);    // Theta angulum āleātōrium optāmus cum distribūtiōne gaussianis 
  double phr = 2.*M_PI*gsl_rng_uniform (r); // Phi angulum āleātōrium optāmus cum distribūtiōne ūnifōrmī
  
  Rphi = {{cos(phr), -sin(phr), 0.}, {sin(phr), cos(phr), 0.}, {0., 0., 1.}}; // Mātrix rotātiōnis circum axem z      
  Rth  = {{cos(thr), 0., -sin(thr)}, {0., 1., 0.}, {sin(thr), 0., cos(thr)}}; // Mātrix rotātiōnis circum axem y

  pt = R.t() * Rphi * Rth * R * pt; // Trimomentō rotātiōnem plēnam applicāmus

  std::vector<double> ps = {pt(0), pt(1), pt(2)};
  
  return ps;

}

//----------------------------------------------------------------------------------//
//                                  Main program                                    //
//                              Programma principālis                               //
//----------------------------------------------------------------------------------//

int main(int argc, const char* argv[]){

  cout << std::setprecision(8);// << scientific << showpoint;

  cout << endl;
  cout << "Mātrix Migrātionis pro analysis neutrinōrum atmosphaericōrum" << endl;
  cout << endl;

  //---------------------------------------------------------------------//
  //                          Parametra in [GeV]                         //
  //---------------------------------------------------------------------//

  double me, mm, mt;
  double mp, mn, mS0, mSp, mSm, mL, mXi0, mXi, mOm, mLc, mScpp, mScp;
  double mpip, mpi0, mK, mKL, mK0, mD0, mDp, mDsp;

  //Leptons
  
  me = 0.511e-3; mm = 105.66e-3; mt = 1.77686;

  //Mesons
  
  mpip = 139.5706e-3; mpi0 = 134.9770e-3;
  mK   = 493.677e-3;  mK0  = 497.611e-3;
  mD0  = 1864.83e-3;  mDp  = 1869.65e-3; mDsp = 1968.34e-3;

  //Baryons

  mp    = 938.2721e-3; mn    = 939.56541e-3;
  mL    = 1115.683e-3; mLc   = 2286.46e-3;
  mS0   = 1192.642e-3; mSp   = 1189.37e-3; mSm   = 1197.449e-3;
  mXi   = 1321.71e-3;  mXi0  = 1314.86e-3;
  mOm   = 1672.45e-3;
  mScpp = 2453.97e-3;  mScp  = 2452.9e-3;
    
  //Kinetic energies thresholds

  double Kthe, Kthm, Ktht, Kthph, KthM, KthB, Kthp;

  Kthe = 0.010; Kthm  = 0.005; Ktht = 1.0; Kthph = 0.010; //Leptons + photons
  KthM = 0.030;               // Mesons
  Kthp = 0.030; KthB = 0.030; // Baryons

  //Angle and momentum resolutions

  double slTh, spTh, spK, seK, smK;

  slTh =  2.*M_PI/180.;
  spTh = 10.*M_PI/180.;

  //  In percentage
  
  spK = 0.10; 
  seK = 0.05;
  smK = 0.05;

  // PDG numbering scheme

  const int ph = 22;

  const int e = 11; const int ne = 12;
  const int m = 13; const int nm = 14;
  const int t = 15; const int nt = 16;

  const int p   = 2212; const int n    = 2112; const int L   = 3122;
  const int Sp  = 3222; const int S0   = 3212; const int Sm  = 3112;
  const int Xi0 = 3322; const int Xi   = 3312; const int Om  = 3334;
  const int Lcp = 4122; const int Scpp = 4222; const int Scp = 4212; 

  const int pip = 211; const int pi0 = 111;
  const int K   = 321; const int K0  = 311;
  const int Dp  = 411; const int D0  = 421; const int Dsp = 431; 

  const int Ar = 1000180400; 

  //-----------------------------------------------------------------------------//
  //                         Tabulāria initūs exitūsque                          //
  //-----------------------------------------------------------------------------//

  double ml;

  // Masses to agree with generator output (Genie-NuWro give total energy -- Fluka gives kinetic energy)
  double meg, mmg, mtg;
  double mpipg, mpi0g, mKg, mK0g;
  double mpg, mng, mS0g, mSpg, mSmg, mXi0g, mXig, mOmg, mLg;
  double mLcg, mScppg, mScpg, mD0g, mDpg, mDspg;

  string inf, onf, gens, EvFf, xscf, nus;

  int inuf;

  cout << endl << "Quaeso, initialem neutrinum optā:" << endl << endl;
  cout << "12 -> nu_e,    " << '\t' << "-12 --> nubar_e," << endl;
  cout << "14 -> nu_mu,   " << '\t' << "-14 --> nubar_mu," << endl;
  cout << "16 -> nu_tau,  " << '\t' << "-16 --> nubar_tau." << endl << endl; 
  cin >> inuf;
  //inuf = 12;

  cout << endl;

  int fer;

  if(inuf > 0) fer = 1; else fer = -1;

  switch(fer*inuf){

  case ne:
    if(inuf > 0) {inf = "nu_e";     onf = "e";   }
    else         {inf = "nubar_e";  onf = "ebar";}
    ml = me;
    break;

  case nm:
    if(inuf > 0) {inf = "nu_mu";     onf = "mu";}
    else         {inf = "nubar_mu";  onf = "mubar";}
    ml = mm;
    break;
    
  case nt:
    if(inuf > 0) {inf = "nu_tau";     onf = "t";}
    else         {inf = "nubar_tau";  onf = "tbar";}
    ml = mt;
    break;

  default:
    break;

  }

  int gen;

  cout << endl << "Quaeso, generatrum optā:" << endl << endl;
  cout << "0 --> GiBUU," << '\t' << "1 --> NuWro," << '\t' << "2 --> Genie," << '\t' << "3 --> Fluka," << endl << endl;
  //gen = 3;
  cin >> gen;

  int Ntot = 1e5;

  switch(gen){

  case 0:
    gens = "GiBUU";
    Ntot = 1e3;
    meg = me; mmg = mm;  mtg = mt;  mpg = mp;
    mpipg = mpip; mpi0g = mpi0;
    mS0g = mS0; mSpg = mSp; mSmg = mSm; mXi0g = mXi0; mXig = mXi;
    mOmg = mOm; mLg = mL;
    mLcg = mLc; mScppg = mScpp; mScpg = mScp;
    mD0g = mD0; mDpg  = mDp; mDspg = mDsp;
    break;

  case 1:
    gens = "NuWro";
    Ntot = 1e5;
    meg = me; mmg = mm;  mtg = mt;  mpg = mp;
    mpipg = mpip; mpi0g = mpi0;
    mD0g = mD0; mDpg  = mDp; mDspg = mDsp;
    mS0g = mS0; mSpg = mSp; mSmg = mSm; mXi0g = mXi0; mXig = mXi;
    mOmg = mOm; mLg = mL;
    mLcg  = mLc; mScppg = mScpp; mScpg = mScp;
    break;
    
  case 2:
    gens = "Genie";
    Ntot = 1e5;
    meg = me;  mmg = mm;  mtg = mt;
    mpipg = mpip; mpi0g = mpi0; mKg = mK; mK0g = mK0;
    mD0g = mD0; mDpg  = mDp; mDspg = mDsp;
    mpg = mp;
    mS0g = mS0; mSpg = mSp; mSmg = mSm;
    mXi0g = mXi0; mXig = mXi;
    mOmg = mOm; mLg = mL;
    mLcg  = mLc; mScppg = mScpp; mScpg = mScp;
    break;

  case 3:
    gens = "Fluka";
    Ntot = 1e4;
    meg  = 0.; mmg = 0.;  mtg = 0.;  
    mpipg = 0.; mpi0g = 0.; mKg = 0.; mK0g = 0.;
    mD0g = 0.; mDpg  = 0.; mDspg = 0.;
    mpg  = 0.;
    mS0g = 0.; mSpg = 0.; mSmg = 0.;
    mXi0g = 0.; mXig = 0.;
    mOmg = 0.; mLg = 0.;
    mLcg  = 0.; mScppg = 0.; mScpg = 0.;
    break;

  default:
    break;

  }

  cout << endl << "Optāvisti " << inf << " ē generatrō " << gens << endl;

  EvFf = "./" + gens + "/Events/data/Evs_" + inf + ".txt";

  char const * evs_nu_c = EvFf.c_str();

  ifstream evs{evs_nu_c, ios::in};
  
  if (!evs.is_open()){
    cerr << endl << evs_nu_c << " tabulārium aperire non possum." << endl;
    return EXIT_FAILURE;
  }

  ofstream mm0p, mm1p, mm2p;
  string   MM0ps, MM1ps, MM2ps;

  MM0ps = "../MM/MM-" + onf + "-1l0p-" + gens + ".txt";
  MM1ps = "../MM/MM-" + onf + "-1l1p-" + gens + ".txt";
  MM2ps = "../MM/MM-" + onf + "-1l2p-" + gens + ".txt";

  mm0p.open(MM0ps.c_str(), ios::out);
  mm1p.open(MM1ps.c_str(), ios::out);
  mm2p.open(MM2ps.c_str(), ios::out);

  mm0p << scientific << showpoint << std::setprecision(6);
  mm1p << scientific << showpoint << std::setprecision(6);
  mm2p << scientific << showpoint << std::setprecision(6);

  //------------------------------------------------------------------------------------// 
  //             Energiam veram reficientem angulumque reficientem dividimus            //
  //------------------------------------------------------------------------------------//

  auto start = std::chrono::system_clock::now();
      
  std::time_t start_time = std::chrono::system_clock::to_time_t(start);

  cout << endl << endl << "Programma incepit in: " << std::ctime(&start_time) << endl << endl;

  std::vector<double> xsec(NEtru + 1);
  
  double Enut[NEtru + 1];

  double Enur[NErec + 1];
  for(int b = 0; b <= NErec; b++) Enur[b] = pow(10., -2. + (2. + 2.)*b/NErec);

  // double Enur[NErec + 1] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.,
  //   1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8,
  //   5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, 6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8,
  //   9., 9.2, 9.4, 9.6, 9.8,
  //   10., 12., 14., 16., 18., 20., 22., 24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 46., 48.,
  //   50., 52., 54., 56., 58., 60., 62., 64., 66., 68., 70., 72., 74., 76., 78., 80., 82., 84., 86., 88.,
  //   90., 92., 94., 96., 98., 100.};

  double Thr[Ntrec + 1];

  for(int t = 0; t <= Ntrec; t++) Thr[t] = M_PI*t/Ntrec;

  //---------------------------------------------------------------------//
  //                 Mātrīcem migrātiōnis computāmus                     //
  //---------------------------------------------------------------------//

  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_rng_env_setup();

  T = gsl_rng_ranlxs2;
  r = gsl_rng_alloc (T);

  std::vector<std::vector<double>> Evs0p0pi, Evs1p0pi, Evs2p0pi, Evsoth;

  string rec;
  double x;

  double Kth;

  int Nl = 0; int Np = 0; int Npip = 0; int NM = 0; int NB = 0; int Nph = 0; int NAr = 0;

  std::vector<std::vector<double>> Nu, Event;

  gsl_histogram2d * h0p = gsl_histogram2d_alloc(NErec, Ntrec); gsl_histogram2d_set_ranges(h0p, Enur, NErec + 1, Thr, Ntrec + 1);
  gsl_histogram2d * h1p = gsl_histogram2d_alloc(NErec, Ntrec); gsl_histogram2d_set_ranges(h1p, Enur, NErec + 1, Thr, Ntrec + 1);
  gsl_histogram2d * h2p = gsl_histogram2d_alloc(NErec, Ntrec); gsl_histogram2d_set_ranges(h2p, Enur, NErec + 1, Thr, Ntrec + 1);

  double Enr0p, pxr0p, pyr0p, pzr0p, thr0p;
  double Enr1p, pxr1p, pyr1p, pzr1p, thr1p;
  double Enr2p, pxr2p, pyr2p, pzr2p, thr2p;

  cube M0p(NEtru + 1, NErec, Ntrec, fill::zeros);
  cube M1p(NEtru + 1, NErec, Ntrec, fill::zeros);
  cube M2p(NEtru + 1, NErec, Ntrec, fill::zeros);

  mat R(3, 3, fill::eye), Rphi(3, 3, fill::eye), Rth(3, 3, fill::eye);

  int Nev = 0;

  int n0p = 0, n1p = 0, n2p = 0;

  int ie = 0;

  if(fer*inuf == nm) {ie = 1; Enut[0] = 0.1;}
  
  for(int in = ie; in < NEtru; in++){

    Nev = 0; n0p = 0, n1p = 0, n2p = 0;
    
    while(std::getline(evs, rec) && Nev < Ntot){
      
      istringstream is(rec);
      
      std::vector<double> row((std::istream_iterator<double>(is)), std::istream_iterator<double>());
      
      //We write down initial neutrino properties
      if (row.size() == 2 || row.size() == 3) Enut[in] = row[1]; //cout << row[1] << endl;}
      
      //We verify what are the visible particles
      else if (row.size() == 5) {

	//switch((int)row[0]){
	switch((int)fabs(row[0])){
	
	case e:// electrons
	  if(row[1] - meg >= Kthe) {Event.push_back(row); Nl++;}
	  break;
	  
	case m:// muons
	  if(row[1] - mmg >= Kthm) {Event.push_back(row); Nl++;}
	  break;
	
	case t:// taus
	  if(row[1] - mtg >= Ktht) {Event.push_back(row); Nl++;}
	  break;
	
	case p:// protons
	  if(row[1] - mpg > Kthp) {Event.push_back(row); Np++;}
	  break;

	case pip:// charged pions
	  if(row[1] - mpipg >= KthM) {Event.push_back(row); Npip++;}
	  break;
	  
	case pi0:// neutral pions
	  if(row[1] - mpi0g >= KthM) NM++;
	  break;
	  
	case K:// charged Kaons
	  if(row[1] - mKg >= KthM)  NM++;
	  break;
	  
	case K0:// neutral Kaons
	  if(row[1] - mK0g >= KthM)  NM++;
	  break;

	// case D0:// D0
	//   if(row[1] - mD0g >= KthM)  NM++;
	//   break;

	// case Dp:// D+
	//   if(row[1] - mDpg >= KthM)  NM++;
	//   break;

	// case Dsp:// Ds+
	//   if(row[1] - mDspg >= KthM)  NM++;
	//   break;

	// case L:// Lambdas
	//   if(row[1] - mLg >= KthB)  NB++;
	//   break;

	// case S0:// Sigma0
	//   if(row[1] - mS0g >= KthB)  NB++;
	//   break;

	// case Sp:// Sigma+
	//   if(row[1] - mSpg >= KthB)  NB++;
	//   break;

	// case Sm:// Sigma-
	//   if(row[1] - mSmg >= KthB)  NB++;
	//   break;

	// case Xi0:// Xi0
	//   if(row[1] - mXi0g >= KthB)  NB++;
	//   break;

	// case Xi:// Xi
	//   if(row[1] - mXig >= KthB)  NB++;
	//   break;

	// case Om:// Omega-
	//   if(row[1] - mOmg >= KthB)  NB++;
	//   break;

	// case Lcp:// Lambda_c+
	//   if(row[1] - mLcg >= KthB)  NB++;
	//   break;
	  
	// case Scpp:// Sigma_c++
	//   if(row[1] - mScppg >= KthB)  NB++;
	//   break;

	// case Scp:// Sigma_c+
	//   if(row[1] - mScpg >= KthB)  NB++;
	//   break;
	  
	case ph:// photons
	  if(row[1] >= Kthph) Nph++;
	  break;

	case Ar:// Argon nucleus
	  NAr++;
	  break;
	  
	default:
	  break;
	  
	}
	
      }
      
      // We clasify the events according the number of visible protons and mesons
      else if(rec.length() == 0) {

	//cout << Np << '\t' << Nl << '\t' << Npip << endl;

	//if(in == 10 && Nev <= 100){

	//cout << Enut[in] << '\t' << Nl << '\t' << Np << '\t'
	//     << Npip << '\t' << NM << '\t' << NB << '\t' << NAr << '\t' << Nph << endl;
	
	  // for ( const auto &row : Event ){
	  //   for ( double x : row ) std::cout << x << '\t';
	  //   std::cout << std::endl;
	  // }
	  // std::cout << std::endl;

	//}

	// Keep the event if there is only one visible lepton and no meson or heavy baryon
	// Ēventum sī ūnum tantum lepton vīsibilis erit neque meson vel baryon eriut tenemus
	if(Nl == 1 && Npip == 0 && NM == 0 && Nph == 0 && NB == 0){
	
	  if(Np >= 3) Np = 3;
	  
	  switch(Np){
	    
	  case 0:// CC1l0p0pi
	    for(const auto &row: Event) Evs0p0pi.push_back(row);
	    n0p++;
	    break;
	    
	  case 1:// CC1l1p0pi
	    for(const auto &row: Event) Evs1p0pi.push_back(row);
	    n1p++;
	    break;
	    
	  case 2:// CC1l2p0pi
	    for(const auto &row: Event) Evs2p0pi.push_back(row);
	    n2p++;
	    break;
	    
	  case 3:// CC1l3+p0pi
	    for(const auto &row: Event) Evsoth.push_back(row);
	    break;
	    
	  default:
	    break;
	    
	  } 
	  
	}

	Event.erase(Event.begin(), Event.end());
	
	Npip = 0; NM = 0; NB = 0; Nph = 0; NAr = 0;
	
        Nl = 0; Np = 0; 
	
	Nev++; 
	
      }
      
    }

    cout << in << '\t' << Enut[in] << endl;
    /*
    if(in == 44){

      int c = 0;
    
      cout << "1l-0p Events:" << endl << endl;
      
      for ( const auto &row : Evs0p0pi ){

    	for ( double x : row ) std::cout << x << '\t';
    	std::cout << std::endl;

    	c++;

    	if (c == 150) break;
	
      }
      
      std::cout << std::endl;

    }
      
      
      cout << "1l-1p Events:" << endl << endl;
      
      for ( const auto &row : Evs1p0pi ){
	
    	for ( double x : row ) std::cout << x << '\t';
    	std::cout << std::endl;
 
      }

    }
    
      std::cout << std::endl;
      
      cout << "1l-2p Events:" << endl << endl;
      
      for ( const auto &row : Evs2p0pi ){
	
	for ( double x : row ) std::cout << x << '\t';
	std::cout << std::endl;
      }

      std::cout << std::endl;
    }
    
      
      cout << "1l-3+p Events:" << endl << endl;
      
      for ( const auto &row : Evsoth )
      {
      for ( double x : row ) std::cout << x << '\t';
      std::cout << std::endl;
      } 
    */
    // We compute next the reconstructed neutrino energy and angle
    
    //+++++++++++++++++++++++++++++++++++++++//
    //              1l0p events              //
    //+++++++++++++++++++++++++++++++++++++++//

    int c0p = 0, c1p = 0, c2p = 0;
    
    Enr0p = 0.; pxr0p = 0.; pyr0p = 0.; pzr0p = 0.;
    
    for(const auto &row0p : Evs0p0pi){

      double unp = gsl_ran_gaussian(r, seK);

      Enr0p = sqrt(ml*ml + (1. + unp) * (row0p[2]*row0p[2] + row0p[3]*row0p[3] + row0p[4]*row0p[4]));

      std::vector<double>  pl = {row0p[2], row0p[3], row0p[4]};      
      std::vector<double> pls = angsmp(pl, slTh, r);

      thr0p = acos(pls[2]/sqrt(pls[0]*pls[0] + pls[1]*pls[1] + pls[2]*pls[2]));
      
      gsl_histogram2d_increment(h0p, Enr0p, thr0p);
      
    }

    
    
    //+++++++++++++++++++++++++++++++++++++++//
    //              1l1p events              //
    //+++++++++++++++++++++++++++++++++++++++//
    
    Enr1p = 0.; pxr1p = 0.; pyr1p = 0.; pzr1p = 0.;
    
    Np = 0;
    
    for(const auto &row1p : Evs1p0pi){
      
      if((int)fabs(row1p[0]) == (int)e || (int)fabs(row1p[0]) == (int)m || (int)fabs(row1p[0]) == (int)t) {

	double unpl1p = gsl_ran_gaussian(r, seK);

	Enr1p = sqrt(ml*ml + (1. + unpl1p) * (row1p[2]*row1p[2] + row1p[3]*row1p[3] + row1p[4]*row1p[4]));
	
	std::vector<double> pl1p  = {row1p[2], row1p[3], row1p[4]};
	std::vector<double> pls1p = angsmp(pl1p, slTh, r);

	pxr1p = (1. + unpl1p) * pls1p[0]; pyr1p = (1. + unpl1p) * pls1p[1]; pzr1p = (1. + unpl1p) * pls1p[2];

     }
      
      if((int)row1p[0] == (int)p) {

	double unpp = gsl_ran_gaussian(r, spK);

	Enr1p += sqrt(mp*mp + (1. + unpp) * (row1p[2]*row1p[2] + row1p[3]*row1p[3] + row1p[4]*row1p[4])) - mp;

	std::vector<double> pp1p  = {row1p[2], row1p[3], row1p[4]};
	std::vector<double> pps1p = angsmp(pp1p, spTh, r);
	
	pxr1p += (1. + unpp) * pps1p[0]; pyr1p += (1. + unpp) * pps1p[1]; pzr1p += (1. + unpp) * pps1p[2];
	
	Np++;
	
      }
      
      if(Np == 1) {
	
	thr1p = acos(pzr1p/sqrt(pxr1p*pxr1p + pyr1p*pyr1p + pzr1p*pzr1p));
	
	gsl_histogram2d_increment(h1p, Enr1p, thr1p);
	
	Np = 0;
	
      }

      c1p++;
      
    }

    //+++++++++++++++++++++++++++++++++++++++//
    //              1l2p events              //
    //+++++++++++++++++++++++++++++++++++++++//
    
    Enr2p = 0.; pxr2p = 0.; pyr2p = 0.; pzr2p = 0.;
    
    Np = 0;
    
    for(const auto &row2p : Evs2p0pi){

      
      if((int)fabs(row2p[0]) == (int)e || (int)fabs(row2p[0]) == (int)m || (int)fabs(row2p[0]) == (int)t) {

	double unpl2p = gsl_ran_gaussian(r, seK);

	Enr2p = sqrt(ml*ml + (1. + unpl2p) * (row2p[2]*row2p[2] + row2p[3]*row2p[3] + row2p[4]*row2p[4]));
	
	std::vector<double> pl2p  = {row2p[2], row2p[3], row2p[4]};
	std::vector<double> pls2p = angsmp(pl2p, slTh, r);	

	pxr2p = (1. + unpl2p) * pls2p[0]; pyr2p = (1. + unpl2p) * pls2p[1]; pzr2p = (1. + unpl2p) * pls2p[2];
	
      }
      
      if((int)row2p[0] == (int)p) {

	double unpp = gsl_ran_gaussian(r, spK);

	Enr2p += sqrt(mp*mp + (1. + unpp) * (row2p[2]*row2p[2] + row2p[3]*row2p[3] + row2p[4]*row2p[4])) - mp;
	
	std::vector<double> pp2p  = {row2p[2], row2p[3], row2p[4]};
	std::vector<double> pps2p = angsmp(pp2p, spTh, r);
	
	pxr2p += (1. + unpp) * pps2p[0]; pyr2p += (1. + unpp) * pps2p[1]; pzr2p += (1. + unpp) * pps2p[2];
	
	Np++;
	
      }
      
      if(Np == 2) {
	
	thr2p = acos(pzr2p/sqrt(pxr2p*pxr2p + pyr2p*pyr2p + pzr2p*pzr2p));
	
	gsl_histogram2d_increment(h2p, Enr2p, thr2p);
	
	Np = 0;

	c2p++;
	
      }
      
    }

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    //          Writing down the element of the migration matrix              //
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
    
    for(int i = 0; i < NErec; i++){
      
      for(int j = 0; j < Ntrec; j++){

	M0p(in, i, j) = gsl_histogram2d_get(h0p, i, j)/Ntot;
	M1p(in, i, j) = gsl_histogram2d_get(h1p, i, j)/Ntot;
	M2p(in, i, j) = gsl_histogram2d_get(h2p, i, j)/Ntot;

	// cout << Enur[i] << '\t' << Thr[j] << '\t' << Enut[in] << '\t'
	//      << M0p(in, i, j) << '\t'
	//      << M1p(in, i, j) << '\t'
	//      << M2p(in, i, j) << endl;
	//      << gsl_histogram2d_get(h0p, i, j) <<'\t'
	// << gsl_histogram2d_get(h1p, i, j) <<'\t'
	// << gsl_histogram2d_get(h2p, i, j) << endl;

      }
      
    }

    gsl_histogram2d_reset(h0p);
    gsl_histogram2d_reset(h1p);
    gsl_histogram2d_reset(h2p);

    Evs0p0pi.erase(Evs0p0pi.begin(), Evs0p0pi.end());
    Evs1p0pi.erase(Evs1p0pi.begin(), Evs1p0pi.end());
    Evs2p0pi.erase(Evs2p0pi.begin(), Evs2p0pi.end());
    
    Evsoth.erase(Evsoth.begin(), Evsoth.end());

    Nev = 0;

    //cout << Enut[in] << '\t' << xsec[in] << endl;

  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  //                                  Writing the Migration Matrix                                 //
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
  
  for(int i = 0; i < NErec; i++){
    
    for(int j = 0; j < Ntrec; j++){
      
      for(int in = 0; in < NEtru; in++){

	mm0p << Enur[i] << '\t' << Thr[j] << '\t' << Enut[in] << '\t' <<  M0p(in, i, j) << endl;
	mm1p << Enur[i] << '\t' << Thr[j] << '\t' << Enut[in] << '\t' <<  M1p(in, i, j) << endl;
	mm2p << Enur[i] << '\t' << Thr[j] << '\t' << Enut[in] << '\t' <<  M2p(in, i, j) << endl;
		
	
      }
      
    }
    
  }
  
  /*    
  for(int in = 0; in < NEtru; in++){

    for(int i = 0; i < NErec; i++){

      for(int j = 0; j < Ntrec; j++){	
	
	mm0p << Thr[j] << '\t' << Enur[i] << '\t' << Enut[in] << '\t' <<  M0p(in, i, j) << endl;
	mm1p << Thr[j] << '\t' << Enur[i] << '\t' << Enut[in] << '\t' <<  M1p(in, i, j) << endl;
	mm2p << Thr[j] << '\t' << Enur[i] << '\t' << Enut[in] << '\t' <<  M2p(in, i, j) << endl;
		
	
      }
      
    }
    
  }
  */
  auto end = std::chrono::system_clock::now();
      
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::time_t end_time = std::chrono::system_clock::to_time_t(end);

  cout << endl
       << "Mātrix Migrātionis ē " << gens << " facta est. Tempus elapsum: " << elapsed_seconds.count() << "s. Programma ā "
       << std::ctime(&end_time) << " perfecit." << endl;
  
  evs.close();
  
  return 0;
 
}
