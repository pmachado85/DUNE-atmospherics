//*******************************************************
// Calculo el numero de eventos. Le paso el parametro de
// masas solar, t23, dcp y t12 al programa.
// Hago un scan sobre t12 y la masa solar.
// Hacemos bines en czenit. 
// Usamos la migration matrix correcta.
// Uso la migration matrix de NuWro
// El numero de target está dividido por la masa de Ar
//*******************************************************
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <time.h>
#include <algorithm>
#include <vector>
#include <iomanip>
#include "./LibProb/atm.fluxes.hh"
#include "./LibMath/errmsg.hh"
#include "./LibMath/fitter-P.hh"
#include "./LibMath/progress.hh"
#include "./LibMath/constants.hh"
#include "./LibMath/misc.hh"
#include "./LibMath/interp.hh"
#include "./LibProb/prb.earth.hh"
#include "./LibMath/progress.hh"


using namespace std;
using namespace Probs::Atmos;
using namespace Tools;
using namespace Probs;


enum FLV {FLV_elc, FLV_mu, FLV_tau, NUM_FLV};
enum CH {CH_p, CH_ap, NUM_CH};
enum FPS {F1L0P, F1L1P, F1L2P, F1L3P, NPRO};


const double GeV = 1.e3;        //From GeV to MeV
const double Me = 0.5;          //Electron mass in MeV
const double Mm = 105.65;       //Muon mass in MeV
const double DUNEm = 4.e10;     //In grams. We assume a mass of 40 kton.
const double NA = 6.022e23;     //Avogadro number
const double eV = 1.e9;         //From GeV to eV
const double Arm = 39.948;   //Argon mass 


static double RADIAN = M_PI/180.;


const int NErecM = 100, NTrecM = 200, NEtrueM = 109;
double dltte1 = 0.05, dltte2 = 0.2, dltte3 = 2.;
double Ermin = 0.01, Ermax = 1e3, dltlr = (log10(Ermax) - log10(Ermin))/double(NErecM);
float MigMat[NUM_CH][NUM_FLV][NErecM][NTrecM][NEtrueM];
float Xsec[NUM_CH][NUM_FLV][NEtrueM];
float ErecE[NErecM], EtrueM[NEtrueM+1], TrecM[NTrecM+1];
void ReadMig(){
  
  ifstream infile;
  string line;

//  char files[NUM_CH][NUM_FLV][100] = {"./MMNuwro/MMNuE1L0P.dat",
//				      "./MMNuwro/MMNuMu1L0P.dat",
//				      "./MMNuwro/MMNuTau1L0P.dat",
//				      "./MMNuwro/MMNuEBar1L0P.dat",
//				      "./MMNuwro/MMNuMuBar1L0P.dat",
//				      "./MMNuwro/MMNuTauBar1L0P.dat"};


  char files[NUM_CH][NUM_FLV][100] = {"./MMNuwro/MMNuE1L1P.dat",
				      "./MMNuwro/MMNuMu1L1P.dat",
				      "./MMNuwro/MMNuTau1L1P.dat",
				      "./MMNuwro/MMNuEBar1L1P.dat",
				      "./MMNuwro/MMNuMuBar1L1P.dat",
				      "./MMNuwro/MMNuTauBar1L1P.dat"};



  

  
  
//  char files[NUM_CH][NUM_FLV][100] = {"./MMGenie/GMM_NuE_1L1P.dat",
//				      "./MMGenie/GMM_NuMu_1L1P.dat",
//				      "./MMGenie/GMM_NuTau_1L1P.dat",
//				      "./MMGenie/GMM_NuEBar_1L1P.dat",
//				      "./MMGenie/GMM_NuMuBar_1L1P.dat",
//				      "./MMGenie/GMM_NuTauBar_1L1P.dat"};


  char fxsec[NUM_CH][NUM_FLV][100] = {"./MMNuwro/XSECNuE.dat",
				      "./MMNuwro/XSECNuMu.dat",
				      "./MMNuwro/XSECNuTau.dat",
				      "./MMNuwro/XSECNuEBar.dat",
				      "./MMNuwro/XSECNuMuBar.dat",
				      "./MMNuwro/XSECNuTauBar.dat"};


  
//  char fxsec[NUM_CH][NUM_FLV][100] = {"./MMGenie/GXSECNuE.dat",
//				      "./MMGenie/GXSECNuMu.dat",
//				      "./MMGenie/GXSECNuTau.dat",
//				      "./MMGenie/GXSECNuEBar.dat",
//				      "./MMGenie/GXSECNuMuBar.dat",
//				      "./MMGenie/GXSECNuTauBar.dat"};

  
  for(int e=0; e <= NEtrueM; e++){
    
    if(e <= 18)
      EtrueM[e] = 0.1 + double(e) * dltte1;
    else if(e > 18 && e <= 63)
      EtrueM[e] = 1. + double(e - 18) * dltte2;
    else
      EtrueM[e] = 10. + double(e - 63) * dltte3;

  }

  
  for(int t=0; t <= NTrecM; t++)
    TrecM[t] = M_PI*double(t)/double(NTrecM);

  for(int e=0; e < NErecM; e++)
    ErecE[e] = exp10(log10(Ermin) + double(e) * dltlr);
    
  double del1, del2, del3, mig;
  //int index = 0;

  //Initialize to zero the migration matri
  for(int s=0; s < NUM_CH; s++)
    for(int fl=0; fl < NUM_FLV; fl++)
      for(int er=0; er < NErecM; er++)
	for(int tr=0; tr < NTrecM; tr++)
	  for(int et=0; et < NEtrueM; et++)
	    MigMat[s][fl][er][tr][et] = 0.; 
  

  for(int s=0; s < NUM_CH; s++)
    for(int fl=0; fl < NUM_FLV; fl++)
      for(int et=0; et < NEtrueM; et++)
	Xsec[s][fl][et] = 0.;

  for(int s=0; s < NUM_CH; s++)
    for(int fl=0; fl < NUM_FLV; fl++){
      
      infile.open(fxsec[s][fl]);
      for(int t=0; t < NEtrueM; t++){
	getline(infile,line);
	sscanf(line.c_str(), "%lf    %lf\n", &del1, &mig);
	Xsec[s][fl][t] = mig;
      }
      infile.close();
    }
  
  
  //Read and combine the migration matrix using the reconstructed bins in energy
  double sum;
  for(int s=0; s < NUM_CH; s++)
    for(int fl=0; fl < NUM_FLV; fl++){

      infile.open(files[s][fl]);
      printf("%s\n",files[s][fl]);
      
      for(int et=0; et < NEtrueM; et++){
	sum = 0.;
	
	for(int er=0; er < NErecM; er++)
	  for(int zr=0; zr < NTrecM; zr++){
	    getline(infile,line);
	    sscanf(line.c_str(), "%lf    %lf     %lf     %lf\n", &del1,  &del2,  &del3,  &mig);
	    MigMat[s][fl][er][zr][et] += fabs(mig);                                  //Divide by the size of the angle bin
	    sum += fabs(mig);
	  }
	
	
	if(sum > 1){
	  printf("%f   %d\n",sum, et);
	  printf("Error al calcular la migration matrix\n");
	  //exit(1);
	}

	Xsec[s][fl][et] *= sum;
      }
      infile.close();
    }
  
}



double xsec(int sg, int flv, double ene){
  
  int bet;

  if(ene <= 1)
    bet = int((ene - 0.1)/dltte1);
  else if(ene > 1 && ene <= 10)
    bet = int((ene - 1)/dltte2) + 18;
  else
    bet = int((ene - 10)/dltte3) + 63;
  
  return Xsec[sg][flv][bet];
  
  
}




//Tabulamos los bines que vamos a usar
const int NTzen = 100;                                             //Number of points in true zenit
static double CTmax = 1., CTmin = -1., Cint[NTzen];                //Limits of true zenit
static double dltzi = (CTmax - CTmin) / NTzen;
const int Nczen = 40;                                              //Number of bins in reconstructed zenit
static double czenmax = 1., czenmin = -1, czen[Nczen+1];           //Limits of reconstructed zenit
static double dltczen = (czenmax - czenmin)/Nczen;  
const  int Nene = 60;                                               //Number of bins in reconstructed energy
static double Emin = 0.1, Emax = 100., Ene[Nene+1];                 //Limits of reconstructed energy
static double dltle = (log10(Emax) - log10(Emin))/Nene;
const int Ninte = 1000;                                              //Number of points in true energy
static double Eint[Ninte], ETmax = 100., ETmin = 0.1;               //Limits of true energy
static double dltlei = (log10(ETmax) - log10(ETmin)) / Ninte;


void Tabbins(){

  //Bins of true zenit
  for(int tz=0; tz < NTzen; tz++)
    Cint[tz] = CTmin + dltzi * double(tz);

  //Bins of reconstructed zenit
  for(int c=0; c <= Nczen; c++)
    czen[c] = czenmin + double(c)*dltczen;

  //Bins of reconstructed energy
  for(int e=0; e<= Nene; e++)
    Ene[e] = exp10(log10(Emin) + dltle * double(e));

  //Bins of true energy
  for(int ie=0; ie< Ninte; ie++)
    Eint[ie] = exp10(log10(ETmin) + double(ie) * dltlei);

}




struct vecc{
  double x = 0.;
  double y = 0.;
  double z = 0.;
};

double BinZ(double zenr, double phir, double zent){


  vecc rot, rotdb;

  rot.x = sin(phir) * sin(zenr);
  rot.y = cos(phir) * sin(zenr);
  rot.z = cos(zenr);

  rotdb.x = rot.x;
  rotdb.y = rot.y * cos(zent) + rot.z * sin(zent);
  rotdb.z = -rot.y * sin(zent) + rot.z * cos(zent);

  return rotdb.z;
  
}






static float Eps[NUM_CH][NUM_FLV][NErecM][NEtrueM][NTzen][Nczen];
void TabEps(){

  //double zen[Nczen+1], czint[Nczen][180], dltzen[Nczen];
  

  int Nzenr = 200;
  double ZENrmax = M_PI, ZENrmin = 0., zenr, dltzenr = (ZENrmax - ZENrmin)/Nzenr;                   //Normalized to one
  double zent;
  
  int Nphr = 100;
  double PHrmax = 2*M_PI, PHrmin = 0., dltphr = (PHrmax - PHrmin)/Nphr;                           //Normalized to one
  double phir;

  for(int s=0; s < NUM_CH; s++)
    for(int fl=0; fl < NUM_FLV; fl++)
      for(int er=0; er < NErecM; er++)
	for(int et=0; et < NEtrueM; et++)
	  for(int tz=0; tz < NTzen; tz++)
	    for(int cz=0; cz < Nczen; cz++)
	      Eps[s][fl][er][et][tz][cz] = 0.;
	      
  //progress();
  for(int zr=0; zr < Nzenr; zr++){
    zenr = ZENrmin + double(zr) * dltzenr;
    
    //int bzr = int((zenr - TrecM[0])/(TrecM[1] - TrecM[0]));
    
    for(int phr=0; phr < Nphr; phr++){
      phir = PHrmin + dltphr * double(phr);
      
      for(int tz=0; tz < NTzen; tz++){
	zent = acos(Cint[tz]);
	
	int  bin = int((BinZ(zenr, phir, zent) - czenmin)/dltczen);

	for(int s=0; s < NUM_CH; s++)
	  for(int fl=0; fl < NUM_FLV; fl++)
	    for(int er=0; er < NErecM; er++)
	      for(int et=0; et < NEtrueM; et++)
		Eps[s][fl][er][et][tz][bin] += MigMat[s][fl][er][zr][et] * 1./double(Nphr);
      }
    }
    //progress(double(zr+1)/double(Nzenr));
  }  
  




//for(int er=0; er < NErecM; er++)
//  for(int et=0; et < NEtrueM; et++)
//    for(int tz=0; tz < NTzen; tz++)
//	for(int cz=0; cz < Nczen; cz++){
//	  printf("%f\n", Eps[1][1][er][et][tz][cz]);
//	}
  
  
}



double Effic(const int flv, const int sg, const int Er, const int zr, const int zt, const double Enet){

    
  int Et;     // Reduzco en una unidad para corregir el bin seleccionado
  int EEr = Er + 20;


  if(Enet <= 1)
    Et = int((Enet - 0.1)/dltte1);
  else if(Enet > 1 && Enet <= 10)
    Et = int((Enet - 1)/dltte2) + 18;
  else
    Et = int((Enet - 10)/dltte3) + 63;


  
  return Eps[sg][flv][EEr][Et][zt][zr];
  
}



const int Ngamma = 40;
static double PROBTBL[NTzen][Ninte][NUM_SG][NUM_FLV][NUM_FLV] = {};
static double Event[Ngamma][NUM_CH][NUM_FLV-1][NUM_FLV][Nene][Nczen];
int main(){

  ReadMig();        
  Tabbins();
  //TabEps();
  exit(1);
  
  printf("End of tabulating things\n");

  double Ntargt = NA * DUNEm;
  double Time = 10 * 365 * 24 * 3600;

  double flx[NUM_SG][NUM_FLA];
  Fluxes atm_flx("./Data/flx.kamioka.dat");
  

  printf("Tabulando probabilidades\n");
  //progress();


  ifstream infile;
  string line;

  int Ndens = 11;
  double del1, del2;
  double dens1, dens2, dens3;
  double dlt3d = 22., dlt2d = 10.;  
  string name1, name2;
  string str[3];
  string nexo[3];
  

  for(int ds2 = 0; ds2 < Ndens; ds2++)
    for(int ds3 = 0; ds3 < Ndens; ds3++){


  //for(int ds2 = 0; ds2 < 1; ds2++)
  //for(int ds3 = 0; ds3 < 1; ds3++){

      dens1 = 6.0;
      dens2 = 0.0 + dlt2d * double(ds2);
      dens3 = 0.0 + dlt3d * double(ds3);

      ostringstream strs[3];
      strs[0] << (dens1/10.);
      strs[1] << (dens2/10.);
      strs[2] << (dens3/10.);

      str[0] = strs[0].str();
      str[1] = strs[1].str();
      str[2] = strs[2].str();


      if(fmod(dens1,10.) == 0.)
	nexo[0] = ".0_";
      else
	nexo[0] = "_";

      if(fmod(dens2,10.) == 0.)
	nexo[1] = ".0_";
      else
	nexo[1] = "_";

      if(fmod(dens3,10.) == 0.)
	nexo[2] = ".0";
      else
	nexo[2] = "";

  

      name1 = "./DensityScan/OscProbs_" + str[0] + nexo[0] + str[1] + nexo[1] + str[2] + nexo[2] + ".dat";

      name2 = "./DensityScan/EventDens_" + str[0] + "_" + str[1] + "_" + str[2] + ".bin";


      infile.open(name1.c_str());
      printf("%s\n", name1.c_str());
      
      if(!infile.is_open()){
	printf("Error al abrir la densidad:  %s\n", name1.c_str());
	exit(1);
      }
            
      for(int ie=0; ie< Ninte; ie++)
	for(int iz=0; iz < NTzen; iz++){
	  
	  getline(infile,line);
	  sscanf(line.c_str(), "%lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf    %lf\n", &del1, &del2,
		 &PROBTBL[iz][ie][SG_n][FLV_elc][FLV_elc],
		 &PROBTBL[iz][ie][SG_n][FLV_mu][FLV_elc],
		 &PROBTBL[iz][ie][SG_n][FLV_mu][FLV_mu],
		 &PROBTBL[iz][ie][SG_n][FLV_elc][FLV_mu],
		 &PROBTBL[iz][ie][SG_a][FLV_elc][FLV_elc],
		 &PROBTBL[iz][ie][SG_a][FLV_mu][FLV_elc],
		 &PROBTBL[iz][ie][SG_a][FLV_mu][FLV_mu],
		 &PROBTBL[iz][ie][SG_a][FLV_elc][FLV_mu]);

	  PROBTBL[iz][ie][SG_a][FLV_elc][FLV_tau] = 1 - PROBTBL[iz][ie][SG_a][FLV_elc][FLV_elc] - PROBTBL[iz][ie][SG_a][FLV_elc][FLV_mu];
	  PROBTBL[iz][ie][SG_a][FLV_mu][FLV_tau] = 1 - PROBTBL[iz][ie][SG_a][FLV_mu][FLV_elc] - PROBTBL[iz][ie][SG_a][FLV_mu][FLV_mu];
	  PROBTBL[iz][ie][SG_n][FLV_elc][FLV_tau] = 1 - PROBTBL[iz][ie][SG_a][FLV_elc][FLV_elc] - PROBTBL[iz][ie][SG_a][FLV_elc][FLV_mu];
	  PROBTBL[iz][ie][SG_n][FLV_mu][FLV_tau] = 1 - PROBTBL[iz][ie][SG_a][FLV_mu][FLV_elc] - PROBTBL[iz][ie][SG_a][FLV_mu][FLV_mu];

	  
	}
      
      infile.close();
           
      //Systematics
      double Gmin = -0.1, Gmax = 0.1, dltg = (Gmax - Gmin)/double(Ngamma);
      double Tilt[Ngamma][Ninte], Gamma[Ngamma];
      
      for(int g=0; g < Ngamma; g++){
	Gamma[g] = Gmin + dltg * double(g);
	for(int ie=0; ie< Ninte; ie++)
	  Tilt[g][ie] = pow(Eint[ie],Gamma[g]);
      }
    
 
      FILE *outfile;
      

      outfile = fopen(name2.c_str(), "w");

      fwrite(Gamma, sizeof(double), Ngamma, outfile);

      fwrite(Ene, sizeof(double), Nene+1, outfile);

      fwrite(czen, sizeof(double), Nczen+1, outfile);
  

      Param prm;
      prm.the12 = 33.82*RADIAN;
      prm.the13 = 8.61*RADIAN;  
      prm.the23 = 49.6 * RADIAN;
      prm.dmqSL = 7.39*1e-5;
      prm.dmqAT = 2.525*1e-3;
      prm.dltCP = 0.;
  
  
      printf("Start computing the number of events\n");
      //progress();
      for(int m=0; m <= 0; m++){
	//dm = DmMin + double(m) * dltdm;
   
	for(int ss=0; ss <= 0; ss++){
	  //sinq12 = Sin12Min + double(ss) * dltsin12;
    
	  for(int s=0; s <= 0; s++){
	    //sinq  = SinMin + double(s) * dltsin;
	    //prm.the23 = asin(sqrt(sinq));
	
	
	    for(int mm=0; mm <= 0; mm++){
	      //dm = DmMin + double(m) * dltdm;
	      //dmAT = DmATMin + double(mm) * dltdmAT;

	      for(int g=0; g < Ngamma; g++)
		for(int s=0; s < NUM_CH; s++)
		  for(int i=FLV_elc; i < NUM_FLV-1; i++)
		    for(int f=FLV_elc; f < NUM_FLV; f++)
		      for(int e=0; e < Nene; e++)
			for(int z=0; z < Nczen; z++)
			  Event[g][s][i][f][e][z] = 0.;
	    
	    
	      for(int iz=0; iz < NTzen; iz++)
		for(int ie=0; ie< Ninte; ie++){
		
		  atm_flx.get(Eint[ie], Cint[iz], flx);
		
		  //Anti neutrinos
		  for(int g=0; g < Ngamma; g++)
		    for(int eb=0; eb < Nene; eb++)   
		      for(int zb=0; zb < Nczen; zb++){
		    
			Event[g][SG_a][FLV_elc][FLV_elc][eb][zb] += xsec(SG_a, FLV_elc, Eint[ie]) * flx[SG_a][FLV_elc] * PROBTBL[iz][ie][SG_a][FLV_elc][FLV_elc] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_elc, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];
			
			Event[g][SG_a][FLV_elc][FLV_mu][eb][zb] += xsec(SG_a, FLV_mu, Eint[ie]) * flx[SG_a][FLV_elc] * PROBTBL[iz][ie][SG_a][FLV_elc][FLV_mu] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_mu, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

			Event[g][SG_a][FLV_elc][FLV_tau][eb][zb] += xsec(SG_a, FLV_tau, Eint[ie]) * flx[SG_a][FLV_elc] * PROBTBL[iz][ie][SG_a][FLV_elc][FLV_tau] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_tau, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];


			Event[g][SG_a][FLV_mu][FLV_elc][eb][zb] += xsec(SG_a, FLV_elc, Eint[ie]) * flx[SG_a][FLV_mu] * PROBTBL[iz][ie][SG_a][FLV_mu][FLV_elc] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_elc, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

		      
			Event[g][SG_a][FLV_mu][FLV_mu][eb][zb] += xsec(SG_a, FLV_mu, Eint[ie]) * flx[SG_a][FLV_mu] * PROBTBL[iz][ie][SG_a][FLV_mu][FLV_mu] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_mu, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

			Event[g][SG_a][FLV_mu][FLV_tau][eb][zb] += xsec(SG_a, FLV_tau, Eint[ie]) * flx[SG_a][FLV_mu] * PROBTBL[iz][ie][SG_a][FLV_mu][FLV_tau] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_tau, SG_a, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

			
		      }
	    
	    
		  //Neutrinos
		  for(int g=0; g < Ngamma; g++)
		    for(int eb=0; eb < Nene; eb++)
		      for(int zb=0; zb < Nczen; zb++){
		    
		      
			Event[g][SG_n][FLV_elc][FLV_elc][eb][zb] += xsec(SG_n, FLV_elc, Eint[ie]) * flx[SG_n][FLV_elc] * PROBTBL[iz][ie][SG_n][FLV_elc][FLV_elc] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_elc, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];
			
			
			Event[g][SG_n][FLV_elc][FLV_mu][eb][zb] += xsec(SG_n, FLV_mu, Eint[ie]) * flx[SG_n][FLV_elc] * PROBTBL[iz][ie][SG_n][FLV_elc][FLV_mu] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_mu, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

			Event[g][SG_n][FLV_elc][FLV_tau][eb][zb] += xsec(SG_n, FLV_tau, Eint[ie]) * flx[SG_n][FLV_elc] * PROBTBL[iz][ie][SG_n][FLV_elc][FLV_tau] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_tau, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];


			Event[g][SG_n][FLV_mu][FLV_elc][eb][zb] += xsec(SG_n, FLV_elc, Eint[ie]) * flx[SG_n][FLV_mu] * PROBTBL[iz][ie][SG_n][FLV_mu][FLV_elc] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_elc, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

		      
			Event[g][SG_n][FLV_mu][FLV_mu][eb][zb] += xsec(SG_n, FLV_mu, Eint[ie]) * flx[SG_n][FLV_mu] * PROBTBL[iz][ie][SG_n][FLV_mu][FLV_mu] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_mu, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];

			Event[g][SG_n][FLV_mu][FLV_tau][eb][zb] += xsec(SG_n, FLV_tau, Eint[ie]) * flx[SG_n][FLV_mu] * PROBTBL[iz][ie][SG_n][FLV_mu][FLV_tau] * dltzi * M_LN10 * Eint[ie] *  dltlei * Effic(FLV_tau, SG_n, eb, zb, iz, Eint[ie]) * Tilt[g][ie];
		      
		      }
		}

	      
	      //Normalizo con el numero de target y el tiempo
	      for(int g=0; g < Ngamma; g++)
		for(int s=0; s < NUM_SG; s++)
		  for(int i=FLV_elc; i < NUM_FLV-1; i++)
		    for(int f=FLV_elc; f < NUM_FLV; f++)
		      for(int eb=0; eb < Nene; eb++)
			for(int zb=0; zb < Nczen; zb++)
			  Event[g][s][i][f][eb][zb] *= Ntargt * Time;


	      
	      /*
	      for(int g=0; g < Ngamma; g++)
		for(int eb=0; eb < Nene; eb++)
		  for(int zb=0; zb < Nczen; zb++)
		    outfile << std::setprecision(5)
			    << Event[g][SG_n][FLV_elc][FLV_elc][eb][zb] << "   " << Event[g][SG_n][FLV_mu][FLV_elc][eb][zb] << "   "
			    << Event[g][SG_n][FLV_elc][FLV_mu][eb][zb] << "   " << Event[g][SG_n][FLV_mu][FLV_mu][eb][zb] << "   "
			    << Event[g][SG_n][FLV_elc][FLV_tau][eb][zb] << "   " << Event[g][SG_n][FLV_mu][FLV_tau][eb][zb] << "   "
			    << Event[g][SG_a][FLV_elc][FLV_elc][eb][zb] << "   " << Event[g][SG_a][FLV_mu][FLV_elc][eb][zb] << "   "
			    << Event[g][SG_a][FLV_elc][FLV_mu][eb][zb] << "   " << Event[g][SG_a][FLV_mu][FLV_mu][eb][zb] << "   "
			    << Event[g][SG_a][FLV_elc][FLV_tau][eb][zb] << "   " << Event[g][SG_a][FLV_mu][FLV_tau][eb][zb] << "\n";
	      */

      	      
	      fwrite(Event, sizeof(double), Ngamma * NUM_SG * (NUM_FLV-1) * NUM_FLV * Nene * Nczen, outfile);	    
	    	    	      
	    }
	  }
	}
      }
      
      fclose(outfile);

    }
  
  
  return 0;

}
