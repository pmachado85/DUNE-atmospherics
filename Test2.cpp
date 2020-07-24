//--------------------------------------------------------------------------------------------------------------------------
// Kevin J. Kelly: kkelly12@fnal.gov
//--------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <math.h>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <vector> 
using namespace std;

const complex<double> I(0,1);
const int nspecies = 3;
const double Ye = 0.5;

const double RE = 6371.0;

int gsl_linalg_exponential_ss(const gsl_matrix * A, gsl_matrix * eA, gsl_mode_t mode);

struct OscParams{
  double T12;
  double T13;
  double T23;
  double d;
  double Dm12;
  double Dm13;
};

struct MDProfile{
  std::vector<double> radii;
  std::vector<double> densities;
  int nlayers;
};

void DotMat3x3(complex<double> test[3][3], complex<double> First[3][3], complex<double> Last[3][3])
{

  for (int j=0; j <= 2; j++){
    for (int k=0; k <= 2; k++){
      test[j][k]=0;
      for (int m=0;m<=2;m++)
        test[j][k]+=First[j][m]*Last[m][k];
    }
  }
  return;
}

void UMatFn(complex<double> test2[3][3], double T12, double T13, double T23, double d, int MatType){

  complex<double> test[3][3] = {{0.}};
  test[0][0] = cos(T12)*cos(T13);
  test[0][1] = sin(T12)*cos(T13);
  test[0][2] = sin(T13)*exp(-I*d);

  test[1][0] = -sin(T12)*cos(T23)-cos(T12)*sin(T23)*sin(T13)*exp(I*d);
  test[1][1] = cos(T12)*cos(T23)-sin(T12)*sin(T23)*sin(T13)*exp(I*d);
  test[1][2] = sin(T23)*cos(T13);

  test[2][0] = -exp(I*d)*cos(T12)*cos(T23)*sin(T13)+sin(T12)*sin(T23);
  test[2][1] = -cos(T12)*sin(T23)-sin(T12)*cos(T23)*sin(T13)*exp(I*d);
  test[2][2] = cos(T13)*cos(T23);

  for (int h =0; h < 3; h++){
    for (int w = 0; w < 3; w++){
      if (MatType == 0)
	test2[h][w] = test[h][w];
      else if (MatType == 1)
	test2[h][w] = conj(test[w][h]);
      else if (MatType == 2)
	test2[h][w] = conj(test[h][w]);
      else if (MatType == 3)
	test2[h][w] = test[w][h];
    }
  }
  
  return;
}

void KMatFn(complex<double> test[3][3], double z, double L, double Dm12, double Dm13){

  for (int j = 0; j < 3; j++){
    for (int k = 0; k < 3; k++){
      test[j][k] = 0.0;}}
  test[1][1] = 2*1.267*Dm12*L/z;
  test[2][2] = 2*1.267*Dm13*L/z;

  return;
}

void VMassFn(complex<double> Result[3][3], double L, double dens, bool Mode, double T12, double T13, double T23, double d){
  const double A0 = 3.845e-04*Ye;
  complex<double> VFl[3][3] = {{0.}};
  if (Mode)
    VFl[0][0] = L*dens*A0;
  else
    VFl[0][0] = -L*dens*A0;

  complex<double> UMat[3][3] = {{0.}};
  UMatFn(UMat, T12, T13, T23, d, 0);
  complex<double> UDagMat[3][3] = {{0.}};
  UMatFn(UDagMat, T12, T13, T23, d, 1);
  complex<double> UStarMat[3][3] = {{0.}};
  UMatFn(UStarMat, T12, T13, T23, d, 2);
  complex<double> UTranMat[3][3] = {{0.}};
  UMatFn(UTranMat, T12, T13, T23, d, 3);
  
  complex<double> IntMat[3][3] = {{0.}};

  if (Mode){
    DotMat3x3(IntMat, VFl, UMat);
    DotMat3x3(Result, UDagMat, IntMat);
  }
  else{
    DotMat3x3(IntMat, VFl, UStarMat);
    DotMat3x3(Result, UTranMat, IntMat);
  }

  return;
}

void ExpHFn(complex<double> RedExpHMat[3][3], double z, double L, double dens, bool Mode, double T12, double T13, double T23, double d, double Dm12, double Dm13){
  complex<double> KinMat[3][3] = {{0.}};
  KMatFn(KinMat, z, L, Dm12, Dm13);
  complex<double> VMassMat[3][3] = {{0.}};
  VMassFn(VMassMat, L, dens, Mode, T12, T13, T23, d);

  complex<double> HMat[3][3] = {{0.}};
  for (int h = 0; h < 3; h++){
    for (int w = 0; w < 3; w++){
      HMat[h][w] = -I*(KinMat[h][w] + VMassMat[h][w]);
    }
  }

  double THam[36];

  gsl_matrix_view VTHam;
  gsl_matrix_view VExpTHam;
  double BigExpHMat[36];

  for (int j=0;j<=2;j++){
    for (int k=0; k<=2; k++){
      THam[12*j+2*k] = HMat[j][k].real();
      THam[12*j+2*k+1] = -HMat[j][k].imag();
      THam[12*j+2*k+6] = HMat[j][k].imag();
      THam[12*j+2*k+7] = HMat[j][k].real();
    }
  }

  VTHam = gsl_matrix_view_array(THam, 6, 6);
  VExpTHam = gsl_matrix_view_array(BigExpHMat, 6, 6);

  gsl_linalg_exponential_ss(&VTHam.matrix, &VExpTHam.matrix, 0);
  
  for (int j = 0; j <= 2; j++){
    for (int k = 0; k <= 2; k++)
      RedExpHMat[j][k] = BigExpHMat[12*j+2*k] - I*BigExpHMat[12*j+2*k+1];
  }

  return;
}

double *GetProb3Nu(std::vector<double> energy, const int length, double ctz, MDProfile MDP, bool NuNuBar, OscParams OP){
  //P is an array of the oscillation probability as a function of energy, to be returned.
  double *P = new double[(nspecies)*(nspecies)*(length+1)];
  double T12(OP.T12), T13(OP.T13), T23(OP.T23), Dm12(OP.Dm12), Dm13(OP.Dm13), d(OP.d);
  std::vector<double> radii = MDP.radii; std::vector<double> densities = MDP.densities; int nlayers = MDP.nlayers;
  double *RegionLens = new double[2*nlayers];
  double *RegionDens = new double[2*nlayers];

  double AtmosHeight = 15.0; //Assumed height of atmosphere

  if (ctz < -1.0 || ctz > 1.0){
    for (int k = 0; k < (nspecies)*(nspecies)*(length+1); k++){
      P[k] = 0.0;
    }
    return P;
  }
  else{
    if (ctz >= 0.0){ //Neutrinos that do not pass through the earth. Assume they travel through 15 km (check) of atmosphere with zero density.
      if (ctz == 0){
	      RegionLens[0] = pow( 2.0*RE*AtmosHeight - AtmosHeight*AtmosHeight, 0.5);
	      RegionDens[0] = 0.0;
      }
      else{
	      RegionLens[0] = RE*ctz*(pow(1.0 + (2.0*AtmosHeight/(RE*ctz*ctz)) + (AtmosHeight*AtmosHeight)/(RE*RE*ctz*ctz), 0.5) - 1.0);
	      RegionDens[0] = 0.0;
      }
      
      for (int kr = 1; kr < 2*nlayers; kr++){
	      RegionLens[kr] = 0.0; RegionDens[kr] = 0.0;
      }
    }
    else{
      for (int kr=0; kr < nlayers; kr++){
	      double RPathMin = RE*pow(1.0 - ctz*ctz, 0.5);
	      if (RPathMin <= radii[kr]){
	        RegionLens[kr] = -radii[kr]*ctz;
	        if (kr > 0){
	          RegionLens[kr-1] -= -radii[kr]*ctz;
	        }
	        RegionDens[kr] = densities[kr];
	      }
	      else{
	        RegionLens[kr] = 0.0; RegionDens[kr] = 0.0;
	      }
      }

      for (int kr = nlayers; kr < 2*nlayers; kr++){
	      RegionLens[kr] = RegionLens[(2*nlayers) - kr - 1];
	      RegionDens[kr] = RegionDens[(2*nlayers) - kr - 1];
      }
    }
    complex<double> UMat[3][3] = {{0.}};
    UMatFn(UMat, T12, T13, T23, d, 0);
    complex<double> UDagMat[3][3] = {{0.}};
    UMatFn(UDagMat, T12, T13, T23, d, 1);
    complex<double> UStarMat[3][3] = {{0.}};
    UMatFn(UStarMat, T12, T13, T23, d, 2);
    complex<double> UTranMat[3][3] = {{0.}};
    UMatFn(UTranMat, T12, T13, T23, d, 3);

    complex<double> NewExpHMat[3][3] = {{0.}};
    complex<double> OldExpHMat[3][3] = {{0.}};
    complex<double> MatProd[3][3] = {{0.}};
    
    complex<double> IntMat[3][3] = {{0.}};
    complex<double> ProbMatNu[3][3] = {{0.}};
    for (int m = 0; m <= length; m++){

      for (int h = 0; h < 3; h++){
	      for (int w = 0; w < 3; w++){
	        OldExpHMat[h][w] = 0.0;
	      }
	      OldExpHMat[h][h] = 1.0;
      }

      for (int k = 0; k < 2*nlayers; k++){
	      ExpHFn(NewExpHMat, energy[m], RegionLens[k], RegionDens[k], NuNuBar, T12, T13, T23, d, Dm12, Dm13);
	      DotMat3x3(MatProd, NewExpHMat, OldExpHMat);
	      for (int h = 0; h < 3; h++){
	        for (int w = 0; w < 3; w++){
	          OldExpHMat[h][w] = MatProd[h][w];
	        }
	      }
      }
      
      if (NuNuBar){
	      DotMat3x3(IntMat, OldExpHMat, UDagMat);
	      DotMat3x3(ProbMatNu, UMat, IntMat);
      }
      else if (!NuNuBar){
	      DotMat3x3(IntMat, OldExpHMat, UTranMat);
	      DotMat3x3(ProbMatNu, UStarMat, IntMat);
      }
      for (int ft = 0; ft < nspecies; ft++){
	      for (int it = 0; it < nspecies; it++){
	        P[ft*nspecies*(length+1) + it*(length+1) + m] = abs(ProbMatNu[ft][it])*abs(ProbMatNu[ft][it]);
	      }
      }
    }

    return P;
  }
}

double *AverageCTZ(std::vector<double> energy, const int length, double ctzmin, double ctzmax, const int NAvg, MDProfile MDP, bool NuNuBar, OscParams OP){
  double ctz(0.0);
  double AvgPF(0.0);
  double *AveragedProbs = new double[(nspecies)*(nspecies)*(length+1)];
  for (int m = 0; m < (nspecies)*(nspecies)*(length+1); m++)
    AveragedProbs[m] = 0.0;

  for (int k = 0; k <= NAvg; k++){
    ctz = ctzmin + (ctzmax-ctzmin)/(double(NAvg))*k;
    double *ThisProb;
    ThisProb = GetProb3Nu(energy, length, ctz, MDP, NuNuBar, OP);
    if (k == 0 || k == NAvg)
      AvgPF = 0.5/(double(NAvg));
    else
      AvgPF = 1.0/(double(NAvg));
    for (int m = 0; m < (nspecies)*(nspecies)*(length+1); m++){
      AveragedProbs[m] += AvgPF*ThisProb[m];
    }
    delete ThisProb;
  }
  return AveragedProbs;
}

double *AverageEnu(std::vector<double> EFine, const int lFine, const int NAvgE, double ctzmin, double ctzmax, const int NAvgCTZ, MDProfile MDP, bool NuNuBar, OscParams OP){
  int lAvg = lFine/NAvgE;
  double AvgPF(0.0);

  double *AvgProbs = new double[(nspecies)*(nspecies)*(lAvg)];
  for (int m = 0; m < (nspecies)*(nspecies)*(lAvg); m++)
    AvgProbs[m] = 0.0;
  
  double *FineProb;
  FineProb = AverageCTZ(EFine, lFine, ctzmin, ctzmax, NAvgCTZ, MDP, NuNuBar, OP);

  for (int k = 0; k < lAvg; k++){

    for (int efinei = k*NAvgE; efinei <= (k+1)*NAvgE; efinei++){
      if (efinei == k*NAvgE || efinei == (k+1)*NAvgE)
        AvgPF = 0.5/(double(NAvgE));
      else
        AvgPF = 1.0/(double(NAvgE));

      for (int ft = 0; ft < nspecies; ft++){
        for (int it = 0; it < nspecies; it++){
          AvgProbs[ft*nspecies*(lAvg) + it*(lAvg) + k] += AvgPF*FineProb[ft*nspecies*(lFine+1) + it*(lFine+1) + efinei];
        }
      }
    }
  }
  return AvgProbs;
}

//const int lenAvg = 1000;
//const int NCTZ = 100
int main(int argc, char* argv[])
{
    double T12_Phys(0.5905), T13_Phys(0.150265), T23_Phys(0.843529), Dm12_Phys(7.39e-05), Dm13_Phys(2.523e-03), d_Phys(3.875);
    OscParams OP;
    OP.T12 = T12_Phys; OP.T13 = T13_Phys; OP.T23 = T23_Phys; OP.d = d_Phys; OP.Dm12 = Dm12_Phys; OP.Dm13 = Dm13_Phys;

    double lEMin = log10(0.1);
    double lEMax = log10(100.0);
    const int lenAvg = 1000;
    const int NAvgE = 4;
    const int len = lenAvg*NAvgE;
    std::vector<double> energy(len+1,0.0);
    for (int k = 0; k <= len; k++){
      energy[k] = pow(10.0, lEMin + (lEMax-lEMin)/(double(len))*k);
    }
    std::vector<double> energyBinned(lenAvg+1, 0.0);
    for (int k = 0; k <= lenAvg; k++){
      energyBinned[k] = pow(10.0, lEMin + (lEMax-lEMin)/(double(lenAvg))*k);
    }
    
    const int nlayers(3);
    std::vector<double> radii(nlayers,0.0);
    radii[0] = 6731.0; radii[1] = 5700.0; radii[2] = 3400.0;
    std::vector<double> densities(nlayers,0.0);
    //densities[0] = 3.0; 
    densities[0] = atof(argv[1]);
    densities[1] = atof(argv[2]); 
    densities[2] = atof(argv[3]);
    MDProfile MDP;
    MDP.radii = radii; MDP.densities = densities; MDP.nlayers = nlayers;

    double ctzMin(-1.0), ctzMax(1.0);
    const int NCTZ = 100;
    const int NAvgCTZ(4);
    std::vector<double> ctzVec(NCTZ+1, 0.0);
    for (int k = 0; k <= NCTZ; k++)
      ctzVec[k] = ctzMin + (ctzMax-ctzMin)/NCTZ*k;

    double PMuE[NCTZ][lenAvg] = {{0.}};
    double PEE[NCTZ][lenAvg] = {{0.}};
    double PMuMu[NCTZ][lenAvg] = {{0.}};
    double PEMu[NCTZ][lenAvg] = {{0.}};

    double PMuEBar[NCTZ][lenAvg] = {{0.}};
    double PEEBar[NCTZ][lenAvg] = {{0.}};
    double PMuMuBar[NCTZ][lenAvg] = {{0.}};
    double PEMuBar[NCTZ][lenAvg] = {{0.}};
    
    double ctz(-1.0), ctz2(-0.99);
    for (int ctzi = 0; ctzi < NCTZ; ctzi++){
      //std::cout << ctzi << endl;
      double *ProbAvg, *ProbAvgBar;
      ctz = ctzVec[ctzi];
      ctz2 = ctzVec[ctzi+1];
      ProbAvg = AverageEnu(energy, len, NAvgE, ctz, ctz2, NAvgCTZ, MDP, true, OP);
      ProbAvgBar = AverageEnu(energy, len, NAvgE, ctz, ctz2, NAvgCTZ, MDP, false, OP);
      for (int kE = 0; kE < lenAvg; kE++){
        PMuE[ctzi][kE] = ProbAvg[0*nspecies*lenAvg + 1*lenAvg + kE];
        PEE[ctzi][kE] = ProbAvg[0*nspecies*lenAvg +0*lenAvg + kE];
        PMuMu[ctzi][kE] = ProbAvg[1*nspecies*lenAvg + 1*lenAvg + kE];
        PEMu[ctzi][kE] = ProbAvg[1*nspecies*lenAvg + 0*lenAvg + kE];

        PMuEBar[ctzi][kE] = ProbAvgBar[0*nspecies*lenAvg + 1*lenAvg + kE];
        PEEBar[ctzi][kE] = ProbAvgBar[0*nspecies*lenAvg +0*lenAvg + kE];
        PMuMuBar[ctzi][kE] = ProbAvgBar[1*nspecies*lenAvg + 1*lenAvg + kE];
        PEMuBar[ctzi][kE] = ProbAvgBar[1*nspecies*lenAvg + 0*lenAvg + kE];
      }
      delete ProbAvg;
      delete ProbAvgBar;
    }

    for (int ei = 0; ei < lenAvg; ei++){
      for (int ctzi = 0; ctzi < NCTZ; ctzi++){
        std::cout << 0.5*(energyBinned[ei] + energyBinned[ei+1]) << " " << 0.5*(ctzVec[ctzi] + ctzVec[ctzi+1]) << " " << PEE[ctzi][ei] << " " << PMuE[ctzi][ei] << " " << PMuMu[ctzi][ei] << " " << PEMu[ctzi][ei] << " " << PEEBar[ctzi][ei] << " " << PMuEBar[ctzi][ei] << " " << PMuMuBar[ctzi][ei] << " " << PEMuBar[ctzi][ei] << endl;
      }
    }

    return 0;
}