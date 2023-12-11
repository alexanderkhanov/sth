// control parameters

double Afactor = 1; // =1 -> absorbed into f-matrix
double beta = 1;
int order = 1; // 1 = normal, 2 = inverted

// complex number calculations
struct ComplexNumber {
  ComplexNumber(double aA, double aB): a(aA), b(aB) {}
  ComplexNumber operator*(const ComplexNumber& other) {
    return ComplexNumber(a*other.a - b*other.b, a*other.b + b*other.a);
  }
  double a,b;
};

void eval(double* par, double* res)
{
  // input: par, model parameters
  // output: res, predictions

  // f-matrix
  int jpar = 0;
  double f12 = par[jpar++];
  double f13 = par[jpar++];
  double f23 = par[jpar++];

  MX_matrix f(0);
  f.Re(0,1) = f12;
  f.Re(0,2) = f13;
  f.Re(1,0) =-f12;
  f.Re(1,2) = f23;
  f.Re(2,0) =-f13;
  f.Re(2,1) =-f23;

#include "lepton_masses.h"

  MX_matrix m(0);
  m.Re(0,0) = me;
  m.Re(1,1) = mmu;
  m.Re(2,2) = mtau;
  const double mij[N*N] = {
    me,me,me,
    me,mmu,mmu,
    me,mmu,mtau
  };

  // y-matrix
  MX_matrix y;
  for (int i = 0; i<N; ++i) {
    for (int j = 0; j<N; ++j) {
      y.Re(i,j) = par[jpar++]*mij[i*N+j];
    }
  }
  for (int i = 0; i<N; ++i) {
    for (int j = 0; j<N; ++j) {
      y.Im(i,j) = par[jpar++]*mij[i*N+j];
    }
  }

  double fb = -1./sqrt(2)/sin(beta);

  MX_matrix N = f*m*m + m*m*f.T() + (f*m*y + y.T()*m*f.T())*fb;
  N *= Afactor;
  MX_matrix A = N*N.dagger();

  MX_matrix VN;
  MX_matrix DN = A.eigen(VN,order);

  double tan12_2 = (VN.Re(0,1)*VN.Re(0,1) + VN.Im(0,1)*VN.Im(0,1))
		  /(VN.Re(0,0)*VN.Re(0,0) + VN.Im(0,0)*VN.Im(0,0));
  double tan23_2 = (VN.Re(1,2)*VN.Re(1,2) + VN.Im(1,2)*VN.Im(1,2))
		  /(VN.Re(2,2)*VN.Re(2,2) + VN.Im(2,2)*VN.Im(2,2));
  double sin13_2 = (VN.Re(0,2)*VN.Re(0,2) + VN.Im(0,2)*VN.Im(0,2));

  double sin12_2 = 1/(1+1/tan12_2);
  double sin23_2 = 1/(1+1/tan23_2);

  ComplexNumber uprod =
    ComplexNumber(VN.Re(0,0),-VN.Im(0,0)) *
    ComplexNumber(VN.Re(0,2), VN.Im(0,2)) *
    ComplexNumber(VN.Re(2,0), VN.Im(2,0)) *
    ComplexNumber(VN.Re(2,2),-VN.Im(2,2));
  double s12 = sqrt(sin12_2);
  double s13 = sqrt(sin13_2);
  double s23 = sqrt(sin23_2);
  double c12 = sqrt(1-sin12_2);
  double c13 = sqrt(1-sin13_2);
  double c23 = sqrt(1-sin23_2);
  double cprod = c12*c13*c13*c23*s13;
  double sprod = s12*s23;
  double cosdel = (uprod.a/cprod + c12*c23*s13)/sprod;
  double sindel = (uprod.b/cprod)/sprod;
  double delta = -atan2(sindel,cosdel);
  if (delta<0) delta += 2*M_PI; // 0..2*pi

  double m12 = DN.Re(0,0);
  double m22 = DN.Re(1,1);
  double m32 = DN.Re(2,2);

  res[0] = m22 - m12;
  res[1] = m32 - m12;
  res[2] = sin12_2;
  res[3] = sin23_2;
  res[4] = sin13_2;
  res[5] = delta;
}
