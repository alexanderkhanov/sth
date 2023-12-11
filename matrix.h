#ifndef _matrix_h_
#define _matrix_h_

#include <cmath>
#include <algorithm>

#include <TMatrixD.h>
#include <TMatrixDSymEigen.h>

inline double sq(double x) { return x*x; }

// 3x3 matrix manipuation

// matrices C=A+iB are stored in a 18-dim array
// a11,a12,a13,..,a33,b11,b12,b13,..,b33

namespace MX {

const int N = 3;
const int M = 18;

class MX_matrix
{
 public:

  // uninitialized matrix
  MX_matrix() {}

  // scalar matrix
  MX_matrix(double a) {
    for (int i = 0; i<M; ++i) data[i] = 0;
    data[0] = data[4] = data[8] = a;
  }

  // scalar matrix (a+b*i)
  MX_matrix(double a, double b) {
    for (int i = 0; i<M; ++i) data[i] = 0;
    data[0] = data[4] = data[8] = a;
    data[9] = data[13] = data[17] = b;
  }

  // matrix from array (can't make it a constructor: what is MX_matrix(0) ?)
  void setData(double* array) {
    for (int i = 0; i<M; ++i) data[i] = array[i];
  }

  double getData(int i) { return data[i]; }

  // define Hermitian matrix:
  //
  // | a[0]              0.5*(a[3]+i*a[4])   0.5*(a[5]+i*a[6]) |
  // | 0.5*(a[3]-i*a[4]) a[1]                0.5*(a[7]+i*a[8]) |
  // | 0.5*(a[5]-i*a[6]) 0.5*(a[7]-i*a[8])   a[2]              |
  //
  // produces random Gaussian matrix if inputs = Gaus(0,1)

  void setHermitian(double* array) {
    data[0] = array[0];
    data[1] = array[3]*0.5;
    data[2] = array[5]*0.5;
    data[3] = array[3]*0.5;
    data[4] = array[1];
    data[5] = array[7]*0.5;
    data[6] = array[5]*0.5;
    data[7] = array[7]*0.5;
    data[8] = array[2];
    data[9] = 0;
    data[10] = array[4]*0.5;
    data[11] = array[6]*0.5;
    data[12] =-array[4]*0.5;
    data[13] = 0;
    data[14] = array[8]*0.5;
    data[15] =-array[6]*0.5;
    data[16] =-array[8]*0.5;
    data[17] = 0;
  }

  // define complex symmetric matrix:
  //
  // | (a[0]+i*a[1])/2  (a[6]+i*a[7])/4    (a[8]+i*a[9])/4   |
  // | (a[6]+i*a[7])/4  (a[2]+i*a[3])/2    (a[10]+i*a[11])/4 |
  // | (a[8]+i*a[9])/4  (a[10]+i*a[11])/4  (a[4]+i*a[5])/2   |
  //
  // produces random c.s. matrix if inputs = Gaus(0,1)

  void setSymmetric(double* array) {
    data[0] = array[0];
    data[1] = array[6]*0.5;
    data[2] = array[8]*0.5;
    data[3] = array[6]*0.5;
    data[4] = array[2];
    data[5] = array[10]*0.5;
    data[6] = array[8]*0.5;
    data[7] = array[10]*0.5;
    data[8] = array[4];

    data[9] = array[1];
    data[10] = array[7]*0.5;
    data[11] = array[9]*0.5;
    data[12] =-array[7]*0.5;
    data[13] = array[3];
    data[14] = array[11]*0.5;
    data[15] =-array[9]*0.5;
    data[16] =-array[11]*0.5;
    data[17] = array[5];
  }

  // define unitary matrix
  // origin: J.Phys.A: Math.Gen.27 (1994) 4235
  // produces random unitary matrix if inputs = Uniform(0,1)

  void setUnitary(double* r)
  {
    // angles
    const double twopi = 2*M_PI;
    double alpha = r[0]*twopi;

    double psi12 = r[1]*twopi;
    double psi13 = r[2]*twopi;
    double psi23 = r[3]*twopi;

    double ksi12 = r[4]*twopi;
    double ksi13 = r[5]*twopi;

    double phi12 = asin(pow(r[6],0.5));
    double phi13 = asin(pow(r[7],0.5));
    double phi23 = asin(pow(r[8],0.25));

    double cospsi12 = cos(psi12), sinpsi12 = sin(psi12);
    double cospsi13 = cos(psi13), sinpsi13 = sin(psi13);
    double cospsi23 = cos(psi23), sinpsi23 = sin(psi23);

    double cosksi12 = cos(ksi12), sinksi12 = sin(ksi12);
    double cosksi13 = cos(ksi13), sinksi13 = sin(ksi13);
    double cosksi23 = 1,          sinksi23 = 0;

    double cosphi12 = cos(phi12), sinphi12 = sin(phi12);
    double cosphi13 = cos(phi13), sinphi13 = sin(phi13);
    double cosphi23 = cos(phi23), sinphi23 = sin(phi23);

    // elementary rotations
    MX_matrix E12(0); // i=0, j=1
    E12.Re(2,2) = 1;
    E12.Re(0,0) = cosphi12*cospsi12; E12.Im(0,0) = cosphi12*sinpsi12;
    E12.Re(0,1) = sinphi12*cosksi12; E12.Im(0,1) = sinphi12*sinksi12;
    E12.Re(1,0) =-sinphi12*cosksi12; E12.Im(1,0) = sinphi12*sinksi12;
    E12.Re(1,1) = cosphi12*cospsi12; E12.Im(1,1) =-cosphi12*sinpsi12;

    MX_matrix E13(0); // i=0, j=2
    E13.Re(1,1) = 1;
    E13.Re(0,0) = cosphi13*cospsi13; E13.Im(0,0) = cosphi13*sinpsi13;
    E13.Re(0,2) = sinphi13*cosksi13; E13.Im(0,2) = sinphi13*sinksi13;
    E13.Re(2,0) =-sinphi13*cosksi13; E13.Im(2,0) = sinphi13*sinksi13;
    E13.Re(2,2) = cosphi13*cospsi13; E13.Im(2,2) =-cosphi13*sinpsi13;

    MX_matrix E23(0); // i=1, j=2
    E23.Re(0,0) = 1;
    E23.Re(1,1) = cosphi23*cospsi23; E23.Im(1,1) = cosphi23*sinpsi23;
    E23.Re(1,2) = sinphi23*cosksi23; E23.Im(1,2) = sinphi23*sinksi23;
    E23.Re(2,1) =-sinphi23*cosksi23; E23.Im(2,1) = sinphi23*sinksi23;
    E23.Re(2,2) = cosphi23*cospsi23; E23.Im(2,2) =-cosphi23*sinpsi23;

    // composite rotations
    MX_matrix E1 = E12;
    MX_matrix E2 = E23*E13;

    // phase
    MX_matrix E0(cos(alpha),sin(alpha));

    // complete transformation
    *this = E0*E1*E2;
  }

  // assignment operator
  MX_matrix& operator=(const MX_matrix& other) {
    for (int i = 0; i<M; ++i) data[i] = other.data[i];
    return *this;
  }

  // access matrix elements
  double& Re(int i, int j) { return data[N*i+j]; }
  double& Im(int i, int j) { return data[N*N+N*i+j]; }

  // add two matrices
  MX_matrix operator+(const MX_matrix& other) {
    MX_matrix result;
    for (int i = 0; i<M; ++i) result.data[i] = data[i] + other.data[i];
    return result;
  }

  // multiply matrices elementwise
  MX_matrix& hadamard(const MX_matrix& other) {
    const int N2 = N*N;
    for (int i = 0; i<N2; ++i) {
      double a = data[i], b = data[i+N2], c = other.data[i], d = other.data[i+N2];
      data[i] = a*c - b*d;
      data[i+N2] = a*d + b*c;
    }
    return *this;
  }

  // multiply by a real number
  MX_matrix operator*(double a) {
    MX_matrix result;
    for (int i = 0; i<M; ++i) result.data[i] = data[i]*a;
    return result;
  }

  // multiply itself by a real number
  MX_matrix& operator*=(double a) {
    for (int i = 0; i<M; ++i) data[i] *= a;
    return *this;
  }

  // multiply two matrices
  MX_matrix operator*(const MX_matrix& other) {
    MX_matrix result;

    double* a = data;
    const double* b = other.data;
    double* c = result.data;

    c[0] = a[0]*b[0] - a[9]*b[9] + a[1]*b[3] - a[10]*b[12] + a[2]*b[6] - a[11]*b[15];
    c[9] = a[0]*b[9] + a[9]*b[0] + a[1]*b[12] + a[10]*b[3] + a[2]*b[15] + a[11]*b[6];
    c[1] = a[0]*b[1] - a[9]*b[10] + a[1]*b[4] - a[10]*b[13] + a[2]*b[7] - a[11]*b[16];
    c[10] = a[0]*b[10] + a[9]*b[1] + a[1]*b[13] + a[10]*b[4] + a[2]*b[16] + a[11]*b[7];
    c[2] = a[0]*b[2] - a[9]*b[11] + a[1]*b[5] - a[10]*b[14] + a[2]*b[8] - a[11]*b[17];
    c[11] = a[0]*b[11] + a[9]*b[2] + a[1]*b[14] + a[10]*b[5] + a[2]*b[17] + a[11]*b[8];
    c[3] = a[3]*b[0] - a[12]*b[9] + a[4]*b[3] - a[13]*b[12] + a[5]*b[6] - a[14]*b[15];
    c[12] = a[3]*b[9] + a[12]*b[0] + a[4]*b[12] + a[13]*b[3] + a[5]*b[15] + a[14]*b[6];
    c[4] = a[3]*b[1] - a[12]*b[10] + a[4]*b[4] - a[13]*b[13] + a[5]*b[7] - a[14]*b[16];
    c[13] = a[3]*b[10] + a[12]*b[1] + a[4]*b[13] + a[13]*b[4] + a[5]*b[16] + a[14]*b[7];
    c[5] = a[3]*b[2] - a[12]*b[11] + a[4]*b[5] - a[13]*b[14] + a[5]*b[8] - a[14]*b[17];
    c[14] = a[3]*b[11] + a[12]*b[2] + a[4]*b[14] + a[13]*b[5] + a[5]*b[17] + a[14]*b[8];
    c[6] = a[6]*b[0] - a[15]*b[9] + a[7]*b[3] - a[16]*b[12] + a[8]*b[6] - a[17]*b[15];
    c[15] = a[6]*b[9] + a[15]*b[0] + a[7]*b[12] + a[16]*b[3] + a[8]*b[15] + a[17]*b[6];
    c[7] = a[6]*b[1] - a[15]*b[10] + a[7]*b[4] - a[16]*b[13] + a[8]*b[7] - a[17]*b[16];
    c[16] = a[6]*b[10] + a[15]*b[1] + a[7]*b[13] + a[16]*b[4] + a[8]*b[16] + a[17]*b[7];
    c[8] = a[6]*b[2] - a[15]*b[11] + a[7]*b[5] - a[16]*b[14] + a[8]*b[8] - a[17]*b[17];
    c[17] = a[6]*b[11] + a[15]*b[2] + a[7]*b[14] + a[16]*b[5] + a[8]*b[17] + a[17]*b[8];

    return result;
  }

  // hermitian conjugate
  MX_matrix dagger() {
    MX_matrix result;

    double* a = data;
    double* b = result.data;

    b[0] = a[0]; b[9] = -a[9];
    b[1] = a[3]; b[10] = -a[12];
    b[2] = a[6]; b[11] = -a[15];
    b[3] = a[1]; b[12] = -a[10];
    b[4] = a[4]; b[13] = -a[13];
    b[5] = a[7]; b[14] = -a[16];
    b[6] = a[2]; b[15] = -a[11];
    b[7] = a[5]; b[16] = -a[14];
    b[8] = a[8]; b[17] = -a[17];

    return result;
  }

  // transpose
  MX_matrix T() {
    MX_matrix result;

    double* a = data;
    double* b = result.data;

    b[0] = a[0]; b[9] = a[9];
    b[1] = a[3]; b[10] = a[12];
    b[2] = a[6]; b[11] = a[15];
    b[3] = a[1]; b[12] = a[10];
    b[4] = a[4]; b[13] = a[13];
    b[5] = a[7]; b[14] = a[16];
    b[6] = a[2]; b[15] = a[11];
    b[7] = a[5]; b[16] = a[14];
    b[8] = a[8]; b[17] = a[17];

    return result;
  }

  // calculate eigen vectors and eigen values of a Hermitian matrix

  MX_matrix eigen() {
    MX_matrix result(0);
    eigenvalues(result);
    return result;
  }

  MX_matrix eigen(MX_matrix& result1, int order = 0) {
    MX_matrix result(0);
    eigenvalues(result);
    eigenvectors(result,result1);
    // result is m1<m2<m3
    if (order==1) { // swap 1<->3
      std::swap(result.data[0],result.data[8]);

      std::swap(result1.data[0],result1.data[2]);
      std::swap(result1.data[3],result1.data[5]);
      std::swap(result1.data[6],result1.data[8]);
      std::swap(result1.data[9],result1.data[11]);
      std::swap(result1.data[12],result1.data[14]);
      std::swap(result1.data[15],result1.data[17]);
    } else if (order==2) { // swap 1<->2
      std::swap(result.data[0],result.data[4]);

      std::swap(result1.data[0],result1.data[1]);
      std::swap(result1.data[3],result1.data[4]);
      std::swap(result1.data[6],result1.data[7]);
      std::swap(result1.data[9],result1.data[10]);
      std::swap(result1.data[12],result1.data[13]);
      std::swap(result1.data[15],result1.data[16]);
    }
    return result;
  }

  void eigenvalues(MX_matrix& result) {

    // the code assumes (but does not verify) that the matrix is hermitian

    // origin: arXiv:physics/0610206

    double* a = data;

    double c2 = -a[0] - a[4] - a[8];

    double a12a2 = sq(a[1]) + sq(a[10]);
    double a13a2 = sq(a[2]) + sq(a[11]);
    double a23a2 = sq(a[5]) + sq(a[14]);

    double a11a22 = a[0]*a[4];
    double a11a33 = a[0]*a[8];
    double a22a33 = a[4]*a[8];

    double c1 = a11a22 + a11a33 + a22a33 - a12a2 - a13a2 - a23a2;

    double a11a22a33 = a11a22*a[8];

    double a13sa12a23re = (a[2]*a[1] + a[11]*a[10])*a[5] + (a[11]*a[1] - a[2]*a[10])*a[14];

    double c0 = a[0]*a23a2 + a[4]*a13a2 + a[8]*a12a2 - a11a22a33 - 2*a13sa12a23re;

    // eigenvalues are roots of x^3 + c2*x^2 + c1*x + c0 = 0

    double roots[N];

    if (c0==0) {

      roots[2] = 0;

    } else {

      // If the roots are largely different, Cardano's formula
      // will only give numerically stable results for the largest root,
      // but we need the smallest root. That one can be found by solving
      // the equation for 1/x:
      // (1/x)^3 + c1/c0*(1/x)^2 + c2/c0*(1/x) + 1/c0 = 0

      double d0 = 1/c0;
      double d1 = c2*d0;
      double d2 = c1*d0;

      double p = d2*d2 - 3*d1;
      double q = -13.5*d0 - d2*d2*d2 + 4.5*d2*d1;
      double phi = atan(sqrt(27*(0.25*d1*d1*(p-d1) + d0*(q+6.75*d0)))/q)/3;
      double xinv = (2*cos(phi)*sqrt(p) - d2)/3;

      roots[2] = 1/xinv;
    }

    // if x1 is a root of x^3+c2*x^2+c1*x+c0=0,
    // then the two remaining roots are the roots of x^2+b*x+c=0,
    // where b=c2+x1, c=b*x1+c1

    double b = c2+roots[2], c = b*roots[2]+c1;
    double d = sqrt(b*b-4*c);
    if (b>0) { roots[1] = (-b-d)/2; roots[0] = c/roots[1]; }
    else     { roots[0] = (-b+d)/2; roots[1] = c/roots[0]; }

    result.data[0] = roots[0];
    result.data[4] = roots[1];
    result.data[8] = roots[2];
  }

  void eigenvectors(MX_matrix& result, MX_matrix& result1) {

    // calculate eigenvectors: v_i = [(A1 - l_i*e1) x (A2 - l_i*e2)]^*

    double roots[N];
    roots[0] = result.data[0];
    roots[1] = result.data[4];
    roots[2] = result.data[8];

    double* a = data;

    //complex<double> c11(a[0],a[9]);
    //complex<double> c12(a[1],a[10]);
    //complex<double> c21(a[3],a[12]);
    //complex<double> c22(a[4],a[13]);
    //complex<double> c31(a[6],a[15]);
    //complex<double> c32(a[7],a[16]);
    //double norm_a1 = norm(c11) + norm(c21) + norm(c31);
    //double norm_a2 = norm(c12) + norm(c22) + norm(c32);

    double norm_c11 = sq(a[0]) + sq(a[9]);
    double norm_c21 = sq(a[3]) + sq(a[12]);
    double norm_c31 = sq(a[6]) + sq(a[15]);

    double norm_c12 = sq(a[1]) + sq(a[10]);
    double norm_c22 = sq(a[4]) + sq(a[13]);
    double norm_c32 = sq(a[7]) + sq(a[16]);

    double norm_a1 = norm_c11 + norm_c21 + norm_c31;
    double norm_a2 = norm_c12 + norm_c22 + norm_c32;

    double c21c32re = a[3]*a[7] - a[12]*a[16];
    double c21c32im = a[3]*a[16] + a[12]*a[7];

    double c31c12re = a[6]*a[1] - a[15]*a[10];
    double c31c12im = a[6]*a[10] + a[15]*a[1];

    double c12c21re = a[1]*a[3] - a[10]*a[12];
    double c12c21im = a[1]*a[12] + a[10]*a[3];

    for (int j = 0; j<N; ++j) {

      //complex<double> v_1 = c21*c32 - c31*(c22 - roots[j]);
      //complex<double> v_2 = c31*c12 - c32*(c11 - roots[j]);
      //complex<double> v_3 = (c11 - roots[j])*(c22 - roots[j]) - c12*c21;

      double c11r = a[0] - roots[j];
      double c22r = a[4] - roots[j];

      double vxre = c21c32re - a[6]*c22r;
      double vxim = c21c32im - a[15]*c22r;
      double vyre = c31c12re - a[7]*c11r;
      double vyim = c31c12im - a[16]*c11r;
      double vzre = c11r*c22r - c12c21re;
      double vzim =           - c12c21im;

      //double norm_v = norm(v_1) + norm(v_2) + norm(v_3);
      double norm_v = vxre*vxre + vxim*vxim + vyre*vyre + vyim*vyim + vzre*vzre + vzim*vzim;

      // check if (A1 - l_i*e1) and (A2 - l_i*e2) are parallel
      const double eps = 64*std::numeric_limits<double>::epsilon();
      if (sqrt(norm_v/(norm_a1*norm_a2))<eps) {
	double norm1 = norm_c12;
	double norm2 = c22r*c22r;
	double norm3 = norm_c32;
	double mure, muim;
	if (norm1>norm2 && norm1>norm3) {
	  //mu = (c11 - roots[j])/c12;
	  mure = c11r*a[1]/norm_c12;
	  muim = -c11r*a[10]/norm_c12;
	} else if (norm2>norm1 && norm2>norm3) {
	  //mu = c21/(c22 - roots[j]);
	  mure = a[3]/c22r;
	  muim = a[12]/c22r;
	} else {
	  //mu = c31/c32;
	  mure = (a[6]*a[7] + a[15]*a[16])/norm_c32;
	  muim = (a[15]*a[7] - a[6]*a[16])/norm_c32;
	}
	//v_1 = 1;
	//v_2 = -mu;
	//v_3 = 0;
	//norm_v = 1 + norm(mu);
	vxre = 1; vxim = 0;
	vyre = -mure; vyim = -muim;
	vzre = vzim = 0;
	norm_v = 1 + mure*mure + muim*muim;
      }
      double vabs = sqrt(norm_v);
      result1.data[j]         = vxre/vabs;
      result1.data[N*N+j]     =-vxim/vabs;
      result1.data[N+j]       = vyre/vabs;
      result1.data[N*N+N+j]   =-vyim/vabs;
      result1.data[N*2+j]     = vzre/vabs;
      result1.data[N*N+N*2+j] =-vzim/vabs;
    }
  }

  // find inverse matrix (use root)
  MX_matrix inverse() {
    MX_matrix result(0); // set all elements to 0

    // method: (A+iB)*(C+iD) = I is equivalent to
    //
    //  | A -B | | C -D |   | I  0 |
    //  |      |*|      | = |      |
    //  | B  A | | D  C |   | 0  I |

    // construct the 2Nx2N matrix
    TMatrixD S(2*N,2*N);
    for (int i = 0; i<N; ++i) {
      for (int j = 0; j<N; ++j) {
	S(i,j) = data[N*i+j];
	S(i+N,j+N) = data[N*i+j];
	S(i,j+N) = -data[N*N+N*i+j];
	S(i+N,j) = data[N*N+N*i+j];
      }
    }

    S.InvertFast();

    // fill out the result
    for (int i = 0; i<N; ++i) {
      for (int j = 0; j<N; ++j) {
	result.Re(i,j) = S(i,j);
	result.Im(i,j) = S(i+N,j);
      }
    }

    return result;
  }

  // distance between two matrices
  double squared_distance_to(const MX_matrix& other) const {
    double d2 = 0;
    for (int i = 0; i<M; ++i) {
      double di = data[i] - other.data[i];
      d2 += di*di;
    }
    return d2;
  }

 private:
  double data[M];
};

} // namespace

#endif // _matrix_h_
