#include <TRandom3.h>

#include <TMinuit.h>

#include <iostream>
#include <fstream>
using namespace std;

#include "matrix.h"
using namespace MX;

bool debug = false;

void print(MX_matrix& matrix)
{
  for (int i = 0; i<N; ++i) {
    for (int j = 0; j<N; ++j) {
      cout << TString::Format("%11.4g + i*%9.4g",matrix.Re(i,j),matrix.Im(i,j));
    }
    cout << endl;
  }
}

TRandom3 rnd;

// model description

#include "eval.h"

// operations

#include "minimizer.h"
#include "output.h"

int main()
{
  ofstream fout("par.txt");

  const int natt = 1000000000;
  int nall = 0;
  const int ntarget = 1000;
  int ngood = 0;
  for (int iatt = 0; iatt<natt; ++iatt) {

    cerr << iatt << "/" << ngood << '\r'; cerr.flush();

    const int npar = 21;
    double par[21] = {0}; // (f12,f13,f23), (Re yij, Im yij)
    double res[6];

    Afactor = 1;
    //beta = 1;
    //beta = 0.0897582;
    //beta = 1.56673;
    beta = rnd.Uniform(0.0897582,1.56673);
    //aij_attractor_value = 3;
    aij_attractor_value = rnd.Uniform(0.5,9.5);
    if (debug) cout << "Attempt " << iatt << ", beta = " << beta << ", aij attractor = " << aij_attractor_value << endl;
    order = 1; // normal

    par[0] = pow(10.,rnd.Uniform(-11.,-9.)); // f12
    par[1] = pow(10.,rnd.Uniform(-11.,-9.)); // f13
    par[2] = pow(10.,rnd.Uniform(-11.,-9.)); // f23

    eval(par,res);

    if (debug) {
      cout << "dm1 = " << res [0] << endl;
      cout << "dm2 = " << res [1] << endl;
      cout << "sin12_2 = " << res [2] << endl;
      cout << "sin23_2 = " << res [3] << endl;
      cout << "sin13_2 = " << res [4] << endl;
    }

    double err[npar];
    minimizer(npar, par, par, err);

    if (debug) print_min = true;
    double r = chi2(par);
    print_min = false;

    ++nall;
    if (r<10.) {
      output(fout,npar,par);
      if (++ngood>=ntarget) break;
    }
  }

  fout.close();

  if (debug) cout << "Successful runs: " << ngood << "/" << nall << endl;
}
