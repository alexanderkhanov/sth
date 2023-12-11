#include "resval.h"

bool print_min = false;

double aij_attractor_value = -1; // <0 -> no preferred aij value
double aij_attractor_weight = 0.5;

double chi2(double* par, int j0 = 0, int j1 = nres)
{
  double res[nres];
  eval(par, res);

  if (print_min) {
    cout << "-> chi2 evaluation" << endl;
    for (int ires = 0; ires<nres; ++ires) {
      double dres = res[ires] - resval[ires];
      if (res_is_phi[ires]) dres -= floor(dres/(2*M_PI))*2*M_PI;
      cout << ires << ": (" << res[ires] << " - " << resval[ires] << ") / " << reserr[ires] << " = " << dres/reserr[ires] << endl;
    }
  }

  double sum = 0;
  for (int ires = j0; ires<j1; ++ires) {
    double dres = res[ires] - resval[ires];
    if (res_is_phi[ires]) dres -= floor(dres/(2*M_PI))*2*M_PI;
    double dsum = dres/reserr[ires];
    sum += dsum*dsum;
  }
  if (print_min) cout << "chi2[" << j0 << "," << j1 << "] = " << sum << endl;

  if (aij_attractor_value>=0) {
    for (int ipar = 3; ipar<21; ++ipar) {
      sum += aij_attractor_weight/18 * pow(fabs(par[ipar])-aij_attractor_value,2);
    }
  }

  return sum;
}

int xlog0 = 0, xlog1 = 0;

void xlog(int& np, double* deriv, double& f, double par[], int flag)
{
  f = chi2(par,xlog0,xlog1);
}

void minimizer(int npar, double* parin, double* parout, double* errout)
{
  // input: parin = initial parameter values
  // output: parout = final parameter values, errout = their errors

  TMinuit* minuit = new TMinuit(npar);
  if (!print_min) minuit->SetPrintLevel(-1);
  minuit->SetFCN(xlog);

  xlog0 = 0;
  xlog1 = 6;

  for (int ipar = 0; ipar<npar; ++ipar) {
    TString sp = "par"; sp += ipar;
    double step = 0.1;
    if (parin[ipar]!=0) step = 0.1*fabs(parin[ipar]);
    double parmin = 0, parmax = 0;
    if (ipar>=3) { parmin = -10; parmax = 10; }
    minuit->DefineParameter(ipar,sp,parin[ipar],step,parmin,parmax);
  }

  //for (int ipar = 12; ipar<npar; ++ipar) { minuit->FixParameter(ipar); }

  double arglist[2] = {100000,0.1}; // max number of steps, tolerance
  int ierflag;
  minuit->mnexcm("MIGRAD",arglist,2,ierflag);

  for (int ipar = 0; ipar<npar; ++ipar) {
    minuit->GetParameter(ipar,parout[ipar],errout[ipar]);
  }

  delete minuit;
}
