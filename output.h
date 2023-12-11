void output(ofstream& fout, int npar, double* par)
{
  for (int ipar = 0; ipar<npar; ++ipar) {
    if (ipar) fout << " ";
    fout << setprecision(numeric_limits<double>::max_digits10) << par[ipar];
  }
  fout << endl;
}
