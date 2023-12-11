const int nres = 6;

double resval[nres] = {
  //7.56e-23,-2.49e-21, 0.3219, 0.5955, 0.021385, 1.455*M_PI // inverted
  7.56e-23, 2.55e-21, 0.3219, 0.4310, 0.021385, 1.455*M_PI // normal
};

double reserr[nres] = {
  //0.19e-23, 0.04e-21, 0.0170, 0.0175, 0.000835, 0.245*M_PI // inverted
  0.19e-23, 0.04e-21, 0.0170, 0.019, 0.000835, 0.245*M_PI // normal
};

int res_is_phi[nres] = {
  0,0,0,0,0,1
};
