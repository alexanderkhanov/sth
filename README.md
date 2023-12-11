# sth
Probabilistic appraisal of underdetermined models:

https://arxiv.org/abs/1612.07787

Compile and run (need ROOT installed):

gmake all

./d_ana

The model is described in eval.h, it takes npar inputs par
and calculates nres measurables res.

Starting from random initial set of parameters par, the model
tries to converge to a point where res are close to values defined
in resval.h .
