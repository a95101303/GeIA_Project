#include "math.h"

Double_t tanh_fit(Double_t *x, Double_t *par)
{
   double result;

   result = 0.5*par[0]*(tanh(par[3]*(x[0]-par[2]))) + par[1];

   return result;
}

