#include <stdio.h>
#include <stdlib.h>
#include <Rcpp.h> 

using namespace Rcpp;

//RcppExport SEXP getProposal (SEXP currentState, SEXP psi, SEXP min, SEXP max);

// The objects declared as arguments do not need to be declared again below. Warning: Can not use 'inf' to pass an infinite number from R to C. C only works with finite numbers, better to the user to put an high number entry, like 100 or even 1000.
RcppExport SEXP getProposal (SEXP _currentState, SEXP _psi, SEXP _min, SEXP _max) 
{
  double currentState = as<double>(_currentState)
  double psi = as<double>(_psi)
  double max = as<double>(_min)
  double max = as<double>(_max)
  double random; // This object stores the random value [runif(1)].
  double prop; // The proposal.

  // The next three lines are generating a random number from 0 to 1 (inclusive). See http://pubs.opengroup.org/onlinepubs/9699919799/functions/drand48.html and http://www.thinkage.ca/english/gcos/expl/c/lib/erand4.html.
  unsigned short work [3];
  double erand48();
  random = erand48(work);

  prop = currentState + ( (random - 0.5) * psi ); // The very same equation of the original function.

  if (prop < min)
    prop = min + (min - prop);
  if (prop > max)
    prop = max - (prop - max);

  return (prop);
}
