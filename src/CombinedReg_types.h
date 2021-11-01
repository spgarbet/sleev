#ifndef TYPES_H
#define TYPES_H

#include <RcppEigen.h>

using namespace Eigen;

// Define our typedefs
typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;
typedef Map<MatrixXd> MapMatd;
typedef Map<MatrixXi> MapMati;

// Explained here: https://stackoverflow.com/a/42886223/11279006
#endif