#pragma once

#include <Eigen/Sparse>

#include "mesh.h"

// Ax + c = b
// where A is the matrix and c is the
// product of the coefficients with known values of x
// returns a pair of A and c
std::pair<Eigen::SparseMatrix<float>, Eigen::VectorXf>
getDiffusionMatrices(const CartesianMesh& mesh);