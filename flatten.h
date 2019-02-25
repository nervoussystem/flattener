#pragma once
#include <Eigen/Dense>
#include <vector>

void flatten(Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & V_flat, std::vector<int> & flattenIndices);
