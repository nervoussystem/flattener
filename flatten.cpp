#include "flatten.h"
#include "arap2.h"

void flatten(Eigen::MatrixXd & V, Eigen::MatrixXi & T, Eigen::MatrixXd & V_flat, std::vector<int> & flattenIndices) {
	assert(flattenIndices.size() > 0);
	igl::ARAPData arap_data;
	//arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;

	arap_data.max_iter = 100;
	Eigen::VectorXi b[3];// (flattenIndices.size());
	//Eigen::MatrixXd bc(flattenIndices.size(),3);
	Eigen::VectorXd bc[3];
	b[0].resize(1);
	b[0](0) = flattenIndices[0];
	b[1].resize(1);
	b[1](0) = flattenIndices[0];
	b[2].resize(flattenIndices.size());
	bc[0].resize(1);
	bc[0](0) = V(flattenIndices[0], 0);
	bc[1].resize(1);
	bc[1](0) = V(flattenIndices[0], 1);
	bc[2].resize(flattenIndices.size());
	for (int i = 0; i < flattenIndices.size(); ++i) {
		//b(index, 2) = 1;
		b[2](i) = flattenIndices[i];
		bc[2](i) = 0;
		//bc.row(i) = V.row(flattenIndices[i]);
		//bc(i, 2) = 0;
	}
	int index0 = flattenIndices[0];
	//b(index0, 0) = 1;
	//b(index0, 1) = 1;

	std::cout << "precompute" << std::endl;
	igl::arap_precomputation2(V, T, V.cols(), &b[0], arap_data);

	//V_flat = V;

	igl::arap_solve2(&bc[0], arap_data, V_flat);
}
