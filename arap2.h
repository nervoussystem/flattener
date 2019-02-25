#pragma once
#include "igl/arap.h"
namespace igl {
	// Compute necessary information to start using an ARAP deformation
	//
	// Inputs:
	//   V  #V by dim list of mesh positions
	//   F  #F by simplex-size list of triangle|tet indices into V
	//   dim  dimension being used at solve time. For deformation usually dim =
	//     V.cols(), for surface parameterization V.cols() = 3 and dim = 2
	//   b  #b list of "boundary" fixed vertex indices into V
	// Outputs:
	//   data  struct containing necessary precomputation
	template <
		typename DerivedV,
		typename DerivedF,
		typename Derivedb>
		bool arap_precomputation2(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			const int dim,
			Eigen::PlainObjectBase<Derivedb> * b,
			ARAPData & data);
	// Inputs:
	//   bc  #b by dim list of boundary conditions
	//   data  struct containing necessary precomputation and parameters
	//   U  #V by dim initial guess
	template <
		typename Derivedbc,
		typename DerivedU>
		bool arap_solve2(
			const Eigen::PlainObjectBase<Derivedbc> * bc,
			ARAPData & data,
			Eigen::PlainObjectBase<DerivedU> & U);
};

#include "arap2.cpp"