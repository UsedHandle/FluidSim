#include "fvsolver.h"

#include <array>
#include <vector>

std::pair<Eigen::SparseMatrix<float>, Eigen::VectorXf>
getDiffusionMatrices(const CartesianMesh& mesh) {
	const size_t size = mesh.faceIDs.size();

	Eigen::SparseMatrix<float> finalMat(size, size);
	Eigen::VectorXf finalKnowns = Eigen::VectorXf::Zero(size);

	typedef Eigen::Triplet<float> Tri;
	std::vector<Tri> triList;
	triList.reserve(size * 5);

	for (size_t i = 0; i < size; ++i) {
		std::array<float, 4> faceCoeff = { 1.f, 1.f, 1.f, 1.f };
		float A_p = 0.f;
		for (size_t facei = 0; facei < 4; ++facei) {
			const auto& face = mesh.faceIDs[i].faces[facei];
			if (std::get<bool>(face)) {
				triList.push_back(Tri(i, std::get<size_t>(face), faceCoeff[facei]));
				A_p -= faceCoeff[facei];
			}
			else {
				// dirchlet BC where boundary is 1
				if (i < 20) {
					finalKnowns(i) += 2.f * 1.f;
				}

				// dirchlet BC where boundary is zero
				A_p -= 2.0f;
			}
		}

		triList.push_back(Tri(i, i, A_p));
	}

	finalMat.setFromTriplets(triList.begin(), triList.end());
	return std::make_pair(finalMat, finalKnowns);
}