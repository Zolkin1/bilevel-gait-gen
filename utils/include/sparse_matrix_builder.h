//
// Created by zolkin on 1/3/24.
//

#ifndef BILEVEL_GAIT_GEN_SPARSE_MATRIX_BUILDER_H
#define BILEVEL_GAIT_GEN_SPARSE_MATRIX_BUILDER_H

#include <Eigen/SparseCore>

namespace utils {
    class SparseMatrixBuilder {
        using triplet_t = Eigen::Triplet<double>;
        using vector_t = Eigen::VectorXd;
        using matrix_t = Eigen::MatrixXd;

    public:
        SparseMatrixBuilder();

        void Reserve(int num_nz);

        void SetDiagonalMatrix(double val, int row_start, int col_start, int num_diag);

        void SetMatrix(const matrix_t& M, int row_start, int col_start);

        void SetVectorDiagonally(const vector_t& vec, int row_start, int col_start);

        const std::vector<triplet_t>& GetTriplet() const;

    protected:
    private:
        std::vector<triplet_t> triplet_;
    };
}


#endif //BILEVEL_GAIT_GEN_SPARSE_MATRIX_BUILDER_H
