//
// Created by zolkin on 2/20/24.
//
#ifndef BILEVEL_GAIT_GEN_QPPARTIALS_H
#define BILEVEL_GAIT_GEN_QPPARTIALS_H

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace mpc {
    using matrix_t = Eigen::MatrixXd;
    using vector_t = Eigen::VectorXd;
    using sp_matrix_t = Eigen::SparseMatrix<double>;

    struct QPPartials {
        sp_matrix_t dA;
        sp_matrix_t dP;
        sp_matrix_t dG;
        vector_t dl;
        vector_t du;
        vector_t dq;
        vector_t db;
        vector_t dh;

        void SetZero() {
            dP.setZero();
            dA.setZero();
            dG.setZero();
            dq.setZero();
            dl.setZero();
            du.setZero();
            db.setZero();
            dh.setZero();
        }
    };
} // mpc

#endif //BILEVEL_GAIT_GEN_QPPARTIALS_H
