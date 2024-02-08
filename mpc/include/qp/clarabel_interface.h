//
// Created by zolkin on 2/7/24.
//

#ifndef BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H
#define BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H

#include <Eigen/SparseCore>

#include <Clarabel>

#include "qp_interface.h"
#include "timer.h"

namespace mpc {
    template<typename T>
    void print_array(Eigen::Map<Eigen::VectorX<T>> &vec)
    {
        printf("[");
        size_t n = vec.size();
        for (size_t i = 0; i < n; i++)
        {
            printf("%.10f", vec.data()[i]);
            if (i < n - 1)
            {
                printf(", ");
            }
        }
        printf("]\n");
    }

    template<typename T>
    void print_solution(clarabel::DefaultSolution<T> &solution)
    {
        printf("Solution (x)\t = ");
        print_array(solution.x);
        printf("Multipliers (z)\t = ");
        print_array(solution.z);
        printf("Slacks (s)\t = ");
        print_array(solution.s);
    }
    using namespace clarabel;

    using sp_matrix_t = Eigen::SparseMatrix<double>;

    class ClarabelInterface : public QPInterface {

    public:
        ClarabelInterface(const QPData& data, bool verbose);

        void SetupQP(QPData& data, const vector_t& warm_start) override;

        vector_t Solve(const QPData& data) override;

        void ConfigureForInitialRun() const override;

        void ConfigureForRealTime(double run_time_iters) const override;

        vector_t GetDualSolution() const override;

        void SetupDerivativeCalcs(vector_t& dx, vector_t& dy_l, vector_t& dy_u);

        void CalcDerivativeWrtMats(Eigen::SparseMatrix<double>& dP, Eigen::SparseMatrix<double>& dA);

        void CalcDerivativeWrtVecs(vector_t& dq, vector_t dl, vector_t du);

        void Computedx(const Eigen::SparseMatrix<double>& P, const vector_t& q, const vector_t& xstar);

        vector_t Getdx() const;

        SolveQuality GetSolveQuality() const override;

    protected:
    private:
        void ModifyConstraintsForSolver(QPData& data);

        std::unique_ptr<DefaultSolver<double>> solver_;
        DefaultSettings<double> settings_;
        std::vector<SupportedConeT<double>> cones_;

        SolveQuality solve_quality_;

        vector_t dual_;
    };
}


#endif //BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H
