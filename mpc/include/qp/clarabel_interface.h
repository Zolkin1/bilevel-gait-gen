//
// Created by zolkin on 2/7/24.
//

#ifndef BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H
#define BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H

#include <fstream>

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

        ClarabelInterface(const ClarabelInterface& other);

        ClarabelInterface& operator=(const ClarabelInterface& other);

        void SetupQP(QPData& data, const vector_t& warm_start) override;

        vector_t Solve(const QPData& data) override;

        void ConfigureForInitialRun();

        void ConfigureForRealTime(double run_time_iters);

        vector_t GetDualSolution() const override;

        // TODO: Clean up these functions
        void SetupDerivativeCalcs(vector_t& dx, vector_t& dy_l, vector_t& dy_u,
                                  const QPData& data);

        void CalcDerivativeWrtMats(matrix_t& dP, matrix_t& dA, matrix_t& dG);

        void CalcDerivativeWrtVecs(vector_t& dq, vector_t& db, vector_t& dh);

        vector_t Computedx(const Eigen::SparseMatrix<double>& P, const vector_t& q,
                           const vector_t& xstar);

        vector_t Getdx() const;

        SolveQuality GetSolveQuality() const override;

        void SetVerbosity(bool verbose);

    protected:
    private:
        void ModifyConstraintsForSolver(QPData& data);

        std::unique_ptr<DefaultSolver<double>> solver_;
        DefaultSettings<double> settings_;
        std::vector<SupportedConeT<double>> cones_;

        SolveQuality solve_quality_;

        int num_equality_constraints_;
        int num_inequality_constraints_;

        vector_t dual_;
        vector_t primal_;
        vector_t slacks_;

        vector_t d_;
        vector_t dx_;
        vector_t lam_;
        vector_t nu_;

//        std::ofstream qp_file;
    };
}


#endif //BILEVEL_GAIT_GEN_CLARABEL_INTERFACE_H
