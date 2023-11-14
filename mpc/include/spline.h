//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SPLINE_H
#define BILEVEL_GAIT_GEN_SPLINE_H

#include <vector>
#include <array>

#include <Eigen/Core>

namespace mpc {
    using vector_t = Eigen::VectorXd;

    /**
     * Spline class.
     * Each continuous parameterization is given by a spline, which is composed of a set of polynomials.
     * Each polynomial has a start and an end time and a set of defining variables.
     */
    class Spline {
    public:
        static int constexpr POLY_ORDER = 4;    // actually a cubic

        Spline(int num_polys, const std::vector<double>& times, bool start_on_poly);

        double ValueAt(double time) const;

        double DerivWrtTime(double time) const;

        double DerivWrtVars(double time) const;

        void SetPolyVars(double time);

        /**
         * Sets the polynomial variables for a given polynomial in the spline.
         * @param poly_num Index of the polynomial to modify
         * @param vars new polynomial variables
         */
        void SetPolyVars(int poly_num, const std::array<double, POLY_ORDER>& vars);

        void SetPolyVar(int poly_num, int var_idx);

        /**
         *
         * @return The total number of scalars used to describe the spline.
         */
        int GetTotalPolyVars() const;

        /**
         *
         * @return The number of polynomial segments used to describe the spline.
         */
        int GetTotalPoly() const;

        /**
         *
         * @param time the time at which to evaluate the spline for its variables
         * @return The scalar variables used to describe this portion of the spline.
         */
        Eigen::Vector<double, POLY_ORDER> GetPolyVars(double time) const;

        /**
         *
         * @return a vector of all the scalar variables used to describe the spline.
         */
        vector_t GetAllPolyVars() const;

        /**
         * Gets the coefficients for the polynomial variables.
         * @param time to query the spline at
         * @return Vector of coefficients of the polynomial variables
         */
        Eigen::Vector<double, POLY_ORDER> GetPolyVarsLin(double time) const;

        /**
         * Given a time, gets the index for the associated polynomial.
         * @param time to query the spline at
         * @return index of the associated polynomial
         */
        int GetPolyIdx(double time) const;
    protected:
    private:

        // Time should be [0, DeltaT]
        static double EvalPoly(const std::array<double, POLY_ORDER>& poly_vals, double time, double DeltaT);

        static std::array<double, POLY_ORDER> constexpr ZERO_POLY = {0, 0, 0, 0};

        std::vector<std::array<double, POLY_ORDER>> poly_vars_;   // All the variables that define each polynomial
        // The poly vals are: x0, x0dot, x1, x1dot. DeltaT is determined by the switching times

        std::vector<double> poly_times_;    // The start and end of each polynomial

        int total_poly_;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SPLINE_H
