//
// Copyright (c) 2023 Zachary Olkin. All rights reserved.
//

#ifndef BILEVEL_GAIT_GEN_SPLINE_H
#define BILEVEL_GAIT_GEN_SPLINE_H

#include <vector>
#include <array>

namespace mpc {
    /**
     * Spline class.
     * Each continuous parameterization is given by a spline, which is composed of a set of polynomials.
     * Each polynomial has a start and an end time and a set of defining variables.
     */
    class Spline {
    public:
        Spline(int num_polys, const std::vector<double>& times, bool start_on_poly);

        double ValueAt(double time) const;

        double DerivWrtTime(double time) const;

        double DerivWrtVars(double time) const;

        void SetPolyVars(double time);

        void SetPolyVars(int poly_num);

        void SetPolyVar(int poly_num, int var_idx);
    protected:
    private:
        static int constexpr POLY_ORDER = 4;    // actually a cubic

        // Time should be [0, DeltaT]
        static double EvalPoly(const std::array<double, POLY_ORDER>& poly_vals, double time, double DeltaT);

        static std::array<double, POLY_ORDER> constexpr ZERO_POLY = {0, 0, 0};

        std::vector<std::array<double, POLY_ORDER>> poly_vars_;   // All the variables that define each polynomial
        // The poly vals are: x0, x0dot, x1, x1dot. DeltaT is determined by the switching times

        std::vector<double> poly_times_;    // The start and end of each polynomial
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SPLINE_H
