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

    // 11/15/23: Starting a spline rework
    /*
     * Note that the spline should always be continuous. In my case, it will always be a set of nth order polynomials
     * followed by a constant and so forth.
     *
     * The constant is parameterized by a single variable.
     *
     * The polynomials are each, say, 4 variables. But due to continuity, for each pair, there will only be 6 free
     * variables: start:[x0, x0dot] connecting:[x1, x1dot] end:[x2, x2dot] So we have removed two variables.
     *
     * Let's examine a stance phase followed by a swing phase. Stance has one: x0. Swing (for 2 polynomials) has 6
     * nominally, but we have an initial constraint, so it actually has 4 variables. So there are a total of 5.
     *
     * Now consider stance, swing, stance. stance (1) + swing (2) + stance (1). Stance has been reduced to 2 (only the
     * middle position and derivative). If the swing had 3 polynomials then we get an extra 2 variables for a total of 8.
     *
     * If we end on a swing then that last end gets a position and derivative.
     *
     * So each spline will be made of sections of constants and polynomials.
     *
     * Note: assumes that the polynomials in between the constants are made up of at least 2 polynomials
     */

    /**
     * Spline class.
     * Each continuous parameterization is given by a spline, which is composed of a set of polynomials.
     * Each polynomial has a start and an end time and a set of defining variables.
     */
    class Spline {
    public:
        static int constexpr POLY_ORDER = 4;    // actually a cubic

        /**
         *
         * Note that the times is the times between constant and non-constant polynomials.
         * Need to convert this to times for each individual polynomial. Will do this by
         * evenly splitting the space.
         *
         * @param num_polys number of polynomials between each constant polynomial
         * @param times times to switch between constant and non-constant polynomials
         * @param start_on_poly false if we start on a constant polynomial
         */
        Spline(int num_polys, const std::vector<double>& times, bool start_on_poly);

        double ValueAt(double time) const;

        double DerivWrtTime(double time) const;

        double DerivWrtVars(double time) const;

        /**
         * Sets the variables at a given index.
         * @param poly_time Index of the polynomial switching time
         * @param vars new polynomial variables
         */
        void SetPolyVars(int poly_time, const std::vector<double>& vars);

        void SetAllSplineVars(const std::vector<std::vector<double>>& vars);

        /**
         * Note that this function returns the minimum number of variables used to describe
         * the spline. i.e. all the constant polynomial sections are collapsed to one variable.
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
        std::array<double, Spline::POLY_ORDER> GetPolyVars(double time) const;

        /**
         * Note that this will return one number per constant, as opposed to the internal representation using 2.
         * @return a vector of all the scalar variables used to describe the spline.
         */
        vector_t GetAllPolyVars() const;

        /**
         * Gets the coefficients for the polynomial variables.
         * @param time to query the spline at
         * @return Vector of coefficients of the polynomial variables
         */
        vector_t GetPolyVarsLin(double time) const;

        /**
         * Given a time, gets the index for the associated polynomial.
         * @param time to query the spline at
         * @return index of the associated polynomial
         */
        int GetPolyIdx(double time) const;

        std::pair<int, int> GetVarsIndexEnd(double time) const;

        int GetNumPolyVars(int idx) const;

        int GetNumPolyTimes() const;

        void UpdatePolyVar(int idx, const std::vector<double>& vars);

        const std::vector<std::vector<double>>& GetPolyVars() const;

        const std::vector<double>& GetPolyTimes() const;

        int GetNumConstant() const;

        /**
         * Sets all position variables to positions without touching the derivatives.
         * @param position
         */
        void SetAllPositions(double position);

        int GetNumNonConstantValParams() const;

        double GetEndTime() const;

        // Time gives the dt that this segment takes up, not the absolute time
        // Note that this time is the time in between constant and non-contant,
        // NOT the time a polynomial takes up.
        void AddPoly(double time);

        void RemoveUnused(double time);

    protected:
    private:
        bool IsConstantPoly(int idx) const;

        // Time should be [0, DeltaT]
        static double EvalPoly(const std::array<double, POLY_ORDER>& vars, double time, double DeltaT);

        /*
         * Holds all the information to re-construct the spline values.
         * This has the same length as poly_times_. Each element of the out vector holds a vector of the
         * variables that describe the polynomials. Note that constant polynomials still have two values to
         * describe them even though they only need one.
         */
        std::vector<std::vector<double>> poly_vars_;   // All the variables that define each polynomial

        std::vector<double> poly_times_;    // The start and end of each polynomial

        int total_poly_;

        bool start_on_constant_;

        int num_constant_;

        int num_polys_;

        int num_all_poly_vars;
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_SPLINE_H
