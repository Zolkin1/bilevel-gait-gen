//
// Created by zolkin on 1/30/24.
//

#ifndef BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H
#define BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H

#include "spline/spline_node.h"

namespace mpc {
    enum TimeType {
        LiftOff = 0,
        TouchDown = 1,
        Inter = 2
    };

    class SplineTimes {
    public:
        SplineTimes(double time, TimeType type);
        SplineTimes();

        double GetTime() const;
        TimeType GetType() const;

        void SetTime(double time);
    protected:
    private:
        double time_;
        TimeType type_;
    };

    using time_v = std::vector<SplineTimes>;

    class EndEffectorSplines {
        using node_v = std::vector<SplineNode>;
        using vector_4t = Eigen::Vector4d;
    public:
        enum SplineType {
            Force = 0,
            Position = 1
        };

        EndEffectorSplines(int num_contacts,
                           const std::vector<double>& times,
                           bool start_in_contact,
                           int num_force_polys);

        EndEffectorSplines& operator=(const EndEffectorSplines& ee_spline);


        double ValueAt(SplineType type, int coord, double time) const;

        vector_t GetPolyVarsLin(SplineType type, int coord, double time) const;

        std::pair<int, int> GetVarsIdx(SplineType type, int coord, double time) const;

        bool IsForceMutable(double time) const;

        void AddPoly(double additional_time);

        void RemovePoly(double start_time);

        double ComputePartialWrtTime(SplineType type, int coord, double time, int time_idx) const;     // Note: time_idx is only of contact times

        vector_t ComputeCoefPartialWrtTime(SplineType type, int coord, double time, int time_idx, double dtwdth) const; // Note: time_idx is only of contact times

        vector_t ComputeCoefPartialWrtTime(SplineType type, int coord, double time, int time_idx) const; // Note: time_idx is only of contact times

        bool IsInContact(double time) const;

        // ------------------ Setters ------------------ //
        void SetVars(SplineType type, int coord, int node_idx, const vector_2t& vars);

        void SetContactTimes(time_v& contact_times);

        // ------------------ Getters ------------------ //
        NodeType GetNodeType(SplineType type, int coord, int node_idx) const;

        int GetNumNodes() const;        // Nodes (includes intermediate nodes)

        std::vector<int> GetMutableNodes(SplineType type, int coord) const;    // Should return compressed view nodes (i.e. no double counts)

        std::vector<double> GetTimes() const;

        vector_t GetSplineAsQPVec(SplineType type, int coord) const;

        double GetEndTime() const;

        double GetStartTime() const;

        int GetTotalPolyVars(SplineType type, int coord) const;

        int GetNumContacts() const;

        std::vector<double> GetContactTimeValues() const;

        time_v GetContactTimes() const;

        double GetNextTouchDownTime(double time) const;

        void SetToTouchdown(double time);

        // Function list:
        // - Spline Value
        // - Vars coefficients -- Compressed view
        // - Vars index -- Compressed view
        // - Num vars effecting -- Compressed view
        // - Setting vars
        // - Total vars
        // - Num nodes in the spline
        // - Num of each type

    protected:
    private:
        int GetLowerNodeIdx(SplineType type, int coord, double time) const;
        int GetUpperNodeIdx(SplineType type, int coord, double time) const;

        int ConvertContactNodeToSplineNode(int contact_idx) const;

        double Getx0Coef(double time, double deltat) const;
        double Getx1Coef(double time, double deltat) const;
        double Getx0dotCoef(double time, double deltat) const;
        double Getx1dotCoef(double time, double deltat) const;

        double Getx0CoefPartial(double time, double DeltaT, double dtdth, double dDtdth, bool wrt_t1) const;
        double Getx1CoefPartial(double time, double DeltaT, double dtdth, double dDtdth, bool wrt_t1) const;
        double Getx0dotCoefPartial(double time, double DeltaT, double dtdth, double dDtdth, bool wrt_t1) const;
        double Getx1dotCoefPartial(double time, double DeltaT, double dtdth, double dDtdth, bool wrt_t1) const;

        node_v& SelectSpline(SplineType type, int coord);

        const node_v& SelectSpline(SplineType type, int coord) const;


        static void UnsupportedSpline() ;

        std::array<node_v, 3> forces_;
        std::array<node_v, 3> positions_;
        time_v times_;
        int num_force_polys_;
        int spline_stride_;
        std::vector<NodeType> force_type_pattern_;
        std::vector<NodeType> position_type_pattern_;
        std::vector<NodeType> z_position_type_pattern_;
        std::vector<TimeType> time_type_pattern_;
        static constexpr int POS_VARS = 3;

        static constexpr double FORCE_MULT = 100; // 100 seems to work best (?)
    };
} // mpc


#endif //BILEVEL_GAIT_GEN_END_EFFECTOR_SPLINES_H
