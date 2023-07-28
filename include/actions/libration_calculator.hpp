#include "jres_finder.hpp"
#include "mestel_spiral.hpp"
#include <optional>
#include <assert.h>
#include <cmath>


namespace actions {

class JacobiSpecialFunctions {
    const int N_intervals;
    double phi(double u, double m) const;
    public:
        JacobiSpecialFunctions(int N_intervals_): N_intervals(N_intervals_) {}
        JacobiSpecialFunctions(const JacobiSpecialFunctions& old): 
        N_intervals(old.N_intervals) {}
        double cn(double u, double m) const;
        double sn(double u, double m) const;
        double cnInv(double x, double m) const;
        double K(double m) const;
};


class LibrationCalculator {
    // Notation from Monari et al 2017;
    const int N_R;
    const int N_phi;
    const double J_f;
    const double J_s_res;
    const double F;
    const double G;

    std::optional<double> k;
    std::optional<double> C;


    double calculateEp(double J_s, double theta_s) const;
    double calculatek(double J_s, double theta_s) const;

    const JacobiSpecialFunctions jfuncs;

    public:
        const double g;
        double omega_0() const;
        // Gives both separatrix curve J_s values as function of theta_s
        std::array<double, 2> separatrices(double theta_s) const;

        double J_a;
        std::array<double, 2> trajectory(double t) const;
     


        LibrationCalculator(double J_f_, const JresFinder& finder,
                            const potential::LibrationFuncs& libfuncs,
                            double delta_J_s, int N_jacobi_intervals);
        // LibrationCalculator(double J_f_, const JresFinder& finder,
        //                     std::function<double(double, double)> F_,
        //                     std::function<double(double, double)> g_,
        //                     double delta_J_s, int N_jacobi_intervals);
        LibrationCalculator(const LibrationCalculator& old);
        
        void setICs(double J_s, double theta_s);
};



}