#include "mestel_spiral.hpp"
#include <iostream>

namespace potential {

MestelSpiralEnvelope::MestelSpiralEnvelope(double v_c_, double q, int global_,
                             int m_, double Omega_p_, double alpha_,
                             double f_):
                             v_c(v_c_), global(global_), m(m_), 
                             Omega_p(Omega_p_), f(f_), alpha(alpha_),
                             R_0(v_c_/Omega_p_), 
                Sigma_0(q*pow(v_c_, 2)/(2*std::numbers::pi*Units::G*R_0)) {}

MestelSpiralEnvelope::MestelSpiralEnvelope(const MestelSpiralEnvelope& old):
                             v_c(old.v_c), R_0(old.R_0), Omega_p(old.Omega_p),
                             m(old.m), alpha(old.alpha), Sigma_0(old.Sigma_0),
                             f(old.f), global(old.global) {}

double MestelSpiralEnvelope::s(double R) const {
    return m/(sqrt(2)*v_c)*(Omega_p*R - v_c);
}

double MestelSpiralEnvelope::A(double R) const {
    if(global) {
        return 1;
    }
    double s_ = s(R);
    if(s_ > -1 && s_ < 1) {
        return cos(std::numbers::pi/2*s_);
    }
    else {
        return 0;
    }
}

double MestelSpiralEnvelope::dAdR(double R) const {
    if(global) {
        return 0;
    }
    double s_ = s(R);
    if(s_ > -1 && s_ < 1) {
        double prefactor = -m*std::numbers::pi*Omega_p/(2*sqrt(2)*v_c);
        return prefactor*sin(std::numbers::pi/2*s_);
    }
    else {
        return 0;
    }
}

std::function<double(double, double, double)>
getMestelSpiralDensity(const MestelSpiralEnvelope& envelope) {
    return [=] (double R, double phi, double t) {
        double prefactor = envelope.A(R)*envelope.f*envelope.Sigma_0;
        double phase = envelope.m*(1/tan(envelope.alpha)
                    *(1 - envelope.R_0/R) + phi - envelope.Omega_p*t);
        return prefactor*cos(phase);
    };
}

MestelSpiralPotential::
MestelSpiralPotential(const MestelSpiralEnvelope& envelope_):
                      envelope(envelope_) {}

MestelSpiralPotential::
MestelSpiralPotential(const MestelSpiralPotential& old):
                      envelope(old.envelope) {}

double MestelSpiralPotential::potential(double R, double phi, double t) const {
    double prefactor = envelope.A(R)*2*std::numbers::pi*Units::G*envelope.f
                       *envelope.Sigma_0*sin(envelope.alpha)/envelope.m;
    double phase = envelope.m*(1/tan(envelope.alpha)
                   *(1 - envelope.R_0/R) + phi - envelope.Omega_p*t);
    return prefactor*cos(phase);    
}

std::array<double, 2> MestelSpiralPotential::force(double R, double phi, double t) const {
    double prefactor = 2*std::numbers::pi*Units::G*envelope.f*envelope.Sigma_0
                       *sin(envelope.alpha);
    double phase = envelope.m*(1/tan(envelope.alpha)
                   *(1 - envelope.R_0/R) + phi - envelope.Omega_p*t);
    double R_terms = envelope.A(R)*envelope.R_0/(pow(R, 2)
                     *tan(envelope.alpha))*sin(phase)
                   - envelope.dAdR(R)*cos(phase);
    double phi_term = envelope.A(R)*envelope.Omega_p*sin(phase);
    std::array<double, 2> output = {prefactor*R_terms, prefactor*phi_term};
    return output;
}

MestelSpiralLibrationFuncs::
MestelSpiralLibrationFuncs(const MestelSpiralEnvelope& envelope_):
                           envelope(envelope_) {}

MestelSpiralLibrationFuncs::
MestelSpiralLibrationFuncs(const MestelSpiralLibrationFuncs& old):
    envelope(old.envelope) {}


std::complex<double> 
MestelSpiralLibrationFuncs::clm(double J_f, double J_s, int l) const {
    using namespace std::complex_literals;
    double J_phi = envelope.m*J_s;
    double J_R = l*J_s + J_f;
    double R_g = Rg(J_phi);
    double sgn_l = 2*(l >= 0) - 1;
    
    std::complex<double> temp = 1i*(double)envelope.m/tan(envelope.alpha)
                                 *(1 - envelope.R_0/R_g);


    std::complex<double> phase = std::exp(1i*(double)envelope.m/tan(envelope.alpha)
                                 *(1 - envelope.R_0/R_g));
    

    if(l == 0) {
        double prefactor = envelope.A(R_g)*2*std::numbers::pi*Units::G
                          *envelope.f*envelope.Sigma_0*envelope.R_0
                          *sin(envelope.alpha)/envelope.m;
        return prefactor*phase;
    }
    else if(pow(l, 2) == 1) {
        double prefactor = std::numbers::pi*Units::G*envelope.f
                           *envelope.Sigma_0*envelope.R_0
                           *e(J_R, J_phi);
        std::complex<double> terms = envelope.A(R_g)*sin(envelope.alpha)
                                     *sgn_l*gamma(J_phi)
                     - envelope.dAdR(R_g)*R_g*sin(envelope.alpha)/envelope.m
                     + 1i*envelope.A(R_g)*envelope.R_0/R_g*cos(envelope.alpha);
        return prefactor*terms*phase;
    }
    else {
        return 0;
    }
}


double MestelSpiralLibrationFuncs::F(double J_f, double J_s, int l) const {
    return -std::abs(clm(J_f, J_s, l));
}

double MestelSpiralLibrationFuncs::g(double J_f, double J_s, int l) const {
    return std::arg(clm(J_f, J_s, l));
}

double MestelSpiralLibrationFuncs::Rg(double J_phi)  const {
    return J_phi/envelope.v_c;
}

double MestelSpiralLibrationFuncs::Omega(double R) const {
    return envelope.v_c/R;
}

double MestelSpiralLibrationFuncs::kappa(double R) const {
    return sqrt(2)*Omega(R);
}

double MestelSpiralLibrationFuncs::e(double J_R, double J_phi) const {
    double R_g = Rg(J_phi);
    return sqrt(2*J_R/(kappa(R_g)*pow(R_g, 2)));
}

double MestelSpiralLibrationFuncs::gamma(double J_phi) const {
    double R_g = Rg(J_phi);
    return 2*Omega(R_g)/pow(kappa(R_g), 2);
}

MestelObjects::MestelObjects(const MestelObjects& old):
    mestel(old.mestel), envelope(old.envelope), libfuncs(old.libfuncs),
    spiral(old.spiral), R_res(old.R_res), res_finders(old.res_finders),
    N_phi(old.N_phi) {}

MestelObjects::MestelObjects(double v_c, double R_0, double Omega_p, double q,
                  int global, double alpha, double f, int N_phi_,
                  int N_table_intervals, int N_tau, int N_bisect):
                  N_phi(N_phi_), mestel(getMestel(v_c, R_0)), 
                  envelope(v_c, q, global, N_phi, Omega_p,
                  alpha, f), libfuncs(envelope), spiral(envelope),
                  R_res({{R_0*(1 - sqrt(2)/N_phi)/(1 + sqrt(2)/N_phi),
                          R_0/(1 + sqrt(2)/N_phi), R_0}}),
                  res_finders(std::move(constructFinders(N_phi, Omega_p, 
                  N_table_intervals, N_tau, N_bisect))) {}

MestelObjects::MestelObjects(const MestelSpiralEnvelope& envelope_,
                             int N_table_intervals, int N_tau, int N_bisect):
                             N_phi(envelope_.m), 
                             mestel(getMestel(envelope_.v_c, envelope_.R_0)),
                             envelope(envelope_), libfuncs(envelope_),
                             spiral(envelope_),
                             R_res({{envelope_.R_0*(1 - sqrt(2)/N_phi)/(1 + sqrt(2)/N_phi),
                             envelope_.R_0/(1 + sqrt(2)/N_phi), envelope.R_0}}),
                             res_finders(std::move(constructFinders(N_phi, envelope_.Omega_p, 
                             N_table_intervals, N_tau, N_bisect))) {}

std::array<std::optional<actions::JresFinder>, 3> 
MestelObjects::constructFinders(int N_phi, double Omega_p,
                                int N_table_intervals, int N_tau,
                                int N_bisect) const {
    std::vector<std::function<std::optional<actions::JresFinder>()>>
    builder_functions(3);
    for(int i = 0; i < 3; ++i) {
        int N_R = i - 1;
        builder_functions[i] = [=, this]() {
            std::optional<actions::JresFinder> output;
            output.emplace(mestel, N_R, N_phi, Omega_p,
                           N_table_intervals, N_tau, N_bisect);
            return output;    
        };
    }
    std::vector<std::optional<actions::JresFinder>>
    finders = std::move(multithreading::executeInParallelOpt(builder_functions));
    std::array<std::optional<actions::JresFinder>, 3> output;
    for(int i = 0; i < 3; ++i) {
        output[i].emplace(finders[i].value());
    }
    return output;
}

MestelSpiralEnvelope getStandardEnvelope(int global) {
    double v_c = 220*Units::kms;
    double R_0 = 8*Units::kpc;
    int N_phi = 2;
    double q = 0.5;
    double alpha = std::numbers::pi/6;
    double f = 0.01;
    double Omega_p = v_c*(1 + sqrt(2)/N_phi)/R_0;

    return MestelSpiralEnvelope(v_c, q, global,
                    N_phi, Omega_p, alpha, f);
}


std::function<double(double, double, double)>
getStandardSpiralDensity() {
    std::function<double(double, double, double)>
    output = getMestelSpiralDensity(getStandardEnvelope(0));
    return output;
}

std::function<double(double, double, double)>
getGlobalSpiralDensity() {
    std::function<double(double, double, double)>
    output = getMestelSpiralDensity(getStandardEnvelope(1));
    return output;
}

MestelObjects getStandardSpiral() {
    // Half-mass mestel disc with same v_c as solar neighbourhood.
    // Hosts an m = 2 spiral with OLR at 8 kpc (solar neighbourhood).
    // Spiral has peak surface density fraction (of disc) = 0.01.
    // This is modulated by a cos(s) envelope between the ILR and OLR
    
    // Jresfinder paramaters
    int N_table_intervals = 100;
    int N_tau = 1000;
    int N_bisect = 10;

    MestelObjects output(getStandardEnvelope(0),
                         N_table_intervals, N_tau, N_bisect);
    return output;
}

MestelObjects getGlobalSpiral() {
    // Half-mass mestel disc with same v_c as solar neighbourhood.
    // Hosts an m = 2 spiral with OLR at 8 kpc (solar neighbourhood).
    // Spiral has peak surface density fraction (of disc) = 0.01.
    // This is modulated by a cos(s) envelope between the ILR and OLR
    
    // Jresfinder paramaters
    int N_table_intervals = 100;
    int N_tau = 1000;
    int N_bisect = 10;

    MestelObjects output(getStandardEnvelope(1),
                         N_table_intervals, N_tau, N_bisect);
    return output;
}


}