#include "bfe.hpp"

using namespace std::complex_literals;
using special_functions::Gamma;
using special_functions::GammaHalf;
using special_functions::getPochhammerInt;
using special_functions::getPochhammerHalfInt;

namespace basis_functions {

namespace {

double P(int k, int l, int n) {
    double output = (2*k + l + 2*n + 0.5)*GammaHalf(2*k + l + n)
                    /(Gamma(2*k + n + 1)*pow(Gamma(l + 1), 2)*Gamma(n + 1))
                    *GammaHalf(l + n);
    assert(output > 0);
    output = sqrt(output);
    assert(std::isfinite(output));
    return output;
}

double S(int k, int l, int n) {
    double output = (2*k + l + 2*n + 0.5)*Gamma(2*k + n + 1)
                         *GammaHalf(2*k + l + n)
                         /(GammaHalf(l + n)*Gamma(n + 1));
    assert(output != 0);
    output = sqrt(output);
    output = output*Gamma(k + 1)/(std::numbers::pi*Gamma(2*k + 1)
                                  *GammaHalf(k));
    assert(output != 0);
    assert(std::isfinite(output));
    return output;
}

LooongDouble alpha_Ka(int k, int l, int n, int i, int j) {
    // std::cout << k << ", " << l << ", " << n << ", " << i << ", " << j << std::endl;
    LooongDouble output = getPochhammerInt(-k, i)*getPochhammerHalfInt(l, i)
                    *getPochhammerHalfInt(2*k + l + n, j)
                    /(getPochhammerInt(l + 1, i)*getPochhammerInt(1, i)
                      *getPochhammerInt(l + i + 1, j)*getPochhammerHalfInt(l, j) 
                      *getPochhammerInt(1, j))
                    *getPochhammerHalfInt(i + l, j)*getPochhammerInt(-n, j);
    assert(std::isfinite(output.convert_to<double>()));
    return output;
}

LooongDouble beta_Ka(int k, int l, int n, int j) {
    LooongDouble output = getPochhammerHalfInt(2*k + l + n, j)*getPochhammerInt(k + 1, j)
                    *getPochhammerInt(-n, j)
                    /(getPochhammerInt(2*k + 1, j)*getPochhammerHalfInt(k, j)
                      *getPochhammerInt(1, j));
    double check = output.convert_to<double>();
    if(!std::isfinite(check)) {
        std::cout << k << ", " << l << ", " << n << ", " << j << ", " << output << std::endl;
        std::cout << getPochhammerHalfInt(2*k + l + n, j) << ", " 
                  << getPochhammerInt(k + 1, j) << ", " 
                  << getPochhammerInt(-n, j) << ", " 
                  << getPochhammerInt(2*k + 1, j) << ", " << getPochhammerHalfInt(k, j)
                  << ", " << getPochhammerInt(1, j) << std::endl;
    }
    assert(std::isfinite(check));
    return output;
}

}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

namespace {

double U(int k, int n, int l, double R_norm) {
    // Missing R_Ka dependence, which will be reapplied later
    LooongDouble output = 0;
    for(int i = 0; i <= k; ++i) {
        for(int j = 0; j <= n; ++j) {
            LooongDouble term = getAlphaKa(k, l, n, i, j)*pow(R_norm, 2*i + 2*j + l);
            output += term;
        }
    }
    double output_double = output.convert_to<double>();
    output_double *= -sqrt(Units::G)*getP(k, l, n);
    // if(R > 0) {
    //     assert(output != 0);
    // }
    return output_double;

}

double UPrime(int k, int n, int l, double R_norm) {
    // Again need to reapply R_Ka dependence
    LooongDouble output = 0;
    for(int i = 0; i <= k; ++i) {
        for(int j = 0; j <= n; ++j) {
            if(i == 0 && j == 0 && l == 0) {
                continue;
            }
            LooongDouble term = getAlphaKa(k, l, n, i, j)
                          *(2*i + 2*j + l)*pow(R_norm, 2*i + 2*j + l - 1);
            output += term;
        }
    }
    double output_double = output.convert_to<double>();
    output_double *= -sqrt(Units::G)*getP(k, l, n);
    return output_double;
}

double D(int k, int n, int l, double R_norm) {
    // As above
    LooongDouble output = 0;
    for(int j = 0; j <= n; ++j) {
        LooongDouble term = getBetaKa(k, l, n, j)*pow(1 - pow(R_norm, 2), j);
        // std:: cout << "D, " << term << std::endl;
        output += term;
    }
    double output_double = output.convert_to<double>();
    output_double *= pow(-1, n)/(sqrt(Units::G))
             *getS(k, l, n)*pow(1 - pow(R_norm, 2), k - 0.5)*pow(R_norm, l);
    // if(R > 0 && R < R_Ka) {
    //     assert(output != 0);
    // }
    return output_double;
}



}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
namespace {
utility::vector4d<double> 
calculateUUpDValues(const std::function<double(int, int, int, double)>& which) {
    std::array<int, 4> shape = {k_max + 1, n_max + 1, l_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    // Bundling up to prevent under-utilisation of cores
    std::vector<std::function<std::vector<double>()>>
    calculation_functions(shape[0]*shape[1]*shape[2]);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                calculation_functions[i*shape[1]*shape[2]
                                      + j*shape[2] + k]
                    = [i, j, k, N_R_tabulated, &which] {
                    std::vector<double> output(N_R_tabulated);
                    for(int l = 0; l < N_R_tabulated; ++l) {
                        // Calculate central value in R bin
                        double R_norm = (double)(l + 0.5)/N_R_tabulated;
                        output[l] = which(i, j, k, R_norm);
                    }
                    return output;
                };
            }
        }
    }
    std::vector<std::vector<double>> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector4d<double> output = utility::makeShape<double>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    output[i][j][k][l] = 
                    flat[i*shape[1]*shape[2]
                         + j*shape[2] + k][l];                    
                }
            }
        }
    }
    return output;    
} 
}

utility::vector4d<double> calculateUValues() {
    return calculateUUpDValues(U);
}

utility::vector4d<double> calculateUPrimeValues() {
    return calculateUUpDValues(UPrime);
}

utility::vector4d<double> calculateDValues() {
    return calculateUUpDValues(D);
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector5d<LooongDouble>& 
alpha_Ka_values() {
    static utility::vector5d<LooongDouble> x;
    if(x.empty()) {x = getAlphaKaValues();};
    return x;
}
utility::vector4d<LooongDouble>& 
beta_Ka_values() {
    static utility::vector4d<LooongDouble> x;
    if(x.empty()) {x = getBetaKaValues();};
    return x;
}
utility::vector3d<double>& 
P_values() {
    static utility::vector3d<double> x;
    if(x.empty()) {x = getPValues();};
    return x;
}
utility::vector3d<double>& 
S_values() {
    static utility::vector3d<double> x;
    if(x.empty()) {x = getSValues();};
    return x;
}


utility::vector4d<double>& U_values() {
    static utility::vector4d<double> x;
    if(x.empty()) {x = getUValues();};
    return x;
}
utility::vector4d<double>& UPrime_values() {
    static utility::vector4d<double> x;
    if(x.empty()) {x = getUPrimeValues();};
    return x;
}
utility::vector4d<double>& D_values() {
    static utility::vector4d<double> x;
    if(x.empty()) {x = getDValues();};
    return x;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector4d<double> getUValues() {
    std::array<int, 4> shape = {k_max + 1, l_max + 1, n_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::string path = "../cache/basis_functions/u_values.csv";
    if(utility::fileExists(path)) {
        std::vector<double> flat = utility::readCsv(path);
        if(flat.size() == N_values) {
            return utility::reshape(flat, shape);
        }
    }
    utility::vector4d<double> output = calculateUValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector4d<double> getUPrimeValues() {
    std::array<int, 4> shape = {k_max + 1, l_max + 1, n_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::string path = "../cache/basis_functions/u_prime_values.csv";
    if(utility::fileExists(path)) {
        std::vector<double> flat = utility::readCsv(path);
        if(flat.size() == N_values) {
            return utility::reshape(flat, shape);
        }
    }
    utility::vector4d<double> output = calculateUPrimeValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector4d<double> getDValues() {
    std::array<int, 4> shape = {k_max + 1, l_max + 1, n_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::string path = "../cache/basis_functions/d_values.csv";
    if(utility::fileExists(path)) {
        std::vector<double> flat = utility::readCsv(path);
        if(flat.size() == N_values) {
            return utility::reshape(flat, shape);
        }
    }
    utility::vector4d<double> output = calculateDValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

double getU(int k, int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = U_values()[k][n][l][R_bin];
    output /= pow(R_Ka, 0.5);
    return output;
}

double getUPrime(int k, int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = UPrime_values()[k][n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

double getD(int k, int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = D_values()[k][n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector5d<LooongDouble> getAlphaKaValues() {
    std::array<int, 5> shape = {k_max + 1, l_max + 1, n_max + 1,
                                i_max + 1, j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3]*shape[4];
    std::vector<std::function<LooongDouble()>> calculation_functions(N_values);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    for(int m = 0; m < shape[4]; ++m) {
                        calculation_functions[i*shape[1]*shape[2]*shape[3]*shape[4]
                                              + j*shape[2]*shape[3]*shape[4]
                                              + k*shape[3]*shape[4]
                                              + l*shape[4] + m]
                            = [i, j, k, l, m] () {
                            return alpha_Ka(i, j, k, l, m);
                        };
                    }
                }
            }
        }
    }
    std::vector<LooongDouble> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector5d<LooongDouble> output = utility::makeShape<LooongDouble>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    for(int m = 0; m < shape[4]; ++m) {
                        output[i][j][k][l][m] = flat[i*shape[1]*shape[2]*shape[3]*shape[4]
                                              + j*shape[2]*shape[3]*shape[4]
                                              + k*shape[3]*shape[4]
                                              + l*shape[4] + m];
                    }
                }
            }
        }
    }
    return output;
}

utility::vector4d<LooongDouble> getBetaKaValues() {
    std::array<int, 4> shape = {k_max + 1, l_max + 1, n_max + 1,
                                j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::vector<std::function<LooongDouble()>> calculation_functions(N_values);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    calculation_functions[i*shape[1]*shape[2]*shape[3]
                                              + j*shape[2]*shape[3]
                                              + k*shape[3] + l]
                            = [i, j, k, l] () {
                        return beta_Ka(i, j, k, l);
                    };
                }
            }
        }
    }
    std::vector<LooongDouble> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector4d<LooongDouble> output = utility::makeShape<LooongDouble>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                    output[i][j][k][l] = flat[i*shape[1]*shape[2]*shape[3]
                                              + j*shape[2]*shape[3]
                                              + k*shape[3] + l];
                }
            }
        }
    }
    return output;
}
namespace {
utility::vector3d<double> getPSValues(std::function<double(int, int, int)> which) {
    std::array<int, 3> shape = {k_max + 1, l_max + 1, n_max + 1};
    int N_values = shape[0]*shape[1]*shape[2];
    std::vector<std::function<double()>> calculation_functions(N_values);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                calculation_functions[i*shape[1]*shape[2] + j*shape[2] + k]
                            = [i, j, k, &which] () {
                    return which(i, j, k);
                };
            }
        }
    }
    std::vector<double> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector3d<double> output = utility::makeShape<double>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                output[i][j][k] = flat[i*shape[1]*shape[2] + j*shape[2] + k];
            }
        }
    }
    return output;
}
}

utility::vector3d<double> getPValues() {
    return getPSValues(P);
}

utility::vector3d<double> getSValues() {
    return getPSValues(S);
}


LooongDouble getAlphaKa(int k, int l, int n, int i, int j) {
    assert(k >= 0);
    assert(k <= k_max);
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(i >= 0);
    assert(i <= i_max);
    assert(j >= 0);
    assert(j <= j_max);

    return alpha_Ka_values()[k][l][n][i][j];
}

LooongDouble getBetaKa(int k, int l, int n, int j) {
    assert(k >= 0);
    assert(k <= k_max);
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(j >= 0);
    assert(j <= j_max);

    return beta_Ka_values()[k][l][n][j];
}

double getP(int k, int l, int n) {
    assert(k >= 0);
    assert(k <= k_max);
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    return P_values()[k][l][n];
}

double getS(int k, int l, int n) {
    assert(k >= 0);
    assert(k <= k_max);
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    return S_values()[k][l][n];
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
              const std::function<std::complex<double>(double, double)>& density,
              double R_max, int N_R, int N_phi) {
    double dR = R_max/N_R;
    double dphi = 2*std::numbers::pi/N_phi;
    std::complex<double> output = 0;
    for(int i = 0; i < N_R; ++i) {
        for(int j = 0; j < N_phi; ++j) {
            // Avoid possible singularities at R = 0 and R = R_Ka
            double R = (double)(i + 0.5)/N_R*R_max;
            double phi = (double)j/N_phi*2*std::numbers::pi;
            output += R*dR*dphi*pot_conj(R, phi)*density(R, phi);
        }
    }
    return output;
}

std::complex<double>
scalarProduct(const std::function<std::complex<double>(double, double)>& pot_conj,
              const std::vector<std::array<double, 3>>& density, double R_max) {
    // Each particle has array {R, phi, m}, so integral becomes a sum
    int N_particles = density.size();
    std::complex<double> output = 0;
    for(int i = 0; i < N_particles; ++i) {
        assert(density[i][0] < R_max);
        output += pot_conj(density[i][0], density[i][1])*density[i][2];
    }
    return output;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


BFE::BFE(int k_Ka_, double R_Ka_, int N_R_, int N_phi_): 
k_Ka(k_Ka_), R_Ka(R_Ka_), N_R(N_R_), N_phi(N_phi_) {}
BFE::BFE(const BFE& old): k_Ka(old.k_Ka), R_Ka(old.R_Ka),
                          N_R(old.N_R), N_phi(old.N_phi) {}

std::function<std::complex<double>(double, double)> 
BFE::psi(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [*this, n, l] (double R, double phi) {
        double prefactor = getU(k_Ka, n, l, R, R_Ka);
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        return prefactor*phase;
    };
}

std::function<std::complex<double>(double, double)> 
BFE::rho(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [*this, n, l] (double R, double phi) {
        double prefactor = getD(k_Ka, n, l, R, R_Ka);
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        return prefactor*phase;
    };
}

std::function<std::array<std::complex<double>, 2>(double, double)>
BFE::psi_f(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [*this, n, l] (double R, double phi) {
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        std::array<std::complex<double>, 2> output;
        output[0] = -getUPrime(k_Ka, n, l, R, R_Ka)*phase;
        output[1] = -1i*(double)l/R*getU(k_Ka, n, l, R, R_Ka)*phase;
        return output;
    };
}

// // Calculate as LooongDouble, but return as double
// double BFE::U(int n, int l, double R) const {
//     LooongDouble output = 0;
//     for(int i = 0; i <= k_Ka; ++i) {
//         for(int j = 0; j <= n; ++j) {
//             LooongDouble term = getAlphaKa(k_Ka, l, n, i, j)*pow(R/R_Ka, 2*i + 2*j + l);
//             output += term;
//         }
//     }
//     double output_double = output.convert_to<double>();
//     output_double *= -sqrt(Units::G)/sqrt(R_Ka)*getP(k_Ka, l, n);
//     // if(R > 0) {
//     //     assert(output != 0);
//     // }
//     return output_double;
// }

// double BFE::UPrime(int n, int l, double R) const {
//     LooongDouble output = 0;
//     for(int i = 0; i <= k_Ka; ++i) {
//         for(int j = 0; j <= n; ++j) {
//             if(i == 0 && j == 0 && l == 0) {
//                 continue;
//             }
//             LooongDouble term = getAlphaKa(k_Ka, l, n, i, j)
//                           *(2*i + 2*j + l)*pow(R/R_Ka, 2*i + 2*j + l - 1);
//             output += term;
//         }
//     }
//     double output_double = output.convert_to<double>();
//     output_double *= -sqrt(Units::G)/pow(R_Ka, 1.5)*getP(k_Ka, l, n);
//     return output_double;
// }

// double BFE::D(int n, int l, double R) const {
//     LooongDouble output = 0;
//     for(int j = 0; j <= n; ++j) {
//         LooongDouble term = getBetaKa(k_Ka, l, n, j)*pow(1 - pow(R/R_Ka, 2), j);
//         // std:: cout << "D, " << term << std::endl;
//         output += term;
//     }
//     double output_double = output.convert_to<double>();
//     output_double *= pow(-1, n)/(sqrt(Units::G)*pow(R_Ka, 1.5))
//              *getS(k_Ka, l, n)*pow(1 - pow(R/R_Ka, 2), k_Ka - 0.5)*pow(R/R_Ka, l);
//     // if(R > 0 && R < R_Ka) {
//     //     assert(output != 0);
//     // }
//     return output_double;
// }

std::complex<double>
BFE::getCoefficient(int n, int l, const std::function<double(double, double)>&
               density) const {
    return -scalarProduct(utility::conjugateFunction(psi(n, l)),
                            density, R_Ka, N_R, N_phi);
}

std::complex<double>
BFE::getCoefficient(int n, int l, const std::vector<std::array<double, 3>>&
               density) const {
    return -scalarProduct(utility::conjugateFunction(psi(n, l)),
                            density, R_Ka);
}

}