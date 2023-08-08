#include "bfe.hpp"
#include <iostream>
#include <mutex>

std::mutex mtx;

using namespace std::complex_literals;
using special_functions::Gamma;
using special_functions::GammaHalf;
using special_functions::getPochhammerInt;
using special_functions::getPochhammerHalfInt;
namespace mp = boost::multiprecision;

namespace basis_functions {

namespace {
// double only goes up to 1e307 or so
double P(int k, int l, int n) {
    double output = (LooongDouble(2*k + l + 2*n + 0.5)
                          *LooongDouble(GammaHalf(2*k + l + n))
                    /(LooongDouble(Gamma(2*k + n + 1))
                      *LooongDouble(pow(Gamma(l + 1), 2))
                      *LooongDouble(Gamma(n + 1)))
                    *LooongDouble(GammaHalf(l + n))).convert_to<double>();
    // mtx.lock();
    // std::cout << k << ", " << l << ", " << n << ", " << output << std::endl;
    // mtx.unlock();
    assert(output > 0);
    output = sqrt(output);
    assert(std::isfinite(output));
    return output;
}

double S(int k, int l, int n) {
    LooongDouble output = LooongDouble(2*k + l + 2*n + 0.5)
                         *LooongDouble(Gamma(2*k + n + 1))
                         *LooongDouble(GammaHalf(2*k + l + n))
                         /(LooongDouble(GammaHalf(l + n))
                           *LooongDouble(Gamma(n + 1)));
    assert(output.convert_to<double>() != 0);
    output = sqrt(output);
    output = output*LooongDouble(Gamma(k + 1))
                /(LooongDouble(std::numbers::pi)
                  *LooongDouble(Gamma(2*k + 1))
                    *LooongDouble(GammaHalf(k)));
    double output_double = output.convert_to<double>();
    assert(output_double > 0);
    if(!std::isfinite(output_double)) {
        mtx.lock();
        std::cout << "S: " << k << ", " << l << ", " << n << ", " 
                  << 2*k + l + 2*n + 0.5 << ", " << Gamma(2*k + n + 1) << ", " 
                  << GammaHalf(2*k + l + n) << ", " << GammaHalf(l + n) << ", "
                  << Gamma(n + 1) << ", " << output_double << std::endl;
        mtx.unlock();
    }
    assert(std::isfinite(output_double));
    return output_double;
}

LooongDouble alpha_Ka(int k, int l, int n, int i, int j) {
    LooongDouble output = getPochhammerInt(-k, i)*getPochhammerHalfInt(l, i)
                    *getPochhammerHalfInt(2*k + l + n, j)
                    /(getPochhammerInt(l + 1, i)*getPochhammerInt(1, i)
                      *getPochhammerInt(l + i + 1, j)*getPochhammerHalfInt(l, j) 
                      *getPochhammerInt(1, j))
                    *getPochhammerHalfInt(i + l, j)*getPochhammerInt(-n, j);
    double output_double = output.convert_to<double>();
    // if(output_double == 0) {
    //     std::cout << k << ", " << l << ", " << n << ", " << i << ", " << j << ", "
    //               << getPochhammerInt(-k, i) << ", " << getPochhammerHalfInt(l, i) 
    //               << ", " << getPochhammerHalfInt(2*k + l + n, j) << ", " 
    //               << getPochhammerInt(l + 1, i) << ", "
    //               << getPochhammerInt(1, i) << ", " << getPochhammerInt(l + i + 1, j)
    //               << ", " << getPochhammerHalfInt(l, j) << ", " 
    //               << ", " << getPochhammerInt(1, j) << ", "
    //               << getPochhammerHalfInt(i + l, j) << ", " 
    //               << getPochhammerInt(-n, j) << std::endl;
    // }
    assert(std::isfinite(output_double));
    // assert(output_double != 0);
    return output;
}

LooongDouble beta_Ka(int k, int l, int n, int j) {
    LooongDouble output = getPochhammerHalfInt(2*k + l + n, j)*getPochhammerInt(k + 1, j)
                    *getPochhammerInt(-n, j)
                    /(getPochhammerInt(2*k + 1, j)*getPochhammerHalfInt(k, j)
                      *getPochhammerInt(1, j));
    double check = output.convert_to<double>();
    // if(!std::isfinite(check)) {
    //     std::cout << k << ", " << l << ", " << n << ", " << j << ", " << output << std::endl;
    //     std::cout << getPochhammerHalfInt(2*k + l + n, j) << ", " 
    //               << getPochhammerInt(k + 1, j) << ", " 
    //               << getPochhammerInt(-n, j) << ", " 
    //               << getPochhammerInt(2*k + 1, j) << ", " << getPochhammerHalfInt(k, j)
    //               << ", " << getPochhammerInt(1, j) << std::endl;
    // }
    assert(std::isfinite(check));
    // assert(check != 0);
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
            // Powers of R_norm must also be computed to high precision
            LooongDouble term = getAlphaKa(l, n, i, j)
                        *mp::pow(LooongDouble(R_norm), 2*i + 2*j);
            output += term;
        }
    }
    double output_double = output.convert_to<double>();
    output_double *= -sqrt(Units::G)*getP(l, n)*pow(R_norm, l);
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
            LooongDouble term = getAlphaKa(l, n, i, j)
                        *(2*i + 2*j + l)*mp::pow(LooongDouble(R_norm), 2*i + 2*j);
            output += term;
        }
    }
    double output_double = output.convert_to<double>();
    output_double *= -sqrt(Units::G)*getP(l, n)*pow(R_norm, l - 1);
    return output_double;
}

double D(int k, int n, int l, double R_norm) {
    // As above
    LooongDouble output = 0;
    for(int j = 0; j <= n; ++j) {
        LooongDouble term = getBetaKa(l, n, j)
            *mp::pow(LooongDouble(1) - mp::pow(LooongDouble(R_norm), 2), j);
        // std:: cout << "D, " << term << std::endl;
        output += term;
    }
    double output_double = output.convert_to<double>();
    output_double *= pow(-1, n)/(sqrt(Units::G))
             *getS(l, n)*pow(1 - pow(R_norm, 2), k - 0.5)*pow(R_norm, l);
    // if(R > 0 && R < R_Ka) {
    //     assert(output != 0);
    // }
    return output_double;
}



}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
namespace {
utility::vector3d<double> 
calculateUUpDValues(const std::function<double(int, int, int, double)>& which) {
    std::array<int, 3> shape = {n_max + 1, l_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2];
    // Bundling up to prevent under-utilisation of cores
    std::vector<std::function<std::vector<double>()>>
    calculation_functions(shape[0]*shape[1]);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            calculation_functions[i*shape[1] + j]
                = [i, j, N_R_tabulated, &which] {
                // std::cout << "UUpD: " << i << ", " << j << std::endl;
                std::vector<double> output(N_R_tabulated);
                for(int k = 0; k < N_R_tabulated; ++k) {
                    // Calculate central value in R bin
                    double R_norm = (double)(k + 0.5)/N_R_tabulated;
                    output[k] = which(k_Ka, i, j, R_norm);
                }
                return output;
            };
        }
    }
    std::vector<std::vector<double>> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector3d<double> output = utility::makeShape<double>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                output[i][j][k] = 
                flat[i*shape[1] + j][k];
            }
        }
    }
    return output;    
} 
}

utility::vector3d<double> calculateUValues() {
    return calculateUUpDValues(U);
}

utility::vector3d<double> calculateUPrimeValues() {
    return calculateUUpDValues(UPrime);
}

utility::vector3d<double> calculateDValues() {
    return calculateUUpDValues(D);
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector4d<LooongDouble>& 
alpha_Ka_values() {
    static utility::vector4d<LooongDouble> x;
    if(x.empty()) {x = getAlphaKaValues();};
    return x;
}
utility::vector3d<LooongDouble>& 
beta_Ka_values() {
    static utility::vector3d<LooongDouble> x;
    if(x.empty()) {x = getBetaKaValues();};
    return x;
}
utility::vector2d<double>& 
P_values() {
    static utility::vector2d<double> x;
    if(x.empty()) {x = getPValues();};
    return x;
}
utility::vector2d<double>& 
S_values() {
    static utility::vector2d<double> x;
    if(x.empty()) {x = getSValues();};
    return x;
}


utility::vector3d<double>& U_values() {
    static utility::vector3d<double> x;
    if(x.empty()) {x = getUValues();};
    return x;
}
utility::vector3d<double>& UPrime_values() {
    static utility::vector3d<double> x;
    if(x.empty()) {x = getUPrimeValues();};
    return x;
}
utility::vector3d<double>& D_values() {
    static utility::vector3d<double> x;
    if(x.empty()) {x = getDValues();};
    return x;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

namespace {
std::unique_ptr<utility::vector3d<double>> readUUprimeDValues(std::string path) {
    std::unique_ptr<utility::vector3d<double>> output(nullptr);
    std::array<int, 3> shape = {l_max + 1, n_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2];
    if(utility::fileExists(path)) {
        std::vector<double> flat = utility::readCsv(path);
        if(flat.size() == N_values) {
            std::unique_ptr<utility::vector3d<double>>
            // Pointers are annoying :-)
            output = std::make_unique<utility::vector3d<double>>
                        (utility::reshape(flat, shape));
        }
    }
    return std::move(output);
}
}

utility::vector3d<double> getUValues() {
    std::string path = "../cache/basis_functions/u_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::unique_ptr<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return *values;
    }
    std::cout << "U values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateUValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector3d<double> getUPrimeValues() {
    std::string path = "../cache/basis_functions/u_prime_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::unique_ptr<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return *values;
    }
    std::cout << "U prime values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateUPrimeValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector3d<double> getDValues() {
    std::string path = "../cache/basis_functions/d_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::unique_ptr<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return *values;
    }
    std::cout << "D values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateDValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

double getU(int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = U_values()[n][l][R_bin];
    output /= pow(R_Ka, 0.5);
    return output;
}

double getUPrime(int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = UPrime_values()[n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

double getD(int n, int l, double R, double R_Ka) {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = D_values()[n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector4d<LooongDouble> getAlphaKaValues() {
    std::array<int, 4> shape = {l_max + 1, n_max + 1,
                                i_max + 1, j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::vector<std::function<std::vector<LooongDouble>()>> 
    calculation_functions(N_values/(shape[2]*shape[3]));
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            calculation_functions[i*shape[1] + j]
                = [i, j, shape] () {
                std::vector<LooongDouble> output(shape[2]*shape[3]);
                for(int k = 0; k < shape[2]; ++k) {
                    for(int l = 0; l < shape[3]; ++l) {
                        output[k*shape[3] + l] = alpha_Ka(k_Ka, i, j, k, l);
                    }
                }
                return output;
            };
        }
    }
    std::vector<std::vector<LooongDouble>> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector4d<LooongDouble> output = utility::makeShape<LooongDouble>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                for(int l = 0; l < shape[3]; ++l) {
                        output[i][j][k][l] = flat[i*shape[1] + j][k*shape[3] + l];
                }
            }
        }
    }
    return output;
}

utility::vector3d<LooongDouble> getBetaKaValues() {
    std::array<int, 3> shape = {l_max + 1, n_max + 1,
                                j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2];
    std::vector<std::function<std::vector<LooongDouble>()>> 
    calculation_functions(N_values/shape[2]);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            calculation_functions[i*shape[1] + j]
                = [i, j, shape] () {
                std::vector<LooongDouble> output(shape[2]);
                for(int k = 0; k < shape[2]; ++k) {
                        output[k] = beta_Ka(k_Ka, i, j, k);
                }
                return output;
            };
        }
    }
    std::vector<std::vector<LooongDouble>> 
    flat = multithreading::executeInParallel(calculation_functions);
    utility::vector3d<LooongDouble> output = utility::makeShape<LooongDouble>(shape);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            for(int k = 0; k < shape[2]; ++k) {
                    output[i][j][k] = flat[i*shape[1] + j][k];
            }
        }
    }
    return output;
}
namespace {
utility::vector2d<double> getPSValues(std::function<double(int, int, int)> which) {
    std::array<int, 2> shape = {l_max + 1, n_max + 1};
    int N_values = shape[0]*shape[1];
    std::vector<std::function<std::vector<double>()>> 
    calculation_functions(N_values/shape[1]);
    for(int i = 0; i < shape[0]; ++i) {
        calculation_functions[i]
                    = [i, shape, &which] () {
            std::vector<double> output(shape[1]);
            for(int j = 0; j < shape[1]; ++j) {
                output[j] = which(k_Ka, i, j);
            }
            return output;
        };
    }
    utility::vector2d<double>
    output = multithreading::executeInParallel(calculation_functions);
    return output;
}
}

utility::vector2d<double> getPValues() {
    return getPSValues(P);
}

utility::vector2d<double> getSValues() {
    return getPSValues(S);
}


LooongDouble getAlphaKa(int l, int n, int i, int j) {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(i >= 0);
    assert(i <= i_max);
    assert(j >= 0);
    assert(j <= j_max);

    assert(alpha_Ka_values()[l][n][i][j].convert_to<double>() != 0);

    return alpha_Ka_values()[l][n][i][j];
}

LooongDouble getBetaKa(int l, int n, int j) {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(j >= 0);
    assert(j <= j_max);

    assert(beta_Ka_values()[l][n][j].convert_to<double>() != 0);

    return beta_Ka_values()[l][n][j];
}

double getP(int l, int n) {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    if(P_values()[l][n] <= 0) {
        mtx.lock();
        std::cout << l << ", " << n << ", " << P_values()[l][n] << std::endl;
        mtx.unlock();
    }

    assert(P_values()[l][n] > 0);

    return P_values()[l][n];
}

double getS(int l, int n) {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    if(S_values()[l][n] <= 0) {
        mtx.lock();
        std::cout << l << ", " << n << ", " << S_values()[l][n] << std::endl;
        mtx.unlock();
    }

    assert(S_values()[l][n] > 0);

    return S_values()[l][n];
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


BFE::BFE(double R_Ka_, int N_R_, int N_phi_): 
R_Ka(R_Ka_), N_R(N_R_), N_phi(N_phi_) {}
BFE::BFE(const BFE& old): R_Ka(old.R_Ka),
                          N_R(old.N_R), N_phi(old.N_phi) {}

std::function<std::complex<double>(double, double)> 
BFE::psi(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [*this, n, l] (double R, double phi) {
        double prefactor = getU(n, l, R, R_Ka);
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        return prefactor*phase;
    };
}

std::function<std::complex<double>(double, double)> 
BFE::rho(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [*this, n, l] (double R, double phi) {
        double prefactor = getD(n, l, R, R_Ka);
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
        output[0] = -getUPrime(n, l, R, R_Ka)*phase;
        output[1] = -1i*(double)l/R*getU(n, l, R, R_Ka)*phase;
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