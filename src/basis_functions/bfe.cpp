#include "bfe.hpp"
#include <iostream>
#include <mutex>

std::mutex mtx;

using namespace std::complex_literals;
namespace mp = boost::multiprecision;

namespace basis_functions {

// double only goes up to 1e307 or so
double BFETables::P(int k, int l, int n) const {
    double output = (LooongDouble(2*k + l + 2*n + 0.5)
                          *LooongDouble(gtables.value().GammaHalf(2*k + l + n))
                    /(LooongDouble(gtables.value().Gamma(2*k + n + 1))
                      *LooongDouble(pow(gtables.value().Gamma(l + 1), 2))
                      *LooongDouble(gtables.value().Gamma(n + 1)))
                    *LooongDouble(gtables.value().GammaHalf(l + n))).convert_to<double>();
    // mtx.lock();
    // std::cout << k << ", " << l << ", " << n << ", " << output << std::endl;
    // mtx.unlock();
    assert(output > 0);
    output = sqrt(output);
    assert(std::isfinite(output));
    return output;
}

double BFETables::S(int k, int l, int n) const {
    LooongDouble output = LooongDouble(2*k + l + 2*n + 0.5)
                         *LooongDouble(gtables.value().Gamma(2*k + n + 1))
                         *LooongDouble(gtables.value().GammaHalf(2*k + l + n))
                         /(LooongDouble(gtables.value().GammaHalf(l + n))
                           *LooongDouble(gtables.value().Gamma(n + 1)));
    assert(output.convert_to<double>() != 0);
    output = sqrt(output);
    output = output*LooongDouble(gtables.value().Gamma(k + 1))
                /(LooongDouble(std::numbers::pi)
                  *LooongDouble(gtables.value().Gamma(2*k + 1))
                    *LooongDouble(gtables.value().GammaHalf(k)));
    double output_double = output.convert_to<double>();
    assert(output_double > 0);
    if(!std::isfinite(output_double)) {
        mtx.lock();
        std::cout << "S: " << k << ", " << l << ", " << n << ", " 
                  << 2*k + l + 2*n + 0.5 << ", " << gtables.value().Gamma(2*k + n + 1) << ", " 
                  << gtables.value().GammaHalf(2*k + l + n) << ", " << gtables.value().GammaHalf(l + n) << ", "
                  << gtables.value().Gamma(n + 1) << ", " << output_double << std::endl;
        mtx.unlock();
    }
    assert(std::isfinite(output_double));
    return output_double;
}

LooongDouble BFETables::alpha_Ka(int k, int l, int n, int i, int j) const {
    LooongDouble output = ptables.value().getPochhammerInt(-k, i)*ptables.value().getPochhammerHalfInt(l, i)
                    *ptables.value().getPochhammerHalfInt(2*k + l + n, j)
                    /(ptables.value().getPochhammerInt(l + 1, i)*ptables.value().getPochhammerInt(1, i)
                      *ptables.value().getPochhammerInt(l + i + 1, j)*ptables.value().getPochhammerHalfInt(l, j) 
                      *ptables.value().getPochhammerInt(1, j))
                    *ptables.value().getPochhammerHalfInt(i + l, j)*ptables.value().getPochhammerInt(-n, j);
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

LooongDouble BFETables::beta_Ka(int k, int l, int n, int j) const {
    LooongDouble output = ptables.value().getPochhammerHalfInt(2*k + l + n, j)*ptables.value().getPochhammerInt(k + 1, j)
                    *ptables.value().getPochhammerInt(-n, j)
                    /(ptables.value().getPochhammerInt(2*k + 1, j)*ptables.value().getPochhammerHalfInt(k, j)
                      *ptables.value().getPochhammerInt(1, j));
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


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

BFETables::BFETables(): U_values(getUValues()), 
                        UPrime_values(getUPrimeValues()),
                        D_values(getDValues()) {
    // Clear out intermediate tables from memory
    gtables.reset();
    ptables.reset();
    alpha_Ka_values.reset();
    beta_Ka_values.reset();
    P_values.reset();
    S_values.reset();
}

BFETables::BFETables(const BFETables& old):
    U_values(old.U_values), UPrime_values(old.UPrime_values),
    D_values(old.D_values) {}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


double BFETables::U(int k, int n, int l, double R_norm) const {
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

double BFETables::UPrime(int k, int n, int l, double R_norm) const {
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

double BFETables::D(int k, int n, int l, double R_norm) const {
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




////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

void BFETables::ensureSubTablesExist() {
    if(!gtables) {
        std::cout << "Calculating Gammas" << std::endl;
        gtables.emplace();
    }
    if(!ptables) {
        std::cout << "Calculating Pochhammers" << std::endl;
        ptables.emplace();
    }
    if(!alpha_Ka_values) {
        std::cout << "Calculating alpha_Ka" << std::endl;
        alpha_Ka_values.emplace(getAlphaKaValues());
    }
    if(!beta_Ka_values) {
        std::cout << "Calculating beta_Ka" << std::endl;
        beta_Ka_values.emplace(getBetaKaValues());
    }
    if(!P_values) {
        std::cout << "Calculating P" << std::endl;
        P_values.emplace(getPValues());
    }
    if(!S_values) {
        std::cout << "Calculating S" << std::endl;
        S_values.emplace(getSValues());
    }
}

utility::vector3d<double> 
BFETables::calculateUUpDValues(
        std::function<double(int, int, int, double)> which) {
    ensureSubTablesExist();
    std::cout << "Calculating U/U_prime/D" << std::endl;
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

utility::vector3d<double> BFETables::calculateUValues() {
    std::function<double(int, int, int, double)>
    f = [this] (int k, int n, int l, double R_norm) {
        return U(k, n, l, R_norm);
    };
    return calculateUUpDValues(f);
}

utility::vector3d<double> BFETables::calculateUPrimeValues() {
    std::function<double(int, int, int, double)>
    f = [this] (int k, int n, int l, double R_norm) {
        return UPrime(k, n, l, R_norm);
    };
    return calculateUUpDValues(f);
}

utility::vector3d<double> BFETables::calculateDValues() {
    std::function<double(int, int, int, double)>
    f = [this] (int k, int n, int l, double R_norm) {
        return D(k, n, l, R_norm);
    };
    return calculateUUpDValues(f);
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

std::optional<utility::vector3d<double>> 
BFETables::readUUprimeDValues(std::string path) {
    std::optional<utility::vector3d<double>> output;
    std::array<int, 3> shape = {n_max + 1, l_max + 1, N_R_tabulated};
    int N_values = shape[0]*shape[1]*shape[2];
    if(utility::fileExists(path)) {
        std::vector<double> flat = utility::readCsv(path);
        if(flat.size() == N_values) {
            output.emplace(utility::reshape(flat, shape));
        }
    }
    return output;
}

utility::vector3d<double> BFETables::getUValues() {
    std::string path = "../cache/basis_functions/u_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::optional<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return values.value();
    }
    std::cout << "U values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateUValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector3d<double> BFETables::getUPrimeValues() {
    std::string path = "../cache/basis_functions/u_prime_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::optional<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return values.value();
    }
    std::cout << "U prime values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateUPrimeValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}

utility::vector3d<double> BFETables::getDValues() {
    std::string path = "../cache/basis_functions/d_values_k=" 
                       + std::to_string(k_Ka) + ".csv";
    std::optional<utility::vector3d<double>>
    values = readUUprimeDValues(path);
    if(values) {
        return values.value();
    }
    std::cout << "D values not cached; caching..." << std::endl;
    utility::vector3d<double> output = calculateDValues();
    utility::writeCsv(path, utility::flatten(output));
    return output;
}



////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

double BFETables::getU(int n, int l, double R, double R_Ka) const {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = U_values[n][l][R_bin];
    output /= pow(R_Ka, 0.5);
    return output;
}

double BFETables::getUPrime(int n, int l, double R, double R_Ka) const {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = UPrime_values[n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

double BFETables::getD(int n, int l, double R, double R_Ka) const {
    int R_bin = R/R_Ka*N_R_tabulated;
    double output = D_values[n][l][R_bin];
    output /= pow(R_Ka, 1.5);
    return output;
}

std::function<double(double, double)>
BFETables::getUFunction(int n, int l) const {
    // std::vector<double> values = U_values[n][l];
    return [this, n, l, N_R_tabulated=N_R_tabulated] (double R, double R_Ka) {
        int R_bin = R/R_Ka*N_R_tabulated;
        double output = U_values[n][l][R_bin];
        output /= pow(R_Ka, 0.5);
        return output;
    };
}

std::function<double(double, double)>
BFETables::getUPrimeFunction(int n, int l) const {
    // std::vector<double> values = UPrime_values[n][l];
    return [this, n, l, N_R_tabulated=N_R_tabulated] (double R, double R_Ka) {
        int R_bin = R/R_Ka*N_R_tabulated;
        double output = UPrime_values[n][l][R_bin];
        output /= pow(R_Ka, 1.5);
        return output;
    };
}

std::function<double(double, double)>
BFETables::getDFunction(int n, int l) const {
    // std::vector<double> values = D_values[n][l];
    return [this, n, l, N_R_tabulated=N_R_tabulated] (double R, double R_Ka) {
        int R_bin = R/R_Ka*N_R_tabulated;
        double output = D_values[n][l][R_bin];
        output /= pow(R_Ka, 1.5);
        return output;
    };
}


////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

utility::vector4d<LooongDouble> BFETables::getAlphaKaValues() const {
    std::array<int, 4> shape = {l_max + 1, n_max + 1,
                                i_max + 1, j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2]*shape[3];
    std::vector<std::function<std::vector<LooongDouble>()>> 
    calculation_functions(N_values/(shape[2]*shape[3]));
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            calculation_functions[i*shape[1] + j]
                = [i, j, shape, this] () {
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

utility::vector3d<LooongDouble> BFETables::getBetaKaValues() const {
    std::array<int, 3> shape = {l_max + 1, n_max + 1,
                                j_max + 1};
    int N_values = shape[0]*shape[1]*shape[2];
    std::vector<std::function<std::vector<LooongDouble>()>> 
    calculation_functions(N_values/shape[2]);
    for(int i = 0; i < shape[0]; ++i) {
        for(int j = 0; j < shape[1]; ++j) {
            calculation_functions[i*shape[1] + j]
                = [i, j, shape, this] () {
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

utility::vector2d<double> 
BFETables::getPSValues(std::function<double(int, int, int)> which) const {
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

utility::vector2d<double> BFETables::getPValues() const {
    std::function<double(int, int, int)>
    f = [this] (int k, int l, int n) {
        return P(k, l, n);
    };
    return getPSValues(f);
}

utility::vector2d<double> BFETables::getSValues() const {
    std::function<double(int, int, int)>
    f = [this] (int k, int l, int n) {
        return S(k, l, n);
    };
    return getPSValues(f);
}


LooongDouble BFETables::getAlphaKa(int l, int n, int i, int j) const {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(i >= 0);
    assert(i <= i_max);
    assert(j >= 0);
    assert(j <= j_max);

    assert(alpha_Ka_values.value()[l][n][i][j].convert_to<double>() != 0);

    return alpha_Ka_values.value()[l][n][i][j];
}

LooongDouble BFETables::getBetaKa(int l, int n, int j) const {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);
    assert(j >= 0);
    assert(j <= j_max);

    assert(beta_Ka_values.value()[l][n][j].convert_to<double>() != 0);

    return beta_Ka_values.value()[l][n][j];
}

double BFETables::getP(int l, int n) const {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    if(P_values.value()[l][n] <= 0) {
        mtx.lock();
        std::cout << l << ", " << n << ", " << P_values.value()[l][n] << std::endl;
        mtx.unlock();
    }

    assert(P_values.value()[l][n] > 0);

    return P_values.value()[l][n];
}

double BFETables::getS(int l, int n) const {
    assert(l >= 0);
    assert(l <= l_max);
    assert(n >= 0);
    assert(n <= n_max);

    if(S_values.value()[l][n] <= 0) {
        mtx.lock();
        std::cout << l << ", " << n << ", " << S_values.value()[l][n] << std::endl;
        mtx.unlock();
    }

    assert(S_values.value()[l][n] > 0);

    return S_values.value()[l][n];
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
    int id = std::rand();
    int N_particles = density.size();
    std::complex<double> output = 0;
    for(int i = 0; i < N_particles; ++i) {
        // In simulation, particle may leave basis functions domain.
        // We'll say then that it doesn't contribute to any scalar product.
        // assert(density[i][0] < R_max);
        if(density[i][0] < R_max) {
            output += pot_conj(density[i][0], density[i][1])*density[i][2];
        }
        // for debugging
        else {
            mtx.lock();
            std::cout << "Position: " << density[i][0] << ", " 
                      << R_max << std::endl;
            mtx.unlock();
        }
    }
    return output;
}

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////


// BFE::BFE(double R_Ka_, int N_R_, int N_phi_): 
// R_Ka(R_Ka_), N_R(N_R_), N_phi(N_phi_) {}
// BFE::BFE(const BFE& old): R_Ka(old.R_Ka),
//                           N_R(old.N_R), N_phi(old.N_phi),
//                           tables(old.tables) {
//     mtx.lock();
//     std::cout << "Copied BFE" << std::endl;
//     mtx.unlock();
// }

BFE::BFE(double R_Ka_, int N_R_, int N_phi_):
R_Ka(R_Ka_), N_R(N_R_), N_phi(N_phi_), tables(getTables()) {}
BFE::BFE(const BFE& old): R_Ka(old.R_Ka), N_R(old.N_R),
                          N_phi(old.N_phi), tables(old.tables) {}

std::vector<std::shared_ptr<const BFETables>>
BFE::getTables() const {
    std::vector<std::shared_ptr<const BFETables>> output(N_cores);
    output[0] = std::make_shared<const BFETables>();
    for(int i = 1; i < N_cores; ++i) {
        output[i] = std::make_shared<const BFETables>(*output[0]);
    }
    return output;
}

std::shared_ptr<const BFETables> BFE::accessTables() const {
    mtx.lock();
    int cpu = sched_getcpu();
    std::cout << "cpu: " << cpu << std::endl;
    std::shared_ptr<const BFETables> output = tables[cpu];
    mtx.unlock();
    return output;
}

std::function<std::complex<double>(double, double)> 
BFE::psi(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [this, n, l, R_Ka=R_Ka] (double R, double phi) {
        double prefactor = accessTables()->getU(n, l, R, R_Ka);
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        return prefactor*phase;
    };
}

std::function<std::complex<double>(double, double)> 
BFE::rho(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [this, n, l, R_Ka=R_Ka] (double R, double phi) {
        double prefactor = accessTables()->getD(n, l, R, R_Ka);
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        return prefactor*phase;
    };
}

std::function<std::array<std::complex<double>, 2>(double, double)>
BFE::psi_f(int n, int l) const {
    assert(n >= 0);
    assert(l >= 0);
    return [this, n, l, R_Ka=R_Ka] (double R, double phi) {
        utility::SimpleTimer timer;
        timer.start();
        std::complex<double> phase = std::exp(1i*(double)l*phi);
        std::array<std::complex<double>, 2> output;
        output[0] = -accessTables()->getUPrime(n, l, R, R_Ka)*phase;
        output[1] = -1i*(double)l/R*accessTables()->getU(n, l, R, R_Ka)*phase;
        timer.stop();
        double duration = std::chrono::duration_cast
                            <std::chrono::nanoseconds>(timer.getDuration()).count();
        mtx.lock();
        utility::debug_print("psi_f function: " + std::to_string(duration)
                             + "ns", 0);
        mtx.unlock();
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