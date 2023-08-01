#include "bfe.hpp"
#include "add_functions.hpp"
#include "comp_funcs.hpp"
#include "flatten.hpp"
#include "execute_in_parallel.hpp"

namespace basis_functions {

class PotentialFromDensity {

    const std::array<int, 2> nl_max;
    const std::vector<std::vector<std::complex<double>>> coefficients;


    template <typename DensityType>
    std::vector<std::vector<std::complex<double>>> 
    getCoefficients(const DensityType& density, const BFE& expansion) const;

    template <typename DataType>
    std::function<DataType(double, double)>
    getTruncFunction(std::function<std::function<DataType(double, double)>
                        (int, int)> BFE_member_function,
                     const BFE& expansion) const;

    std::function<double(double, double)> 
    getTruncDensity(const BFE& expansion) const;
    std::function<double(double, double)> 
    getTruncPotential(const BFE& expansion) const;
    std::function<std::array<double, 2>(double, double)> 
    getTruncForce(const BFE& expansion) const;

    public:
        const std::function<double(double, double)> trunc_density;
        const std::function<double(double, double)> trunc_potential;
        const std::function<std::array<double, 2>(double, double)> trunc_force;

        template <typename DensityType>
        PotentialFromDensity(const DensityType& density, const BFE& expansion, 
                             int n_max, int l_max);
        
        // For diagnostic purposes
        std::vector<std::vector<std::complex<double>>> getCoefficients() const;
};


}