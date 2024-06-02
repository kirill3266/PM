//
// Created by kirill3266 on 24.05.24.
//

#ifndef PM_UTILITY_H
#define PM_UTILITY_H

#include <fftw3.h>
#include <complex>
#include <vector>

namespace Utility {

    double logPower(const std::complex<double>& t_x, double scale);

    std::vector<double> logPower(const std::vector<std::complex<double>>& t_x, double scale);

    std::vector<double> magnitude(const std::vector<std::complex<double>>& t_x);

    std::vector<std::complex<double>> toComplex(const std::vector<double>& t_x);

    void goertzel_cmplx(std::vector<std::complex<double>> &in, std::vector<double> f,
                        std::vector<std::complex<double>> &out);

    std::vector<std::complex<double>> FFT(const std::vector<std::complex<double>> &t_data, int n_threads);

    std::vector<std::complex<double>> winFFT(const std::vector<std::complex<double>> &t_data, int n_threads);

    std::vector<std::complex<double>> iFFT(const std::vector<std::complex<double>> &t_data, int n_threads);

    void unwrap(std::vector<double>& phi, double lev, double mar);

    std::vector<double> conv(const std::vector<double> &t_x, const std::vector<double> &t_y);

    std::vector<std::complex<double>> filter(const std::vector<std::complex<double>> &t_x, std::vector<std::complex<double>> t_y,int t_len, int n_threads);

}
#endif //PM_UTILITY_H
