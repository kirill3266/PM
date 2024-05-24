//
// Created by kirill3266 on 24.05.24.
//

#ifndef PM_UTILITY_H
#define PM_UTILITY_H

#include <fftw3.h>
#include <complex>
#include <vector>

namespace Utility {

    double logPower(fftw_complex in, double scale);

    void goertzel_cmplx(std::vector<std::complex<double>> &in, std::vector<double> f,
                        std::vector<std::complex<double>> &out);

    double logPower(double I, double Q, double scale);

    std::vector<double> FFT(const std::vector<double> &t_data, const int t_Fs);

    std::vector<double> FFT(std::vector<std::complex<double>> &t_data, int t_Fs, const int n_threads);

    void unwrap(std::vector<double>& phi, double lev, double mar);

}
#endif //PM_UTILITY_H
