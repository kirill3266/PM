//
// Created by kirill3266 on 24.05.24.
//

#include "Utility.h"
#include <fftw3.h>
#include <algorithm>
#include <iostream>
#include <algorithm>

void Utility::goertzel_cmplx(std::vector<std::complex<double>> &in, std::vector<double> f,
                             std::vector<std::complex<double>> &out) {

        int m, p;
        std::complex<double> w;
        double alpha;
        std::complex<double> v[3];

        for (p = 0; p < f.size(); p++) {
                w.real(cos(2 * std::numbers::pi * f[p] / static_cast<double>(in.size())));
                w.imag(sin(2 * std::numbers::pi * f[p] / static_cast<double>(in.size())));

                alpha = 2.0 * std::real(w);
                v[0] = {0.0, 0.0};
                v[1] = {0.0, 0.0};
                v[2] = {0.0, 0.0};

                for (m = 0; m < in.size(); m++) {
                        v[2].real(std::real(v[1]));
                        v[1].real(std::real(v[0]));
                        v[0].real(std::real(in[m]) + alpha * std::real(v[1]) - std::real(v[2]));

                        v[2].imag(std::imag(v[1]));
                        v[1].imag(std::imag(v[0]));
                        v[0].imag(std::imag(in[m]) + alpha * std::imag(v[1]) - std::imag(v[2]));
                }
                out[p].real(w.real() * v[0].real() - w.imag() * v[0].imag() - std::real(v[1]));
                out[p].imag(w.real() * v[0].imag() + w.imag() * v[0].real() - std::imag(v[1]));
        }
}

// Вычисление спектральной плотности энергии и перевод её в дб
double Utility::logPower(const std::complex<double> &t_x, double scale) {
        double re = t_x.real() * scale;
        double im = t_x.imag() * scale;
        double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
        return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

std::vector<double> Utility::logPower(const std::vector<std::complex<double>> &t_x, double scale) {
        std::vector<double> res(t_x.size());
        double re, im, magsq;
        for (std::size_t i = 0; i < t_x.size(); ++i) {
                re = t_x[i].real() * scale;
                im = t_x[i].imag() * scale;
                magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
                res[i] = std::log(magsq) * 10.0 / std::log(10.0); // перевод спектральной плотности энергии в дБ
        }
        return res;
}

std::vector<double> Utility::magnitude(const std::vector<std::complex<double>> &t_x) {
        std::vector<double> res(t_x.size());
        for (std::size_t i = 0; i < t_x.size(); ++i) {
                res[i] = std::hypot(t_x[i].real(), t_x[i].imag());
        }
        return res;
}

std::vector<std::complex<double>> Utility::toComplex(const std::vector<double> &t_x) {
        std::vector<std::complex<double>> res(t_x.size());
        for (std::size_t i = 0; i < t_x.size(); ++i) {
                res[i] = std::complex<double>{t_x[i], 0.0};
        }
        return res;
}

std::vector<std::complex<double>>
Utility::FFT(const std::vector<std::complex<double>> &t_data, const int n_threads) {
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<std::complex<double>> window(fftSize), res(fftSize);
        if (!fftw_init_threads()) throw std::runtime_error("Error! Unable to initialize threads!");
        in = fftw_alloc_complex(fftSize);
        out = fftw_alloc_complex(fftSize);
        fftw_plan_with_nthreads(n_threads);
        p = fftw_plan_dft_1d(fftSize, in, out, FFTW_FORWARD, FFTW_MEASURE);

        // Заполнение входных параметров БПФ и окном
        for (int i = 0; i < fftSize; ++i) {
                in[i][0] = t_data[i].real();
                in[i][1] = t_data[i].imag();
        }

        fftw_execute(p);

        for (int i = 0; i < fftSize; i++) {
                res[i].real(out[i][0]);
                res[i].imag(out[i][1]);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup_threads();
        return res;
}

// TODO add choosing win type
std::vector<std::complex<double>>
Utility::winFFT(const std::vector<std::complex<double>> &t_data, const int n_threads) {
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<double> window(fftSize);
        std::vector<std::complex<double>> res(fftSize);
        if (!fftw_init_threads()) throw std::runtime_error("Error! Unable to initialize threads!");
        in = fftw_alloc_complex(fftSize);
        out = fftw_alloc_complex(fftSize);
        fftw_plan_with_nthreads(n_threads);
        p = fftw_plan_dft_1d(fftSize, in, out, FFTW_FORWARD, FFTW_MEASURE);
        // Заполнение окна Блэкмана при alpha = 0.16
        for (int i = 0; i < fftSize; ++i) {
                window[i] = 0.42 + 0.5 * std::cos((2 * std::numbers::pi * i) / (fftSize - 1)) +
                            0.08 * std::cos((4 * std::numbers::pi * i) / (fftSize - 1));
        }

        // Заполнение входных параметров БПФ и окном
        for (int i = 0; i < fftSize; ++i) {
                in[i][0] = t_data[i].real() * window[i];
                in[i][1] = t_data[i].imag() * window[i];
        }

        fftw_execute(p);

        for (int i = 0; i < fftSize; i++) {
                res[i].real(out[i][0]);
                res[i].imag(out[i][1]);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup_threads();
        return res;
}

std::vector<std::complex<double>> Utility::iFFT(const std::vector<std::complex<double>> &t_data, int n_threads) {
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<std::complex<double>> res(fftSize);
        if (!fftw_init_threads()) throw std::runtime_error("Error! Unable to initialize threads!");
        in = fftw_alloc_complex(fftSize);
        out = fftw_alloc_complex(fftSize);
        fftw_plan_with_nthreads(n_threads);
        p = fftw_plan_dft_1d(fftSize, in, out, FFTW_BACKWARD, FFTW_MEASURE);

        // Заполнение входных параметров БПФ и окном
        for (int i = 0; i < fftSize; ++i) {
                in[i][0] = t_data[i].real();
                in[i][1] = t_data[i].imag();
        }

        fftw_execute(p);

        for (int i = 0; i < fftSize; i++) {
                res[i].real(out[i][0]);
                res[i].imag(out[i][1]);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup_threads();
        return res;
}

// phi - phase angle, lev - unwrap to which value, mar = threshold,
void Utility::unwrap(std::vector<double> &phi, double lev, double mar) {
        double a[2] = {0.0, 0.0};
        double d;
        double th;
        int k;
        int flag = 1;

        th = mar * lev;
        while (flag) {
                flag = 0;
                a[0] = a[1] = 0.0;
                for (k = 0; k < phi.size() - 1; k++) {
                        d = phi[k + 1] - phi[k];
                        if (d > th) {
                                a[0] -= lev;
                                flag = 1;
                        }
                        if (d < -th) {
                                a[0] += lev;
                                flag = 1;
                        }
                        phi[k] += a[1];
                        a[1] = a[0];
                }
                phi[phi.size() - 1] += a[1];
        }
}

std::vector<double> Utility::conv(const std::vector<double> &t_x, const std::vector<double> &t_y) {
        if ((t_x.empty()) && (t_y.empty())) {
                return {};
        }

        std::vector<double> a;
        std::vector<double> b;
        if (t_x.size() < t_y.size()) {
                a = t_x;
                b = t_y;
        } else {
                a = t_y;
                b = t_x;
        }

        std::vector<double> result(a.size() + b.size() - 1, 0);
        for (size_t k = 0; k < a.size(); k++) {
                for (size_t l = 0; l < b.size(); l++) {
                        result[l + k] += a[k] * b[l];
                }
        }
        return result;
}

// метод перекрытия с накоплением, t_len должен быть больше минимального размера входных векторов и предпочтительнее степенью двойки
std::vector<std::complex<double>>
Utility::filter(const std::vector<std::complex<double>> &t_x, std::vector<std::complex<double>> t_y,
                const int t_block_len, const int n_threads) {
        if (t_block_len < static_cast<int>(t_y.size()))
                throw std::invalid_argument("Error! Param t_block_len must be more or equal than kih filter size!");
        int M = static_cast<int>(t_y.size());
        int overlap = M - 1;
        int step_size = static_cast<int>(t_block_len) - overlap;
        int pos = 0;
        std::vector<std::complex<double>> res, kih(t_block_len), block(t_block_len);
        t_y.resize(t_block_len, 0);
        // FFT for t_y
        kih = FFT(t_y, n_threads);

        while (pos + t_block_len <= t_x.size()) {
                std::copy(t_x.begin() + pos, t_x.begin() + pos + t_block_len, block.begin());
                block = FFT(block, n_threads);
                for (int i = 0; i < t_block_len; ++i)
                        block[i] *= kih[i];
                block = iFFT(block, n_threads);
                res.resize(pos+step_size);
                std::copy(block.begin() + M, block.end(), res.begin() + pos);
                pos += step_size;
        }

        return res;
}
