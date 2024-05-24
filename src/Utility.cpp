//
// Created by kirill3266 on 24.05.24.
//

#include "Utility.h"
#include <algorithm>

void Utility::goertzel_cmplx(std::vector<std::complex<double>>& in, std::vector<double> f, std::vector<std::complex<double>>& out)
{

        int m, p;
        std::complex<double> w;
        double alpha;
        std::complex<double> v[3];

        for(p = 0; p < f.size(); p++)
        {
                w.real(cos(2 * std::numbers::pi * f[p] / static_cast<double>(in.size())));
                w.imag(sin(2 * std::numbers::pi * f[p] / static_cast<double>(in.size())));

                alpha = 2.0 * std::real(w);
                v[0] = {0.0,0.0};
                v[1] = {0.0,0.0};
                v[2] = {0.0,0.0};

                for(m = 0; m < in.size(); m++)
                {
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
double Utility::logPower(fftw_complex in, double scale) {
        double re = in[0] * scale;
        double im = in[1] * scale;
        double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
        return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

// Вычисление спектральной плотности энергии и перевод её в дб
double Utility::logPower(double I,double Q, double scale) {
        double re = I * scale;
        double im = Q * scale;
        double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
        return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

std::vector<double> Utility::FFT(const std::vector<double> &t_data, const int t_Fs) {
        // Объявление переменных
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<double> window(fftSize);
        std::vector<double> pwr(fftSize / 2);
        in = fftw_alloc_complex(fftSize);
        out = fftw_alloc_complex(fftSize);
        p = fftw_plan_dft_1d(fftSize, in, out, FFTW_FORWARD, FFTW_MEASURE);

        // Заполнение окна Блэкмана при alpha = 0.16
        for (int i = 0; i < fftSize; ++i) {
                window[i] = 0.42 + 0.5 * std::cos((2 * M_PI * i) / (fftSize - 1)) +
                            0.08 * std::cos((4 * M_PI * i) / (fftSize - 1));
        }

        // Заполнение входных параметров БПФ
        for (int i = 0; i < fftSize; ++i) {
                in[i][0] = t_data[i] * window[i];
        }

        fftw_execute(p);

        // Запись только положительной части симметричного спектра
        for (int i = 0; i < fftSize / 2; i++) {
                pwr[i] = logPower(out[i], 1.0 / fftSize);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        return pwr;
}

std::vector<double> Utility::FFT(std::vector<std::complex<double>> &t_data, int t_Fs, const int n_threads) {
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<double> window(fftSize), res(fftSize);
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
                res[i] = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup_threads();
        return res;
}

// phi - phase angle, lev - unwrap to which value, mar = threshold,
void Utility::unwrap(std::vector<double>& phi, double lev, double mar)
{
        double a[2] = {0.0, 0.0};
        double d;
        double th;
        int k;
        int flag = 1;

        th = mar*lev;
        while(flag)
        {
                flag = 0;
                a[0] = a[1] = 0.0;
                for(k = 0; k<phi.size()-1; k++)
                {
                        d = phi[k+1] - phi[k];
                        if( d > th)
                        {
                                a[0] -= lev;
                                flag = 1;
                        }
                        if( d < -th)
                        {
                                a[0] += lev;
                                flag = 1;
                        }
                        phi[k]+=a[1];
                        a[1] = a[0];
                }
                phi[phi.size()-1]+=a[1];
        }
}
