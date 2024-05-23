#include <iostream>
#include <vector>
#include <cmath>
#include <matplot/matplot.h>
#include <numeric>
#include <fftw3.h>

const int Fs = 1000; // 1 КГц частота дискретизации
const int Fc = 50; // 50 Гц Несущая частота
const double Fc1 = 10; // 5 Гц Частота модулирующего сигнала
const double M = 0.698; // Индекс модуляции
const double Amplitude = 127;


double diff(double x, double y) {
        double h = 1e-10;
        return (x - y) / h;
}

void goertzel_cmplx(std::vector<std::complex<double>>& in, std::vector<double> f, std::vector<std::complex<double>>& out)
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
double logPower(fftw_complex in, double scale) {
        double re = in[0] * scale;
        double im = in[1] * scale;
        double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
        return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

// Вычисление спектральной плотности энергии и перевод её в дб
double logPower(double I,double Q, double scale) {
        double re = I * scale;
        double im = Q * scale;
        double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
        return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

std::vector<double> FFT(std::vector<std::complex<double>> &t_data, int t_Fs, const int n_threads) {
        fftw_complex *in, *out;
        fftw_plan p;
        int fftSize = static_cast<int>(t_data.size());
        std::vector<double> window(fftSize), res(fftSize / 2);
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

        for (int i = 0; i < fftSize / 2; i++) {
                res[i] = std::sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);
        }

        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
        fftw_cleanup_threads();
        return res;
}

int main() {
        using namespace matplot;
        std::vector<double> t(Fs); // Временные отсчёты
        std::vector<double> Sm(Fs); // Модулирующие отсчёты
        std::vector<double> Sn(Fs); // Несущие отсчёты
        std::vector<double> I(Fs); // Синфазная составляющая
        std::vector<double> Q(Fs); // Квадратурная составляющая
        std::vector<double> S(Fs); // Сигнальные отсчёты
        std::vector<std::complex<double>> complexS(Fs);
        std::vector<double> diffI(Fs), diffQ(Fs), diffS(Fs);
        std::vector<double> TheorS(Fs); // Теоретический выходной сигнал
        std::vector<double> Demod(Fs); // Демодулированный сигнал
        std::vector<double> Unwrap(Fs); // Unwrapped демодулированный сигнал
        std::complex<double> j{0.0,1.0};
        std::complex<double> cmp;
        std::vector<double> F;
        std::vector<double> spectrum;

        // Генерация временных отсчётов
        std::iota(t.begin(), t.end(), 0); // [0,1,2....]
        for (auto &i: t) i /= Fs;

        // Генерация сигнала
        for (int i = 0; i < Fs; ++i) {
                Sm[i] = std::sin(2 * std::numbers::pi * Fc1 * t[i]);
                Sn[i] = Amplitude * std::sin(2 * std::numbers::pi * Fc * t[i]);
                cmp = std::exp(j*M*Sm[i]);
                I[i] = Amplitude*std::real(cmp); // or I[i] = Amplitude * std::cos(M * Sm[i]);
                Q[i] = Amplitude * std::imag(cmp); // or Q[i] = -Amplitude * std::sin(M * Sm[i]);
                S[i] = I[i] * std::cos(2 * std::numbers::pi * Fc * t[i]) +
                       Q[i] * std::sin(2 * std::numbers::pi * Fc * t[i]);
                TheorS[i] = Amplitude * std::cos(2 * std::numbers::pi * Fc * t[i] + M * Sm[i]);
        }

        // Построение исходного сигнала
        auto f = figure(true);
        f->size(1920, 1080);

        auto ax = axes(f);
        auto ax1 = subplot(f, 4, 1, 0);
        plot(ax1, t, Sm);
        title(ax1, "График модулирующего сигнала");
        auto ax2 = subplot(f, 4, 1, 1);
        plot(ax2, t, Sn);
        title(ax2, "График несущего сигнала");
        auto ax3 = subplot(f, 4, 1, 2);
        plot(ax3, t, TheorS);
        title(ax3, "График теоретического сигнала");
        auto ax4 = subplot(f, 4, 1, 3);
        plot(ax4, t, S);
        title(ax4, "График итогового сигнала");
        show(f);
        save(f, "../img/modulation.svg");


        // Демодуляция
        for (int i = 0; i < Fs; ++i) {

                I[i] = std::cos(2 * std::numbers::pi * Fc * t[i]) * S[i] -
                       0.5 * std::sin(2 * 2 * std::numbers::pi * Fc * t[i] + M * Sm[i]);
                Q[i] = -std::sin(2 * std::numbers::pi * Fc * t[i]) * S[i] -
                       0.5 * std::cos(2 * 2 * std::numbers::pi * Fc * t[i] + M * Sm[i]);
                Demod[i] = std::atan(-1 * I[i] / Q[i]);
                complexS[i] = S[i]*std::exp(j*2.0*std::numbers::pi*static_cast<double>(Fc)*t[i]);
        }

        spectrum = FFT(complexS,Fs,6);
        F.resize(spectrum.size());
        for (int i = 0; i < F.size(); ++i) {
                F[i] = static_cast<double>(i) * static_cast<double>(Fs) / (static_cast<double>(F.size() * 2));
        }

        auto f1 = figure(true);
        f1->size(1920,1080);

        auto axSpectrum = axes(f1);

        plot(axSpectrum,F,spectrum);
        title(axSpectrum, "Положительный спектр принятого на демодулятор сигнала");
        show(f1);

        std::vector<std::complex<double>> out(2);
        goertzel_cmplx(complexS,std::vector<double>{10,20},out);
        std::cout << logPower(out[0].real(),out[0].imag(),1) << " " << logPower(out[1].real(),out[1].imag(),1) << std::endl;

        save(f1,"../img/spectrum.svg");

//        for (int i = 0; i < Fs; ++i) {
//                if (i == 0) {
//                        diffI[i] = diff(I[i], 0);
//                        diffQ[i] = diff(Q[i], 0);
//                }
//                diffI[i] = diff(I[i], I[i - 1]);
//                diffQ[i] = diff(Q[i], Q[i - 1]);
//                Unwrap[i] = (diffI[i] * Q[i] - diffQ[i] * I[i]) / (I[i] * I[i] + Q[i] * Q[i]);
//        }

//        // Построение демодулированного сигнала
//        auto f1 = figure(true);
//        f1->size(1920, 1080);
//
//        auto axDemod = axes(f1);
//        auto axDemod1 = subplot(f1, 3, 1, 0);
//        plot(axDemod1, t, Sm);
//        title(axDemod1, "График модулирующего сигнала");
//        auto axDemod2 = subplot(f1, 3, 1, 1);
//        plot(axDemod2, t, Demod);
//        title(axDemod2, "График демодулированного сигнала");
//        auto axDemod3 = subplot(f1, 3, 1, 2);
//        plot(axDemod3, t, Unwrap);
//        title(axDemod3, "График unwrapped сигнала");
//        show(f1);
//        save(f1, "../img/demodulation.svg");

        return 0;
}