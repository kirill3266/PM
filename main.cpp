#include "Utility.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <matplot/matplot.h>

const int Fs = 1000; // 1 КГц частота дискретизации
const int Fc = 50; // 50 Гц Несущая частота
const double Fc1 = 30; // 5 Гц Частота модулирующего сигнала
const double M = 2.0; // Индекс модуляции
const double Amplitude = 127; // Амплитуда

int main() {
        using namespace matplot;
        std::vector<double> t(Fs); // Временные отсчёты
        std::vector<double> Sm(Fs); // Модулирующие отсчёты
        std::vector<double> Sn(Fs); // Несущие отсчёты
        std::vector<double> I(Fs); // Синфазная составляющая
        std::vector<double> Q(Fs); // Квадратурная составляющая
        std::vector<double> S(Fs); // Сигнальные отсчёты
        std::vector<std::complex<double>> complexS(Fs);
        std::vector<double> TheorS(Fs); // Теоретический выходной сигнал
        std::vector<double> Demod(Fs); // Демодулированный сигнал
        std::vector<double> demodI(Fs); // Демодулированная синфазная составляющая
        std::vector<double> demodQ(Fs); // Демодулированная Квадратурная составляющая
        std::vector<double> Unwrap(Fs); // Unwrapped сигнал
        std::complex<double> j{0.0, 1.0}; // комплексная единица
        std::complex<double> cmp; // временная комплексная перемення
        std::vector<double> F; // Вектор частот спектра
        std::vector<double> spectrum; // Спектр сигнала

        // Генерация временных отсчётов
        std::iota(t.begin(), t.end(), 0); // [0,1,2....]
        for (auto &i: t) i /= Fs;

        // Генерация сигнала
        for (int i = 0; i < Fs; ++i) {
                Sm[i] = std::sin(2 * std::numbers::pi * Fc1 * t[i]);
                Sn[i] = Amplitude * std::sin(2 * std::numbers::pi * Fc * t[i]);
                cmp = std::exp(j * M * Sm[i]);
                I[i] = Amplitude * std::real(cmp); // or I[i] = Amplitude * std::sin(M * Sm[i]);
                Q[i] = Amplitude * std::imag(cmp); // or Q[i] = -Amplitude * std::cos(M * Sm[i]);
                S[i] = I[i] * std::sin(2 * std::numbers::pi * Fc * t[i]) +
                       Q[i] * std::cos(2 * std::numbers::pi * Fc * t[i]);
                TheorS[i] = Amplitude * std::sin(2 * std::numbers::pi * Fc * t[i] + M * Sm[i]);
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

        // TODO добавить фильтрацию
        // Демодуляция
        for (int i = 0; i < Fs; ++i) {
                demodI[i] = std::cos(2 * std::numbers::pi * Fc * t[i]) * S[i];
                demodQ[i] = -std::sin(2 * std::numbers::pi * Fc * t[i]) * S[i];
                complexS[i] =
                        S[i] * std::exp(j * 2.0 * std::numbers::pi * static_cast<double>(Fc) * t[i]);
                Demod[i] = std::arg(complexS[i]);
                Unwrap[i] = std::arg(complexS[i]);
        }

        Utility::unwrap(Unwrap, 2 * std::numbers::pi, 0.8);

        // Построение исходного сигнала
        auto f1 = figure(true);
        f1->size(1920, 1080);

        auto axDemod = axes(f1);
        auto axDemod1 = subplot(f1, 3, 1, 0);
        plot(axDemod1, t, Sm);
        title(axDemod1, "График модулирующего сигнала");
        auto axDemod2 = subplot(f1, 3, 1, 1);
        plot(axDemod2, t, Demod);
        title(axDemod2, "График демодулированного сигнала");
        auto axDemod3 = subplot(f1, 3, 1, 2);
        plot(axDemod3, t, Unwrap);
        title(axDemod3, "График unwrapped сигнала");
        show(f1);
        save(f1, "../img/demodulation.svg");


        spectrum = Utility::FFT(Demod, Fs);
        F.resize(spectrum.size());
        // Move spectrum
        for (int i = 0; i < F.size(); ++i) {
                F[i] = static_cast<double>(i) * static_cast<double>(Fs) / (static_cast<double>(F.size()));
        }

//        spectrum = Utility::FFT(complexS, Fs, 6);
//        F.resize(spectrum.size());
//        // Move spectrum
//        for (int i = 0; i < F.size(); ++i) {
//                if (i < Fs / 2)
//                        F[i] = static_cast<double>(i) * static_cast<double>(Fs) / (static_cast<double>(F.size()));
//                else
//                        F[i] = -static_cast<double>(F.size() - i) * static_cast<double>(Fs) /
//                               (static_cast<double>(F.size()));
//        }

        // Герцель
        std::vector<std::complex<double>> out(2);
        Utility::goertzel_cmplx(complexS, std::vector<double>{10, 30}, out);
        std::cout << Utility::logPower(out[0].real(), out[0].imag(), 1) << " "
                  << Utility::logPower(out[1].real(), out[1].imag(), 1)
                  << std::endl;

        auto f2 = figure(true);
        f2->size(1920, 1080);

        auto axSpectrum = axes(f2);

        plot(axSpectrum, F, spectrum);
        title(axSpectrum, "Положительный спектр принятого на демодулятор сигнала");
        show(f2);

        save(f2, "../img/spectrum.svg");

        return 0;
}