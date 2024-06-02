#include "Utility.h"
#include "kih5.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <matplot/matplot.h>

const int Fs = 1000; // 200 КГц частота дискретизации
const int N = 5 * Fs; // число генерируемых отсчётов
const int Fc = 50; // 66,(6) Гц Несущая частота
const double Fc1 = 5; // 100 Гц Частота модулирующего сигнала
const double M = 5; // Индекс модуляции
const double Amplitude = 127; // Амплитуда

int main() {
        using namespace matplot;
        std::vector<double> t(N), Sm(N), Sn(N), I(N), Q(N), S(N), TheorS(N); // Модулированные сигналы
        std::vector<double> Demod(N + kih.size() - 1), Unwrap(N + kih.size() - 1);
        std::vector<double> unfilteredS(N), filteredS(N);
        std::vector<std::complex<double>> cmpS(N), demodCmpS(N), filteredCmpS(N), cmpKih(kih.size());
        std::vector<double> spectrumCmpS, spectrumDemodCmpS, spectrumKih, spectrumFilteredCmpS; // Спектр сигнала
        std::complex<double> j{0.0, 1.0}; // комплексная единица
        std::vector<double> F, demodF, filteredF, kihF; // Вектор частот спектра

        // Генерация временных отсчётов
        std::iota(t.begin(), t.end(), 0); // [0,1,2....]
        for (auto &i: t) i /= Fs;

        // Генерация сигнала
        for (int i = 0; i < N; ++i) {
                Sm[i] = std::sin(2 * std::numbers::pi * Fc1 * t[i]); // Модулирующий сигнал
                Sn[i] = Amplitude * std::sin(2 * std::numbers::pi * Fc * t[i]); // Несущий сигнал
                cmpS[i] = std::exp(j * M * Sm[i]);
                I[i] = Amplitude * std::cos(M * Sm[i]); // or I[i] = Amplitude * std::real(cmpS[i]);
                Q[i] = Amplitude * std::sin(M * Sm[i]); // or Q[i] = Amplitude * std::imag(cmpS[i]);
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

        // Демодуляция
        for (int i = 0; i < N; ++i) {
                demodCmpS[i] = std::exp(j * 2.0 * std::numbers::pi * static_cast<double>(Fc) * t[i]) * S[i];
                unfilteredS[i] = S[i] * std::cos(2.0 * std::numbers::pi * static_cast<double>(Fc) * t[i]) +
                                 S[i] * std::sin(2.0 * std::numbers::pi * static_cast<double>(Fc) * t[i]);
//                demodCmpS[i].real(std::cos(2 * std::numbers::pi * Fc * t[i]) * S[i]);
//                demodCmpS[i].imag(std::sin(2 * std::numbers::pi * Fc * t[i]) * S[i]);
        }

        // Фильтрация
//        filteredS = Utility::conv(S, kih);
//        filteredCmpS = Utility::toComplex(filteredS);
        filteredCmpS = Utility::filter(demodCmpS, Utility::toComplex(kih), 512, 6);
        // TODO добавлять нули в начале а не в конце
        filteredCmpS.resize(demodCmpS.size());
        for (int i = 0; i < filteredCmpS.size(); ++i) {
                Demod[i] = std::arg(filteredCmpS[i]);
        }
        std::copy(Demod.begin(), Demod.end(), Unwrap.begin());
        Utility::unwrap(Unwrap, 2 * std::numbers::pi, 0.8);

        // Построение исходного сигнала
        auto f1 = figure(true);
        f1->size(1920, 1080);

        auto axDemod = axes(f1);
        auto axDemod1 = subplot(f1, 4, 1, 0);
        plot(axDemod1, t, Sm);
        title(axDemod1, "График модулирующего сигнала");
        auto axDemod2 = subplot(f1, 4, 1, 1);
        plot(axDemod2, t, unfilteredS);
        title(axDemod2, "График принятого сигнала до фильтрации");
        auto axDemod3 = subplot(f1, 4, 1, 2);
        plot(axDemod3, t, Demod);
        title(axDemod3, "График принятого сигнала после фильтрации");
        auto axDemod4 = subplot(f1, 4, 1, 3);
        plot(axDemod4, t, Unwrap);
        title(axDemod4, "График unwrapped сигнала после фильтрации");
        show(f1);
        save(f1, "../img/demodulation.svg");

        spectrumCmpS = Utility::magnitude(Utility::winFFT(cmpS, 6));
        F.resize(spectrumCmpS.size());
        // Move spectrum
        for (int i = 0; i < F.size(); ++i) {
                if (i < F.size() / 2)
                        F[i] = static_cast<double>(i) * static_cast<double>(Fs) / (static_cast<double>(F.size()));
                else
                        F[i] = -static_cast<double>(F.size() - i) * static_cast<double>(Fs) /
                               (static_cast<double>(F.size()));
        }

        spectrumDemodCmpS = Utility::magnitude(Utility::winFFT(demodCmpS, 6));
        demodF.resize(spectrumDemodCmpS.size());
        for (int i = 0; i < demodF.size(); ++i) {
                if (i < demodF.size() / 2)
                        demodF[i] = static_cast<double>(i) * static_cast<double>(Fs) /
                                    static_cast<double>(demodF.size());
                else
                        demodF[i] = -static_cast<double>(demodF.size() - i) * static_cast<double>(Fs) /
                                    (static_cast<double>(demodF.size()));
        }

        for (int i = 0; i < kih.size(); ++i) {
                cmpKih[i] = std::complex<double>{kih[i], 0};
        }

        spectrumKih = Utility::magnitude(Utility::FFT(cmpKih, 6));
        kihF.resize(spectrumKih.size());
        for (int i = 0; i < kihF.size(); ++i) {
                if (i < kihF.size() / 2)
                        kihF[i] = static_cast<double>(i) * static_cast<double>(Fs) /
                                  static_cast<double>(kihF.size());
                else
                        kihF[i] = -static_cast<double>(kih.size() - i) * static_cast<double>(Fs) /
                                  (static_cast<double>(kihF.size()));
        }

        spectrumFilteredCmpS = Utility::magnitude(Utility::FFT(filteredCmpS, 6));
        filteredF.resize(spectrumFilteredCmpS.size());
        // Move spectrum
        for (int i = 0; i < filteredF.size(); ++i) {
                if (i < filteredF.size() / 2)
                        filteredF[i] = static_cast<double>(i) * static_cast<double>(Fs) /
                                       (static_cast<double>(filteredF.size()));
                else
                        filteredF[i] = -static_cast<double>(filteredF.size() - i) * static_cast<double>(Fs) /
                                       static_cast<double>(filteredF.size());
        }

        auto f2 = figure(true);
        f2->size(1920, 1080);

        auto axSpectrum = axes(f2);

        auto axSpectrum1 = subplot(f2, 4, 1, 0);
        plot(axSpectrum1, F, spectrumCmpS);
        title(axSpectrum1, "Cпектр исходного сигнала");

        auto axSpectrum2 = subplot(f2, 4, 1, 1);
        plot(axSpectrum2, demodF, spectrumDemodCmpS);
        title(axSpectrum2, "Cпектр принятого сигнала");

        auto axSpectrum3 = subplot(f2, 4, 1, 2);
        plot(axSpectrum3, kihF, spectrumKih);
        title(axSpectrum3, "Cпектр фильтра");

        auto axSpectrum4 = subplot(f2, 4, 1, 3);
        plot(axSpectrum4, filteredF, spectrumFilteredCmpS);
        title(axSpectrum4, "Cпектр фильтрованного сигнала");

        show(f2);
        save(f2, "../img/spectrum.svg");

        // Герцель
        std::vector<std::complex<double>> out(7);
        Utility::goertzel_cmplx(filteredCmpS, std::vector<double>{30, 35, 40, 45, 50, 55, 60}, out);
        std::cout << Utility::logPower(out[0], 1.0 / Fs) << " "
                  << Utility::logPower(out[1], 1.0 / Fs) << " "
                  << Utility::logPower(out[2], 1.0 / Fs) << " "
                  << Utility::logPower(out[3], 1.0 / Fs) << " "
                  << Utility::logPower(out[4], 1.0 / Fs) << " "
                  << Utility::logPower(out[5], 1.0 / Fs) << " "
                  << Utility::logPower(out[6], 1.0 / Fs) << " "
                  << std::endl;

        return 0;
}