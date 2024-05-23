#include <iostream>
#include <vector>
#include <cmath>
#include <matplot/matplot.h>
#include <numeric>

const int Fs = 1000; // 1 КГц частота дискретизации
const int Fc = 50; // 50 Гц Несущая частота
const double Fc1 = 5; // 5 Гц Частота модулирующего сигнала
const double M = 0.698; // Индекс модуляции
const double Amplitude = 127;

void unwrap_array(const std::vector<double> &in, std::vector<double> out) {
        out[0] = in[0];
        for (int i = 1; i < in.size(); i++) {
                double d = in[i] - in[i - 1];
                d = d > std::numbers::pi ? d - 2 * std::numbers::pi : (d < -std::numbers::pi ? d + 2 * std::numbers::pi
                                                                                             : d);
                out[i] = out[i - 1] + d;
        }
}

int main() {
        using namespace matplot;
        std::vector<double> t(Fs); // Временные отсчёты
        std::vector<double> Sm(Fs); // Модулирующие отсчёты
        std::vector<double> Sn(Fs); // Несущие отсчёты
        std::vector<double> I(Fs); // Синфазная составляющая
        std::vector<double> Q(Fs); // Квадратурная составляющая
        std::vector<double> S(Fs); // Сигнальные отсчёты
        std::vector<double> TheorS(Fs); // Теоретический выходной сигнал

        // Генерация временных отсчётов
        std::iota(t.begin(), t.end(), 0); // [0,1,2....]
        for (auto &i: t) i /= Fs;

        // Генерация сигнала
        for (int i = 0; i < Fs; ++i) {
                Sm[i] = std::sin(2 * std::numbers::pi * Fc1 * t[i]);
                Sn[i] = Amplitude * std::sin(2 * std::numbers::pi * Fc * t[i]);
                I[i] = Amplitude * std::cos(M * Sm[i]);
                Q[i] = -Amplitude * std::sin(M * Sm[i]);
                S[i] = I[i] * std::cos(2 * std::numbers::pi * Fc * t[i]) + Q[i] * std::sin(2 * std::numbers::pi * Fc * t[i]);
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
        save(f, "../img/subplot2.svg");

        return 0;
}