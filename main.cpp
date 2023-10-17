#include <iostream>
#include <vector>
#include <cmath>
#include <matplot/matplot.h>
#include "fftw3.h"

const int Fs = 192000; // 192000 Гц частота дискретизации
const int Fc = 66667; // Несущая частота
const double Fc1 = 100; // Частота модулирующего сигнала
const int N = 50000; // Количество отсчётов

// Вычисление спектральной плотности энергии и перевод её в дб
double logPower(fftw_complex in, double scale) {
    double re = in[0] * scale;
    double im = in[1] * scale;
    double magsq = re * re + im * im; // квадрат амплитудного спектра или спектральная плотность
    return (std::log(magsq) * 10.0 / std::log(10.0)); // перевод спектральной плотности энергии в дБ
}

double Goertzel(std::vector<double> t_data, double t_freq) {
    int k, i;
    double Sn, Sn1 = 0, Sn2 = 0;

    k = std::ceil((static_cast<double>(t_data.size()) * t_freq) / Fs);

    for (i = 0; i < t_data.size(); i++) {
        Sn = 2 * std::cos((2.0 * M_PI * k) / static_cast<double>(t_data.size())) * Sn1 - Sn2 + t_data[i];
        Sn2 = Sn1;
        Sn1 = Sn;
    }
    return std::sqrt(
            Sn1 * Sn1 + Sn2 * Sn2 - 2 * std::cos((2.0 * M_PI * k) / static_cast<double>(t_data.size())) * Sn1 * Sn2);
}

std::vector<double> Filter(const std::vector<double> &in) {
    const int filter_length = 1024; // Длина фильтра
    const double Fp = 500; // Частота полосы пропускания
    const double Fd = 550; // Частота полосы затухания

    double H[filter_length] = {0}; // Импульсная характеристика фильтра
    double H_id[filter_length] = {0}; // Идеальная импульсная характеристика
    double W[filter_length] = {0}; // Весовая функция

    // Расчёт импульсной характеристики фильтра
    double Ir = (Fp + Fd) / (2 * Fs);

    for (int i = 0; i < filter_length; ++i) {
        if (!i) H_id[i] = 2 * M_PI * Ir;
        else H_id[i] = std::sin(2 * M_PI * Ir * i) / (M_PI * i);

        // Весовая функция Блекмена
        W[i] = 0.42 + 0.5 * std::cos((2 * M_PI * i) / (filter_length - 1)) +
               0.08 * std::cos((4 * M_PI * i) / (filter_length - 1));
        H[i] = H_id[i] * W[i];
    }

    // Нормировка импульсной характеристики
    double sum = 0;
    for (double i: H) sum += i;
    for (double &i: H) i /= sum;

    std::vector<double> res(in.size());

    // Фильтрация (линейная свёртка)
    for (int i = 0; i < in.size(); ++i) {
        res[i] = 0;
        for (int j = 0; j < filter_length - 1; ++j)
            if (i - j >= 0)
                res[i] += H[j] * in[i - j];
    }
    return res;
}

int main() {
    // Объявление переменных
    fftw_complex *in, *out;
    fftw_plan p;
    in = fftw_alloc_complex(N);
    out = fftw_alloc_complex(N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);
    std::vector<double> t(N); // Временные отсчёты
    std::vector<double> x(N / 2); // Частотные отсчёты
    std::vector<double> S(N); // Сигнальные отсчёты
    std::vector<double> tmp(N / 2); // Сигнальные отсчёты
    std::vector<double> S1(N); // Отсчёты модулирующего сигнала
    std::vector<double> S2(N);
    std::vector<double> S3(N);
    std::vector<double> demod(N); // Сигнальные отсчёты демодулированного сигнала

    // Генерация временных отсчётов
    std::iota(t.begin(), t.end(), 0); // [0,1,2....]
    for (auto &i: t) i /= Fs;

    // Генерация сигнала
    for (int i = 0; i < N; ++i) {
        S1[i] = std::cos(2 * M_PI * Fc1 * t[i]);
        S[i] = std::cos(2 * M_PI * Fc * t[i] + 0.698 * S1[i]);
    }

    // Построение исходного сигнала
    matplot::figure();
    matplot::title("График сигнала");
    matplot::xlabel("Время, с");
    matplot::ylabel("Амплитуда");
    matplot::plot(t, S);

    // Заполнение входных параметров БПФ
    for (int i = 0; i < N; ++i) {
        in[i][0] = S[i] * (0.53836 - 0.46164 * std::cos(2 * M_PI / (N - 1)));
    }
    std::copy(S.begin(), S.end(), demod.begin());

    // БПФ
    fftw_execute(p);

    // Копирование выходных частотных отсчётов БПФ
    for (int i = 0; i < N; ++i) {
        S[i] = logPower(out[i], 1.0 / N);
    }

    // Значение частотных отсчётов бпф
    for (int i = 0; i < N / 2; ++i) {
        x[i] = static_cast<double>(i) * Fs / N;
    }
    std::copy(S.begin(), S.begin() + S.size() / 2, tmp.begin());

    // Построение спектра исходного сигнала
    matplot::figure();
    matplot::title("Спектр исходного сигнала");
    matplot::xlabel("Частота, гц");
    matplot::ylabel("Амплитуда, дб");
    matplot::plot(x, tmp);

    // Демодулирование путём умножения на несущий сигнал
    for (int i = 0; i < N; ++i) {
        S2[i] = demod[i] * std::cos(2 * M_PI * Fc * t[i]);
        S3[i] = demod[i] * std::sin(2 * M_PI * Fc * t[i]);
    }

    S2 = Filter(S2);
    S3 = Filter(S3);

    std::vector<double> res(N);
    for (int i = 0; i < N; ++i) {
        res[i] = std::atan(S3[i] / S2[i]) / 0.698;
    }

    // Построение модулирующего сигнала
    matplot::figure();
    matplot::title("Модулирующий сигнал");
    matplot::xlabel("Частота, гц");
    matplot::ylabel("Амплитуда");
    matplot::plot(t, res);

    // Заполнение входных параметров БПФ с применением окна
    for (int i = 0; i < N; ++i) {
        in[i][0] = res[i] * (0.53836 - 0.46164 * std::cos(2 * M_PI / (N - 1)));
    }
    std::copy(res.begin(), res.end(), demod.begin());

    // БПФ
    fftw_execute(p);

    // Копирование выходных частотных отсчётов БПФ
    for (int i = 0; i < N; ++i) {
        res[i] = logPower(out[i], 1.0 / N);
    }
    std::copy(res.begin(), res.begin() + res.size() / 2, tmp.begin());

    // Построение спектра модулирующего сигнала
    matplot::figure();
    matplot::title("Спектр модулирующего сигнала");
    matplot::xlabel("Частота, гц");
    matplot::ylabel("Амплитуда, дб");
    matplot::plot(x, tmp);

    std::cout << "Значения алгоритма Герцеля на частоте 100 гц: " << Goertzel(demod, 100) << std::endl;
    std::cout << "Значения алгоритма Герцеля на частоте 312.5 гц: " << Goertzel(demod, 312) << std::endl;
    std::cout << "Значения БПФ на частоте 100 гц: " << res[26] << std::endl; // 100*N(50000)/Fs(192000)
    std::cout << "Значения БПФ на частоте 312.5 гц: " << res[81] << std::endl; // 312.5*N(50000)/Fs(192000)

    matplot::show();
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    return 0;
}