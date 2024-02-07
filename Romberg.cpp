#include "Romberg.h"

#include <random>
#include <numeric>
#include <vector>
#include <iostream>


Romberg::Romberg(double _T, double _sig, double _r, double _S0, double _K, int _N, int _R, int _M_R, int _P_R):
    T(_T), sig(_sig), r(_r), S0(_S0), K(_K), N(_N), R(_R), M_R(_M_R), P_R(_P_R), R_M(R * M_R), Npas(_P_R, 0.0),
    erreul(_P_R, 0.0), liceul (_P_R, 0.0), errmil(_P_R, 0.0), licmil(_P_R, 0.0)
{}

void Romberg::calculate3()
{
    // Choix de l'algorithme de génération de nombres aléatoires
    std::random_device rd;
    std::mt19937 gen(time(0));  // Générateur de nombres aléatoires Mersenne Twister
    std::normal_distribution<double> dist(0.0, 1.0);  // Distribution normale

    for (int i = 0; i < P_R; ++i) {
        double dt = T / N;
        double sigdt = sig * std::sqrt(dt);
        double rdt = r * dt;
        double der = rdt - std::pow(sigdt, 2) / 2;

        double payeul = 0.0;
        double careul = 0.0;
        double paymil = 0.0;
        double carmil = 0.0;

        for (int j = 0; j < R; ++j) {
            std::vector<double> S(M_R, S0);
            std::vector<double> Se(M_R, S0);
            std::vector<double> S2e(M_R, S0);
            std::vector<double> Sm(M_R, S0);
            std::vector<double> S2m(M_R, S0);

            for (int k = 0; k < N; ++k) {
                std::vector<double> g1(M_R), g2(M_R), g(M_R);
                std::default_random_engine gen;
                std::normal_distribution<double> dist(0.0, 1.0);


                for (int m = 0; m < M_R; ++m) {
                    g1[m] = dist(gen) / std::sqrt(2);
                    g2[m] = dist(gen) / std::sqrt(2);
                    g[m] = g1[m] + g2[m];


                    // formule de récurrence
                    S[m] = S[m] * std::exp(sigdt * g[m] + der); // inchangé
                    Se[m] = Se[m] * (1 + sigdt * g[m] + rdt);   // inchangé
                    S2e[m] = S2e[m] * (1 + sigdt * g1[m] + rdt / 2) * (1 + sigdt * g2[m] + rdt / 2);
                    Sm[m] = Sm[m] * (1 + sigdt * g[m] + std::pow(sigdt * g[m], 2) / 2 + der);  // inchangé
                    S2m[m] = S2m[m] * (1 + sigdt * g1[m] + std::pow(sigdt * g1[m], 2) / 2 + der / 2) * (1 + sigdt * g2[m] + std::pow(sigdt * g2[m], 2) / 2 + der / 2);
                }
            }
            std::vector<double> paymc_vec(M_R, 0.0);

            std::vector<double> payeul_vec(M_R, 0.0);
            std::vector<double> payeul2_vec(M_R, 0.0);

            std::vector<double> paymil_vec(M_R, 0.0);
            std::vector<double> paymil2_vec(M_R, 0.0);

            std::vector<double> eul_square_vec(M_R, 0.0);
            std::vector<double> mil_square_vec(M_R, 0.0);

            for (int m = 0; m < M_R; ++m) {
                paymc_vec.at(m) = std::max(double(0.0), K - S[m]);
                payeul2_vec.at(m) = std::max(double(0.0), K - S2e[m]);
                payeul_vec.at(m) = std::max(double(0.0), K - Se[m]);
                paymil2_vec.at(m) = std::max(double(0.0), K - S2m[m]);
                paymil_vec.at(m) = std::max(double(0.0), K - Sm[m]);

                eul_square_vec.at(m) = pow(2 * payeul2_vec.at(m) - payeul_vec.at(m) - paymc_vec.at(m), 2);
                mil_square_vec.at(m) = pow(2 * paymil2_vec.at(m) - paymil_vec.at(m) - paymc_vec.at(m), 2);


            }

            double paymc_acc(std::accumulate(paymc_vec.begin(), paymc_vec.end(), 0.0));
            double payeul_acc(std::accumulate(payeul_vec.begin(), payeul_vec.end(), 0.0));
            double payeul2_acc(std::accumulate(payeul2_vec.begin(), payeul2_vec.end(), 0.0));
            double paymil_acc(std::accumulate(paymil_vec.begin(), paymil_vec.end(), 0.0));
            double paymil2_acc(std::accumulate(paymil2_vec.begin(), paymil2_vec.end(), 0.0));

            payeul += (2 * payeul2_acc - payeul_acc - paymc_acc) / M_R;
            paymil += (2 * paymil2_acc - paymil_acc - paymc_acc) / M_R;

            careul += std::accumulate(eul_square_vec.begin(), eul_square_vec.end(), 0.0) / M_R;
            carmil += std::accumulate(mil_square_vec.begin(), mil_square_vec.end(), 0.0) / M_R;
        }

        erreul[i] = std::exp(-r * T) * payeul / R;
        liceul[i] = 1.96 * std::exp(-r * T) * std::sqrt((careul / R - std::pow(payeul / R, 2)) / R_M);

        errmil[i] = std::exp(-r * T) * paymil / R;
        licmil[i] = 1.96 * std::exp(-r * T) * std::sqrt((carmil / R - std::pow(paymil / R, 2)) / R_M);

        Npas[i] = N;
        N = N * 2;

    }
}


void Romberg::displayResults3()
{   // Display the results

    std::cout << "Romberg Euler" << std::endl;
    for (int i = 0; i < P_R; ++i) {
        std::cout << erreul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Euler" << std::endl;
    for (int i = 0; i < P_R; ++i) {
        std::cout << liceul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Romberg Milshtein" << std::endl;
    for (int i = 0; i < P_R; ++i) {
        std::cout << errmil[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Milshtein" << std::endl;
    for (int i = 0; i < P_R; ++i) {
        std::cout << licmil[i] << " ";
    }
    std::cout << std::endl;
}

