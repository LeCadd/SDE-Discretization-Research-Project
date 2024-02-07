#include "Cvg_forte.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>  // Pour time(0)
#include <random>


Cvg_forte::Cvg_forte(double _T, double _sig, double _r, double _S0, int _N, int _M, int _P):
    T(_T), sig(_sig), r(_r), S0(_S0), N(_N), M(_M), P(_P), Npas(_P, 0),
    erreul(_P, 0.0), liceul(_P, 0.0), errmil(_P, 0.0), licmil(_P, 0.0)
{}

void Cvg_forte::calculate() {

    std::mt19937 gen(time(0));  // Générateur de nombres aléatoires Mersenne Twister
    std::normal_distribution<double> dist(0.0, 1.0);  // Distribution normale

    for (int i = 0; i < P; ++i) {
        // Paramètres utiles pour la discrétisation avec N pas
        double dt = T / N;
        double sigdt = sig * sqrt(dt);
        double rdt = r * dt;
        double der = rdt - pow(sigdt, 2) / 2;

        std::vector<double> maxerreul(M, 0.0);
        std::vector<double> maxerrmil(M, 0.0);

        for (int j = 0; j < M; ++j) {
            // Initialisation des M traj de l'EDS
            double S = S0;
            // Initialisation des M traj du schéma d'Euler
            double Se = S0;
            // Initialisation des M traj du schéma de Milstein
            double Sm = S0;

            // Initialisation de l'écart maximum EDS-Euler en valeur absolue
            double maxErrEul = 0.0;
            // Initialisation de l'écart maximum EDS-Milstein en valeur absolue
            double maxErrMil = 0.0;

            for (int k = 0; k < N; ++k) {
                // Génération d'un vecteur de gaussiennes centrées réduites
                // double g = rand() / (RAND_MAX + 1.0); La fonction rand() me donnait des résultats incohérents avec ceux trouvés via Python
                double g = dist(gen);  // Utilisation de la distribution normale

                // Évolution du sous-jacent S, du schéma d'Euler Se et du schéma de Milstein Sm 
                S = S * exp(sigdt * g + der);
                Se = Se * (1 + sigdt * g + rdt);
                Sm = Sm * (1 + sigdt * g + pow(sigdt * g, 2) / 2 + der);

                maxErrEul = fmax(maxErrEul, fabs(S - Se));
                maxErrMil = fmax(maxErrMil, fabs(S - Sm));
            }

            // Stockage de l'écart maximum pour chaque simulation
            maxerreul[j] = maxErrEul;
            maxerrmil[j] = maxErrMil;
        }

        // Calcul des statistiques
        double someul = 0.0;
        double careul = 0.0;
        double sommil = 0.0;
        double carmil = 0.0;

        for (int k = 0; k < M; ++k) {
            someul += pow(maxerreul[k], 2);
            careul += pow(maxerreul[k], 4);
            sommil += pow(maxerrmil[k], 2);
            carmil += pow(maxerrmil[k], 4);
        }

        erreul[i] = someul / M;
        liceul[i] = 1.96 * sqrt((careul / M - pow(someul / M, 2)) / M);
        errmil[i] = sommil / M;
        licmil[i] = 1.96 * sqrt((carmil / M - pow(sommil / M, 2)) / M);
        Npas[i] = N;
        N = 2 * N;
        }
    }

    void Cvg_forte::displayResults() {
        // Affichage des vecteurs d'erreurs et des demi-largeurs des IC à 95%
        std::cout << "Erreur euler" << std::endl;
        for (int i = 0; i < P; ++i) {
            std::cout << erreul[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "LIC Euler" << std::endl;
        for (int i = 0; i < P; ++i) {
            std::cout << liceul[i] << " ";
        }
        std::cout << std::endl;

        // Affichage des vecteurs d'erreurs et des demi-largeurs des IC à 95%
        std::cout << "Erreur Milstein" << std::endl;
        for (int i = 0; i < P; ++i) {
            std::cout << errmil[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "LIC Milstein" << std::endl;
        for (int i = 0; i < P; ++i) {
            std::cout << licmil[i] << " ";
        }
        std::cout << std::endl;
    }



