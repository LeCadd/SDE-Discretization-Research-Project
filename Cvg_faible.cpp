#include "Cvg_faible.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>
#include <numeric>


double normalCDF(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}


Cvg_faible::Cvg_faible(double _T, double _sig, double _r, double _S0, double _K, int _N, int _M, int _P):
    T(_T), sig(_sig), r(_r), S0(_S0), K(_K), BS(0.0), N(_N), M(_M), P(_P), Npas(_P, 0.0),
    erreul(_P, 0.0), liceul(_P, 0.0), conterreul(_P, 0.0), contliceul(_P, 0.0),
    errmil(_P, 0.0), licmil(_P, 0.0), conterrmil(_P, 0.0), contlicmil(_P, 0.0),
    careul_vec(_P, 0.0), puteul(_P, 0.0), contputeul(_P, 0.0),
    carmil_vec(_P, 0.0), putmil(_P, 0.0), contputmil(P, 0.0),
    putmc(_P, 0.0)
{}

void Cvg_faible::calculateBS()
{
    double d1 = (log(S0 / K) + r * T) / (sig * sqrt(T)) + sig * sqrt(T) / 2;
    double d2 = d1 - sig * sqrt(T);
    BS = K * exp(-r * T) * normalCDF(-d2) - S0 * normalCDF(-d1); // formule fermée du prix du put européen dans le modèle de Black-Scholes
    std::cout << "BS: " << BS << std::endl;
}


void Cvg_faible::calculate2()
{
    // Choix de l'algorithme de génération de nombres aléatoires
    std::mt19937 gen(time(0));  // Générateur de nombres aléatoires Mersenne Twister
    std::normal_distribution<double> dist(0.0, 1.0);  // Distribution normale


    for (int i = 0; i < P; ++i) {

        std::vector<double> S(M, S0);
        std::vector<double> Se(M, S0);
        std::vector<double> Sm(M, S0);

        // paramètres du problème
        double dt = T / N;
        double sigdt = sig * sqrt(dt);
        double rdt = r * dt;
        double der = rdt - pow(sigdt, 2) / 2;

        std::vector<double> maxerreul(M, 0.0);
        std::vector<double> maxerrmil(M, 0.0);

        // Initialisation et calcul des vecteurs contputeul, contcareul et contputmil, contcarmil
        std::vector<double> contputeul_vec(M, 0.0);
        std::vector<double> contputmil_vec(M, 0.0);
        std::vector<double> contcareul_vec(M, 0.0);
        std::vector<double> contcarmil_vec(M, 0.0);

        double contcareul = 0.0;
        double contcarmil = 0.0;
        double careul = 0.0;
        double carmil = 0.0;

        std::vector<double> g(M, 0.0);

        for (int k = 0; k < N; ++k) { // boucle sur les pas de temps

            for (int m = 0; m < M; ++m) {
                g[m] = dist(gen); // génération d'un vecteur de gaussiennes centrées réduites
                S[m] = S[m] * exp(sigdt * g[m] + der);
                Se[m] = Se[m] * (1 + sigdt * g[m] + rdt);
                Sm[m] = Sm[m] * (1 + sigdt * g[m] + pow(sigdt * g[m], 2) / 2 + der);
            }
        }

        // Initialisation des vecteurs paymc, payeul et paymil
        std::vector<double> paymc(M, 0.0);
        std::vector<double> payeul(M, 0.0);
        std::vector<double> paymil(M, 0.0);

        // Calcul des vecteurs paymc, payeul et paymil
        for (int m = 0; m < M; ++m) {
            paymc.at(m) = std::max(double(0.0), K - S[m]);
            payeul.at(m) = std::max(double(0.0), K - Se[m]);
            paymil.at(m) = std::max(double(0.0), K - Sm[m]);
        }

        puteul.at(i) = std::accumulate(payeul.begin(), payeul.end(), 0.0);
        putmil.at(i) = std::accumulate(paymil.begin(), paymil.end(), 0.0);
        putmc.at(i) = std::accumulate(paymc.begin(), paymc.end(), 0.0);

        contputeul.at(i) = puteul.at(i) - putmc.at(i);
        contputmil.at(i) = putmil.at(i) - putmc.at(i);

        for (int j = 0; j < M; ++j) {
            careul_vec.at(i) += payeul[j] * payeul[j];
            carmil_vec.at(i) += paymil[j] * paymil[j];
            contcareul_vec[j] = pow(payeul[j] - paymc[j], 2);
            contcarmil_vec[j] = pow(paymil[j] - paymc[j], 2);
        }

        contcareul = std::accumulate(contcareul_vec.begin(), contcareul_vec.end(), double(0.0));
        contcarmil = std::accumulate(contcarmil_vec.begin(), contcarmil_vec.end(), double(0.0));

        erreul[i] = exp(-r * T) * puteul.at(i) / M - BS;
        liceul[i] = 1.96 * exp(-r * T) * sqrt((careul_vec.at(i) / M - pow(puteul.at(i) / M, 2)) / M);
        conterreul[i] = exp(-r * T) * contputeul.at(i) / M;
        contliceul[i] = 1.96 * exp(-r * T) * sqrt((contcareul / M - pow(contputeul.at(i) / M, 2)) / M);

        errmil[i] = exp(-r * T) * putmil.at(i) / M - BS;
        licmil[i] = 1.96 * exp(-r * T) * sqrt((carmil_vec.at(i) / M - pow(putmil.at(i) / M, 2)) / M);
        conterrmil[i] = exp(-r * T) * contputmil.at(i) / M;
        contlicmil[i] = 1.96 * exp(-r * T) * sqrt((contcarmil / M - pow(contputmil.at(i) / M, 2)) / M);

        Npas[i] = N;
        N *= 2;

    } 
}

void Cvg_faible::displayResults2()
{// Display the results
    std::cout << "Vitesse faible des schemas d'Euler et de Milshtein" << std::endl;

    std::cout << "Erreur faible Euler" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << erreul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Euler" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << liceul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Erreur faible Euler VA controle" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << conterreul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Euler VA controle" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << contliceul[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Erreur faible Milstein" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << errmil[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Milstein" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << licmil[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Erreur faible Milstein VA controle" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << conterrmil[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "LIC Milshtein VA controle" << std::endl;
    for (int i = 0; i < P; ++i) {
        std::cout << contlicmil[i] << " ";
    }
    std::cout << std::endl;
}
