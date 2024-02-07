#include <iostream>
#include <algorithm>
#include <numeric>
#include "Cvg_forte.h"
#include "Cvg_faible.h"
#include "Romberg.h"


int main() {

    double T = 1.0;
    double sig = 0.2;
    double r = 0.05;
    double S0 = 100.0;
    double K = 110.0;
    int N = 1;
    int M = 100000;
    int P = 6;
    int R = 10;
    int M_R = 1000000;
    int P_R = 4;

    // Etude de la vitesse forte des schemas d'Euler et de Milshtein
    std::cout << "Convergence forte des schemas d'Euler et Milshtein" << std::endl;
    Cvg_forte vitesse_forte(T, sig, r, S0, N, M, P);
    vitesse_forte.calculate();
    vitesse_forte.displayResults();

    
    // Etude de la vitesse faible des schemas d'Euler et de Milshtein
    std::cout << "Vitesse faible des schemas d'Euler et Milshtein" << std::endl;
    Cvg_faible european_put(T, sig, r, S0, K, N, M, P);
    european_put.calculateBS();  // Calculer la valeur Black-Scholes
    european_put.calculate2();
    european_put.displayResults2();

    // Etude de la vitesse faible des schemas d'Euler et de Milshtein par la méthode d'extrapolation de Romberg
    std::cout << "Etude de l'extrapolation de Romberg" << std::endl;
    Romberg put_romberg(T, sig, r, S0, K, N, R, M_R, P_R);
    put_romberg.calculate3();
    put_romberg.displayResults3();

	return 0;

};