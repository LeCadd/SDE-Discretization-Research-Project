#pragma once
#include <iostream>
#include <vector>

class Cvg_forte {
private:
    double T;    // maturité
    double sig;  // volatilité
    double r;    // taux sans risque
    double S0;   // valeur initiale du sous-jacent
    int N;       // nombre initial de pas de discrétisation
    int M;       // nombre de simulations indépendantes
    int P;       // nombre d'itérations sur la valeur de pas

    std::vector<double> erreul;
    std::vector<double> liceul;
    std::vector<double> errmil;
    std::vector<double> licmil;
    std::vector<int> Npas;

public:
    //Contructor
    Cvg_forte(double _T, double _sig, double _r, double _S0, int _N, int _M, int _P);

    void calculate();

    void displayResults();

};


