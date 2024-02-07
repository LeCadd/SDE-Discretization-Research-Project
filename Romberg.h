#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>

class Romberg
{
private:
    double T;
    double sig;
    double r;
    double S0;
    double K;
    int N;
    int R;
    int M_R;
    int P_R;
    int R_M;

    std::vector<double> erreul;
    std::vector<double> liceul;
    std::vector<double> errmil;
    std::vector<double> licmil;
    std::vector<double> Npas;

public:

    // Constructor
    Romberg(double _T, double _sig, double _r, double _S0, double _K, int _N, int _R, int _M_R, int _P_R);

    void calculate3();

    void displayResults3();

};

