#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>


class Cvg_faible {

private:
    double T;
    double sig;
    double r;
    double S0;
    double K;
    int N;
    int M;
    int P;
    double BS;

    std::vector<double> erreul;
    std::vector<double> liceul;
    std::vector<double> conterreul;
    std::vector<double> contliceul;
    std::vector<double> errmil;
    std::vector<double> licmil;
    std::vector<double> conterrmil;
    std::vector<double> contlicmil;

    // Vectors used for statistic calculus
    std::vector<double> careul_vec;
    std::vector<double> carmil_vec;
    std::vector<double> puteul;
    std::vector<double> putmil;
    std::vector<double> putmc;
    std::vector<double> contputeul;
    std::vector<double> contputmil;


    std::vector<int> Npas;
    std::vector<double> S;
    std::vector<double> Se;
    std::vector<double> Sm;

public:
    // Constructor
    Cvg_faible(double _T, double _sig, double _r, double _S0, double _K, int _N, int _M, int _P);

    void calculateBS();

    void calculate2();

    void displayResults2();


};