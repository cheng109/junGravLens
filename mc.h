/*
 * mc.h
 *
 *  Created on: Apr, 2016
 *      Author: En-Hsin Peng
 */

#ifndef MC_H_
#define MC_H_

#include <vector>
#include <random>
#include <functional>
#include <cmath>
#include "Model.h"

class MC {
public:
    MC(unsigned seed);
    double gauss(double sigma, double mu, double x);
    void makeCgauss();
    double cgauss();
    double random();
    double stepPar(MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, vector<vector<size_t>> &iter, int &j, int &k);
    void stepPar(MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, int &iter);
    double stepPar(vector<vec> &src, double cfac, size_t &iter);
    void checkPoint(string outputFileName, MultModelParam &param, vector<vector<size_t>> &freePar, double cfac, vector<vector<size_t>> &iter, double L);
    void checkPoint(string outputFileName, MultModelParam &param, vector<vector<size_t>> &freePar, double cfac, int iter, double L);
    double load(string fileName, MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, vector<vector<size_t>> &iter);
    double load(string fileName, MultModelParam &param, vector<vector<size_t>> &freePar, double &cfac, int &iter);

private:
    std::vector<double> cgArr;
    std::mt19937 rng_engine;
    //std::default_random_engine rng_engine;
    std::function<double()> rng;
    std::function<double()> normal;
};

#endif /* MC_H_ */
