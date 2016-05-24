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
    double stepPar(MultModelParam &param, double cfac, size_t &iter);

private:
    std::vector<double> cgArr;
    std::mt19937 rng_engine;
    //std::default_random_engine rng_engine;
    std::function<double()> rng;
};

#endif /* MC_H_ */
