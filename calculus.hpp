#ifndef CALCULUS_LIB
#define CALCULUS_LIB

#include <cmath>
#include <functional>
#include <vector>
#include <string>

namespace calculus {
    /*
    Calculates the aproximate value of the derivative of a function at point x, using the derivative definition, where h is the step_size. 
    */
    double derivative(std::function<double(double)> f, double x, double step_size = 0.0001) {
        return (f(x + step_size) - f(x - step_size)) / (2 * step_size);
    }
    /*
    Calculates the aproximate value of the definate integral of a function
    between b and a, using the integral definition and rectangle method, where delta x is the step_size. 
    */
    double integral(std::function<double(double)> f, double a, double b, double step_size = 0.0001){
        double range = b - a;
        double n = range / step_size;

        double res = 0;
        for(int i = 0; i < n; i++) res += f(i * step_size + a) * step_size;

        return res;
    }
}

#endif