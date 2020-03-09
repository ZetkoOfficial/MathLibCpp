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

    /*
    Calculates x-coordinate of the local minimum of the function using gradient descent(factor - step_size, steps - step_count).
    It uses the derivative method(using derivative_step_size) defined in the calculus namespace.
    */
    double minimum(std::function<double(double)> f, double start_x, int step_count = 10000, double step_size = 0.001, double derivative_step_size = 0.0001) {
        for(int i = 0; i < step_count; i++) start_x += -step_size * derivative(f, start_x, derivative_step_size);
        return start_x;
    }
}

#endif