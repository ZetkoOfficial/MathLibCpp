#ifndef COMPLEX_LIB
#define COMPLEX_LIB

#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include "matrixvector.hpp"

namespace complex {
    class c_double {
        double a, b;

        public:

        c_double() = default;
        c_double(double a, double b = 0){ this->a = a; this->b = b; }
        c_double(matrixvector::vector vec) { this->a = vec[0], this->b = vec[1]; };

        //Returns vector representation of this complex number.
        matrixvector::vector to_vector() const { return {a, b}; }

        //Returns the real part of this complex number.
        double re() const { return a; }

        //Returns the imaginary part of this complex number.
        double im() const { return b; }

        //Returns the conjugate of this complex number.
        c_double conjugate() const { return {a, -b}; }

        //Creates complex number from polar form. 
        static c_double from_polar(double r, double angle) { return { r * cos(angle), r * sin(angle) }; }
    };

    //The imaginary constant i.
    c_double I = {0, 1};

    /*
    Returns the squared absolute value of a complex number.
    */
    double abs_sq(const c_double& c) { return c.re() * c.re() + c.im() * c.im(); }

    /*
    Returns the absolute value of a complex number.
    */
    double abs(const c_double& c) { return std::sqrt(abs_sq(c)); }

    /*
    Returns the argument(the angle) of a complex number.
    */
    double arg(const c_double& c) { return std::atan2((long double) c.im(), (long double) c.re()); }

    /*
    Adds two complex numbers.
    */
    c_double operator+(const c_double& a, const c_double& b) {
        return { a.re() + b.re(), a.im() + b.im() };
    }

    /*
    Subtracts two complex numbers.
    */
    c_double operator-(const c_double& a, const c_double& b) {
        return { a.re() - b.re(), a.im() - b.im() };
    }

    /*
    Multiplies two complex numbers.
    */
    c_double operator*(const c_double& a, const c_double& b) { 
        return { a.re() * b.re() - a.im() * b.im(), a.re() * b.im() + a.im() * b.re() };
    }

    /*
    Divides a complex number by a double.
    */
    c_double operator/(const c_double& a, const double b) { 
        return { a.re() / b, a.im() / b };
    }

    /*
    Divides two complex numbers.
    */
    c_double operator/(const c_double& a, const c_double& b) {
        return a * b.conjugate() / abs_sq(b);
    }

    /*
    Writes the complex number to a stream.
    Format: a + bi, a - bi or a.
    */
    std::ostream& operator<<(std::ostream& stream, const c_double& c) {
        if(c.im() > 0) stream << c.re() << " + " << c.im() << "i";
        if(c.im() < 0) stream << c.re() << " - " << -c.im() << "i";
        if(c.im() == 0) stream << c.re();
        
        return stream;
    }
    
    /*
    Uses binary exponentiation to quickly calculate power of complex number.
    */
    c_double pow(const c_double& mat, int power) {
        if(power == 0) return 1;

        if(power % 2 == 0) return pow(mat, power/2) * pow(mat, power/2);
        return mat * pow(mat, power - 1);
    }

    /*
    Calculates the p-th root of a complex number using its polar form. 
    */
    c_double root_p(const c_double& c, int p) { return c_double::from_polar(std::pow(abs(c), 1.0/p), arg(c) / p); }
    
    /*
    Calculates the square root of a complex number.
    */
    c_double sqrt(const c_double& c) { 
        return { std::sqrt( (abs(c) + c.re()) / 2 ), (c.im() > 0 ? 1 : -1) * std::sqrt( (abs(c) - c.re()) / 2 ) };
    }
}

/*
Allows usage of operators without namespace.
*/
using complex::operator*;
using complex::operator+;
using complex::operator/;
using complex::operator<<;

#endif