#include <iostream>
#include <cmath>
#include <iomanip>  // Para setprecision
#include <chrono>
double hbar2 = 7.6199682; // Planck's constant squared (in appropriate units)
double V0 = 10.0; // Potential well depth in eV
double a = 3.0; // Width of the potential well in Angstroms
double m = 1.0; // Mass of the particle (in appropriate units)
double alpha(double x){
    return std::sqrt(2 * m * x / hbar2);

}
double beta(double x) {
    return std::sqrt(2 * m * (V0 - x) / hbar2);
}


double function(double x) {
    double f =beta(x)*cos(alpha(x)*a)-sin(alpha(a)*a)*alpha(x);
    f=alpha(x)*1/tan(alpha(x)*a)+beta(x);
    return f; 
}

double derivative(double x) {
    return cos(x)-1.0/2; 
}
double bisect(double a,double b,int iteration=0) {
    if (iteration>30){std::cout << "Number of iterations exceeded the limit!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();}
    
    double fl = function(a);
    double m = (a+b)/2;
    double fm=function(m);
    if (std::fabs(fm)< 1e-9) {return m;}
    else{
        if (fl*fm<0){
            return bisect(a,m,++iteration);
        }
        else{
            return bisect(m,b,++iteration);
        }
    }
}
double newton(double x,int iteration=0) {
    /*in this method we use the derivative of the function 'couse f(x)
        is aproximately equal to f(a) - (x-a)f'(a) """
    */
    if(iteration>30){std::cout << "Number of iterations exceeded the limit!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();}
    double fx = function(x);
    double dfx = derivative(x);
    if (std::fabs(fx)< 1e-9) {return x;}
    else{
        return newton(x-fx/dfx,++iteration);
    }
}

double false_position(double a, double b, int iteration = 0) {
    if (iteration > 30) {
        std::cout << "Number of iterations exceeded the limit!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    double fa = function(a);
    double fb = function(b);
    double c = a - fa * (b - a) / (fb - fa);
    double fc = function(c);

    if (std::fabs(fc) < 1e-9) {
        return c;
    } else if (fa * fc < 0) {
        return false_position(a, c, ++iteration);
    } else {
        return false_position(c, b, ++iteration);
    }
}

double secant(double a, int iteration = 0) {
    if (iteration > 30) {
        std::cout << "Number of iterations exceeded the limit!" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    double fa = function(a);
    double fb = function(a + 1e-6);
    if(std::fabs(fa) < 1e-9) {
        return a;
    }
    double c = a - fa * 1e-6 / (fb - fa);
    return secant(c, ++iteration);
}

int main() {
    double left, right, x_newton;
    std::cin >> left;
    std::cin >> right;
    std::cin >> x_newton;
    std::cout << std::fixed << std::setprecision(12); 
    auto start = std::chrono::high_resolution_clock::now(); // Start timer
    std::cout<<bisect(left,right)<<std::endl;
    auto end = std::chrono::high_resolution_clock::now(); // End timer
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Bisect time: " << elapsed.count() << " seconds" << std::endl;
// -----------------------------------------------------
    // start = std::chrono::high_resolution_clock::now(); // Start timer
    // std::cout<<newton(x_newton)<<std::endl;
    // end = std::chrono::high_resolution_clock::now(); // End timer
    // elapsed = end - start;
    // std::cout << "Newton time: " << elapsed.count() << " seconds" << std::endl;
// -----------------------------------------------------
    start = std::chrono::high_resolution_clock::now(); // Start timer
    std::cout<<false_position(left,right)<<std::endl;
    end = std::chrono::high_resolution_clock::now(); // End timer
    elapsed = end - start;
    std::cout << "False Position time: " << elapsed.count() << " seconds" << std::endl;
    // -----------------------------------------------------
    start = std::chrono::high_resolution_clock::now(); // Start timer
    std::cout<<secant(left)<<std::endl;
    end = std::chrono::high_resolution_clock::now(); // End timer
    elapsed = end - start;
    std::cout << "Secant time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}