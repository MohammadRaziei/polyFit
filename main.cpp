/// Polynomial Regression (Quadratic Fit) in C++
#include "PolynomialRegression.h"
#include <iostream>
#include <iterator>
#define vectDouble std::vector<double>

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec)
{
    out << "[";
    std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<double>(out, ", "));
    out << *vec.rbegin() << "]";
    return out;
}

int main()
{
    vectDouble x{0, 1, 2, 3, 4};
    vectDouble y{1, 1.8, 1.3, 2.5, 6.3};
    vectDouble x2{0, 1, 2, 3, 4};
    int order = 2;
    try
    {
        std::cout << "x     : " << x << std::endl;
        std::cout << "y     : " << y << std::endl;
        std::cout << "x2    : " << x2 << std::endl;
        std::cout << "order : " << order << std::endl;
        std::cout << std::string(60, '*') << std::endl;
        vectDouble coeffs = PolynomialRegression::fit(x, y, order);
        std::cout << std::endl;
        std::cout << "coeffs = fit(x, y, order) : " << coeffs << std::endl;
        std::cout << "eval(coeffs, x2)          : " << PolynomialRegression::eval(coeffs, x2) << std::endl;
        std::cout << std::endl;
        std::cout << "polyFit(x, y, x2, order) : " << PolynomialRegression::polyFit(x, y, x2, order) << std::endl;
        std::cout << "polyFit(x, y, 5, order)  : " << PolynomialRegression::polyFit(x, y, 5., order) << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
}
