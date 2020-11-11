// Polynomial Regression (Quadratic Fit) in C++
#ifndef POLYNOMIALREGRESSION_H
#define POLYNOMIALREGRESSION_H

/**
 * PURPOSE:
 *
 *  Polynomial Regression aims to fit a non-linear relationship to a set of
 *  points. It approximates this by solving a series of linear equations using
 *  a least-squares approach.
 *
 *  We can model the expected value y as an nth degree polynomial, yielding
 *  the general polynomial regression model:
 *
 *  y = a0 + a1 * x + a2 * x^2 + ... + an * x^n
 *
 * LICENSE:
 *
 * MIT License
 *
 * Copyright (c) 2020 Chris Engelsma
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * @authors Chris Engelsma, Mohammad Raziei
 */
#include <vector>
#include <stdlib.h>
#include <cmath>
#include <stdexcept>

class PolynomialRegression
{
public:
    PolynomialRegression();
    virtual ~PolynomialRegression(){};

    template <class TYPE>
    static std::vector<TYPE> fit(
         const std::vector<TYPE>& x,
         const std::vector<TYPE>& y,
         const int& order);

    template <class TYPE>
    static std::vector<TYPE> eval(
         const std::vector<TYPE>& coeffs,
         const std::vector<TYPE>& x);
    template <class TYPE>
    static TYPE eval(
         const std::vector<TYPE>& coeffs,
         const TYPE& x);

    template <class TYPE, class T>
    static T polyFit(
         const std::vector<TYPE>& x,
         const std::vector<TYPE>& y,
         const T& x2Compute,
         const int& order);
};

template <class TYPE>
std::vector<TYPE> PolynomialRegression::fit(
     const std::vector<TYPE>& x,
     const std::vector<TYPE>& y,
     const int& order)
{
    /// The size of xValues and yValues should be same
    if (x.size() != y.size())
        throw std::runtime_error("The size of x & y arrays are different");

    /// The size of xValues and yValues cannot be 0, should not happen
    if (x.size() == 0 || y.size() == 0)
        throw std::runtime_error("The size of x or y arrays is 0");

    size_t N = x.size();
    int n = order;
    int np1 = n + 1;
    int np2 = n + 2;
    int tnp1 = 2 * n + 1;
    TYPE tmp;

    /// X = vector that stores values of sigma(xi^2n)
    std::vector<TYPE> X(tnp1);
    for (int i = 0; i < tnp1; ++i)
    {
        X[i] = 0;
        for (size_t j = 0; j < N; ++j)
            X[i] += (TYPE)pow(x[j], i);
    }

    /// a = vector to store final coefficients.
    std::vector<TYPE> a(np1);

    /// B = normal augmented matrix that stores the equations.
    std::vector<std::vector<TYPE>> B(np1, std::vector<TYPE>(np2, 0));

    for (int i = 0; i <= n; ++i)
        for (int j = 0; j <= n; ++j)
            B[i][j] = X[i + j];

    /// Y = vector to store values of sigma(xi^n * yi)
    std::vector<TYPE> Y(np1);
    for (int i = 0; i < np1; ++i)
    {
        Y[i] = (TYPE)0;
        for (size_t j = 0; j < N; ++j)
        {
            Y[i] += (TYPE)pow(x[j], i) * y[j];
        }
    }

    /// Load values of Y as last column of B
    for (int i = 0; i <= n; ++i)
        B[i][np1] = Y[i];

    n += 1;
    int nm1 = n - 1;

    /// Pivotisation of the B matrix.
    for (int i = 0; i < n; ++i)
        for (int k = i + 1; k < n; ++k)
            if (B[i][i] < B[k][i])
                for (int j = 0; j <= n; ++j)
                {
                    tmp = B[i][j];
                    B[i][j] = B[k][j];
                    B[k][j] = tmp;
                }

    /// Performs the Gaussian elimination.
    /// (1) Make all elements below the pivot equals to zero
    ///     or eliminate the variable.
    for (int i = 0; i < nm1; ++i)
        for (int k = i + 1; k < n; ++k)
        {
            TYPE t = B[k][i] / B[i][i];
            for (int j = 0; j <= n; ++j)
                B[k][j] -= t * B[i][j]; // (1)
        }

    /// Back substitution.
    /// (1) Set the variable as the rhs of last equation
    /// (2) Subtract all lhs values except the target coefficient.
    /// (3) Divide rhs by coefficient of variable being calculated.
    for (int i = nm1; i >= 0; --i)
    {
        a[i] = B[i][n]; // (1)
        for (int j = 0; j < n; ++j)
            if (j != i)
                a[i] -= B[i][j] * a[j]; // (2)
        a[i] /= B[i][i]; // (3)
    }

    //    std::vector<TYPE> coeffs(a.size());
    //    for (size_t i = 0; i < a.size(); ++i)
    //        coeffs[i] = a[i];

    return a; /// return coefficients
}

template <class TYPE>
TYPE PolynomialRegression::eval(
     const std::vector<TYPE>& coeffs,
     const TYPE& x)
{
    TYPE sum = 0;
    TYPE val = 1;
    size_t coeffSize = coeffs.size();
    for (size_t i = 0; i < coeffSize; ++i)
    {
        sum += coeffs[i] * val;
        val *= x;
    }
    return sum;
}

template <class TYPE>
std::vector<TYPE> PolynomialRegression::eval(
     const std::vector<TYPE>& coeffs,
     const std::vector<TYPE>& x)
{
    size_t xSize = x.size();
    std::vector<TYPE> ret(x.size());
    for (size_t i = 0; i < xSize; ++i)
        ret[i] = eval(coeffs, x[i]);
    return ret;
}

template <class TYPE, class T>
T PolynomialRegression::polyFit(
     const std::vector<TYPE>& x,
     const std::vector<TYPE>& y,
     const T& x2Compute,
     const int& order)
{
    return eval(fit(x, y, order), x2Compute);
}
#endif /// POLYNOMIALREGRESSION_H
