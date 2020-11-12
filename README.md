# polyFit
Same as matlab

## Download:  
download [PolynomialRegression.h](https://github.com/MohammadRaziei/polyFit/raw/master/PolynomialRegression.h)  

## Usage:
see [example](https://github.com/MohammadRaziei/polyFit/raw/master/main.cpp)  

### CODE:  
```c++
int main()
{
    const std::vector<double> x{0, 1, 2, 3, 4};
    const std::vector<double> y{1, 1.8, 1.3, 2.5, 6.3};
    const std::vector<double> x2{0, 1, 2, 3, 4};
    const int order = 2;
    try
    {
        std::cout << "x     : " << x << std::endl;
        std::cout << "y     : " << y << std::endl;
        std::cout << "x2    : " << x2 << std::endl;
        std::cout << "order : " << order << std::endl;
        std::cout << std::string(60, '*') << std::endl;
        const std::vector<double> coeffs = PolynomialRegression::fit(x, y, order);
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
```


### OUTPUT:  
```
x     : [0, 1, 2, 3, 4]
y     : [1, 1.8, 1.3, 2.5, 6.3]
x2    : [0, 1, 2, 3, 4]
order : 2
************************************************************

coeffs = fit(x, y, order) : [1.42, -1.07, 0.55]
eval(coeffs, x2)          : [1.42, 0.9, 1.48, 3.16, 5.94]

polyFit(x, y, x2, order) : [1.42, 0.9, 1.48, 3.16, 5.94]
polyFit(x, y, 5, order)  : 9.82
```
