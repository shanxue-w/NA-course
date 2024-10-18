#ifndef __FUNCTION_HPP__
#define __FUNCTION_HPP__

class Function
{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual double integral(double a, double b) const = 0;
    virtual ~Function() {}
private:
    
};

#endif