#ifndef __POLYINTERPOLATION_HPP__
#define __POLYINTERPOLATION_HPP__

class Function
{
public:
    virtual double operator()(double x) const = 0;
    virtual double derivative(double x) const = 0;
    virtual ~Function() {}
};

#endif