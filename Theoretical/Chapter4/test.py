import numpy as np




def test(p):
    a = 1
    for i in range(1, p+1):
        a *= i
    b = 1
    for i in range(p+1, 2*p+1):
        b *= i
    c = p**p
    return (b-c)/a


if __name__ == '__main__':
    print(test(40))