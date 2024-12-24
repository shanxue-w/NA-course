import numpy as np

# y' = x + y, y(0) = 0, y = e^x - x - 1.

def test(N):
    x = np.linspace(0, 1, N)
    h = x[1] - x[0]
    y_n = 0
    for i in range(N-1):
        y_n = y_n + h * (x[i] + y_n + 0.5*h*(1 + x[i] + h/3 + y_n + 1/3 * h *(x[i] + y_n)))
    return y_n - np.exp(1) + 2


import scipy 

n_lists = np.array([10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])
log_n_lists = np.log(n_lists)
log_errors = np.array([np.log(np.abs(test(n))) for n in n_lists])
# 最小二乘法找到拟合直线
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(log_n_lists, log_errors)
print(slope)
print(intercept)
print(r_value)
print(p_value)
print(std_err)

