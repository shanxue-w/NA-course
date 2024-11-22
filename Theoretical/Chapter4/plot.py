import numpy as np
import matplotlib.pyplot as plt

x_lists = np.linspace(0.0001, 1, 100000)
y_lists = x_lists * np.exp(-x_lists) / (1 - np.exp(-x_lists))

plt.plot(x_lists, y_lists)
plt.title('cond$_f(x) = \\frac{xe^{-x}}{1-e^{-x}}$')
plt.savefig('plot.png')