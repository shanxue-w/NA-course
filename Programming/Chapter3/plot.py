# read file log.txt
# plot the following data, in the form
# 10.09,9.60864e-18
# 10.1,2.75573e-17
# 10.11,7.14766e-17
# 10.12,1.70628e-16
# 10.13,3.79901e-16
# 10.14,7.97108e-16
# 10.15,1.58909e-15
# 10.16,3.02996e-15
# 10.17,5.55554e-15
# 10.18,9.83925e-15
# 10.19,1.68956e-14
# 10.2,2.82187e-14
# 10.21,4.59653e-14
# 10.22,7.3192e-14
# 10.23,1.1416e-13
# 10.24,1.74723e-13
# 10.25,2.62807e-13
# 10.26,3.89019e-13
# 10.27,5.67381e-13
# 10.28,8.16239e-13
# 10.29,1.15936e-12
# 10.3,1.62723e-12
# 10.31,2.25868e-12

import matplotlib.pyplot as plt
import numpy as np

# read file log.txt
with open('log.txt', 'r') as f:
    lines = f.readlines()

# split lines into two lists
x = []
y = []
for line in lines:
    line = line.strip()
    if line:
        x.append(float(line.split(',')[0]))
        y.append(float(line.split(',')[1]))

# plot the data
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.show()