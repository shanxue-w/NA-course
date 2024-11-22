import numpy as np

def B(x, x_lists, i):
    #find x in [x_i, x_{i+1}]
    idx = np.searchsorted(x_lists, x)-1
    if (i == 1):
        return (x-x_lists[idx])/(x_lists[idx+1]-x_lists[idx])
    else:
        return 



if __name__ == '__main__':
    x = 0.6
    x_lists = [0, 0.5, 1, 1.5, 2]
    i = 0
    B(x, x_lists, i)