import numpy as np
import matplotlib.pyplot as plt

def B(i: int, x: float) -> float:
    if (x<=i) and (x>=i-1):
        return (x-i+1)**2 / 2
    elif (x<=i+1) and (x>=i):
        return (x-i+1)*(i+1-x) / 2 + (i+2-x)*(x-i) / 2
    elif (x<=i+2) and (x>=i+1):
        return (i+2-x)**2 / 2
    else:
        return 0
    

# Plot B_1(x)
if __name__ == '__main__':
    x = np.linspace(-1, 6, 1000)
    y = [B(1, xi) for xi in x]
    plt.plot(x, y, label='$B_1(x)$')
    y1 = [B(2, xi) for xi in x]
    plt.plot(x, y1, label='$B_2(x)$')
    y2 = [B(3, xi) for xi in x]
    plt.plot(x, y2, label='$B_3(x)$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('$B_i(x)$')
    plt.legend()
    plt.savefig("Plot_B.png")