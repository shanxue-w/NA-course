import numpy as np

def solve(x0: float, N: int = 100, eps: float = 1e-12) -> float:
    x = x0
    print(x)
    for _ in range(N):
        x = (-x*np.cos(x) + 0.5 + np.sin(x)) / (1 - np.cos(x))
        print(x)
        if abs(x - 0.5 - np.sin(x)) < eps:
            break
    
    return x


if __name__ == '__main__':
    x0 = 1.5
    x = solve(x0)
    print(f'x = {x}, x - 0.5 - sin(x) = {x - 0.5 - np.sin(x)}')

