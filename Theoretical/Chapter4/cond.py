import numpy as np
import scipy as scp

# compute the condition number of matrix A
def cond(A, p):
    A_inv = np.linalg.inv(A)
    norm_A = np.linalg.norm(A, p)
    norm_A_inv = np.linalg.norm(A_inv, p)
    return norm_A * norm_A_inv


if __name__ == "__main__":
    n = 10
    a = 0
    A = 2*np.eye(n);
    A[0, 1] = 0.5; A[1, 2] = 1-a; A[2, 1] = 1-a; A[1,0] = a; A[2, 3] = a;
    for i in range(3, n-1):
        A[i, i-1] = 0.5; A[i, i+1] = 0.5;
    A[n-1, n-2] = 0.5;

    print(cond(A, 2))
