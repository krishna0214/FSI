
import numpy as np

def lu_decomposition(A):
    n = A.shape[0]
    L = np.zeros((n, n))
    U = np.zeros((n, n))  
    for i in range(n):
        L[i, i] = 1.0
        for j in range(i, n):
            U[i, j] = A[i, j] - np.dot(L[i, :i], U[:i, j])
        for j in range(i+1, n):
            L[j, i] = (A[j, i] - np.dot(L[j, :i], U[:i, i])) / U[i, i]
    return L, U
def solve_lu(A, b):
    L, U = lu_decomposition(A)
    n = A.shape[0]
    y = np.zeros((n, 1))
    x = np.zeros((n, 1))

    # Solve Ly = b using forward substitution
    for i in range(n):
        y[i, 0] = b[i, 0] - np.dot(L[i, :i], y[:i, 0])

    # Solve Ux = y using backward substitution
    for i in range(n - 1, -1, -1):
        x[i, 0] = (y[i, 0] - np.dot(U[i, i+1:], x[i+1:, 0])) / U[i, i]
    return x


"""
def solve_linear_system(A, b):
    try:
        # Attempt to solve the system using numpy.linalg.solve()
        X = np.linalg.solve(A, b)
        return X
    except np.linalg.LinAlgError:
        # Handle the case where A is singular (not invertible)
        raise ValueError("Matrix A is singular, and the system has no unique solution.")


"""