import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# trying to solve: u'' - \epsilon u^4 = 0; u(0)=1, u(1)=0


def nonlinear_system(x, v):
    """
    Define the system of first-order ODEs
    v[0] = u'
    v[1] = epsilon*u^4
    Returns derivatives [v', v'']
    """
    epsilon = 0.5

    return np.vstack((v[1], epsilon * v[0] ** 4))


def boundary_conditions(va, vb):
    """
    Define boundary conditions at x=a and x=b
    va: solution at x=a
    vb: solution at x=b
    """
    return np.array([va[0] - 1, vb[0] - 0])  # v(0) = 1  # v(1) = 0


def solve_nonlinear_bvp():
    # Generate initial mesh
    x = np.linspace(0, 1, 200)

    # Initial guess for the solution
    u_guess = np.zeros((2, x.size))
    u_guess[0] = x  # Linear interpolation between boundary values
    u_guess[1] = 1  # Guess for u'

    # Solve BVP
    solution = solve_bvp(nonlinear_system, boundary_conditions, x, u_guess)

    if not solution.success:
        raise ValueError("Failed to find solution")

    # Generate dense output for plotting
    x_plot = np.linspace(0, 1, 100)
    u_plot = solution.sol(x_plot)

    # Exact solution

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(x_plot, u_plot[0], "b-", label="u(x)")
    # plt.plot(x_plot, u_plot[1], 'r--', label="u'(x)")
    plt.grid(True)
    plt.legend()
    plt.xlabel(f"$x$")
    plt.ylabel(f"Solution $u(x)$")
    plt.title("Solution of Nonlinear BVP")
    plt.savefig("exam1/nonlinearbvp")

    return solution, x_plot, u_plot


# Solve and display results
solution, x, u = solve_nonlinear_bvp()
print(f"Solution successful: {solution.success}")
print(f"Number of iterations: {solution.niter}")
print(f"Residual norm: {solution.message}")
