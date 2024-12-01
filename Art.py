import numpy as np
import matplotlib.pyplot as plt

def LennardJonesForce(r, sigma, epsilon):
    """Calculate the Lennard-Jones force."""
    return 24 * epsilon / r**2 * (2 * (sigma / r)**12 - (sigma / r)**6)

def simulate_brownian(N, R, L, TotalTime, dt, epsilon, sigma):
    # Simulation parameters
    Nsteps = int(TotalTime / dt)  # Number of time steps
    x = np.random.rand(2, N) * L  # Initial positions of particles

    # Precompute some values
    R2 = R ** 2  # Square of the radius
    sigma2 = sigma ** 2  # Square of the sigma

    # Store trajectories for plotting
    trajectories = np.zeros((2, N, Nsteps))
    trajectories[:, :, 0] = x

    for i in range(1, Nsteps):
        # Generate random noise (Brownian motion)
        xi = np.random.normal(0, 1, (2, N))

        # Initialize forces
        Fi = np.zeros((2, N))

        # Calculate pairwise forces using efficient vectorization
        for m in range(N):
            # Calculate distances to all other particles
            dx = x[0, m] - x[0]
            dy = x[1, m] - x[1]
            dx -= L * np.round(dx / L)  # Periodic boundary
            dy -= L * np.round(dy / L)  # Periodic boundary

            r = np.sqrt(dx**2 + dy**2)

            # Apply Lennard-Jones potential with cutoff
            mask = (r < 2.5 * sigma) & (r > 0)  # Cutoff distance (2.5*sigma to prevent large forces)

            # Calculate forces only for pairs within cutoff
            Fi[0, m] += np.sum(LennardJonesForce(r[mask], sigma, epsilon) * (dx[mask] / r[mask]))
            Fi[1, m] += np.sum(LennardJonesForce(r[mask], sigma, epsilon) * (dy[mask] / r[mask]))

        # Update positions with noise and forces
        x += (dt / R2) * (xi + Fi)

        # Apply periodic boundary conditions
        x = np.mod(x, L)

        # Store the trajectory
        trajectories[:, :, i] = x

    return trajectories

# Simulation parameters
N = 500  # Number of particles
R = 10e-3  # Particle radius in um (10 nm)
L = 1e3  # Box length in nm (1 um)
TotalTime = 1.0  # Total simulation time in seconds
dt = 0.01  # Time step
epsilon = 0.01  # Epsilon for interaction
sigma = 1e-2  # Sigma for interaction (1 nm)

# Run the simulation
trajectories = simulate_brownian(N, R, L, TotalTime, dt, epsilon, sigma)

# Plot the results
plt.figure(figsize=(8, 8))
for i in range(N):
    plt.plot(trajectories[0, i, :], trajectories[1, i, :], alpha=0.5)
plt.title("Brownian Motion of Particles")
plt.xlabel("X Position (nm)")
plt.ylabel("Y Position (nm)")
plt.xlim(0, L)
plt.ylim(0, L)
plt.grid()
plt.show()
