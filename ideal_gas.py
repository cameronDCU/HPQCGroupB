import numpy as np
import sys
import time

# Constants
L = 10.0  # Box size
T = 1.0  # Arbitrary temperature scale for velocity distribution
dt = 0.01  # Time step
steps = 1000  # Number of simulation steps

# Check if the number of particles is passed as a command-line argument
if len(sys.argv) != 2:
    print("Usage: python ideal_gas.py <num_particles>")
    sys.exit(1)

N = int(sys.argv[1])  # Number of particles

# Initialize positions (uniformly distributed in the box)
pos = np.random.uniform(0, L, size=(N, 3))

# Initialize velocities from Maxwell-Boltzmann distribution manually
v_mag = np.sqrt(-2 * T * np.log(np.random.uniform(size=N)))  # Maxwell-Boltzmann speed
v_dir = np.random.normal(size=(N, 3))  # Random directions
v_dir /= np.linalg.norm(v_dir, axis=1, keepdims=True)  # Normalize
vel = v_mag[:, np.newaxis] * v_dir  # Assign velocities

# Track wall collisions for pressure calculation
wall_collisions = 0

def update_positions(pos, vel, L, dt):
    """Update particle positions and handle wall collisions."""
    global wall_collisions
    pos += vel * dt  # Move particles

    for i in range(len(pos)):
        for j in range(3):  # Check each coordinate (x, y, z)
            if pos[i, j] < 0 or pos[i, j] > L:
                vel[i, j] *= -1  # Reverse velocity
                pos[i, j] = max(0, min(L, pos[i, j]))  # Keep inside box
                wall_collisions += 1  # Count collision

    return pos, vel

def net_velocity(vel):
    """Calculate the net velocity vector (should be close to zero)."""
    return np.sum(vel, axis=0)  # Sum of all velocity components

def check_energy_conservation(vel_initial, vel_final):
    """Check if total kinetic energy is conserved."""
    KE_initial = np.sum(0.5 * np.sum(vel_initial**2, axis=1))  # KE before
    KE_final = np.sum(0.5 * np.sum(vel_final**2, axis=1))  # KE after
    return np.isclose(KE_initial, KE_final)

def calculate_temperature(vel):
    """Estimate temperature from mean velocity squared."""
    v_mean_sq = np.mean(np.sum(vel**2, axis=1))
    return v_mean_sq / 3  # Approximate T

def estimate_pressure(wall_collisions, V, N, T):
    """Estimate pressure using the ideal gas approximation."""
    return (N * T) / V  # Using P ≈ NkT/V (simplified)

# Store initial velocities for energy check
vel_initial = vel.copy()

# Start time measurement
start_time = time.time()

# Run simulation
for _ in range(steps):
    pos, vel = update_positions(pos, vel, L, dt)

# Perform testing
net_vel = net_velocity(vel)
energy_conserved = check_energy_conservation(vel_initial, vel)
estimated_T = calculate_temperature(vel)
estimated_P = estimate_pressure(wall_collisions, L**3, N, estimated_T)

# End time measurement
end_time = time.time()

# Display results
print(f"Net velocity vector: {net_vel}")
print(f"Energy conserved: {energy_conserved}")
print(f"Estimated Temperature: {estimated_T:.3f}")
print(f"Estimated Pressure: {estimated_P:.3f}")
print(f"Execution Time: {end_time - start_time:.6f} seconds")
