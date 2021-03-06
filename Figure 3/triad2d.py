# SCSCSCSCSCSCSCSCSCSC

# CSB Lab | Indian Institute of Science, Bengaluru
# Chinmay K Haritas

# Toggle Triad in 2D with Diffusion

print("\nSpatial GRNs - Toggle Triad in 2D\n")

# Import Block
import helpers
from tt_vars import *
from helpers import hill, p_i
np, pd, plt, tqdm, time, njit, prange, animator, mal, os = helpers.libraries()


# Easy Change Parameters
plotName = "Tristable"
plotParent = os.getcwd() + "/Plots"
snapit = False
saveit = True

# Parameter set details
ds_i, pset_class = 1, "Testing"
# Execution timesteps
T = 3000

# Physical Parameters
# Diffusion (& Anisotropy)
D = 7.4*1e-3

print(f"Starting with \033[95m ({nx}, {ny}) \033[0m lattice for \033[92m {T} \033[0m units. \n")

# Setting up the system

# Computational
dx, dy, dt = 0.5, 0.5, 0.05

# Loading Parameters from Dataset
g, k, coop, thr, fold = helpers.parameter_set(ds_i, pd.read_excel('data.xlsx', sheet_name=pset_class))

# Showing Dataset Properties
print(f"Dataset \033[95m{ds_i}\033[0m from \033[92m {pset_class} \033[0m")

# Starting numerical methods
sol = np.ones((T, nx, ny, nodes))*(0)

# Setting Initial Conditions
sol[0] = helpers.patch(0, nx, 0, ny, helpers.fill_rand)

# Synthesising Functions
@njit
def f(N):
    """Enter population number 'N' and receive the change"""
    # Put the function here
    dN = np.empty(nodes)
    # TODO: Enter equation when needed

    # Equation for morphogen interactions
    dN[0] = g[0]*hill(N[2], fold[2][0], coop[2][0], thr[2][0])*hill(N[1], fold[1][0], coop[1][0],  thr[1][0]) - k[0]*N[0]

    dN[1] = g[1]*hill(N[0], fold[0][1], coop[0][1], thr[0][1])*hill(N[2], fold[2][1], coop[2][1],  thr[2][1]) - k[1]*N[1]

    dN[2] = g[2]*hill(N[1], fold[1][2], coop[1][2], thr[1][2])*hill(N[0], fold[0][2], coop[0][2],  thr[0][2]) - k[2]*N[2]
    return dN

@njit()
def DeltaN(prev):
    """Compute deltas for every point in space based on current values"""
    # Create a new empty instance
    delta_earth = np.zeros((nx, ny, nodes))
    for i in (prange((nx+1)*(ny))):
        # Coordinates from just one direction
        x = i // (nx+1)
        y = i - x*(ny+1)
        # Calculating diffusion contribution
        d_next_x, d_next_y = p_i(np.array([x+1, y+1]))
        d_prev_x, d_prev_y = p_i(np.array([x-1, y-1]))
        diffusion = (( prev[d_next_x][y] + prev[d_prev_x][y] - 2*prev[x][y] )/dx**2) + (( prev[x][d_next_y] + prev[x][d_prev_y] - 2*prev[x][y]  )/dy**2)

        # Cumulated Delta
        delta_earth[x][y] = f(prev[x][y]) + diffusion*D

        
    # Finished Iteration
    return delta_earth

# Computing
for t in tqdm(range(1, T)):
    # Get the deltas
    deltaN = DeltaN(sol[t-1])
    # Make the solution
    sol[t] = sol[t-1] + dt*deltaN

# Plotting
sol_p = helpers.gbyk_normalization(np.moveaxis(np.moveaxis(np.moveaxis(sol, 0, -1), 0, -1),0 ,-1), g/k)

# Animation

# Listing the kinds of plots
# 1. Animation
# helpers.rgb_plot(sol_p)
# helpers.animation(sol_p)
# 2. Animation snapshots
# 3. All 3 states at a given time
# 4. Special State Plot
# 5. Histogram of final states

# Saving it
if(saveit): helpers.save_adv(plotParent, plotName, sol_p, {"sheet":pset_class, "pset_id":ds_i, "diff_coeff":D,}, snapshots = snapit)
# Showing Plots vs Saving them
# User controlled

print('All done, Ciao')
