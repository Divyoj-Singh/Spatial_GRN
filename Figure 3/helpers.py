"""# Helper Functions
Provides
1. Parameter Set Initialisation
2. Functional Definition
3. Plotting tools
"""

# SCSCSCSCSC

# CSB Lab | Indian Institute of Science, Bengaluru
# Chinmay K Haritas

# Toggle Triad in 2D with Diffusion
# Helpers Functions

# Import Block
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from numba import njit, prange
import matplotlib.animation as animator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from tt_vars import *
import os

sns.set_theme(context="notebook", style="ticks", font="arial", font_scale=1.2)

# Returns all libraries
def libraries():
    """ # Libraries
    Libraries for Spatial GRNs
    """
    return np, pd, plt, tqdm, time, njit, prange, animator, make_axes_locatable, os


# Return Paramaters
def parameter_set(i: int, ps: pd.DataFrame):
    """
    # Parameters
    Returns the interactions variables for the the given set id, sheet_name and file"""

    global g
    global k
    global coop
    global thr
    global fold

    ds_i = i
    g = np.array([ps["Prod_of_A"][ds_i], ps["Prod_of_B"][ds_i], ps["Prod_of_C"][ds_i]])
    # Degradation
    k = np.array([ps["Deg_of_A"][ds_i], ps["Deg_of_B"][ds_i], ps["Deg_of_C"][ds_i]])

    # Cooperativity
    coop = np.array(
        [
            [ps["Num_of_AToA"][ds_i], ps["Num_of_AToB"][ds_i], ps["Num_of_AToC"][ds_i]],
            [ps["Num_of_BToA"][ds_i], ps["Num_of_BToB"][ds_i], ps["Num_of_BToC"][ds_i]],
            [ps["Num_of_CToA"][ds_i], ps["Num_of_CToB"][ds_i], ps["Num_of_CToC"][ds_i]],
        ]
    )

    # Threshold
    thr = np.array(
        [
            [ps["Trd_of_AToA"][ds_i], ps["Trd_of_AToB"][ds_i], ps["Trd_of_AToC"][ds_i]],
            [ps["Trd_of_BToA"][ds_i], ps["Trd_of_BToB"][ds_i], ps["Trd_of_BToC"][ds_i]],
            [ps["Trd_of_CToA"][ds_i], ps["Trd_of_CToB"][ds_i], ps["Trd_of_CToC"][ds_i]],
        ]
    )

    # Fold
    fold = np.array(
        [
            [ps["Inh_of_AToA"][ds_i], ps["Inh_of_AToB"][ds_i], ps["Inh_of_AToC"][ds_i]],
            [ps["Inh_of_BToA"][ds_i], ps["Inh_of_BToB"][ds_i], ps["Inh_of_BToC"][ds_i]],
            [ps["Inh_of_CToA"][ds_i], ps["Inh_of_CToB"][ds_i], ps["Inh_of_CToC"][ds_i]],
        ]
    )

    print("Starter of dataset", g[0])
    return g, k, coop, thr, fold


# Boundary Conditions
@njit
def p_i(x):
    """Input a [x,y] and recieve the expected coordinates according to periodic boundary conditions"""
    # X is any point to be accessed in space

    # Neumann
    # Send to first if before first asked
    # Send to last if after last asked
    y = np.copy(x)
    if x[0] == -1: y[0] = 0  # X
    if x[0] == nx + 1: y[0] = nx  # X
    if x[1] == -1: y[1] = 0  # Y
    if x[1] == ny + 1: y[1] = ny  # Y
    return y


# Hill Function
@njit
def hill(N, fold, n, thr):
    """# Shifted Hill Function
    1. Simulates the connection between two nodes.
    2. Below 1, lambda represses and above, it activates.
    """
    n_hill = 1 / (1 + (N / thr) ** n)
    return n_hill + fold * (1 - n_hill)


# IC Patching
def patch(x1, x2, y1, y2, fill_value):
    """Enter interval of patch and a function(x,y) to compute initial conditions at (x,y). """
    sol_1 = np.zeros((nx, ny, nodes))
    for x in range(x1, x2):
        for y in range(y1, y2):
            sol_1[x][y] = fill_value(x, y)
    sol_1_p = np.moveaxis(np.moveaxis(sol_1, 0, -1), 0, -1)
    # print(f"Before {np.shape(sol_1)}, after {np.shape(sol_1_p)}")
    return sol_1


# Special Patchers

# Uniform Fill
def fill_unif(u):
    """Fill space with uniform value 'u' """
    return lambda x, y: np.array([u[0], u[1], u[2]])


# Random Fill
r_max = 15


def fill_rand(x, y):
    """Fill space with random integers in the interval"""
    return np.random.rand(nodes) * r_max


# Custom Fill
def fill_custom(x, y):
    """Fill space with custom function"""
    return [0, 9, 0]

# G By k Normalization
def gbyk_normalization(sol_p, gbyk):
    """Performs g/k normalization for better data representation."""
    for i in range(len(sol_p)):
        sol_p[i] /= gbyk[i]
    return sol_p


# Plotting tools
def animation(sol_p):
    # Artists
    def lilly(s):
        plt.clf()
        m = 1
        probe = s * 100
        # plt.imshow((sol_p[0][probe]), cmap="jet", interpolation="gaussian")
        plt.imshow((sol_p[0][probe]), cmap="jet",)
        plt.title(f"At frame {probe}")
        plt.colorbar()

    anim = animator.FuncAnimation(
        plt.figure(), lilly, range(len(sol_p[0]) // 100), interval=10
    )
    plt.show()


# Saving Tools
def save_adv(savePath, plotName, sol_p, pset_deets, snapshots):
    # Time of saving
    nowString = time.strftime("_%Y_%b_%d__%H_%M_%S")
    folderName = savePath + "/" + plotName + nowString
    """
    # Advanced saving options
    Produces fast reproducible and usable plots with metadata.
    Saves three files.
    1. New Folder with Plot Name (schemes)
    2. Plot
    3. Metadata

        1. Timstamp and Author
        2. Parameter Set Details
        3. Misc

    4. 2D Image Data for reproducing plot
    """

    # Create new folder with timestamp and plot name
    os.mkdir(folderName)

    # Save Plot Image (pref as SVG)
    # TODO: Hardcoded axes variables

    # First the three morphogens individually @last

    # Parallely add the data
    if(snapshots):
        # For repressilator
        for t in range(0, 20):
            for i in range(3):
                plt.clf()
                plt.imshow((sol_p[i][t*7]), cmap="jet", interpolation="gaussian")
                plt.title(f"Morphogen {i+1}")
                plt.colorbar()
                plt.savefig(folderName + "/" + f"Morphogen {i+1}@{t*7}")
    else:
        data = open((folderName+"/data.txt"), "a")
        for i in range(3):
            plt.clf()
            plt.imshow((sol_p[i][-1]), cmap="bwr", interpolation="gaussian")
            # plt.imshow((sol_p[i][-1]), cmap="jet")
            plt.title(f"Morphogen {i+1}")
            plt.colorbar()
            plt.savefig(folderName + "/" + f"Morphogen {i+1}")
            np.savetxt(data, sol_p[i][-1])
        data.close()    


    
    # Synthesize the metadata file and save it
    # Metadata format, look at metadata_template.txt
    # md_temp = f"""SCSCSCSCSCSCSCSCSCSC\n\nAuthor:\nChinmay K Haritas,\n@CSB Lab | Indian Institute of Science, Bengaluru\n\nDate:\n{time.strftime("%Y %b %d - %H:%M:%S")}\n\nPlot Name:\n{plotName}\n\nParameter Set:\nSheet Name: {pset_deets['sheet']}\nIndex: {pset_deets['pset_id']} (Row + 2)\nDiffusion: {pset_deets['diff_coeff']} units\nhttps://indianinstituteofscience-my.sharepoint.com/:x:/g/personal/ushasiroy_iisc_ac_in/EZB1LguReEZOtHy_eycekwUBJjFwGmgtSGyl3wamuqapSQ?e=SBd4Ik\n"""
    md_temp = f"""SCSCSCSCSCSCSCSCSCSC\n\nAuthor:\nChinmay K Haritas,\n@CSB Lab | Indian Institute of Science, Bengaluru\n\nDate:\n{time.strftime("%Y %b %d - %H:%M:%S")}\n\nPlot Name:\n{plotName}\n\nParameter Set:\nSheet Name: {pset_deets['sheet']}\nIndex: {pset_deets['pset_id']} (Row + 2)\nDiffusion: {pset_deets['diff_coeff']} units\n"""
    with open((folderName+"/Meta.txt"), "w") as meta:
        meta.write(md_temp)
