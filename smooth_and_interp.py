import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline, make_lsq_spline, BSpline, CubicSpline


#number of gridpoints before and after jump to skip, and threshold value: magnitude of second-derivative-jump must be at least this much to trigger discontinuity finder
start_gap = 40
end_gap = 100
jump_thresh = 5

#number of gridpoints and total radius
n_gridpoints = 12800
R = 120

#use this global variable to store the known matching location so smoothing can be done on variables with less detectable discontinuities
gap_loc = -1

#smooths discontinuities in second derivative (presumably lower should also work).
def smooth_data(input_file, output_file):
    file = open(input_file, 'r')
    data = np.loadtxt(file, dtype=np.float128)
    
    radii = data[:, 0]
    phi_vals = data[:, 1]
    
    d_start = d_end = -1
    dp_dr = np.gradient(phi_vals, radii)
    d2_phi = np.gradient(dp_dr, radii)
    
    global gap_loc
    
    #find discontinuity (only looks for 1 for now)
    for j in range (3, np.size(phi_vals) - 4):
        if ( abs(d2_phi[j] - d2_phi[j-1]) > jump_thresh * abs(d2_phi[j - 1] - d2_phi[j - 2]) and gap_loc == -1):
            d_start = j - start_gap
            d_end = j + end_gap
            gap_loc = j
            break

    if (gap_loc > -1):
        d_start = gap_loc - start_gap
        d_end = gap_loc + end_gap
        

    if (d_start == -1):
        print("WARNING: discontinuity not found")
    #else:
        #print(d_start)
    
    num_points = np.size(phi_vals) - 1
    
    cubic_range = []
    
    for j in range(0, num_points + 1):
        if (j < d_start or j > d_end):
            cubic_range.append(j)
        
    
    c_spline = CubicSpline (radii[cubic_range], phi_vals[cubic_range])
    phi_cubic = c_spline(radii)
    
    dpc_dr = np.gradient(phi_cubic, radii)
    d2pc_dr = np.gradient(dpc_dr, radii)
    
    smoothed_phi = []
    for j in range(0, num_points + 1):
        smoothed_phi.append(phi_cubic[j])
        
    smoothed_data = np.column_stack((radii, smoothed_phi))
    np.savetxt(output_file, smoothed_data, delimiter='\t', comments='')

#interpolates onto a uniform grid with n gridpoints in [0,R]
def uniform_grid_interp(input_file, output_file, R, n_gridpoints):
    file = open(input_file, 'r')
    data = np.loadtxt(file, dtype=np.float128)
    
    radii = data[:, 0]
    data_vals = data[:, 1]
    
    num_points = np.size(data_vals) - 1
    
    cubic_range = []

    uniform_radii = np.linspace(0, R, n_gridpoints)
    
    for j in range(0, num_points + 1):
        cubic_range.append(j)
        
    
    c_spline = CubicSpline (radii, data_vals)
    uniform_cubic = c_spline(uniform_radii)

    #deriv = np.gradient(uniform_cubic, uniform_radii)
    deriv = np.gradient(uniform_cubic, uniform_radii)

    deriv_2 = np.gradient(deriv, uniform_radii)
    smoothed_data = []
    for j in range(0, n_gridpoints):
        smoothed_data.append(uniform_cubic[j])
        
    smoothed_line = np.column_stack((uniform_radii, smoothed_data))
    np.savetxt(output_file, smoothed_line, delimiter='\t', comments='')


smooth_data("Phi.dat", "Phi_smoothed.dat")

smooth_data("thet.dat", "thet_smoothed.dat")

start_gap = 10
end_gap = 20

smooth_data("A.dat", "A_smoothed.dat")
smooth_data("m.dat", "m_smoothed.dat")

uniform_grid_interp("A_smoothed.dat", "unif_A.dat", R, n_gridpoints)
uniform_grid_interp("Phi_smoothed.dat", "unif_Phi.dat", R, n_gridpoints)
uniform_grid_interp("thet_smoothed.dat", "unif_thet.dat", R, n_gridpoints)
uniform_grid_interp("m_smoothed.dat", "unif_m.dat", R, n_gridpoints)

print ("Data smoothed and interpolated onto uniform grid with n_gridpoints =  ", n_gridpoints, " and R = ", R)
