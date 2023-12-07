def n_calcs(x1, x2, pm_angle):
    '''
    Calculates the number profiles necessary for the calculations we want, starting with a volume calculation and going into count, density, and delta density
    
    Parameters:
    -- x1, x2: The distances in the two dimensions of the projection we are working in
    -- pm_angle: The angle that is used for the cuts. Typically, 360 for the entire range, or the angle of the "slice" around the axis in question
    
    Returns:
    -- n: count density
    -- delta_dens: incrememntal change in n
    -- bin_edge_Mpc: Gives the spatial bins R for plotting 
    '''
    rad_dist = np.sqrt(x1**2 + x2**2)
    bin_edge_Mpc = np.logspace(start=np.log10(0.2), stop=np.log10(40.), num=9)

    bin_edge_Mpc = np.append(np.array([0.]), bin_edge_Mpc)  # The left-most bin includes 0.
    
    bin_area = (bin_edge_Mpc[1:]**2 - bin_edge_Mpc[:-1]**2) * np.pi * pm_angle/90.  # Including the central circle.
    accum_area = bin_edge_Mpc[1:]**2 * np.pi * pm_angle/90.  # Fan-shape areas; len of bin_edge_Mpc - 1.

    count = stats.binned_statistic(rad_dist, None, statistic='count', bins=bin_edge_Mpc)[0] # Len = len of bin_edge_Mpc - 1

    dens = count / bin_area   # Inside those individual bins; len of bin_edge_Mpc - 1.

    accum_count = np.cumsum(count)  # Same length as count_gal = len of bin_edge_Mpc - 1; the most left bin has the value of count_gal's most left bin.

    dens_inside = accum_count / accum_area  # n(<r); len of bin_edge_Mpc - 1

    delta_dens = dens_inside[:-1] - dens[1:]  # len of bin_edge_Mpc - 2
    
    return dens, delta_dens, bin_edge_Mpc

def profile_plotting():
    '''
    Brainstorming a function for making plotting easier
    '''
    
    return None