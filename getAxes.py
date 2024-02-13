import numpy as np
from inertia_tensors import inertia_tensors

def get_axes(r1, r2):
    '''
    This function calculates the major and minor axes of a set of points, regardless of if it is large-scale- or sub-structure
    
    Parameters
    -- r1: distance from central halo to "member" halo, NumPy array of the first spatial dimension
    -- r2: same as r1 but for the second spatial dimension
    
    Returns
    -- (maj_ax, maj_vec): Eigenvalue and eigenvector of the proper major axis
    -- (min_ax, min_vec): Same as above, but for minor axis
    '''
    ## getting the values
    r = np.vstack((r1,r2)).T
    I = inertia_tensors(r)
    
    ## calculating eigenstuff
    evals, evecs = np.linalg.eigh(I)
    evals = evals/np.max(evals) # scaling eigenvalues
    
    ## identifying major/minor axis
    maj_ax = np.sqrt(evals[0,1])
    min_ax = np.sqrt(evals[0,0])
    
    maj_vec = evecs[0,1,:]
    min_vec = evecs[0,0,:]
    
    return (maj_ax, maj_vec), (min_ax, min_vec)

def get_axis_cut(ev, x1, x2, theta, r_max):
    '''
    This function gets "pizza" cut along a major or minor axis in both directions
    
    Parameters:
    -- ev: The eigenvector that it is getting a cut around
    -- x1, x2: Distances in x1, x2 spatial dimension (e.g., x and y) from the central halo to the nearby halo, NumPy array
    
    Returns:
    -- 
    '''
    ## converting to polar coords
    th = np.arctan(x2 / x1)
    
    ## separating into quandrants due to range limits on arctan
    #quadrant 2
    th[np.where((x1<0) & (x2>0))] += np.pi
    #quadrant 3 
    th[np.where((x1<0) & (x2<0))] += np.pi
    #quadrant 4
    th[np.where((x1>0) & (x2<0))] += np.pi
    
    r = np.sqrt(x1**2 + x2**2)
    theta = np.radians(theta)
    
    ## getting eigenvector angle
    ev_theta1 = np.arctan(ev[1]/ev[0])
    theta_lower1 = ev_theta1 - theta
    theta_upper1 = ev_theta1 + theta
    
    ev_theta2 = ev_theta1 + np.pi
    theta_lower2 = ev_theta2 - theta
    theta_upper2 = ev_theta2 + theta
    
    ## getting cuts for within bounds of the vector/axis given
    cut1 = np.where( (r <= r_max) & ( ((th >= theta_lower1) & (th <= theta_upper1)))) 
    cut2 = np.where( (r <= r_max) & ( ((th >= theta_lower2) & (th <= theta_upper2))))
    
    return np.append(cut1, cut2)

def stack_axes(dists_full_x1, dists_full_x2, maj_cut, min_cut):
    '''
    This function takes in all major and minor axis cuts, and returns them all in one. 
    This will include their functions of distance so they can be properly binned .
    
    Parameters
    -- dists: dictionary of all the distances either along the major/minor slices or of all as x1, x2 tuples
    -- major_cuts, minor_cuts: The cuts of the index locations 
    
    Returns
    
    -- stacked_maj, stacked_min: Distances of halos from center along major/minor axis (stacked)
    -- stacked_full: stacked distances of all halos for radial binning
    '''
    
    stacked_full_x1 = np.array([])
    stacked_full_x2 = np.array([])
    
    stacked_maj_x1 = np.array([])
    stacked_maj_x2 = np.array([])
    
    stacked_min_x1 = np.array([])
    stacked_min_x2 = np.array([])
    
    keys = dists_full_x1.keys() # keys should e the same for each dictionary
    
    for key in keys:
        stacked_full_x1 = np.append(stacked_full_x1, dists_full_x1[key])
        stacked_full_x2 = np.append(stacked_full_x2, dists_full_x2[key])
        
        stacked_maj_x1 = np.append(stacked_maj_x1, dists_full_x1[key][maj_cut[key]])
        stacked_maj_x2 = np.append(stacked_maj_x2, dists_full_x2[key][maj_cut[key]])
        
        stacked_min_x1 = np.append(stacked_min_x1, dists_full_x1[key][min_cut[key]])
        stacked_min_x2 = np.append(stacked_min_x2, dists_full_x2[key][min_cut[key]])
    
    return stacked_maj_x1, stacked_maj_x2, stacked_min_x1, stacked_min_x2, stacked_full_x1, stacked_full_x2
    
    
