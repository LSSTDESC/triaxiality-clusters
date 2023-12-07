'''
LSST/DESC - Signatures of Triaxiality, Object Selection & Cuts
Deric Jones

This file is 
'''

## necessary imports
import GCRCatalogs
import numpy as np

def init_list(cat_file):
    '''
    This function gives data required for analysis from the halo catalog.
    
    Parameters:
    -- cat_file: File location of the GCR halo catalog
    
    Returns:
    -- data: all the stuff we need for analysis
    '''
    gc = GCRCatalogs.load_catalog(cat_file)
    data = gc.get_quantities(['halo_id', 'galaxy_id', 'halo_mass', 'stellar_mass', 'is_central', 'position_x', 'position_y', 'position_z', 'redshift'])
    
    return data
    

def select_halos(cat_file, n, min_mass = 0):
    '''
    This function takes the GCRCatalog of halos and selects a certain number of halos, with a minimal mass range. This starts at the largest mass. 
    This returns data necessary for large scale structue and substructure analysis to be stacked
    
    Parameters:
    -- cat_file: File location to be passed into init_list function to get the early halo data 
    -- n: This is an integer number of halos you want to stack over/get data on. 
    -- min_mass: The minimum bound for mass, which will overwrite the number n if it goes below a certain mass range. By default, it is 0.
    
    Returns:
    -- data: all data necessary for analysis
    -- mass_index: indeces for the selected halos from the unique halo cut
    -- central_cut: central galaxies for positions of (unique) halos 
    '''
    ## calling init_list to get all the data
    data = init_list(cat_file)
    
    ## cutting for lss (central galaxies)
    central_cut = np.where(data['is_central'])[0]
    
    data_cent = {}
    for key in data.keys():
        data_cent[key] = data[key][central_cut]
    
    ## organizing by halo mass
    # unique_ids, unique_index = np.unique(data['halo_id'], return_index=True)
    #mass_arr = data['halo_mass'][unique_index]
    mass_sort = np.sort(data_cent['halo_mass'])
    #print(mass_sort[(len(mass_sort)-n):])
    main_index = np.zeros(n)
    for i in range(n):
        main_index[i] = np.where(data_cent['halo_mass'] == mass_sort[len(mass_sort)-(i+1)])[0]
    
    return data, data_cent, main_index

def lss_cyl_cut(data, main_index, plane='xy', dist_lim=10, h=5):
    '''
    This function does a cylindrical cut of nearby halos
    
    Parameters
    -- data: defined with the unique index cut
    -- main_index: this is the index of the halos that will be used
    -- plane: which 2D plane to do a projection of
    -- dist_lim: Distance limit of nearby halos
    -- h: Height of the "3D" direction for the cylindrical cut
    
    Returns
    -- cyl_cut: A cut of the halo indeces in a 2D cyclical cut as a dictionary of the ID they belong to if n>1 
    '''
    ## getting the positions of all central galaxies to unique halos
    xs = data['position_x']
    ys = data['position_y']
    zs = data['position_z']
    
    id_halos = data['halo_id'][main_index.astype(int)]
    
    cyl_cut = {}
    
    ## looping through
    for i in range(len(id_halos)):
        # position of central halo
        central_xs = data['position_x'][np.where(data['halo_id'] == id_halos[i])]
        central_ys = data['position_y'][np.where(data['halo_id'] == id_halos[i])]
        central_zs = data['position_z'][np.where(data['halo_id'] == id_halos[i])]
        
        # establishing a dictionary
        
        if plane == 'xy':
            d = np.sqrt((xs - central_xs)**2 + (ys - central_ys)**2 )
            cyl_cut[id_halos[i]]  = np.where((d <= dist_lim) & ( (zs <= central_zs + h) & (zs >= central_zs - h) ))
        if plane == 'xz':
            d = np.sqrt((xs - central_xs)**2 + (zs - central_zs)**2 )
            cyl_cut[id_halos[i]]  = np.where((d <= dist_lim) & ( (ys <= central_ys + h) & (ys >= central_ys - h) ))
        if plane == 'yz':
            d = np.sqrt((xs - central_ys)**2 + (ys - central_ys)**2 )
            cyl_cut[id_halos[i]]  = np.where((d <= dist_lim) & ( (xs <= central_xs + h) & (xs >= central_xs - h) ))
        
    return cyl_cut
