# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:57:52 2021

@author: dal1858
"""

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%

def cisElementsRadiusGyration(input_file, region_list = False, chromosome = 15):
    
    '''
    This function takes Dip-C data and a list of regions of interest as input,
    and then outputs the radius of gyration of those regions (both mat and pat)
    '''
    
    import pandas as pd 
    
    #read 3dg file into a df
    df = pd.read_csv(input_file, 
                     sep = '\t', 
                     names = ['chr', 'coord', 'x', 'y', 'z'])
    
    #filter for specified allele
    mat_bool = (df['chr'] == f'chr{chromosome}(mat)') | (df['chr'] == f'{chromosome}(mat)')
    pat_bool = (df['chr'] == f'chr{chromosome}(pat)') | (df['chr'] == f'{chromosome}(pat)')
    mat_df = df[mat_bool]
    pat_df = df[pat_bool]
    
    #filter allelic df for regions specified in region_list using the custom filter3DGRegions function 
    if region_list:
        filtered_mat_df = filter3DGRegions(mat_df, region_list)
        filtered_pat_df = filter3DGRegions(pat_df, region_list)
        radius_gyration_mat = radiusGyration(filtered_mat_df)
        radius_gyration_pat = radiusGyration(filtered_pat_df)
    else:
        radius_gyration_mat = radiusGyration(mat_df)
        radius_gyration_pat = radiusGyration(pat_df)
        
    return radius_gyration_mat, radius_gyration_pat

#mat_temp, pat_temp = cisElementsRadiusGyration(my_file, my_regions_list)


#%%

def cisElementsRadiusGyrationAlleleCompare(input_3dg_dir, region_list = False, alpha = 0.4, cell_type = None, hist_bins = 40):
    
    '''
    This function performs cisElementsRadiusGyration on all 3dg files in the 
    specified directory and plots and mat and pat distributions of 
    radius of gyration for the specified loci
    '''
    
    import os 
    
    os.chdir(input_3dg_dir)
    
    coordFiles = [file for file in os.listdir() if '.3dg.txt' in file]
    
    mat_radius_list = []
    pat_radius_list = []
    
    for file in coordFiles: 
        
        try:
            mat_radius, pat_radius = cisElementsRadiusGyration(f'{input_3dg_dir}\\{file}', region_list = region_list)
            mat_radius_list.append(mat_radius * 60) # 1 particle radius ~ 60 nm 
            pat_radius_list.append(pat_radius * 60)
        except KeyError: 
            print(f'Key error. File = {file}')  
    
    import matplotlib.pyplot as plt
    import scipy.stats as stats  
    import statistics 
    
    mat_median = statistics.median(mat_radius_list)
    pat_median = statistics.median(pat_radius_list)
    
    test_statistic, p_value = stats.ks_2samp(mat_radius_list, pat_radius_list, alternative = 'two-sided', mode = 'auto')
    p_value = '{:.2e}'.format(p_value) #convert p-value to 2 digits past decimal 
    
    '''
    #plotting
    if cell_type:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb), {cell_type}'
    else:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb)'
    '''
    fig = plt.figure(figsize = (7, 10))
    plt.subplots_adjust(hspace = 0.3)
    density = fig.add_subplot(2, 1, 1)
    ecdf = fig.add_subplot(2, 1, 2)
    
    density.axvline(mat_median, label = 'Mat median', color = 'darkred', linestyle = '--')
    density.axvline(pat_median, label = 'Pat median', color = 'darkblue', linestyle = '--')
    density.hist(mat_radius_list, hist_bins, color = 'red', alpha = alpha, density = True, label = 'Mat')
    density.hist(pat_radius_list, hist_bins, color = 'blue', alpha = alpha, density = True, label = 'Pat')
    density.legend()
    density.set_ylabel('Probability density')
    density.set_xlabel('Radius of Gyration (nm)')
    #density.set_title(title)
    #density.set_xlim(0, xlim)
    
    import numpy as np
    y = np.arange(0, len(mat_radius_list) ) / len(mat_radius_list)
    
    ecdf.scatter(np.sort(mat_radius_list), y, color = 'red', alpha = 0.4, label = 'Mat')
    ecdf.scatter(np.sort(pat_radius_list), y, color = 'blue', alpha = 0.4, label = 'Pat')
    ecdf.legend()
    ecdf.set_ylabel('Cumulative distribution function')
    ecdf.set_xlabel('Radius of Gyration (nm)')
    #ecdf.set_title(title)
    #ecdf.set_xlim(0, xlim)
    ecdf.text(11, 0.82, f'p-value: {p_value}')
    #return fig 
    
    
    
    
#%%

#Regions 
    
class Locus: 
    def __init__(self, name, chromosome, coord):
        self.name = name 
        self.chromosome = chromosome
        self.coord = coord
    
    def coord_name(self):
        out = '{:.2f}'.format(self.coord/1e6)
        return out

#kcnk9_2 with tc9enh1 generally yields the lowest p-value 
dmr1 = Locus('DMR', 15, 72800000)
dmr2 = Locus('DMR', 15, 72820000)
tc9enh1 = Locus('Tc9enh', 15, 73020000)
tc9enh2 = Locus('Tc9enh', 15, 73040000)
kcnk9_1 = Locus('Kcnk9', 15, 72540000)
kcnk9_2 = Locus('Kcnk9', 15, 72560000)
ago2_1 = Locus('Ago2', 15, 73180000)
ago2_2 = Locus('Ago2', 15, 73200000)
down_contact_1 = Locus('Distal_1', 15, 72200000)
down_contact_2 = Locus('Distal_2', 15, 72080000)
down_contact_3 = Locus('Distal_3', 15, 71920000)



far_up_1 = Locus('Control', 15, 70200000)
far_up_2 = Locus('Control', 15, 70700000)
far_up_3 = Locus('Control', 15, 71200000)
far_down_1 = Locus('Control', 15, 75800000)
far_down_2 = Locus('Control', 15, 70020000)
far_down_3 = Locus('Control', 15, 70200000)
distal_ctcf = Locus('Distal CTCF site', 15, 72100000)

#%%

#additional functions

my_regions_list = [tc9enh1, kcnk9_2, ago2_1, down_contact_1, down_contact_2, down_contact_3]
my_regions_list_1 = [tc9enh1, kcnk9_2, ago2_1]

my_file = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus_large_1\GSM4382149_cortex-p001-cb_001.20k.1.clean..locus.3dg.txt'
my_dir = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus_large_1'
small_dir = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus'

def filter3DGRegions(input_3dg_df, region_list):
    
    '''
    This function takes as input a 3DG file and a list of Locus objects
    defining regions of interest, and outputs a filtered 3DG file containing 
    only the specified regions.
    '''
    
    coord = []
    for region in region_list:
        coord.append(region.coord)
        
    col_names = ['coord']
    regions_df = pd.DataFrame(coord, columns = col_names)
    output_df = pd.merge(input_3dg_df, regions_df, on = ['coord'])
    
    return output_df


def radiusGyration(df_3dg):
    
    '''
    This function takes a 3dg dataframe and calculates the radius of gyration
    of the contained genomic elements
    '''
    
    centroid_x = sum(df_3dg.loc[:]['x'])/len(df_3dg)
    centroid_y = sum(df_3dg.loc[:]['y'])/len(df_3dg)
    centroid_z = sum(df_3dg.loc[:]['z'])/len(df_3dg)

    def dist(x, y, z):
        '''
        Finds the distance between 3d points and the centroid
        '''
        output_dist = ( (x - centroid_x)**2 + (y - centroid_y)**2 + (z - centroid_z)**2 )**(1/2)
        return output_dist
    
    df_3dg['dist_from_centroid'] = df_3dg.apply(lambda row: dist(row['x'], row['y'], row['z']), axis = 1)
    
    radius_gyration = (sum(df_3dg['dist_from_centroid'])/len(df_3dg) )**(1/2)
    
    return radius_gyration











