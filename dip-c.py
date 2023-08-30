# -*- coding: utf-8 -*-
"""
Created on Sat May 30 11:40:10 2020

@author: dal1858
"""

#%%

import pandas as pd 
import os

#%%

def distAlleles(inputFile, chromosome, coord_A, coord_B):
    
    # This function outputs the distance between the given genomic coordinates for both alleles 
    
    import pandas as pd 
    
    df = pd.read_csv(inputFile, 
                     sep = '\t', 
                     names = ['chr', 'coord', 'x', 'y', 'z'])
    
    matBool = (df['chr'] == f'chr{chromosome}(mat)') | (df['chr'] == f'{chromosome}(mat)')
    patBool = (df['chr'] == f'chr{chromosome}(pat)') | (df['chr'] == f'{chromosome}(pat)')
    
    matDF = df[matBool]
    patDF = df[patBool]
    
    matCoordsBool = (matDF['coord'] == coord_A) | (matDF['coord'] == coord_B)
    patCoordsBool = (patDF['coord'] == coord_A) | (patDF['coord'] == coord_B)
    
    matCoords = matDF[matCoordsBool]
    patCoords = patDF[patCoordsBool]
    matCoords = matCoords.reset_index()
    patCoords = patCoords.reset_index()
    
    def dist(coords_df):
        
        euc_dist = ( (coords_df.loc[0]['x'] - coords_df.loc[1]['x'])**2 + (coords_df.loc[0]['y'] - coords_df.loc[1]['y'])**2 + (coords_df.loc[0]['z'] - coords_df.loc[1]['z'])**2)**(1/2)
        return euc_dist 
        
    mat_dist = dist(matCoords)
    pat_dist = dist(patCoords)
    
    return mat_dist, pat_dist

#%%
    
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
far_up_1 = Locus('Control', 15, 70200000)
far_up_2 = Locus('Control', 15, 70700000)
far_up_3 = Locus('Control', 15, 71200000)
far_down_1 = Locus('Control', 15, 75800000)
far_down_2 = Locus('Control', 15, 70020000)
far_down_3 = Locus('Control', 15, 70200000)
down_contact_1 = Locus('Distal_1', 15, 72200000)
down_contact_2 = Locus('Distal_2', 15, 72080000)
down_contact_3 = Locus('Distal_3', 15, 71920000)
distal_ctcf = Locus('Distal CTCF site', 15, 72100000)

h19 = Locus('H19', 7, 142560000)
igf2 = Locus('Igf2', 7, 142660000)
ifitm10 = Locus('Ifitm10', 7, 142340000)
kcnq1ot1 = Locus('Kcnq1ot1', 7, 142340000)

actb_1 = Locus('Actb locus control', 5, 142200000)
actb_2 = Locus('Actb locus control', 5, 142700000)
actb_3 = Locus('Actb locus control', 5, 143000000)
actb_4 = Locus('Actb locus control', 5, 143700000)
actb_5 = Locus('Actb locus control', 5, 143800000)
actb_start = Locus('Actb', 5, 142020000)
actb_end = Locus('Actb', 5, 143980000)

tc9_int_1 = Locus('Trappc9 upstream intron', 15, 72600000)
tc9_int_2 = Locus('Trappc9 upstream intron', 15, 72700000)
tc9_int_3 = Locus('Trappc9 downstream intron', 15, 72900000)
tc9_int_4 = Locus('Trappc9 downstream intron', 15, 73000000)

human_dmr = Locus('DMR', 8, 141100000)
human_kcnk9 = Locus('Kcnk9 promoter', 8, 140700000)
human_kcnk9_intron = Locus('Kcnk9 intron', 8, 140660000)
human_tc9enh = Locus('Tc9enh', 8, 141400000)

locusA = kcnk9_1
locusB = tc9enh1

#myDir = r'W:\Users\dloftus\tan_etal_2019_single_neuron_HiC_mm10\locus_3dg\MOE'

group1 = [ [dmr1, tc9enh1], [kcnk9_1, tc9enh1], [kcnk9_2, tc9enh1], [kcnk9_1, dmr1], [kcnk9_2, dmr1], [kcnk9_1, ago2_1], [kcnk9_1, ago2_2] ]
group2 = [ [dmr2, tc9enh1], [dmr2, ago2_1], [far_up_1, far_up_2], [far_down_1, far_down_2], [far_down_2, far_down_3], [tc9enh1, distal_ctcf] ]
group3 = [ [far_up_1, far_up_3], [far_up_1, distal_ctcf]]
group4 = [ [h19, igf2], [ifitm10, h19], [ifitm10, igf2]]

actb_group = [ [actb_1, actb_2], [actb_2, actb_3], [actb_1, actb_3], [actb_1, actb_5], [actb_4, actb_5], [actb_3, actb_4]  ]

trappc9_intron_group = [ [tc9_int_1, tc9_int_2], [tc9_int_2, tc9_int_3], [tc9_int_3, tc9_int_4], [tc9_int_1, tc9_int_3], [tc9_int_2, tc9_int_4] ]



radial_1 = [dmr1, tc9enh1, kcnk9_1, ago2_1]
radial_2 = [h19, igf2, kcnq1ot1]

#%%
for pair in [ [human_dmr, human_kcnk9_intron], [human_kcnk9_intron, human_tc9enh] ]:
    
    alleleCompare(myDir, pair[0], pair[1], 'GM12878', xlim = 1000)

#%%


def alleleCompare(file_3dg_dir, locusA, locusB, cell_type = None, xlim = 700, hist_bins = 40, alpha = 0.4, particle_radius_nm = 100, hist = False, export_svg = False, svg_name = 'name'):
    
    import os 
    import matplotlib.pyplot as plt
    import scipy.stats as stats  
    import statistics 
    
    os.chdir(file_3dg_dir)
    
    coordFiles = [file for file in os.listdir() if '.3dg.txt' in file]
    
    mat_dist_list = []
    pat_dist_list = []
    


    for file in coordFiles: 
        
        try:
            mat_dist, pat_dist = distAlleles(f'{file_3dg_dir}\\{file}', locusA.chromosome, locusA.coord, locusB.coord)
            mat_dist_list.append(mat_dist * particle_radius_nm) # 1 particle radius ~ 60 nm 
            pat_dist_list.append(pat_dist * particle_radius_nm)
        except KeyError: 
            print(f'Key error. File = {file}')
    
    

    test_statistic, p_value = stats.ks_2samp(mat_dist_list, pat_dist_list, alternative = 'two-sided', mode = 'auto')
    p_value = '{:.2e}'.format(p_value) #convert p-value to 2 digits past decimal 
    
    
    #plotting
    if cell_type:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb), {cell_type}'
    else:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb)'
    
    fig = plt.figure(figsize = (7, 10))
    plt.subplots_adjust(hspace = 0.3)

    ecdf = fig.add_subplot(2, 1, 2)

    if hist:
        mat_median = statistics.median(mat_dist_list)
        pat_median = statistics.median(pat_dist_list)
        
        density = fig.add_subplot(2, 1, 1)
        
        density.axvline(mat_median, label = 'Mat median', color = 'darkred', linestyle = '--')
        density.axvline(pat_median, label = 'Pat median', color = 'darkblue', linestyle = '--')
        density.hist(mat_dist_list, hist_bins, color = 'red', alpha = alpha, density = True, label = 'Mat')
        density.hist(pat_dist_list, hist_bins, color = 'blue', alpha = alpha, density = True, label = 'Pat')
        density.legend()
        density.set_ylabel('Probability density')
        density.set_xlabel('Distance (nm)')
        density.set_title(title)
        density.set_xlim(0, xlim)
    
    import numpy as np
    y = np.arange(0, len(mat_dist_list) ) / len(mat_dist_list)
    
    ecdf.scatter(np.sort(mat_dist_list), y, color = 'red', alpha = 0.4, label = 'Mat')
    ecdf.scatter(np.sort(pat_dist_list), y, color = 'blue', alpha = 0.4, label = 'Pat')
    ecdf.legend()
    ecdf.set_ylabel('Cumulative distribution function')
    ecdf.set_xlabel('Distance (nm)')
    ecdf.set_title(title)
    ecdf.set_xlim(0, xlim)
    ecdf.text(11, 0.82, f'p-value: {p_value}')
    
    if export_svg == True:
        plt.savefig(f"{svg_name}.svg")
        
    #return fig
    
    
    


def radialDistCompare(file_3dg_dir, locusA, cell_type = None, xlim = 2700, hist_bins = 40, alpha = 0.4):
    
    import os 
    
    os.chdir(file_3dg_dir)
    
    coordFiles = [file for file in os.listdir() if '.3dg.txt' in file]
    
    mat_dist_list = []
    pat_dist_list = []
    


    for file in coordFiles: 
        
        try:
            mat_dist, pat_dist = radialDist(f'{file_3dg_dir}\\{file}', locusA.chromosome, locusA.coord)
            mat_dist_list.append(mat_dist * 60) # 1 particle radius ~ 60 nm 
            pat_dist_list.append(pat_dist * 60)
        except KeyError: 
            print(f'Key error. File = {file}')
    
    import matplotlib.pyplot as plt
    import scipy.stats as stats  
    import statistics 
    
    
    mat_median = statistics.median(mat_dist_list)
    pat_median = statistics.median(pat_dist_list)
    
    test_statistic, p_value = stats.ks_2samp(mat_dist_list, pat_dist_list, alternative = 'two-sided', mode = 'auto')
    p_value = '{:.2e}'.format(p_value) #convert p-value to 2 digits past decimal 
    
    
    #plotting
    if cell_type:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb), {cell_type}'
    else:
        title = f'{locusA.name} (Chr{locusA.chromosome}: {locusA.coord_name()} Mb) - {locusB.name} (Chr{locusA.chromosome}: {locusB.coord_name()} Mb)'

    fig = plt.figure(figsize = (7, 10))
    plt.subplots_adjust(hspace = 0.3)
    density = fig.add_subplot(2, 1, 1)
    ecdf = fig.add_subplot(2, 1, 2)
    
    density.axvline(mat_median, label = 'Mat median', color = 'darkred', linestyle = '--')
    density.axvline(pat_median, label = 'Pat median', color = 'darkblue', linestyle = '--')
    density.hist(mat_dist_list, hist_bins, color = 'red', alpha = alpha, density = True, label = 'Mat')
    density.hist(pat_dist_list, hist_bins, color = 'blue', alpha = alpha, density = True, label = 'Pat')
    density.legend()
    density.set_ylabel('Probability density')
    density.set_xlabel('Radial distance (nm)')
    density.set_title(title)
    density.set_xlim(0, xlim)
    
    import numpy as np
    y = np.arange(0, len(mat_dist_list) ) / len(mat_dist_list)
    
    ecdf.scatter(np.sort(mat_dist_list), y, color = 'red', alpha = 0.4, label = 'Mat')
    ecdf.scatter(np.sort(pat_dist_list), y, color = 'blue', alpha = 0.4, label = 'Pat')
    ecdf.legend()
    ecdf.set_ylabel('Cumulative distribution function')
    ecdf.set_xlabel('Radial distance (nm)')
    ecdf.set_title(title)
    ecdf.set_xlim(0, xlim)
    ecdf.text(11, 0.82, f'p-value: {p_value}')
    



#%%



def radialDist(inputFile, chromosome, coord):
    
    import pandas as pd 
    
    df = pd.read_csv(inputFile, 
                     sep = '\t', 
                     names = ['chr', 'coord', 'x', 'y', 'z'])
    
    centroid_x = sum(df.loc[:]['x'])/len(df)
    centroid_y = sum(df.loc[:]['y'])/len(df)
    centroid_z = sum(df.loc[:]['z'])/len(df)
    
    matBool = (df['chr'] == f'chr{chromosome}(mat)')
    patBool = (df['chr'] == f'chr{chromosome}(pat)')
    
    matDF = df[matBool]
    patDF = df[patBool]
    
    matCoordsBool = (matDF['coord'] == coord)
    patCoordsBool = (patDF['coord'] == coord)
    
    matCoords = matDF[matCoordsBool]
    patCoords = patDF[patCoordsBool]
    matCoords = matCoords.reset_index()
    patCoords = patCoords.reset_index()
    
    def dist(coords_df):
        
        euc_dist = ( (coords_df.loc[0]['x'] - centroid_x)**2 + (coords_df.loc[0]['y'] - centroid_y)**2 + (coords_df.loc[0]['z'] - centroid_z)**2)**(1/2)
        return euc_dist 
        
    mat_dist = dist(matCoords)
    pat_dist = dist(patCoords)
    
    
    return mat_dist, pat_dist   
    
    
    
    
    
    
#%%
    
my_file = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\2\cell_info.txt'
cells = pd.read_csv(my_file, sep = '\t', encoding='cp1252')

my_dir = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1'
for cell_type in cells['cell-type cluster'].unique():
    
    try:
        os.chdir(f'{my_dir}\{cell_type}')
    except FileNotFoundError:
        pass 
    
    try:
        alleleCompare(os.getcwd(), kcnk9_2, tc9enh1, cell_type = cell_type)
    except ValueError:
        pass
    

#%%
        
#alleleCompare(os.getcwd(), kcnk9_2, tc9enh1, cell_type = 'All cells', export_svg = True, svg_name = 'kcnk9_tc9enh_cdf')
        
#alleleCompare(os.getcwd(), dmr1, tc9enh2, cell_type = 'All cells', export_svg = True, svg_name = 'dmr_tc9enh_cdf')

#alleleCompare(os.getcwd(), dmr1, kcnk9_2, cell_type = 'All cells', export_svg = False, svg_name = 'dmr_tc9enh_cdf')






alleleCompare(os.getcwd(), actb_start, actb_end, xlim = 1000)

