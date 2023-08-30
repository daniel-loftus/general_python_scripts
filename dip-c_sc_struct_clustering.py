# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:41:46 2021

@author: dal1858
"""

#%%

#For this script I will use the distance matrices of the mat and pat alleles as input to cluster single cell chromosomes 

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%



inputFile = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1\GSM4382150_cortex-p001-cb_002.20k.1.clean.locus.3dg.txt'
df = pd.read_csv(inputFile, 
                 sep = '\t', 
                 names = ['chr', 'coord', 'x', 'y', 'z'])

    
#%%

def distAlleles(inputFile, chromosome = 15, limit = False, filter_locus = False):
    
    # This function outputs the distance matrix for each allele
    
    import pandas as pd 
    from scipy.spatial.distance import pdist
    from scipy.spatial.distance import squareform
    
    df = pd.read_csv(inputFile, 
                     sep = '\t', 
                     names = ['chr', 'coord', 'x', 'y', 'z'])
    
    matBool = (df['chr'] == f'chr{chromosome}(mat)') | (df['chr'] == f'{chromosome}(mat)')
    patBool = (df['chr'] == f'chr{chromosome}(pat)') | (df['chr'] == f'{chromosome}(pat)')
    
    matDF = df[matBool]
    patDF = df[patBool]
    
    if limit:
        matDF = matDF.iloc[4:35, 2:5]
        patDF = patDF.iloc[4:35, 2:5]     
    else:
        matDF = matDF.iloc[:, 2:5]
        patDF = patDF.iloc[:, 2:5]
        
    if filter_locus:
        matDF = matDF[matDF.index.isin([7, 19, 30, 38])]
        patDF = patDF[patDF.index.isin([7, 19, 30, 38])]
    
    mat_dist_matrix = squareform(pdist(matDF))
    pat_dist_matrix = squareform(pdist(patDF))
    
    return mat_dist_matrix, pat_dist_matrix

mat, pat = distAlleles(inputFile, limit = True, filter_locus = False)
plt.imshow(mat)
plt.imshow(pat)


#%%

import os 
file_3dg_dir = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1'
os.chdir(file_3dg_dir)

coordFiles = [file for file in os.listdir() if '.3dg.txt' in file]

#as arrays
mat_arrays = []
pat_arrays = []

coord_count = 31

for file in coordFiles:
    
    mat, pat = distAlleles(f'{file_3dg_dir}\{file}', chromosome = 15, limit = True, filter_locus = False)
    
    mat_arrays.append(mat.flatten())
    pat_arrays.append(pat.flatten())
        
mat_arrays = [array for array in mat_arrays if len(array) == coord_count**2]
pat_arrays = [array for array in pat_arrays if len(array) == coord_count**2]
#%%
mat_arrays = np.vstack(mat_arrays)
pat_arrays = np.vstack(pat_arrays)

mat_len = len(mat_arrays)
pat_len = len(pat_arrays)

combined_arrays = np.concatenate((mat_arrays, pat_arrays))

#%%

# Turn into contact map instead of distance map. Make binary if or if not in contact (< some distance)

def dist_to_contact(dist_matrix, min_dist, max_dist):
    contact_matrix_1 = np.where(dist_matrix < min_dist, min_dist, dist_matrix)
    contact_matrix_2 = np.where(contact_matrix_1 > max_dist, max_dist, contact_matrix_1)
    return contact_matrix_2 

my_min = 0.5
my_max = 3

#mat_arrays = dist_to_contact(mat_arrays, my_min, my_max)
#pat_arrays = dist_to_contact(pat_arrays, my_min, my_max)
    
combined_arrays = dist_to_contact(combined_arrays, my_min, my_max)
#mat_arrays = dist_to_contact(mat_arrays, my_min, my_max)
#%%
#delete columns of zeros, which mess up finding the covariance matrix (no variance)
def remove_stdev_zero(matrix):
    
    maskBool = (np.std(matrix, axis = 0) == 0)
    idx = np.where(maskBool)[0]
    matrix = np.delete(matrix, idx, axis = 1)
    
    return matrix 

mat_arrays = remove_stdev_zero(matrix = mat_arrays)
pat_arrays = remove_stdev_zero(matrix = pat_arrays)

combined_arrays = remove_stdev_zero(matrix = combined_arrays)
#mat_arrays = remove_stdev_zero(matrix = mat_arrays)
    
#%%

#PCA

normed = (combined_arrays - combined_arrays.mean(0)) / combined_arrays.std(0)
#normed = (mat_arrays - mat_arrays.mean(0)) / mat_arrays.std(0)

cov = np.cov(normed)

eig_val, eig_vec = np.linalg.eig(cov)

idx = np.argsort(eig_val)[::-1]
sorted_eig_vectors = eig_vec[:, idx]

plt.scatter(sorted_eig_vectors[:mat_len, 0], sorted_eig_vectors[:mat_len, 1], color = 'red', label = 'Mat')
plt.scatter(sorted_eig_vectors[mat_len:, 0], sorted_eig_vectors[mat_len:, 1], color = 'blue', label = 'Pat')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Mouse Cortex Neonatal Neuron Cluster 1')
plt.legend()
#%%

combined_arrays = np.concatenate((mat_arrays, pat_arrays), axis = 0)

normed = (combined_arrays - combined_arrays.mean(0)) / combined_arrays.std(0)
cov = np.cov(normed)

eig_val, eig_vec = np.linalg.eig(cov)




idx = np.argsort(eig_val, axis=0)[::-1]
sorted_eig_vectors = eig_vec[:, idx]


plt.scatter(eig_vec[:len(mat_arrays) - 1, 0], eig_vec[0:len(mat_arrays) - 1, 1], color = 'red', label = 'Maternal')
plt.scatter(eig_vec[len(mat_arrays) - 1:, 0], eig_vec[len(mat_arrays) - 1:, 1], color = 'blue', label = 'Paternal')

plt.show()
#%%

#plot the first 3 PCs 
fig = plt.figure()
ax = plt.axes(projection='3d')


ax.scatter(eig_vec[:len(mat_arrays) - 1, 0], eig_vec[0:len(mat_arrays) - 1, 1], eig_vec[0:len(mat_arrays) - 1, 2], color = 'red', label = 'Maternal', s = 50)
ax.scatter(eig_vec[len(mat_arrays) - 1:, 0], eig_vec[len(mat_arrays) - 1:, 1], eig_vec[len(mat_arrays) - 1:, 2], color = 'blue', label = 'Paternal', s = 50)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
#%%

#animation (does not currently work)
%matplotlib qt
from matplotlib.animation import FuncAnimation


def init():
    ax.scatter(eig_vec[:len(mat_arrays) - 1, 0], eig_vec[0:len(mat_arrays) - 1, 1], eig_vec[0:len(mat_arrays) - 1, 2], 
                       color = 'red', label = 'Maternal')
    ax.scatter(eig_vec[len(mat_arrays) - 1:, 0], eig_vec[len(mat_arrays) - 1:, 1], eig_vec[len(mat_arrays) - 1:, 2], 
                       color = 'blue', label = 'Paternal')
    return fig,

def animate(i):
    ax.view_init(elev=10., azim=i)
    return fig,

# Animate
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=True)
    
#%%
cumsum = np.cumsum(eig_val[idx]) / np.sum(eig_val[idx])
xint = range(1, len(cumsum) + 1)
plt.plot(xint, cumsum)

plt.xlabel("Number of components")
plt.ylabel("Cumulative explained variance")
plt.xticks(xint)
plt.xlim(1, 5, 1)



#%%

# K-means clustering 


from sklearn.cluster import KMeans


kmeans = KMeans(n_clusters = 2)
kmeans.fit(sorted_eig_vectors)
y_means = kmeans.predict(sorted_eig_vectors)


plt.scatter(sorted_eig_vectors[:, 0], sorted_eig_vectors[:, 1], c=y_means, s=50, cmap='viridis')




#%%

# Hierarchical clustering 

data = mat_arrays

from scipy.cluster.hierarchy import dendrogram, linkage


linked = linkage(data, 'single')

labelList = range(0, len(data))

plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()





#%%

def plot3d(inputFile, chromosome = 15, allele = 'mat'):
    
    #plot the chromosome in 3d
    
    import pandas as pd 
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    
    df = pd.read_csv(inputFile, 
                     sep = '\t', 
                     names = ['chr', 'coord', 'x', 'y', 'z'])
    
    matBool = (df['chr'] == f'chr{chromosome}(mat)') | (df['chr'] == f'{chromosome}(mat)')
    patBool = (df['chr'] == f'chr{chromosome}(pat)') | (df['chr'] == f'{chromosome}(pat)')
    
    matDF = df[matBool]
    patDF = df[patBool]
    
    matDF = matDF.iloc[:, 2:5]
    patDF = patDF.iloc[:, 2:5]
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    if allele == 'mat':
        df = matDF
        title = 'Maternal'
    elif allele == 'pat':
        df = patDF
        title = 'Paternal'
        
    x = df.iloc[:, 0]
    y = df.iloc[:, 1]
    z = df.iloc[:, 2]
    
    
    col = np.arange(len(df))

    N = len(df)
    for i in range(len(df)):
        ax.plot(x[i:i+2], y[i:i+2], z[i:i+2], color=plt.cm.jet(i/len(df)), linewidth = 4)
    
    
    #plot peg13, kcnk9, and tc9enh 
    
    def plotFeature(idx, name = None, color = 'black', size = 100):
        feature_x = df.iloc[idx, 0]
        feature_y = df.iloc[idx, 1]
        feature_z = df.iloc[idx, 2]
        
        ax.scatter(feature_x, feature_y, feature_z, c = color, s = 100, label = name)
    
    plotFeature(19, name = 'Peg13', color = 'blue')
    plotFeature(7, name = 'Kcnk9', color = 'red')
    plotFeature(30, name = 'Tc9enh', color = 'green')
    plotFeature(38, name = 'Ago2', color = 'Pink')
    ax.legend()
    plt.title(title)
    
    return df.iloc[:, 0:3]
    

mat_coords = plot3d(inputFile)
pat_coords = plot3d(inputFile, allele = 'pat')



#%%

#icp 

import icp 

my_icp = icp.icp(mat_coords, pat_coords)
transformation = my_icp[0]
rotation = transformation[0:3, 0:3]
translation = transformation[0:3, 3]
rotated_mat = np.matmul(mat_coords, rotation)
moved_mat = rotated_mat + translation 

ax = plt.axes(projection='3d')
#ax.plot(mat_coords.iloc[:,0], mat_coords.iloc[:,1], mat_coords.iloc[:,2], color = 'red')
#ax.plot(rotated_mat.iloc[:,0], rotated_mat.iloc[:,1], rotated_mat.iloc[:,2], color = 'grey')
ax.plot(moved_mat.iloc[:,0], moved_mat.iloc[:,1], moved_mat.iloc[:,2], color = 'darkred')
ax.plot(pat_coords.iloc[:,0], pat_coords.iloc[:,1], pat_coords.iloc[:,2], color = 'blue')

def icp_coords(coords_1, coords_2):
    
    import icp 
    
    my_icp = icp.icp(coords_1, coords_2, tolerance = 0.000001)
    transformation = my_icp[0]
    rotation = transformation[0:3, 0:3]
    translation = transformation[0:3, 3]
    transformed = np.matmul(coords_1, rotation) + translation

    ax = plt.axes(projection='3d')
    ax.plot(transformed.iloc[:,0], transformed.iloc[:,1], transformed.iloc[:,2], color = 'darkred')
    ax.plot(coords_2.iloc[:,0], coords_2.iloc[:,1], coords_2.iloc[:,2], color = 'blue')

icp_coords(mat_coords, pat_coords)
#%%

file1 = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1\GSM4382150_cortex-p001-cb_002.20k.1.clean.locus.3dg.txt'
file2 = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1\GSM4382156_cortex-p001-cb_008.20k.1.clean.locus.3dg.txt'
file3 = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1\GSM4382158_cortex-p001-cb_010.20k.1.clean.locus.3dg.txt'
file4 = r'W:\Users\dloftus\tan_etal_2020_dipc_mouse_hippocampus\cortex\clean_chr15_3dg\locus\1\Neonatal Neuron 1\GSM4382160_cortex-p001-cb_012.20k.1.clean.locus.3dg.txt'

for file in [file1, file2, file3, file4]:
    
    plot3d(file)
    plot3d(file, allele = 'pat')
    
    

#%%
x = df.iloc[:, 0]
y = df.iloc[:, 1]
z = df.iloc[:, 2]

np.array([x, y, z]).T.reshape(-1, 1, 3)






#%%

#t-sne? 

from sklearn.manifold import TSNE

my_tsne = TSNE(n_components = 3, perplexity = 5, init = 'pca').fit_transform(combined_arrays)
plt.scatter(my_tsne[:,0], my_tsne[:,1])


#%%

#plot t-sne 3d

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.scatter(my_tsne[:,0], my_tsne[:,1], my_tsne[:,2])






#%%

# What about using the n closest points to a given point as input? 




























#%%

# Iterative Closest Points (source = https://github.com/ClayFlannigan/icp/blob/master/icp.py)

