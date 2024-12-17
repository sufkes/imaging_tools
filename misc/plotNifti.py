#!/usr/bin/env python3
import os
import sys

import nibabel as nib
from nibabel.orientations import axcodes2ornt, ornt_transform, apply_orientation

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

from skimage import measure

def reorientToRasPlus(image_nii):
    '''Convert NIFTI image to RAS+ orientation.'''

    # NIFTI images store 3D data in a 3D array. The orientation of this 3D data in real space is determined by the affine transform, and varies from image to image. Here, we convert the orientation to RAS+, in which the x-axis (axis 0) is left->right; the y-axis (axis 1) is posterior->anterior; and the z-axis (axis 2) is inferior->superior. We do this to simplify the code later on.
    current_orientation = nib.orientations.io_orientation(image_nii.affine)
    target_orientation = np.array([[0,1],[1,1],[2,1]])
    transform = nib.orientations.ornt_transform(current_orientation, target_orientation)
    data_reoriented = nib.orientations.apply_orientation(image_nii.get_fdata(), transform)
    affine_reoriented = np.dot(image_nii.affine, nib.orientations.inv_ornt_aff(transform, image_nii.shape))
    image_nii_reoriented = nib.Nifti1Image(data_reoriented, affine_reoriented)
    return image_nii_reoriented

class Image(object):
    '''A NIFTI image and associated size/shape information.'''
    def __init__(self, in_path):
        image_nii = nib.load(in_path)
        
        self.image_nii_original = image_nii # NIFTI object
        image_nii_reoriented = reorientToRasPlus(image_nii)
        self.image_nii = image_nii_reoriented
        self.data = self.image_nii.get_fdata() # 3D NumPy array.
        self.image_shape_num_voxels = self.data.shape # numbers of voxels in x, y, z directions (assumed to correspond to right-left, anterior-posterior, inferior-superior directions).
        self.voxel_shape_mm = self.image_nii.header.get_zooms() # widths of voxels in x, y, z and directions, in mm.

        ## If coronal, axial, and sagittal views are displayed in the same image, we need to set their relative sizes in order to ensure that a pixel represents the same length/width in each view.
        self.image_shape_mm = np.array(self.image_shape_num_voxels)*self.voxel_shape_mm # widths of the 3D image in x, y, z directions.

    def getClusterCentroids(self):
        if not hasattr(self, 'cluster_centroids'):
            self.setClusterCentroids()
        return self.cluster_centroids
       
    def setClusterCentroids(self):
        '''For label images, identity distinct clusters of label voxels, and return a list of their centroid coordinates.'''
        # Change the value of each label cluster to a unique value. E.g. if the input lesion image has 3 lesions all labeled 1, change the value of voxels in the second and third lesion cluster to 2 and 3, respectively.
        data_binarized = (self.data > 0).astype(int)
        label_data_numbered_clusters = measure.label(data_binarized, connectivity=3) # assign unique integer values to each cluster of labels.
        cluster_values = [n for n in np.unique(label_data_numbered_clusters) if n!= 0] # E.g. [1, 2, 3] if there are three lesion clusters in the input image.
        num_clusters = len(cluster_values)

        # Compute stats for each cluster (we only need the centroid). Returns a list of length num_clusters.
        region_props_unsorted = measure.regionprops(label_data_numbered_clusters)
        region_props = sorted(region_props_unsorted, key=lambda x: x['area_filled']) # Sort the clusters by size.
        
        centroids = [np.round(cluster_region_props.centroid).astype(int) for cluster_region_props in region_props] # Round to integer values to be used as array indices.
        self.cluster_centroids = centroids

class SubplotElement(object):
    '''A layer within a subplot.'''
    def __init__(self, image, opacity=1, color_map='gray', clipping_range_min=None, clipping_range_max=None):
        self.image = image # Image object (typical case) or None (e.g. if adding some other feature to the subplot like a legend, line etc.).
        self.opacity = opacity # Opacity of plotted image.
        self.color_map = color_map # Matplotlib color map to use.
        self.clipping_range_min = clipping_range_min # values less than or equal to this will not be plotted.
        self.clipping_range_max = clipping_range_max # values greater than or equal to this will not be plotted.

    def getData(self, plane_axis, plane_index):
        '''Return a 2D array of image data selected from the 3D image. E.g. if plane_axis=0 and plane_index=5, return self.image.data[5,:,:]'''
        plane = np.take(self.image.data, indices=plane_index, axis=plane_axis)
        return plane

class Subplot(object):
    '''One cell in the plot grid'''
    def __init__(self):
        self.subplot_element_list = [] # list of image layers for this subplot.
        self.num_layers = 0
        self.shape_mm = (0,0) # shape of the 2D data shown in this subplot in real space.

    def setSliceAndShapeParameters(self, image, view, plane):
        '''Set the axis and index which will be used to select a 2D image from the 3D data. image is an instance of the Image class; view must be "sagittal", "coronal", or "axial"; plane must be "middle" (to select the middle plane of the image) or an integer (to specify as specific image slice).'''
        if view not in ['sagittal', 'coronal', 'axial']:
            msg = "'view' must be 'sagittal', 'coronal', or 'axial'."
            raise Exception(msg)
        else:
            self.view = view        

        self.plane = plane # the index of the image slice to be plotted (along the axis to be determined based on <view>)
        # Determine aspect ratio of the subplot, and the size of the subplot relative to other subplots, using one of the images that are shown in the subplot.
        if self.view == 'sagittal':
            plane_axis = 0
            horizontal_axis = 1
            vertical_axis = 2
        elif self.view == 'coronal':
            plane_axis = 1
            horizontal_axis = 0
            vertical_axis = 2
        elif self.view == 'axial':
            plane_axis = 2
            horizontal_axis = 0
            vertical_axis = 1
        self.plane_axis = plane_axis
        self.shape_mm = (image.image_shape_mm[vertical_axis], image.image_shape_mm[horizontal_axis])
        self.aspect_ratio = image.voxel_shape_mm[vertical_axis] / image.voxel_shape_mm[horizontal_axis]
            
        if self.plane == 'middle':
            self.plane_index = image.image_shape_num_voxels[plane_axis] // 2
        else:
            self.plane_index = int(plane)            

    def addLayer(self, image, view='axial', plane='middle', opacity=1, color_map='gray', clipping_range_min=None, clipping_range_max=None, show_slice_number=True):
        '''Add an element to the subplot'''
        if (not hasattr(self, view)) or (not hasattr(self, plane)):
            # Set the view and plane index for this subplot based on the first subplot element added to it.
            self.setSliceAndShapeParameters(image, view, plane)
            self.show_slice_number = show_slice_number
        elif (view != self.view) or (plane != self.plane):
            msg = 'All images in a given subplot must have the same view.'
            raise Exception(msg)
        
        # Add the subplot element to the list for this subplot.
        subplot_element = SubplotElement(image, opacity=opacity, color_map=color_map, clipping_range_min=clipping_range_min, clipping_range_max=clipping_range_max)
        self.subplot_element_list.append(subplot_element)
        self.num_layers = len(self.subplot_element_list)

class Plot(object):
    '''A grid of subplots.'''
    def __init__(self, shape):
        self.shape = shape
        self.subplot_list = np.array([[Subplot() for i in range(shape[1])] for j in range(shape[0])]) # 2D array of Subplot instances.

    # Make it so that we can select subplots in the plot like an array (e.g. plot[0,1] for the Subplot at row 0, column 1). 
    def __getitem__(self, key):
        return self.subplot_list[key]
    
    def __setitem__(self, key, value):
        self.subplot_list[key] = value

    # Determine the relative heights of the cols and rows based on the shapes and sizes of the subplots.
    def setRealWidthsAndHeights(self):
        subplot_heights = np.zeros(shape=self.subplot_list.shape)
        subplot_widths = np.zeros(shape=self.subplot_list.shape)
        
        for row in range(self.subplot_list.shape[0]):
            for col in range(self.subplot_list.shape[1]):
                subplot = self.subplot_list[row, col]
                shape_mm = subplot.shape_mm
                subplot_heights[row, col] = shape_mm[0]
                subplot_widths[row, col] = shape_mm[1]

        real_column_heights_mm = subplot_heights.max(axis=1).flatten() # list of widths of widest subplot in each column in real life (i.e. if the image was printed to scale)
        real_row_widths_mm = subplot_widths.max(axis=0).flatten() # list of widths of widest subplot in each column (i.e. if the image was printed to scale)

        mm_per_inch = 25.4

        self.real_column_heights_inch = real_column_heights_mm/mm_per_inch
        self.real_row_widths_inch = real_row_widths_mm/mm_per_inch

        self.real_total_height_inch = self.real_column_heights_inch.sum()
        self.real_total_width_inch = self.real_row_widths_inch.sum()

    def setFigsize(self):
        # The PNG file we save has a DPI (dots per inch) and a width & height in inches. The DPI * (height, width) gives us the number of pixels.
        # In the simplest case, the figure is drawn to scale, with the DPI set to something "high enough". In practice, for a given DPI, we may need to adjust the figure size to keep the image size manageable.
        self.figsize = np.array([self.real_total_width_inch, self.real_total_height_inch])
                
    def makeFigure(self, out_path, dpi=50.8):
        ## Determine the relative widths and heights of the subplot columns and rows, respectively.
        # Find the width and height of each subplot.
        self.setRealWidthsAndHeights()
        self.setFigsize()

        # Initialize the pyplot figure and axes.
        fig, axes = plt.subplots(nrows=self.subplot_list.shape[0], ncols=self.subplot_list.shape[1], squeeze=False, figsize=self.figsize, gridspec_kw={'width_ratios': self.real_row_widths_inch, 'height_ratios':self.real_column_heights_inch}, facecolor='black')

        # Add subplots to their positions in the figure grid.
        for row in range(self.subplot_list.shape[0]):
            for col in range(self.subplot_list.shape[1]):
                # Set subplot background to black, even if there is nothing to plot, so that empty subplots will be black.
                axes[row, col].set_facecolor('black')
                
                subplot = self.subplot_list[row, col]

                # If there are no image layers in this subplot, continue to the next (row, column).
                if subplot.num_layers == 0:
                    continue
                
                aspect_ratio = subplot.aspect_ratio
                plane_axis = subplot.plane_axis
                plane_index = subplot.plane_index

                # Add each image layer to the subplot.
                for subplot_element in subplot.subplot_element_list:
                    plane = subplot_element.getData(plane_axis, plane_index) # 2D array of data to be plotted.
                    opacity = subplot_element.opacity
                    color_map = subplot_element.color_map

                    # Set values outside the clipping range to NaN so that they are not plotted. set values at or below the clipping range minimum to NaN so that they Hare not plotted.
                    clipping_range_min = subplot_element.clipping_range_min
                    clipping_range_max = subplot_element.clipping_range_max
                    if not clipping_range_min is None:
                        plane[plane <= clipping_range_min] = np.nan 
                    if not clipping_range_max is None:
                        plane[plane >= clipping_range_max] = np.nan

                    # Plot the image layer.
                    axes[row, col].imshow(np.rot90(plane), cmap=color_map, alpha=opacity, aspect=aspect_ratio)

                # Add a label indicating the plane being plotted (e.g. 'y=10')
                if subplot.show_slice_number:
                    plane_axis_name = {0:'x', 1:'y', 2:'z'}[plane_axis]

                    label_x_pos = (self.real_row_widths_inch[:col].sum() + 0.5 * self.real_row_widths_inch[col]) / self.real_row_widths_inch.sum()  # label x-coordinate relative to the figure.
                    label_y_pos = (self.real_column_heights_inch[row+1:].sum() + 0.05 * self.real_column_heights_inch[row]) / self.real_column_heights_inch.sum()

                    label_text = f'{plane_axis_name} = {plane_index}'
                    
                    # Add text at a specific position.
                    axes[row, col].text(
                        x=label_x_pos,
                        y=label_y_pos,
                        s=label_text,
                        color='white',
                        fontsize='xx-large',
                        ha='center',  # Horizontal alignment
                        va='center',  # Vertical alignment
                        transform=fig.transFigure
                    )                

        # Remove blank space between subplots and at image edges.
        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)

        # Remove the number and ticked horizontal and vertical axes.
        for axis in fig.axes:
            axis.xaxis.set_visible(False)
            axis.yaxis.set_visible(False)

        # Save the figure.
        fig.savefig(out_path, dpi=dpi)

        # Close the plt.figure object to free memory (?)
        plt.close()

# Define a custom color map for label images: just a solid color.
def solidColorMap(color):
    '''Return pyplot color map consisting of a single solid color'''
    color_map = ListedColormap([color])
    return color_map


## Define functions which generate specific types of plots from input images.
def plotThreeView(in_path, out_path, color_map='gray'):
    '''Plot single image in the sagittal, coronal, and axial views, in a 1x3 grid.'''
    image = Image(in_path)
    p = Plot(shape=(1,3))
    p[0,0].addLayer(image, view='sagittal', plane='middle', opacity=1, color_map=color_map)
    p[0,1].addLayer(image, view='coronal', plane='middle', opacity=1, color_map=color_map)
    p[0,2].addLayer(image, view='axial', plane='middle', opacity=1, color_map=color_map)
    p.makeFigure(out_path)

def plotLabelOnBackgroundImageAllLabelSlices(label_path, background_path, out_path, view='axial', max_cols=3, max_rows=None):
    '''Plot label image overlaid on a background image. Every slice containing a nonzero voxel in the label image will be plotted.'''
    label_image = Image(label_path)
    background_image = Image(background_path)
    # Determine the slices in which the label image is nonzero.
    if view == 'sagittal':
        plane_axis = 0
        other_axes = (1,2)
    elif view == 'coronal':
        plane_axis = 1
        other_axes = (0,2)
    elif view == 'axial':
        plane_axis = 2
        other_axes = (0,1)
    nonzero_label_indices_in_plane_axis = np.where((label_image.data > 0).any(axis=other_axes))[0]
    num_planes_with_nonzero_label = len(nonzero_label_indices_in_plane_axis)

    num_cols = min(num_planes_with_nonzero_label, max_cols)
    num_rows = np.ceil(num_planes_with_nonzero_label/num_cols).astype(int)
    
    p = Plot(shape=(num_rows, num_cols))
    label_color_map = solidColorMap('cyan')
    
    plane_index_index = 0
    all_planes_plotted = False
    for row in range(num_rows):
        for col in range(num_cols):
            plane_index = nonzero_label_indices_in_plane_axis[plane_index_index]

            p[row, col].addLayer(background_image, view=view, plane=plane_index)
            p[row, col].addLayer(label_image, view=view, plane=plane_index, color_map=label_color_map, clipping_range_min=0)
            
            plane_index_index += 1
            if plane_index_index == num_planes_with_nonzero_label:
                all_planes_plotted = True
                break
        if all_planes_plotted:
            break
    p.makeFigure(out_path)

def plotAlignedImagePair(image_1_path, image_2_path, out_path):
    '''Plot a pair of images in the sagittal, coronal, and axial views, in a 2x3 grid. Used as a rough check of alignment.'''
    image_1 = Image(image_1_path)
    image_2 = Image(image_2_path)
    p = Plot(shape=(3,3))
    p[0,0].addLayer(image_1, view='sagittal', plane='middle', opacity=1, color_map='gray')
    p[0,1].addLayer(image_1, view='coronal', plane='middle', opacity=1, color_map='gray')
    p[0,2].addLayer(image_1, view='axial', plane='middle', opacity=1, color_map='gray')
    
    p[1,0].addLayer(image_2, view='sagittal', plane='middle', opacity=1, color_map='gray')
    p[1,1].addLayer(image_2, view='coronal', plane='middle', opacity=1, color_map='gray')
    p[1,2].addLayer(image_2, view='axial', plane='middle', opacity=1, color_map='gray')
    p.makeFigure(out_path)    
    
def plotLabelClustersOnBackground(label_path, background_path, out_path):
    '''Plot a label image (e.g. lesion mask) overlaid on a background image (e.g. T1w image). Disconnected clusters of labels are automatically identified. For each label cluster, plot a row of sagittal, coronal and axial views of the label image overlaid on the background image, and centered on the lesion centroid.'''
    label_image = Image(label_path)
    background_image = Image(background_path)

    label_cluster_centroids = label_image.getClusterCentroids() # get a list containing the coordinates of the centroid of each cluster of labels.
    num_clusters = len(label_cluster_centroids)
    
    # Create the plot: sagittal, coronal and axial views centered on lesion cluster for each lesion cluster.
    p = Plot(shape=(num_clusters, 3))
    label_color_map = solidColorMap('cyan')
    for row, centroid in enumerate(label_cluster_centroids):
        p[row, 0].addLayer(background_image, view='sagittal', plane=centroid[0], opacity=1, color_map='gray')
        p[row, 0].addLayer(label_image, view='sagittal', plane=centroid[0], opacity=1, color_map=label_color_map, clipping_range_min=0)

        p[row, 1].addLayer(background_image, view='coronal', plane=centroid[1], opacity=1, color_map='gray')
        p[row, 1].addLayer(label_image, view='coronal', plane=centroid[1], opacity=1, color_map=label_color_map, clipping_range_min=0)

        p[row, 2].addLayer(background_image, view='axial', plane=centroid[2], opacity=1, color_map='gray')
        p[row, 2].addLayer(label_image, view='axial', plane=centroid[2], opacity=1, color_map=label_color_map, clipping_range_min=0)
    p.makeFigure(out_path)
    
def plotDisconnectomeThreeView(background_path, disconnectome_path, lesion_path, out_path):
    background_image = Image(background_path)
    lesion_image = Image(lesion_path)
    disconnectome_image = Image(disconnectome_path)
    
    # Create the plot: sagittal, coronal and axial views centered on lesion cluster for each lesion cluster.
    p = Plot(shape=(1, 3))
    label_color_map = solidColorMap('cyan')
    
    p[0, 0].addLayer(background_image, view='sagittal', plane='middle', opacity=1, color_map='gray')
    p[0, 0].addLayer(disconnectome_image, view='sagittal', plane='middle', opacity=1, color_map='YlOrRd', clipping_range_min=0)
    p[0, 0].addLayer(lesion_image, view='sagittal', plane='middle', opacity=1, color_map=label_color_map, clipping_range_min=0)

    p[0, 1].addLayer(background_image, view='coronal', plane='middle', opacity=1, color_map='gray')
    p[0, 1].addLayer(disconnectome_image, view='coronal', plane='middle', opacity=1, color_map='YlOrRd', clipping_range_min=0)
    p[0, 1].addLayer(lesion_image, view='coronal', plane='middle', opacity=1, color_map=label_color_map, clipping_range_min=0)

    p[0, 2].addLayer(background_image, view='axial', plane='middle', opacity=1, color_map='gray')
    p[0, 2].addLayer(disconnectome_image, view='axial', plane='middle', opacity=1, color_map='YlOrRd', clipping_range_min=0)
    p[0, 2].addLayer(lesion_image, view='axial', plane='middle', opacity=1, color_map=label_color_map, clipping_range_min=0)
    p.makeFigure(out_path)

def plotDisconnectomeAtLesionCentroids(background_path, disconnectome_path, lesion_path, out_path):
    background_image = Image(background_path)
    lesion_image = Image(lesion_path)
    disconnectome_image = Image(disconnectome_path)
    
    lesion_cluster_centroids = lesion_image.getClusterCentroids()
    num_clusters = len(lesion_cluster_centroids)
    
    # Create the plot: sagittal, coronal and axial views centered on lesion cluster for each lesion cluster.
    p = Plot(shape=(num_clusters, 3))
    label_color_map = solidColorMap('cyan')
    for row, centroid in enumerate(lesion_cluster_centroids):
        p[row, 0].addLayer(background_image, view='sagittal', plane=centroid[0], opacity=1, color_map='gray')
        p[row, 0].addLayer(disconnectome_image, view='sagittal', plane=centroid[0], opacity=1, color_map='YlOrRd', clipping_range_min=0)
        p[row, 0].addLayer(lesion_image, view='sagittal', plane=centroid[0], opacity=1, color_map=label_color_map, clipping_range_min=0)

        p[row, 1].addLayer(background_image, view='coronal', plane=centroid[1], opacity=1, color_map='gray')
        p[row, 1].addLayer(disconnectome_image, view='coronal', plane=centroid[1], opacity=1, color_map='YlOrRd', clipping_range_min=0)
        p[row, 1].addLayer(lesion_image, view='coronal', plane=centroid[1], opacity=1, color_map=label_color_map, clipping_range_min=0)

        p[row, 2].addLayer(background_image, view='axial', plane=centroid[2], opacity=1, color_map='gray')
        p[row, 2].addLayer(disconnectome_image, view='axial', plane=centroid[2], opacity=1, color_map='YlOrRd', clipping_range_min=0)
        p[row, 2].addLayer(lesion_image, view='axial', plane=centroid[2], opacity=1, color_map=label_color_map, clipping_range_min=0)
    p.makeFigure(out_path)

# Preview the original brain image (T1 or T2) that the user uploaded. Create a 1x3 grid of images of the middle planes in the sagittal, coronal, and axial views.
plotThreeView('sub-CC00097XX16_ses-33701_run-07_T2w.nii.gz', 'sub-CC00097XX16_ses-33701_run-07_T2w.png')

# Preview the label image overlaid on the brain image. 
plotLabelClustersOnBackground('sub-CC00097XX16_ses-33701_run-07_lesion_mask.nii.gz', 'sub-CC00097XX16_ses-33701_run-07_T2w.nii.gz', 'lesion_on_t2-clusters.png')

# Plot the age-matched template image and the warped brain image beside each other to ensure they are aligned.
plotAlignedImagePair('sub-CC00097XX16_ses-33701_run-07_T2w-43w_template_space.nii.gz', 'week43_T2w.nii.gz', 'alignment_check.png')

# Preview the warped label image overlaid on the age-matched template image. All axial slices containing a nonzero voxel in the label image will be shown.
plotLabelClustersOnBackground('sub-CC00097XX16_ses-33701_run-07_lesion_mask-43w_template_space.nii.gz', 'week43_T2w.nii.gz', 'lesion_on_age_matched_template-clusters.png')

# Plot a label image (e.g. lesion mask) overlaid on a background image (e.g. T1w image). Disconnected clusters of labels are automatically identified. For each label cluster, plot a row of sagittal, coronal and axial views of the label image overlaid on the background image, and centered on the lesion centroid.
plotDisconnectomeAtLesionCentroids('week40_T2w.nii.gz', 'disconnectome-threshold_0.nii.gz', 'sub-CC00097XX16_ses-33701_run-07_lesion_mask-40w_template_space.nii.gz', 'disconnectome_at_lesion_centroids.png')
