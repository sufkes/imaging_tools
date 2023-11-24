#!/usr/bin/env python3

import os
import sys
import argparse
from collections import OrderedDict
from copy import deepcopy
import warnings
import datetime
import pickle

import numpy as np
import pandas as pd
import bct.algorithms as bct
    ## Matrix requirements for the Brain Connectivity Toolbox (BCT), from their website:
    # """ 
    # Which network matrices can I use with the Brain Connectivity Toolbox? 
    # - Most functions do not explicitly check the validity of the input network matrices; it is crucial to manually ensure that these matrices are suitable for their intended use.
    # - The network matrices should be square; rows and columns in these matrices should represent network nodes, matrix entries should represent network links.
    # - The network matrices should not be too small. As a rule of thumb, the toolbox is designed to be used with networks of greater than 20 nodes.
    # - The network matrices should preferably be in double-precision and non-sparse formats. Sparse, single-precision or logical formats may sometimes cause errors.
    # - The network matrices may be binary or weighted, directed or undirected. Each function specifies the network type for which it is suitable.
    # - The network matrices should not contain self-self connections. In other words, all values on the main diagonal of these matrices should be set to 0.
    # - In most cases, the network matrices should not contain negative weights. However, a substantial number of functions can process matrices with positive and negative weights. These functions typically end with sign.m (for signed networks).
    # - In general, randomization functions are designed for non-dense matrices; many randomization functions will be too slow and/or ineffective in dense matrices. However, some randomization functions are specifically designed for dense and weighted matrices.

    ## Generate a binary undirected network matrix using a density threshold.
    #
    # network density = (number of connections)/(number of possible connections)
    #
    #   number of connections = (sum of adjacency matrix 1s)/2 = (sum of node degrees)/2
    #   number of possible connections = n*(n-1)/2
    # 
    # The adjacency matrix is nxn, so we need to find the correlation threshold such that the adjacency matrix filled with a ratio of (n-1)/n*(density).
    # I.e. The adjacency matrix has n zeros which do not correspond to possible connections, so the proportion of nonzero elements in the adjacency matrix will be slightly less than the network density (a factor of (n-1)/2 less).

def thresholdMatrix(matrix, density_threshold):
    num_nodes = matrix.shape[0]
    percentile = ( 1 - density_threshold * float(num_nodes-1)/float(num_nodes) ) * 100.0
    correlation_threshold = np.percentile(matrix, percentile)

    w_adj = matrix.copy()
    w_adj[w_adj < correlation_threshold] = 0

    b_adj = w_adj.copy()
    b_adj[b_adj>0] = 1

    return b_adj, w_adj

def saveMetric(metric, full_out_dir, out_name):
    out_path = os.path.join(full_out_dir, out_name)
    with open(out_path, 'wb') as file_handle:
        np.save(file_handle, metric)    

def getNetworkDensity(adjacency_matrix):
    num_nodes = adjacency_matrix.shape[0]
    num_actual_connections = np.count_nonzero(adjacency_matrix)/2 # do not double count
    num_possible_connections = num_nodes*(num_nodes-1)/2
    network_density = float(num_actual_connections)/float(num_possible_connections)
    return network_density

def applyWeightsRandomly(b_rand, w_adj, seed=489):
    ## Take the random binary network and apply the weights in the real weighted network matrix. This should yield a weighted network with the same number of nodes, same binary degree distribution, and same edge weight distribution as the original weighted network.
    # Get the upper triangular portion of both matrices.
    b_rand_upper = np.triu(b_rand)
    w_adj_upper = np.triu(w_adj)

    # Get flattened versions of the matrices.
    b_rand_upper_flat = b_rand_upper.flatten()
    w_adj_upper_flat = w_adj_upper.flatten()

    # Get a flat array of the nonzero weights in the weighted matrix.
    nonzero_weights = w_adj_upper_flat[w_adj_upper_flat>0]

    # Randomize the weights.
    randomized_weights = np.random.permutation(nonzero_weights)

    # Randomly replace the binary connections in the random matrix with the weights from the real weighted network. 
    assert len(b_rand_upper_flat[b_rand_upper_flat>0]) == len(randomized_weights)
    b_rand_upper_flat[b_rand_upper_flat>0] = randomized_weights
    
    # Build the randomized weighted network.
    w_rand = b_rand_upper_flat.reshape(w_adj.shape)
    w_rand += w_rand.T

    return w_rand

def generateRandomNetworks(b_adj, w_adj, num_random_networks=100, custom_weighted_algorithm=True, debug=False):
    "Generate random networks."
    # For binary matrices, want to preserve number of nodes and degree distribution. Use randmio_und()
    # For weighted matrices, want to preserve number of nodes and weight distribution. Use null_model_und_sign()
    b_rand_list = [] # list of random binary networks
    w_rand_list = [] # list of random weighted networks

    for ii in range(num_random_networks):
        b_rand = bct.randmio_und(b_adj, itr=5)[0] # Need to specify appromixate number of times network is rewired; default is 5 for weighted case, so let's use that.

        # null_model_und_sign Takes a very long time; not sure why but seems like implementing the method in this paper should be easy and much quicker.
        #   Colon-Perez LM, Couret M, Triplett W, Price CC and Mareci TH (2016) Small Worldness in Dense and Weighted Connectomes. Front. Phys. 4:14. doi: 10.3389/fphy.2016.00014
        # According to this paper, I should be able to just distribute the weights in the real network randomly onto each of the connections in the random binary network.
        # - Checks of degree and weight distribution are similar for both methods.
        # - Characteristic paths lengths are similar between methods.
        
        if custom_weighted_algorithm:
            w_rand = applyWeightsRandomly(b_rand, w_adj)
        else:
            w_rand = bct.null_model_und_sign(w_adj)[0] # very slow
            
        b_rand_list.append(b_rand)
        w_rand_list.append(w_rand)
            
        # Sanity check: Ensure the node degree distribution for the random binary network matches the real network.
        if debug:
            # Binary: check degree distribution.
            sorted_degrees_real_b = sorted(list(b_adj.sum(axis=0)))
            sorted_degrees_rand_b = sorted(list(b_rand.sum(axis=0)))
            assert sorted_degrees_real_b == sorted_degrees_rand_b

            # Weighted: check binary degree distribution.                    
            sorted_bin_degrees_real_w = sorted(list(np.count_nonzero(w_adj, axis=0)))
            sorted_bin_degrees_rand_w = sorted(list(np.count_nonzero(w_rand, axis=0)))
            if sorted_bin_degrees_real_w != sorted_bin_degrees_rand_w:
                difference = [x-y for x,y in zip(sorted_bin_degrees_real_w,sorted_bin_degrees_rand_w)]
                print(pd.DataFrame({'real':sorted_bin_degrees_real_w, 'rand':sorted_bin_degrees_rand_w, 'difference':difference}))
            assert sorted_bin_degrees_real_w == sorted_bin_degrees_rand_w
            
            # Weighted: check edge strength distribution.
            sorted_weights_real_w = sorted(list(w_adj[w_adj>0].flatten()))
            sorted_weights_rand_w = sorted(list(w_rand[w_rand>0].flatten()))
            assert sorted_weights_real_w == sorted_weights_rand_w

            # Weighted: compare weighted degree distribution (but this need not be preserved between the real and random network).
            sorted_degrees_real_w = sorted(list(w_adj.sum(axis=0)))
            sorted_degrees_rand_w = sorted(list(w_rand.sum(axis=0)))
            #difference = [x-y for x,y in zip(sorted_degrees_real_w,sorted_degrees_rand_w)]
            #print pd.DataFrame({'real':sorted_degrees_real_w, 'rand':sorted_degrees_rand_w, 'difference':difference})[['real','rand','difference']] 
    return b_rand_list, w_rand_list

def myCharPathBin(b_adj, custom_algorithm=False):
    #distance_binary_floyd, num_edges_shortest_path_binary, _ = bct.distance.distance_wei_floyd(b_adj, transform='inv') # same result as distance_bin
    distance_binary = bct.distance.distance_bin(b_adj)
    assert (np.diagonal(np.abs(b_adj)) == 0).all() # check that diagnoal is zero.
    #print "Disconnected graph: ", (np.triu(distance_binary, k=1) == np.inf).any()
    if custom_algorithm:
        n = distance_binary.shape[0]
        sum_inverses = 2.0*np.triu(1.0/distance_binary, k=1).sum() # exclude diagnoal; double-count
        charpath_binary = float(1)/(sum_inverses/(n*(n-1)))
    else:
        charpath_binary, efficiency_binary, ecc_binary, radius_binary, diameter_binary = bct.distance.charpath(distance_binary, include_diagonal=False, include_infinite=False)
    
    return charpath_binary

def myCharPathWei(w_adj, custom_algorithm=False):
    #w_length = w_adj.copy() # "length" matrix; 1/weight for each weight in the weight matrix.
    #w_length = 1.0 / w_length
    #distance_weighted, num_edges_shortest_path_weighted = bct.distance.distance_wei(w_length)
    distance_weighted_floyd, num_edges_shortest_path_weighted, _ = bct.distance.distance_wei_floyd(w_adj, transform='inv') # same result as distance_wei (after inverting elements as required)
    assert (np.diagonal(np.abs(w_adj)) == 0).all() # check that diagnoal is zero.
    #print "Disconnected graph: ", (np.triu(distance_weighted_floyd, k=1) == np.inf).any()
    #diff = distance_weighted_floyd - distance_weighted
    #if (distance_weighted_floyd != distance_weighted).any():
        #print 'Normal distance and Floyd distance differ (weighted):'
        #print 'Normal:'
        #print distance_weighted[:6,:6]
        #print 
        #print 'Floyd-Warshall:'
        #print distance_weighted_floyd[:6,:6]
        #print
        #print 'Difference:'
        #print diff[:6,:6]
        #print 'biggest difference:'
        #print np.nanmax(diff)
        #print
    if custom_algorithm:
        n = distance_weighted_floyd.shape[0]
        sum_inverses = 2.0*np.triu(1.0/distance_weighted_floyd, k=1).sum() # exclude diagnoal; double-count
        charpath_weighted = float(1)/(sum_inverses/(n*(n-1)))
    else:
        charpath_weighted, efficiency_weighted, ecc_weighted, radius_weighted, diameter_weighted = bct.distance.charpath(distance_weighted_floyd, include_diagonal=False, include_infinite=False)
    return charpath_weighted

def normalizer(adj, rand_list, function, function_kwargs={}, scalar=True):
    """Returns metric value, the value normalized over a list of random networks, and the mean of the metric over those networks.
If 'function' returns a tuple, only the first value is taken and treated as a metric; the rest are discarded."""
    metric = function(adj, **function_kwargs)
    
    rand_sum = 0 # can be a vector or a scala
    for rand in rand_list:
        rand_metric = function(rand, **function_kwargs)
        if type(rand_metric) == tuple:
            rand_metric = rand_metric[0]
        rand_sum += rand_metric

    if scalar: # if metric is a scalar
        rand_mean = float(rand_sum)/float(len(rand_list))
        metric_normalized = float(metric)/float(rand_mean)
        return metric, metric_normalized, rand_mean
    else: # if metric is a vector
        assert len(rand_sum.shape) == 1 # check that this is a vector not a higher dimensional array.
        assert len(metric.shape) == 1 # check that this is a vector not a higher dimensional array.
        rand_mean = rand_sum.sum()/rand_sum.shape[0]/float(len(rand_list))
        metric_mean = metric.mean() # e.g. mean of nodal clustering coefficient gives "clustering coefficient" as defined by Watts and Strogatz.
        metric_mean_normalized = metric_mean/rand_mean
        return metric, metric_mean, metric_mean_normalized, rand_mean # vector, scalar, scalar, scalar.

def computeAuc(auc_dict, out_dir, thresholds, min_thresh=None, max_thresh=None, weighted=False, atlas_string=None):
    # By default, compute AUC over all thresholds in auc_dict
    if min_thresh is not None:
        thresholds = [t for t in thresholds if t >= min_thresh]
    if max_thresh is not None:
        thresholds = [t for t in thresholds if t <= max_thresh]

    # Initialize a dataframe to store the AUC results for all subjects, to be analyzed later.
    auc_df = pd.DataFrame(dtype=np.float64, index=auc_dict.keys())

    bw_string = 'w' if weighted else 'b'

    # Add metrics for each subject at each threshold, and the AUC value.
    for subject, metric_dict in auc_dict.items():
        for metric_name, thresh_to_val_dict in metric_dict.items():
            metric_auc = None
            val_old = None

            for threshold in thresholds:
                val_new = thresh_to_val_dict[threshold]
                thresh_new = threshold
                
                if metric_auc is None: # if first threshold, initialize AUC array and continue to the next threshold.
                    if np.isscalar(val_new): # if global metric
                        metric_auc = 0
                        metric_mean = val_new # for testing only.
                    else: # if nodal metric
                        metric_auc = np.zeros(shape=val_new.shape, dtype=np.float64)
                        metric_mean = np.zeros(shape=val_new.shape, dtype=np.float64)
                else: # add trapezoid area to AUC.
                    trapezoid_area = (val_old + val_new) * (thresh_new - thresh_old) / 2.0
                    metric_auc += trapezoid_area
                    metric_mean += val_new

                # Add value to dataframe.
                if np.isscalar(val_new):
                    col_name = bw_string+'_'+metric_name+'_thresh_'+"{:.2f}".format(thresh_new)
                    auc_df.loc[subject, col_name] = val_new
                else:
                    for node_number, node_value in enumerate(val_new, 1):
                        col_name = bw_string+'_'+metric_name+'_'+atlas_string+'_'+str(int(node_number))+'_thresh_'+'{:.2f}'.format(thresh_new)
                        auc_df.loc[subject, col_name] = node_value
                val_old = val_new
                thresh_old = thresh_new

            metric_int_mean = metric_auc/(thresholds[-1] - thresholds[0]) # this is what I think should be used in literature.
            metric_mean /= float(len(thresholds)) # mean of values taken at regular density intervals should be very close to the mean based on the numerical integration using the trapezoidal rule.
        
            # Add AUC to dataframe:
            # Add value to dataframe.
            if np.isscalar(metric_auc):
                col_name = bw_string+'_'+metric_name+'_auc'
                auc_df.loc[subject, col_name] = metric_auc
                col_name = bw_string+'_'+metric_name+'_intmean'
                auc_df.loc[subject, col_name] = metric_int_mean
            else:
                for node_number, node_value in enumerate(metric_auc, 1):
                    col_name = bw_string+'_'+metric_name+'_'+atlas_string+'_'+str(int(node_number))+'_auc'
                    auc_df.loc[subject, col_name] = node_value
                for node_number, node_value in enumerate(metric_int_mean, 1):
                    col_name = bw_string+'_'+metric_name+'_'+atlas_string+'_'+str(int(node_number))+'_intmean'
                    auc_df.loc[subject, col_name] = node_value

    # Save the dataframe with all of the AUC results.
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    out_name = bw_string+'_auc_'+'{:.4f}-{:.4f}'.format(thresholds[0], thresholds[-1]) + '.csv'
    out_path = os.path.join(out_dir, out_name)
    auc_df.to_csv(out_path)
    

def fisher(r):
    z = np.arctanh(r)
    return z

def graphTheoryMetrics(matrix_paths, legend_path=None, gm_only=False, out_dir=None, atlas_string=None, density_min=None, density_max=None, density_step=0.01, num_random_networks=100, seed=None, verbose=False, check_rand2=False, transform="none", weight_prefix='conn', add_transpose=False, compute_nth_percentile_of_network_density=False):
    # Set seed for NumPy random number generation. Use same seed in rep
    np.random.seed(seed)

    # Initialize a dataframe to store the raw correlation values between each pair of regions that are included, before thresholding.
    conn_df = pd.DataFrame(dtype=np.float64) # store the correlation values after applying transforms
    conn_original_df = pd.DataFrame(dtype=np.float64) # store the correlation values before applying transforms
    
    # Read the correlation matrix for each subject.
    matrix_dict = OrderedDict()
    legend_parsed = False # set True once the ROI legend has been parsed (or inferred from a connectivity matrix).
    for matrix_path in matrix_paths:
        print(f'Processing {os.path.abspath(matrix_path)}')
        subject = os.path.basename(matrix_path).split('.')[0] # Assume file is named <subject ID>.<extension>
        
        matrix = np.loadtxt(matrix_path, dtype=np.float64)

        if not legend_parsed:
            # Read legend.
            if not legend_path is None:
                legend_df = pd.read_csv(legend_path)
            else:                
                legend_df = pd.DataFrame(OrderedDict([('number', range(1, matrix.shape[0]+1))]))

            # Remove rows from legend if they do not have 'type'=='gm'.
            if gm_only:
                legend_df = legend_df.loc[legend_df['type']=='gm', :]
            keep_indices = [ii-1 for ii in legend_df['number']]
                
            legend_df.reset_index(drop=True, inplace=True) # reset index to [0, 1, ..., 63] but keep the 'number' column.
            
            legend_parsed = True

        # If the input matrix is nonsymmetric and needs to be made so.
        if add_transpose:
            matrix += matrix.T
        
        # Remove unwanted rows and columns.
        matrix = matrix[keep_indices, :][:, keep_indices]

        # Set diagonal entries to zero.
        for ii in range(matrix.shape[0]):
            matrix[ii, ii] = 0

        # Store a copy of the correlations matrix before applying any transform (e.g. absolute value).
        matrix_original = matrix.copy()

        # ("none", "abs", "pos", "z", "z_abs", "z_pos")
        if transform == "none":
            print("Not applying transform to input edge strengths.")
            pass
        elif transform == "abs":
            print("Taking absolute values of input edge strengths.")
            matrix = np.abs(matrix)
        elif transform == "pos":
            print("Zeroing negative input edge strengths")
            matrix[matrix<0] = 0
        elif transform == "z":
            print("Fisher z transforming input edge strengths")
            matrix = fisher(matrix)
        elif transform == "z_abs":
            print("Fisher z transforming absolute values of input edge strengths")
            matrix = np.abs(matrix)
            matrix = fisher(matrix)
        elif transform == "z_pos":
            print("Fisher z transforming positive input edge strengths")
            matrix[matrix<0] = 0
            matrix = fisher(matrix)
        else:
            raise Exception("Invalid transform")
        
        matrix_dict[subject] = matrix

        # Save the correlation values to a dataframe before and after applying a transform.
        for ii in range(matrix.shape[0]-1):
            for jj in range(ii+1, matrix.shape[0]):
                conn = matrix[ii, jj]
                conn_original = matrix_original[ii, jj]
                
                ii_label_num = legend_df.loc[ii, 'number']
                jj_label_num = legend_df.loc[jj, 'number']
                column_name = f'{weight_prefix}_{atlas_string}_{str(ii_label_num)}_{str(jj_label_num)}'
                
                conn_df.loc[subject, column_name] = conn
                conn_original_df.loc[subject, column_name] = conn_original
        
    # Create output directory if it does not exist.
    if (out_dir is None) or (not os.path.isdir(out_dir)):
        os.makedirs(out_dir)

    # Save the dataframe with the correlation values before and after transformation.
    out_name = f'{weight_prefix}_{atlas_string}_n{str(matrix.shape[0])}.csv'
    out_path = os.path.join(out_dir, out_name)
    conn_df.to_csv(out_path)

    if (not transform is None):
        out_name = f'raw_{weight_prefix}_{atlas_string}_n{str(matrix.shape[0])}.csv'
        out_path = os.path.join(out_dir, out_name)
        conn_original_df.to_csv(out_path)
        
    ## Compute graph theory metrics.
    b_auc_dict = OrderedDict() # stores binary network metrics for every subject at each threshold
    w_auc_dict = OrderedDict() # stores weithed network metrics for every subject at each threshold
    metrics_dict_base = OrderedDict([('betweenness',OrderedDict()),
                                     ('charpath',OrderedDict()),
                                     ('charpath_rand_mean',OrderedDict()),
                                     ('lambda',OrderedDict()), # normalized characteristic path length
                                     ('clustering_coefficient',OrderedDict()),
                                     ('global_clustering_coefficient',OrderedDict()),
                                     ('gamma',OrderedDict()), # normalized clustering coefficient
                                     ('global_clustering_coefficient_rand_mean',OrderedDict()),
                                     ('degree',OrderedDict()),
                                     ('sigma',OrderedDict()), # small-worldness
                                     ('global_efficiency',OrderedDict()),
                                     ('local_efficiency',OrderedDict()),
                                     ('network_local_efficiency',OrderedDict()),
                                     ('transitivity',OrderedDict()),
                                     ('modularity',OrderedDict()) ])

    ## Define a list of density values to compute AUC over.
    # If density_min not set, set such that average binary degree > 2*log(N)
    #if density_min is None:
    #    N = matrix_dict[list(matrix_dict.keys())[0]].shape[0]
    #    density_min = np.round(2.0*np.log(N)/(N-1.0), decimals=2)

    # Calculate distribution of densities.
    if not compute_nth_percentile_of_network_density is None:
        densities = []
        for subject, matrix in matrix_dict.items():
            density = (np.count_nonzero(matrix)/2) / (matrix.shape[0]*(matrix.shape[0]-1)/2) # number of connections / # possible connections
            densities.append(density)

        densities.sort()
        T_nth_percentile = np.percentile(densities, compute_nth_percentile_of_network_density)
        print(f'{compute_nth_percentile_of_network_density}th percentile of densities = {T_nth_percentile}')
        sys.exit()
        
    # If density_max not set, set to a very high value and stop loop when sigma<1.1 for any subject.
    if density_max is None:
        temp_density_max = 0.80 + density_step
    else:
        temp_density_max = density_max + density_step

    if density_min >= temp_density_max:
        msg = f'Minimum density ({density_min}) must be less than maximum density ({density_max}).'
        raise Exception(msg)
    thresholds = list(np.arange(density_min, temp_density_max, density_step))
    print('Computing metrics over thresholds:', list(thresholds))

    
    for threshold in thresholds:
        print('Threshold: ', threshold)
        for subject, matrix in matrix_dict.items():
            if verbose:
                print('\nSubject: '+subject)


            # If subject not in the metric dictionary, add them.
            if subject not in b_auc_dict:
                b_auc_dict[subject] = deepcopy(metrics_dict_base)
            if subject not in w_auc_dict:
                w_auc_dict[subject] = deepcopy(metrics_dict_base)

            # Get thresholded matrices (weighted and binary versions).
            b_adj, w_adj = thresholdMatrix(matrix, threshold)

            ## Generate random networks with the same number of nodes, degree distribution, and (for weighted networks) edge weight distribution.
            if verbose:
                print('Generating random networks.')
            b_rand_list, w_rand_list = generateRandomNetworks(b_adj, w_adj, num_random_networks=num_random_networks)
            if check_rand2:
                _, w_rand2_list = generateRandomNetworks(b_adj, w_adj, num_random_networks=num_random_networks, custom_weighted_algorithm=False)

            #### Get network topology metrics.
            ## Betweenness centrality
            if verbose:
                print('Calculating betweenness (binary).')
            betweenness_binary = bct.centrality.betweenness_bin(b_adj)
            b_auc_dict[subject]['betweenness'][threshold] = betweenness_binary
            
            if verbose:
                print('Calculating betweenness (weighted).')
            betweenness_weighted = bct.centrality.betweenness_wei(w_adj)
            w_auc_dict[subject]['betweenness'][threshold] = betweenness_weighted

            ## Characteristic path length
            if verbose:
                print('Calculating characteristic path length (binary).')
            #charpath_binary, lambda_binary, charpath_rand_mean_binary = normalizer(b_adj, b_rand_list, myCharPathBin) # Set disconnected nodes distance to zero.
            charpath_binary, lambda_binary, charpath_rand_mean_binary = normalizer(b_adj, b_rand_list, myCharPathBin, function_kwargs={'custom_algorithm':True}) # Use harmonic mean mean approach

            b_auc_dict[subject]['charpath'][threshold] = charpath_binary
            b_auc_dict[subject]['lambda'][threshold] = lambda_binary
            b_auc_dict[subject]['charpath_rand_mean'][threshold] = charpath_rand_mean_binary
            if verbose:
                print('Path length for real network (binary): '+str(charpath_binary))
            if verbose:
                print('Lambda (binary): '+str(lambda_binary))

            if verbose:
                print('Calculating characteristic path length (weighted).')
            #charpath_weighted, lambda_weighted, charpath_rand_mean_weighted = normalizer(w_adj, w_rand_list, myCharPathWei) # Set disconnected nodes distance to zero.
            charpath_weighted, lambda_weighted, charpath_rand_mean_weighted = normalizer(w_adj, w_rand_list, myCharPathWei, function_kwargs={'custom_algorithm':True}) # Use harmonic mean approach
            if check_rand2:
                _, lambda_weighted_rand2, _ = normalizer(w_adj, w_rand2_list, myCharPathWei, function_kwargs={'custom_algorithm':True}) # Harmonic mean

            w_auc_dict[subject]['charpath'][threshold] = charpath_weighted
            w_auc_dict[subject]['lambda'][threshold] = lambda_binary
            w_auc_dict[subject]['charpath_rand_mean'][threshold] = charpath_rand_mean_weighted
            if verbose:
                print('Path length for real network (weighted): '+str(charpath_weighted))
            if verbose:
                print('Lambda (weighted):'+str(lambda_weighted))
            if check_rand2:
                if verbose:
                    print('Lambda (weighted) (alternate randomization method):'+str(lambda_weighted_rand2))
            
            ## Clustering coefficient
            if verbose:
                print('Calculating clustering coefficient (binary).')
            clustering_coefficient_binary, global_clustering_coefficient_binary, gamma_binary, global_clustering_coefficient_rand_mean_binary = normalizer(b_adj, b_rand_list, bct.clustering.clustering_coef_bu, scalar=False)            
            b_auc_dict[subject]['clustering_coefficient'][threshold] = clustering_coefficient_binary
            b_auc_dict[subject]['global_clustering_coefficient'][threshold] = global_clustering_coefficient_binary
            b_auc_dict[subject]['gamma'][threshold] = gamma_binary
            b_auc_dict[subject]['global_clustering_coefficient_rand_mean'][threshold] = global_clustering_coefficient_rand_mean_binary
            if verbose:
                print('Gamma (binary):', gamma_binary)
            if verbose:
                print('Calculating clustering coefficient (weighted).')
            clustering_coefficient_weighted, global_clustering_coefficient_weighted, gamma_weighted, global_clustering_coefficient_rand_mean_weighted = normalizer(w_adj, w_rand_list, bct.clustering.clustering_coef_wu, scalar=False)
            if check_rand2:
                _, _, gamma_weighted_rand2, _ = normalizer(w_adj, w_rand2_list, bct.clustering.clustering_coef_wu, scalar=False)
            w_auc_dict[subject]['clustering_coefficient'][threshold] = clustering_coefficient_weighted
            w_auc_dict[subject]['global_clustering_coefficient'][threshold] = global_clustering_coefficient_weighted
            w_auc_dict[subject]['gamma'][threshold] = gamma_weighted
            w_auc_dict[subject]['global_clustering_coefficient_rand_mean'][threshold] = global_clustering_coefficient_rand_mean_weighted
            if verbose: print('Gamma (weighted):', gamma_weighted)
            if check_rand2:
                if verbose: print('Gamma (weighted) (rand2):', gamma_weighted_rand2)

            ## Small-worldness
            sigma_binary = gamma_binary/lambda_binary
            b_auc_dict[subject]['sigma'][threshold] = sigma_binary
            if verbose: print('Sigma (binary):', sigma_binary)

            sigma_weighted = gamma_weighted/lambda_weighted
            w_auc_dict[subject]['sigma'][threshold] = sigma_weighted
            if verbose: print('Sigma (weighted):', sigma_weighted)

            ## Degree
            if verbose: print('Calculating degree (binary).')
            degree_binary = bct.degree.degrees_und(b_adj)
            b_auc_dict[subject]['degree'][threshold] = degree_binary

            if verbose: print('Calculating degree (weighted).')
            degree_weighted = bct.degree.strengths_und(w_adj)
            w_auc_dict[subject]['degree'][threshold] = degree_weighted
            
            ## Global efficiency
            if verbose: print('Calculating global efficiency (binary).')
            global_efficiency_binary = bct.distance.efficiency_bin(b_adj, local=False)
            b_auc_dict[subject]['global_efficiency'][threshold] = global_efficiency_binary

            if verbose: print('Calculating global efficiency (weighted).')
            global_efficiency_weighted = bct.distance.efficiency_wei(w_adj, local=False)
            w_auc_dict[subject]['global_efficiency'][threshold] = global_efficiency_weighted
            
            ## Local efficiency
            if verbose: print('Calculating local efficiency (binary).')
            local_efficiency_binary = bct.distance.efficiency_bin(b_adj, local=True)
            b_auc_dict[subject]['local_efficiency'][threshold] = local_efficiency_binary
            
            if verbose: print('Calculating local efficiency (weighted).')
            local_efficiency_weighted = bct.distance.efficiency_wei(w_adj, local=True)
            w_auc_dict[subject]['local_efficiency'][threshold] = local_efficiency_weighted

            ## Network local efficiency
            b_auc_dict[subject]['network_local_efficiency'][threshold] = local_efficiency_binary.mean()

            w_auc_dict[subject]['network_local_efficiency'][threshold] = local_efficiency_weighted.mean()
            
            ## Rich club
            #R_binary, Nk_binary, Ek_binary = bct.core.rich_club_bu(b_adj)

            #Rw_weighted = bct.core.rich_club_wu(w_adj)
            
            ## Transitivity
            if verbose: print('Calculating transitivity (binary).')
            transitivity_binary = bct.clustering.transitivity_bu(b_adj)
            b_auc_dict[subject]['transitivity'][threshold] = transitivity_binary

            if verbose: print('Calculating transitivity (weighted).')
            transitivity_weighted = bct.clustering.transitivity_wu(w_adj)
            w_auc_dict[subject]['transitivity'][threshold] = transitivity_weighted

            ## Modularity # NOTE: Algorithm in BCT appears broke. For now, save NaN if it returns error.
            if verbose: print('Calculating modularity (binary).')
            try:
                ci_binary, q_binary = bct.modularity.modularity_und(b_adj)
            except:
                ci_binary, q_binary = (np.nan, np.nan)
            b_auc_dict[subject]['modularity'][threshold] = q_binary

            if verbose: print('Calculating modularity (weighted).')
            try:
                ci_weighted, q_weighted = bct.modularity.modularity_und(w_adj)
            except:
                ci_weighted, q_weighted = (np.nan, np.nan)
            w_auc_dict[subject]['modularity'][threshold] = q_weighted

        #### Check whether this density threshold should be used as an endpoint in the AUC calculation.
        ## Want average (binary) degree > 2*ln(N) (sets lower bound). First in 2007-Bullmore?
        # This is solely dictacted by the density; need density
        ## Want sigma (binary) > 1.1 for all patients (sets upper bound?). This never seems to occur, and I can't find an argment for the value 1.1.
#        sigma_min = np.inf
#        sigma_mean = 0
#        for subject, metric_dict in w_auc_dict.items():
#            sigma = metric_dict['sigma'][threshold]
#            gamma = metric_dict['gamma'][threshold]
#            lambda_b = metric_dict['lambda'][threshold]
#            if sigma < sigma_min:
#                sigma_min = sigma
#            sigma_mean += sigma
#            print 'Small worldness (binary) for subject '+subject+': '+str(sigma)
#            print 'gamma:', gamma
#            print 'lambda:', lambda_b
#        sigma_mean /= float(len(b_auc_dict))
#        print 'Minimum small worldness (binary):', sigma_min
#        print 'Mean small worldness (binary):', sigma_mean

    ## Save the auc_dict Python objects in case I want to change the threshold range without recomputing everything
#    if not os.path.exists(out_dir):
#        os.makedirs(out_dir)

#    b_auc_args = {'auc_dict':b_auc_dict, 'out_dir':out_dir, 'thresholds':thresholds, 'weighted':False, 'atlas_string':atlas_string}
#    b_pickle_name = 'b_auc_args_'+datetime.datetime.now().strftime("%Y%m%d_%H%M%S")+'.pkl'
#    b_pickle_path = os.path.join(out_dir, b_pickle_name)
#    with open(b_pickle_path, 'w') as f:
#        pickle.dump(b_auc_args, f)

#    w_auc_args = {'auc_dict':w_auc_dict, 'out_dir':out_dir, 'thresholds':thresholds, 'weighted':True, 'atlas_string':atlas_string}
#    w_pickle_name = 'w_auc_args_'+datetime.datetime.now().strftime("%Y%m%d_%H%M%S")+'.pkl'
#    w_pickle_path = os.path.join(out_dir, b_pickle_name)
#    with open(w_pickle_path, 'w') as f:
#        pickle.dump(w_auc_args, f)

#    computeAuc(**b_auc_args)
#    computeAuc(**w_auc_args)

    ## Compute the AUC for each metric, and save all metrics (density-specific and AUC) to CSV.
    computeAuc(b_auc_dict, out_dir, thresholds, weighted=False, atlas_string=atlas_string)
    computeAuc(w_auc_dict, out_dir, thresholds, weighted=True, atlas_string=atlas_string)
    
    return
    
if (__name__ == '__main__'):
    # Create argument parser.
    description = """"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    # Define positional arguments.
    parser.add_argument('matrix_paths', help='Paths to correlation matrices to select regions from.', nargs='+')    
    
    # Define optional arguments.
    parser.add_argument('-l', '--legend_path', help='path to legend CSV file (not required).')
    parser.add_argument('-a', '--atlas_string', help='string denoting which atlas was used. Useful if you have connectivity matrices derived from multiple node parcellations, and want to automatically name files and column headers based on the name of the parcellation.', default='roi')
    parser.add_argument('-g', '--gm_only', help="Keep only grey-matter regions. To use this option, you must input a legend CSV containing a column named 'type'; only rows that have type='gm' will be kept. This may be useful if you are using a parcellation which includes non-grey-matter regions, which you want to exclude.", action='store_true')
    parser.add_argument('-o', '--out_dir', help='directory to save matrices to', type=str, default='.')
    parser.add_argument('-n', '--num_random_networks', help='number of random networks to compare real network to when obtaining normalized metrics (e.g. normalized characteristic path and clustering coefficient). Default = 100', type=int, default=100)
    parser.add_argument('-v', '--verbose', help='print lots of information.', action='store_true')
    parser.add_argument('-s', '--seed', type=int, help='seed (int) for random number generation, used in construction of randomized networks. Set to the same value between runs to obtain the same set of randomized ')
    parser.add_argument('-r', '--check_rand2', action='store_true', help='REMOVE THIS TESTING OPTION')
    parser.add_argument('--density_min', type=float, help='Minimum network density for AUC calculation.', default=0.01)
    parser.add_argument('--density_max', type=float, help='Maximum network density for AUC calculation. Default: 0.80', default=0.80)
    parser.add_argument('--density_step', type=float, help='Network density step size for AUC calculation. Default: 0.01', default=0.01)
    parser.add_argument('--transform', type=str, help='Transformation applied to input correlation values. Must be one of ("none", "abs", "pos", "z", "z_abs", "z_pos"). none: no transformation, abs: take absolute value of correlations, pos: set negative corrrelations to zero, z: Fisher z-transform correlations, z_abs: Fisher z-transform absolute values of the correlations, z_pos: set negative correlations to zero, then Fisher z-transform ', choices=["none", "abs", "pos", "z", "z_abs", "z_pos"], default="none")
    parser.add_argument('--weight_prefix', type=str, default='conn', help='network edge weights will be saved to CSV with columns named <weight_prefix>_<atlas_string>_<ROI 1 number>_<ROI 2 number>')
    parser.add_argument('--add_transpose', action='store_true', help='before computing metrics, add the transpose of the adjancy matrix. FSL ProbtrackX generates nonsymmetric matrices, in which the number of streamlines seeded in ROI 1 and terminating in ROI 2 will be different from the number of steamlines seeded in ROI 2 and terminating in ROI 1.')
    parser.add_argument('--compute_nth_percentile_of_network_density', type=float, help='Compute the nth percentile of network densities across all input matrices, and quit without computing metrics. Use this mode to help choose a maximum density threshold. Pass the desired percentile (e.g. 95) in this argument.', default=None)

    # Print help if no args input.
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    # Parse arguments.
    args = parser.parse_args()

    # Run main function.
    graphTheoryMetrics(**vars(args))
