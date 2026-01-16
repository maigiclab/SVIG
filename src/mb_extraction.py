import musical
import numpy as np
import pickle
import logging
import pandas as pd
import argparse
import sys
import os
from pathlib import Path
sys.path.append(os.path.dirname('utils/'))
from sorted_array import SORTED_CAT_ARRAY
from sorted_array import SORTED_CAT_ARRAY_RT

interactive = sys.stdin.isatty()
parser = argparse.ArgumentParser()
parser.add_argument("--matrix-type", default="non_clustered_rfd_to_transpose", help="Matrix type for extraction")
parser.add_argument("--no-iters", type=int, default=1000, help="Number of iterations")
parser.add_argument("--min-n-components", type=int, default=3, help="Minimum number of components")
parser.add_argument("--max-n-components", type=int, default=20, help="Maximum number of components")
parser.add_argument(
    "--normalize-X",
    action="store_true",
    default=False,  
    help="Normalize X before processing"
)
parser.add_argument("--catalogue-matrix", default='../data/processed/SVmatrices/PCAWG/RFD/Ovary_RFD.csv', help="Input matrix catalogue")
parser.add_argument("--cohort-name", default='gel_sel', help="Input cohortname")


args = parser.parse_args()



# Setup logging
logging_folder='../data/processed/sig_extraction/'
Path(logging_folder).mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=logging_folder + 'SV_SigExtraction_mvnmf_RFD_'+args.matrix_type+'.log',
                    level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Start logging
logging.info('Starting the DenovoSig model fitting process.\n')

# Load your data here
X = pd.read_csv(args.catalogue_matrix, 
                sep=',', index_col=0)

if args.matrix_type=='non_clustered_rfd' or args.matrix_type=='non_clustered_normalized':
    colnames_sel = [c for c in SORTED_CAT_ARRAY if 'non-clustered' in c]
    X = X.loc[colnames_sel,]
elif args.matrix_type=='non_clustered_rfd_to_transpose':
    X=X.T
    colnames_sel = [c for c in SORTED_CAT_ARRAY if 'non-clustered' in c]
    X = X.loc[colnames_sel,]
elif args.matrix_type=='non_clustered_replitime' :
    colnames_sel = [c for c in SORTED_CAT_ARRAY_RT if 'non-clustered' in c]
    X = X.loc[[c for c in colnames_sel if c in X.index],]
elif args.matrix_type=='non_clustered_replitime_to_transpose' :
    X = X.T
    colnames_sel = [c for c in SORTED_CAT_ARRAY_RT if 'non-clustered' in c]
    X = X.loc[[c for c in colnames_sel if c in X.index],]
elif args.matrix_type=='channels32' or args.matrix_type=='channels32_normalized':
    n = 4  # group size
    new_rows = []
    X = X.loc[SORTED_CAT_ARRAY,]
    for i in range(0, X.shape[0], n):
        new_row = X.iloc[i:i+n,:].sum(axis=0)
        new_rows.append(new_row)
    # Combine the new columns into a new dataframe
    counts_summed_32 = pd.concat(new_rows, axis=1).T
    counts_summed_32.index = [s.removesuffix("_L_L") for s in X.index[range(0, X.shape[0], n)]]  
    X=counts_summed_32
elif args.matrix_type=='non_clustered_corr_normalized':
    colnames_sel = [c for c in SORTED_CAT_ARRAY if 'non-clustered' in c]
    corr_df = pd.read_csv('/home/dg204/park_dglodzik/svig/non_clustered_shuffle_simple/rfd_counts_non_clustered_shuffle_simple.csv', index_col=0)
    corr_df['cat32'] = ['_'.join(c.rsplit('_', 2)[:-2]) for c in corr_df.index]
    corr_df['count_per_cat32'] = corr_df.groupby('cat32')['count'].transform('sum')
    corr_df['corr'] = corr_df['count'] / corr_df['count_per_cat32']
    corr_df.sort_index(inplace=True)
    X_nonclust = X.loc[colnames_sel,]
    X_nonclust_corr = X_nonclust.div(corr_df.loc[X_nonclust.index, 'corr'].values, axis=0)
    X = X_nonclust_corr


X = X.apply(pd.to_numeric, errors="coerce")
X=X.loc[:,X.sum(axis=0)>30]
X.to_csv(logging_folder + '/X_'+args.matrix_type+'_'+args.cohort_name+'.csv')


# Initialize the DenovoSig model
model = musical.DenovoSig(
    X=X, 
    min_n_components=args.min_n_components,  # Minimum number of signatures to test
    max_n_components=args.max_n_components,  # Maximum number of signatures to test
    init='random',  # Initialization method
    method='nmf',  # mvnmf or nmf
    n_replicates=20,  # Number of mvnmf/nmf replicates to run per n_components
    ncpu=25,  # Number of CPUs to use
    max_iter=args.no_iters,  # Maximum number of iterations for each mvnmf/nmf run
    bootstrap=True,  # Whether or not to bootstrap X for each run
    tol=1e-8,  # Tolerance for claiming convergence of mvnmf/nmf
    verbose=1,  # Verbosity of output
    normalize_X=args.normalize_X  # Whether or not to L1 normalize each sample in X before mvnmf/nmf
)

logging.info('DenovoSig model initialized.\n')

# Fit the model
model.fit()

with open(logging_folder+'SV_SigExtraction_nmf_RFD_'+args.matrix_type+'_'+args.cohort_name+'.pkl', 'wb') as f:
    pickle.dump(model, f, pickle.HIGHEST_PROTOCOL)

logging.info('Model fitting complete.\n')