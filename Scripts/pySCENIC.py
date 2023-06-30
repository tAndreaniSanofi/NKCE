import warnings
import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import sys
import os
import glob
import pickle
from scipy.sparse import csr_matrix
import pyranges as pr
import requests
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3
from dask.diagnostics import ProgressBar
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from multiprocessing import Pool
from typing import Mapping, Optional
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from sklearn import mixture
from tqdm import tqdm
from pyscenic.binarization import binarize, derive_threshold
from pyscenic .rss import regulon_specificity_scores
from pyscenic.diptest import diptst
from pyscenic.export import export2loom 

warnings.simplefilter(action='ignore', category=FutureWarning)
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = 'path_to_your_workdir'

DATA_FOLDER="data"
OUTPUT_FOLDER = "out"
RESOURCES_FOLDER="resources"
DATABASE_FOLDER = "db"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg19-*.mc9nr.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(os.getcwd(), RESOURCES_FOLDER, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(os.getcwd(), RESOURCES_FOLDER, 'TF_list.csv')
SC_EXP_FNAME = os.path.join(os.getcwd(), DATA_FOLDER, "count_matrix.csv")
CELLTYPE_FNAME = os.path.join(os.getcwd(), DATA_FOLDER, "cell_type.csv")
REGULONS_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "motifs.csv")
MODULES_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "modules.csv")
ADJ_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "adj.csv")
AUCELL_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "auc.csv")
BIN_AUCELL_FNAME = os.path.join(os.getcwd(), OUTPUT_FOLDER, "bin_auc.csv")

ex_matrix = pd.read_csv(SC_EXP_FNAME, sep = ',')
ex_matrix.shape
ex_matrix.index = ex_matrix['Unnamed: 0']
ex_matrix.drop('Unnamed: 0', axis=1, inplace=True)

cell_type = pd.read_csv(CELLTYPE_FNAME, index_col= 0)

tf_names = load_tf_names(MM_TFS_FNAME)
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

dbs = [RankingDatabase(fname=db_fnames[ii], name=name(db_fnames[ii])) for ii in range(len(db_fnames))]

adjacencies_grnboost = grnboost2(expression_data= ex_matrix.T, tf_names= tf_names, verbose=True)
adjacencies_grnboost.to_csv(ADJ_FNAME, sep = '\t', index = False, header = True) 

# adjacencies_genie3 = genie3(expression_data= ex_matrix.T, tf_names= tf_names, verbose=True)
# adjacencies_genie3.to_csv(ADJ_FNAME, sep = '\t', index = False, header = True) 

adjacencies = adjacencies_grnboost 
# adjacencies = adjacencies_genie3

modules = list(modules_from_adjacencies(adjacencies, ex_matrix.T))
modules.to_csv(MODULES_FNAME, sep = '\t', index = False, header = True) 

with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df, df.columns)
regulons = [r.rename(r.name.replace('(+)',' ('+str(len(r))+'g)')) for r in regulons] 


# Save the enriched motifs and the discovered regulons to disk.
df.columns = df.columns.droplevel(0)
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)

auc_mtx = aucell(ex_matrix.T, regulons, num_workers=16)
auc_mtx.to_csv(AUCELL_FNAME)

bin_mtx, th= binarize(auc_mtx, num_workers=16)
bin_mtx.to_csv(BIN_AUCELL_FNAME)

rss = regulon_specificity_scores(auc_mtx, celltype)
rss.to_csv("rss.csv")

export2loom(
    ex_mtx = ex_matrix.T,
    regulons = regulons,
    out_fname = "all_loom.loom",
    cell_annotations= cell_type.type.to_dict(),
    tree_structure = "",
    title = "Regulons",
    nomenclature = "Unknown",
    num_workers = 16,
    auc_mtx= auc_mtx ,
    auc_thresholds=th,
    compress = False
)

