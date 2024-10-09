import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import pickle as pk
import sys
from pathlib import Path
from os import path

from sklearn.decomposition import PCA, IncrementalPCA
import openTSNE 
from sklearn.manifold import TSNE
import umap
import phate


def load_matrix(path, chunksize):
    return pd.concat([x for x in pd.read_csv(path, chunksize=chunksize, names=[i for i in range(101,29779)], sep='\t')])

def load_chunks(path, chunksize):
    return pd.read_csv(path, chunksize=chunksize, names=[i for i in range(101,29779)], sep='\t')

def run_incremental_pca(matrix, chunks, batch_size, n_components):
    ipca = IncrementalPCA(n_components=n_components, batch_size=batch_size)
    for chunk in chunks:
        ipca.partial_fit(chunk)

    return ipca

def transform_ipca(ipca, matrix, n_components):
    PC_embeddings = ipca.transform(matrix)
    PCs_df = pd.DataFrame(data=PC_embeddings,
                      index=matrix.index,
                      columns=[f'PC{x}' for x in range(1,n_components+1)])

    return PC_embeddings, PCs_df


def fit_transform_tSNE(matrix):
    tSNE_object = TSNE(perplexity=150,
                random_state=100)
    
    tSNE_embeddings = tSNE_object.fit_transform(matrix)
    
    tSNE_df = pd.DataFrame(data=tSNE_embeddings,
                           index=matrix.index,
                           columns=['tSNE1', 'tSNE2'])

    return tSNE_object, tSNE_embeddings, tSNE_df

def fit_transform_umap(matrix):
    umap_object = umap.UMAP()    
    umap_embeddings = umap_object.fit_transform(matrix)

    umap_df = pd.DataFrame(data=umap_embeddings,
                           index=matrix.index,
                           columns=['umap1', 'umap2'])
    return umap_object, umap_embeddings, umap_df


def fit_transform_PHATE(matrix):
    '''This function performs PHATE embedding on a given matrix, 
    returning the PHATE object, the embeddings, and a DataFrame containing 
    the embeddings labeled as 'PHATE1' and 'PHATE2'.'''
    phate_object = phate.PHATE(knn=50)
    phate_embeddings = phate_object.fit_transform(matrix)
    phate_df = pd.DataFrame(data=phate_embeddings,
                            index=matrix.index,
                            columns=['PHATE1', 'PHATE2'])
    return phate_object, phate_embeddings, phate_df


def main():
    input_matrix = Path(sys.argv[1])
    chunksize = int(sys.argv[2])
    n_components = 20
    output_path = input_matrix.parents[0]
    
    is_transposed = False

    experiment_condition =  f"{input_matrix.stem}"
        

   #load data
    print('Loading Data...')
    working_matrix = load_matrix(input_matrix, chunksize)
    working_chunks = load_chunks(input_matrix, chunksize)
    print(f'Matrix shape is: {working_matrix.shape}')
    
    # Transpose if is_transposed is True
    if is_transposed:
        working_matrix = working_matrix.T
        working_matrix = working_matrix.loc[(working_matrix != 0).any(axis=1)]
        num_chunks = len(working_matrix) // chunksize
        working_chunks = [working_matrix.iloc[i*chunksize:(i+1)*chunksize, :] for i in range(num_chunks)]

        print(f'Transposed Matrix shape is: {working_matrix.shape}')
    
    if path.exists(f'{output_path}/pca_{experiment_condition}_ipca.pkl'):
        print('Loading ipca object...')
        ipca = pk.load(open(f'{output_path}/pca_{experiment_condition}_ipca.pkl', 'rb'))
    else:
        #Running Incremental PCA
        print('Running Incremental PCA...')
        ipca = run_incremental_pca(working_matrix, working_chunks, chunksize, n_components)
        pk.dump(ipca, open(f'{output_path}/pca_{experiment_condition}_ipca.pkl', 'wb'))

    if path.exists(f'{output_path}/pca_{experiment_condition}_embeddings.pkl') and \
    path.exists(f'{output_path}/pca_{experiment_condition}_embeddings.tsv'):
        print('Loading PC embeddings')
        PC_embeddings = pk.load(open(f'{output_path}/pca_{experiment_condition}_embeddings.pkl', 'rb'))
        PC_df = pd.read_csv(f'{output_path}/pca_{experiment_condition}_embeddings.tsv', sep='\t', index_col=[0])
    else:
        print('Tranforming Incremental PCA...')
        PC_embeddings, PC_df = transform_ipca(ipca, working_matrix, n_components)
        pk.dump(PC_embeddings, open(f'{output_path}/pca_{experiment_condition}_embeddings.pkl', 'wb'))
        PC_df.to_csv(f'{output_path}/pca_{experiment_condition}_embeddings.tsv', '\t')

    if path.exists(f'{output_path}/tSNE_{experiment_condition}_object.pkl'):
        pass
    else:
        print('Running tSNE on PCs...')
        tSNE_object, tSNE_embeddings, tSNE_df = fit_transform_tSNE(PC_df)
        pk.dump(tSNE_object, open(f'{output_path}/tSNE_{experiment_condition}_object.pkl', 'wb'))
        pk.dump(tSNE_embeddings, open(f'{output_path}/tSNE_{experiment_condition}_embeddings.pkl', 'wb'))
        tSNE_df.to_csv(f'{output_path}/tSNE_{experiment_condition}_embeddings.tsv', '\t')

    if path.exists(f'{output_path}/umap_{experiment_condition}_object.pkl'):
        pass
    else:
        print('Running umap on PCs...')
        umap_object, umap_embeddings, umap_df = fit_transform_umap(PC_df)
        pk.dump(umap_object, open(f'{output_path}/umap_{experiment_condition}_object.pkl', 'wb'))
        pk.dump(umap_embeddings, open(f'{output_path}/umap_{experiment_condition}_embeddings.pkl', 'wb'))
        umap_df.to_csv(f'{output_path}/umap_{experiment_condition}_embeddings.tsv', '\t')

    if path.exists(f'{output_path}/phate_{experiment_condition}_object.pkl'):
        pass
    else:
        print('Running PHATE on PCs...')
        phate_object, phate_embeddings, phate_df = fit_transform_PHATE(PC_df)
        pk.dump(phate_object, open(f'{output_path}/phate_{experiment_condition}_object.pkl', 'wb'))
        pk.dump(phate_embeddings, open(f'{output_path}/phate_{experiment_condition}_embeddings.pkl', 'wb'))
        phate_df.to_csv(f'{output_path}/phate_{experiment_condition}_embeddings.tsv', '\t')

def test():
    from pprint import pprint
    test_path="/lustre06/project/6065672/shared/covid-19/database/data/intra/2020_2021"
    test_matrix= f'{test_path}/test.mat'
    pprint('loading matrix...')
    matrix = load_matrix(test_matrix, 1000)
    pprint(f'matrix shape: {matrix.shape}')
    pprint(matrix.head)
    pprint(f'laoding matrix into chuncks')
    chuncks = load_chunks(test_matrix, 1000)
    pprint('matrix chunks loaded')

if __name__ == '__main__':
    #test()
    main()
