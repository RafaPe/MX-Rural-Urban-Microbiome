import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import pandas as pd
import numpy as np


# CLR transformation function
def clr_transform(X):
    """
    Apply Centered Log-Ratio (CLR) transformation to compositional data.
    
    Parameters:
    X : array-like, shape (n_samples, n_features)
        Compositional data (must be positive)
    
    Returns:
    X_clr : array, CLR-transformed data
    """
    # Add small constant to avoid log(0) if there are zeros
    
    # Compute geometric mean for each sample (row)
    geom_mean = np.exp(np.mean(np.log(X), axis=1, keepdims=True))
    
    # CLR: log(x_i / geometric_mean)
    X_clr = np.log(X / geom_mean)
    
    return X_clr



# LEER DATOS DE ABUNDANCIA RELATIVA
data = pd.read_csv('../datasets/latinbiota_merge_metaphlan_data.csv', sep = '\t', skiprows=1, index_col = 0)
metadata = pd.read_excel('../data/metadata_LATINBIOTA_MEXICO.xlsx', sheet_name='Data')

data_sp = data[data.index.astype(str).str.contains("g__")]
data_sp = data_sp[~data_sp.index.astype(str).str.contains("t__")]
data_sp = data_sp[~data_sp.index.astype(str).str.contains("s__")]
data_sp = data_sp[~data_sp.index.astype(str).str.contains("unclass")]

new = []

for item in data_sp.index:
    new.append(item.split('|g__')[1].replace('|t__', ' ').replace('_', ' '))

data_sp.index = new

data_sp = data_sp.T

# SE OBTIENE EL VALOR MÍNIMO 
min_value = data_sp[data_sp > 0].min().min()


# AMALGAMACIÓN DE ESPECIES POCO COMUNES
cols_otros = data_sp.columns[(data_sp == 0).sum() > 180]

data_sp["otros"] = data_sp[cols_otros].sum(axis=1)

# TRATAMIENTO DE CEROS
suma_actual = data_sp.sum(axis=1)
faltante = 100 - suma_actual
data_sp["otros"] += faltante


# pseudocount = min_value
pseudocount = min_value / 2
# print(pseudocount)

# Todas las columnas (incluye "otros" si existe)
cols = data_sp.columns

# Contar cuántos ceros hay por fila en las columnas (excepto que "otros" no tenga ceros normalmente)
num_ceros = (data_sp[cols] == 0).sum(axis=1)

# Reemplazar ceros por el pseudocount
data_sp[cols] = data_sp[cols].replace(0, pseudocount)

# Exceso introducido en cada fila
exceso = num_ceros * pseudocount   # Serie fila a fila

# Cantidad que se debe restar a CADA columna de esa fila
correccion = exceso / len(cols)

# Restar corrección a cada columna para mantener suma total
for c in cols:
    data_sp[c] = data_sp[c] - correccion


X = data_sp.drop(columns=['label']).values
y = data_sp['label'].values

# 1. Apply CLR transformation
# X_clr = X
X_clr = clr_transform(X)

