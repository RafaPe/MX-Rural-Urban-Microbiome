import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score, silhouette_score, calinski_harabasz_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.mixture import GaussianMixture
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import make_scorer
from sklearn.cluster import HDBSCAN
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics.cluster import homogeneity_score, adjusted_mutual_info_score
from sklearn.metrics import silhouette_score
from scipy.stats import ks_2samp
import optuna


class CustomClustering(ClusterMixin, BaseEstimator):
    def __init__(self, algorithm='kmeans', scaler_=StandardScaler(), filter_criteria=None, 
                 min_threshold=0, max_threshold=np.inf, **kwargs):
        """
        Custom clustering class that applies filtering, scaling, and clustering.

        Parameters:
        -----------
        algorithm : str, default='kmeans'
            Clustering algorithm to use ('kmeans', 'gaussian_mixture', 'hdbscan').

        scaler_method : sklearn Scaler, default=StandardScaler()
            Scaling method to normalize the data before clustering.

        filter_criteria : str, default=None
            Feature filtering criterion ('abundance' or 'presence').

        min_threshold : int or float, default=0
            Lower bound for feature filtering.

        max_threshold : int or float, default=np.inf
            Upper bound for feature filtering.

        **kwargs : additional parameters
            Additional parameters for clustering algorithms.
        """
        self.algorithm = algorithm
        self.scaler_ = scaler_
        self.filter_criteria = filter_criteria
        self.min_threshold = min_threshold
        self.max_threshold = max_threshold
        self.kwargs = kwargs

    def fit(self, X, y=None):
        """
        Fits the clustering algorithm to the data after applying the filtering criterion.
        """
        self.data_ = X.copy()

        # Feature filtering based on min/max threshold
        if self.filter_criteria == 'abundance':
            column_sums = X.sum()
            threshold_lower = column_sums.mean() * self.min_threshold
            threshold_upper = column_sums.mean() * self.max_threshold
            valid_features = (column_sums >= threshold_lower) & (column_sums <= threshold_upper)
            self.filtered_data_ = X.loc[:, valid_features]

        elif self.filter_criteria == 'presence':
            binary_df = (X > 0).astype(int)
            presence_counts = binary_df.sum(axis=0)
            valid_features = (presence_counts >= self.min_threshold) & (presence_counts <= self.max_threshold)
            self.filtered_data_ = X.loc[:, valid_features]

        elif self.filter_criteria is None:
            self.filtered_data_ = self.data_

        else:
            raise ValueError("Invalid filter_criteria. Choose 'abundance' or 'presence'.")

        # Store selected feature names
        self.selected_features_ = self.filtered_data_.columns
        # print(self.filter_criteria)
        # print(self.min_threshold, self.max_threshold)
        # print(self.filtered_data_)

        # Scaling the data
        self.data_scaled_ = self.scaler_.fit_transform(self.filtered_data_)

        # Clustering
        if self.algorithm == 'kmeans':
            self.model_ = KMeans(**self.kwargs)
        elif self.algorithm == 'gaussian_mixture':
            self.model_ = GaussianMixture(**self.kwargs)
        elif self.algorithm == 'hdbscan':
            self.model_ = HDBSCAN(**self.kwargs)
        else:
            raise ValueError("Invalid algorithm. Choose 'kmeans', 'gaussian_mixture', or 'hdbscan'.")

        # Fit the model
        self.labels_ = self.model_.fit_predict(self.data_scaled_)

        return self

    def predict(self, X):
        """
        Predicts cluster labels for new data.
        """
        if not hasattr(self, "model_"):
            raise ValueError("The model has not been fitted yet.")
        if not hasattr(self, "scaler_"):
            raise ValueError("Scaler has not been fitted yet. Fit the model first.")
        if not hasattr(self, "selected_features_"):
            raise ValueError("No features were selected during fit. Ensure the model is fitted before predicting.")

        # Filter data to match selected features
        X_filtered = X[self.selected_features_]

        # Scale the input data using the fitted scaler
        X_scaled = self.scaler_.fit_transform(X_filtered)

        # Predict cluster labels
        return self.model_.fit_predict(X_scaled) if self.algorithm == 'hdbscan' else self.model_.predict(X_scaled)

def optimization_function(trial, data, labels):
    algorithm = trial.suggest_categorical('algorithm', ['kmeans', 'gaussian_mixture'])
    filter_criteria = trial.suggest_categorical('filter_criteria', ['presence', 'abundance'])
    filter_threshold = trial.suggest_float('filter_threshold', 0, 120, step=5) if filter_criteria == 'presence' else trial.suggest_float('filter_threshold', 0, 2.7, step=0.05)

    # Convert to integer only if needed
    if filter_criteria == 'presence':
        filter_threshold = int(filter_threshold)

    k = trial.suggest_int('k', 2, 6)

    if algorithm == 'kmeans':
        model = CustomClustering(
            algorithm=algorithm,
            filter_criteria=filter_criteria,
            filter_threshold=filter_threshold,
            n_clusters=k
        )
    else:
        model = CustomClustering(
            algorithm=algorithm,
            filter_criteria=filter_criteria,
            filter_threshold=filter_threshold,
            n_components=k,
            covariance_type='diag',
        )

   

    skf = StratifiedKFold(n_splits=5)

    scores = []

    for train_index, test_index in skf.split(data, labels):
        X_train, X_test = data.iloc[train_index], data.iloc[test_index]
        y_train = [labels[i] for i in train_index]
        y_test = [labels[i] for i in test_index]
        
        model.fit(X_train, y_train)
        
        y_pred = model.predict(X_test)

        score = adjusted_mutual_info_score(y_test, y_pred)
        scores.append(score)

    mean_score = sum(scores) / len(scores)

    return mean_score

def optimization_function2(trial, data):
    algorithm = trial.suggest_categorical('algorithm', ['kmeans', 'gaussian_mixture', 'hdbscan'])
    # filter_criteria = None
    # filter_threshold = trial.suggest_float('filter_threshold', 0, 120, step=5) if filter_criteria == 'presence' else trial.suggest_float('filter_threshold', 0, 2.7, step=0.05)

    # Convert to integer only if needed
    # if filter_criteria == 'presence':
    #     filter_threshold = int(filter_threshold)
    selected_features = [trial.suggest_categorical(name, [True, False]) for name in list(data.columns)]

        # List with names of selected features
    selected_feature_names = [name for name, selected in zip(list(data.columns), selected_features) if selected]

    if algorithm == 'hdbscan':
        k = None
    else:
        k = trial.suggest_int('k', 2, 10)

    if algorithm == 'kmeans':
        model = CustomClustering(
            algorithm=algorithm,
            n_clusters=k
        )
    elif algorithm == 'gaussian_mixture':
        model = CustomClustering(
            algorithm=algorithm,
            n_components=k,
            covariance_type='diag',
        )
    elif algorithm == 'hdbscan':
        model = CustomClustering(
            algorithm=algorithm
        )

    model.fit(data[selected_feature_names])
    try:
        #! MISSING NORMALIZATION
        scaler = StandardScaler()
        score = silhouette_score(scaler.fit_transform(data[selected_feature_names]), model.predict(data[selected_feature_names]))
        return score
    except:
        return float('-inf')

   

    # skf = StratifiedKFold(n_splits=5)

    # scores = []

    # for train_index, test_index in skf.split(data, labels):
    #     X_train, X_test = data.iloc[train_index], data.iloc[test_index]
    #     y_train = [labels[i] for i in train_index]
    #     y_test = [labels[i] for i in test_index]
        
    #     model.fit(X_train, y_train)
        
    #     y_pred = model.predict(X_test)

    #     score = adjusted_mutual_info_score(y_test, y_pred)
    #     scores.append(score)

    # mean_score = sum(scores) / len(scores)

    # return mean_score

def get_clustering_metrics(data, model, cluster_labels, true_labels, scaler_):
    """
    Compute clustering evaluation metrics for the provided data and cluster assignments.

    Args:
        data (pd.DataFrame): A DataFrame of shape (n_samples, n_features) with the input data.
        model (object): A fitted clustering model that has a `selected_features_` attribute 
                        containing the indices of the selected features.
        cluster_labels (np.ndarray): A 1D array of shape (n_samples,) with the predicted cluster labels.
        true_labels (np.ndarray): A 1D array of shape (n_samples,) with the true labels of the data.

    Returns:
        Tuple[float, float, float, float]: A tuple containing:
            - Silhouette score (float): Measures how similar each sample is to its own cluster compared to other clusters.
            - Davies-Bouldin score (float): Measures the average similarity ratio of each cluster with other clusters.
            - Calinski-Harabasz score (float): Measures the ratio of the sum of between-cluster dispersion and within-cluster dispersion.
            - Homogeneity score (float): Measures the extent to which clusters contain only members of a single class.
    """
    # scaler = StandardScaler()
    # data = data[model.selected_features_]
    data = scaler_.fit_transform(data)
    print(f'Silhoutte score: {silhouette_score(data, cluster_labels)}')
    print(f'Davies-Bouldin score: {davies_bouldin_score(data, cluster_labels)}')
    print(f'Calinksi-Harabasz score: {calinski_harabasz_score(data, cluster_labels)}')
    print(f'Homogenity score: {homogeneity_score(true_labels, cluster_labels)}')

    return silhouette_score(data, cluster_labels), davies_bouldin_score(data, cluster_labels), calinski_harabasz_score(data, cluster_labels), homogeneity_score(true_labels, cluster_labels)

def plot_clustering(data:pd.DataFrame, labels:list, clusters, scaler_):
    """
    Visualize clustering results in a 2D PCA plot with additional metadata information.

    Parameters:
    -----------
    clusters_data : pd.DataFrame
        The scaled and filtered dataset used for clustering. Rows represent samples, and columns represent features.
    
    metadata : pd.DataFrame
        Metadata associated with the samples, containing at least the following columns:
        - 'Lane': Unique identifier matching the sample indices in `clusters_data`.
        - 'Lifestyle': Lifestyle category (e.g., rural or urban) for each sample.
        - 'BMI': Body Mass Index (BMI) value for each sample, used for categorization.

    labels : list
        Cluster labels for each sample, as obtained from the clustering algorithm.
    """
    # scaler = StandardScaler()
    data_scaled = scaler_.fit_transform(data)
    pca = PCA(n_components=2)
    components = pca.fit_transform(data.to_numpy())
    components = pd.DataFrame(components)
    components['labels'] = labels
    components.columns = ['x', 'y', 'labels']
    components.index = data.index
    components['clusters'] = clusters
    # components['BMI'] = bmi
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=components, x='x', y='y', hue= 'clusters', style = 'labels', palette="deep", alpha=0.8)
    plt.legend(title='Category', loc='upper left', bbox_to_anchor=(1, 1))

def get_urban_rural_counts(cluster: pd.DataFrame):
    """Returns the count of urban and rural individuals in a cluster."""
    urban_count = (cluster['Lifestyle'] == 'Urban').sum()
    rural_count = (cluster['Lifestyle'] == 'Rural').sum()
    return urban_count, rural_count

def get_health_state(cluster: pd.DataFrame):
    """Returns the count of healthy and diseased individuals in a cluster."""
    healthy_count = (cluster['Health state'] == 'Healthy').sum()
    diseased_count = (cluster['Health state'] == 'Diseased').sum()
    return healthy_count, diseased_count

def get_gender_counts(cluster: pd.DataFrame):
    """Returns the count of male and female individuals in a cluster."""
    male_count = (cluster['Gender'] == 'Male').sum()
    female_count = (cluster['Gender'] == 'Female').sum()
    return male_count, female_count

def cluster_analysis(data: pd.DataFrame, clusters: np.array, metadata: pd.DataFrame):
    """
    Analyzes clusters based on urban/rural lifestyle, health state, gender, age, and BMI.
    Generates visualizations for comparison.
    """
    
    # Assign clusters to data
    data['cluster'] = clusters
    
    # Create metadata subsets for each cluster
    cluster_data = {
        f'Cluster {i+1}': metadata[metadata['Lane'].isin(data[data['cluster'] == i].index)]
        for i in range(2)
    }
    
    # Store aggregated counts
    lifestyle_counts = {}
    health_counts = {}
    gender_counts = {}
    age_data = {}
    bmi_data = {}
    
    for cluster_name, cluster_df in cluster_data.items():
        lifestyle_counts[cluster_name] = dict(zip(['Urban', 'Rural'], get_urban_rural_counts(cluster_df)))
        health_counts[cluster_name] = dict(zip(['Healthy', 'Unhealthy'], get_health_state(cluster_df)))
        gender_counts[cluster_name] = dict(zip(['Male', 'Female'], get_gender_counts(cluster_df)))
        age_data[cluster_name] = cluster_df['Age'].dropna().tolist()
        bmi_data[cluster_name] = cluster_df['BMI'].dropna().tolist()
    
    # Compute relative distributions
    relative_health = {
        cluster: {state: count / sum(states.values()) for state, count in states.items()}
        for cluster, states in health_counts.items()
    }
    relative_gender = {
        cluster: {gender: count / sum(genders.values()) for gender, count in genders.items()}
        for cluster, genders in gender_counts.items()
    }
    
    # Convert data to DataFrames for visualization
    df_lifestyle = pd.DataFrame(lifestyle_counts)
    df_health = pd.DataFrame(relative_health)
    df_gender = pd.DataFrame(relative_gender)
    
    # Create figure for visualizations
    fig, axes = plt.subplots(5, len(df_lifestyle.columns), figsize=(12, 9))
    
    # Lifestyle Pie Charts
    for i, cluster in enumerate(df_lifestyle.columns):
        axes[0, i].pie(
            df_lifestyle[cluster], labels=df_lifestyle.index, autopct="%1.1f%%", colors=["#4CAF50", "#FFC107"]
        )
        axes[0, i].set_title(cluster)
    
    # Age Distributions
    for i, (cluster, ages) in enumerate(age_data.items()):
        sns.kdeplot(ages, ax=axes[1, i], color="skyblue", fill=True, alpha=0.7)
        axes[1, i].set_xlabel("Age")
        axes[1, i].set_xlim(0, 90)
        axes[1, i].axvline(np.mean(ages), color="red", linestyle="--", linewidth=2)
    
    # BMI Distributions
    for i, (cluster, bmi) in enumerate(bmi_data.items()):
        sns.kdeplot(bmi, ax=axes[2, i], color="skyblue", fill=True, alpha=0.7)
        axes[2, i].set_xlabel("BMI")
        axes[2, i].set_xlim(10, 45)
        axes[2, i].axvline(np.mean(bmi), color="red", linestyle="--", linewidth=2)
    
    # Health State Bar Charts
    for i, cluster in enumerate(df_health.columns):
        axes[3, i].bar(df_health.index, df_health[cluster], color=["#4CAF50", "#FFA500"], alpha=0.8, edgecolor="black")
        axes[3, i].set_ylim(0, 1)
    
    # Gender Distribution Bar Charts
    for i, cluster in enumerate(df_gender.columns):
        axes[4, i].bar(df_gender.index, df_gender[cluster], color=["#1E90FF", "#FF6347"], alpha=0.8, edgecolor="black")
        axes[4, i].set_ylim(0, 1)
    
    # Adjust layout and show plot
    plt.tight_layout()
    plt.show()


def get_feature_importance(data: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Computes feature importance by comparing the distributions of features (OTUs) between urban and rural samples.
    
    Parameters:
    data (pd.DataFrame): Feature abundance data with sample IDs as index.
    metadata (pd.DataFrame): Metadata containing sample IDs and lifestyle classification ('Urban' or 'Rural').
    
    Returns:
    pd.DataFrame: DataFrame containing OTU names, KS statistic, p-values, and species names, sorted by KS statistic.
    """
    
    # Extract sample IDs based on lifestyle
    rural_ids = metadata.loc[metadata['Lifestyle'] == 'Rural', 'Lane']
    urban_ids = metadata.loc[metadata['Lifestyle'] == 'Urban', 'Lane']
    
    # Filter data to keep only relevant samples
    rural_features = data.loc[data.index.intersection(rural_ids)]
    urban_features = data.loc[data.index.intersection(urban_ids)]
    
    # Initialize lists to store results
    otus, kstest, p_values = [], [], []
    
    for otu in data.columns:
        rural_dist = rural_features[otu].dropna()
        urban_dist = urban_features[otu].dropna()
        
        if rural_dist.empty or urban_dist.empty:
            continue  # Skip if there is no data for comparison
        
        ks_stat, p_value = ks_2samp(rural_dist, urban_dist)
        
        otus.append(otu)
        kstest.append(ks_stat)
        p_values.append(p_value)
    
    # Create DataFrame with KS test results
    feature_importance = pd.DataFrame({
        'OTU': otus,
        'KS statistic': kstest,
        'p-value': p_values
    })
    
    # Load taxonomy table and create species names
    otus_df = pd.read_csv('../data/taxonomy_table_otus.csv', index_col=0)
    otus_df['SpeciesFull'] = otus_df['Genus'] + ' ' + otus_df['Species']
    
    # Map OTU to species
    feature_importance['Specie'] = feature_importance['OTU'].map(otus_df['SpeciesFull'])
    
    return feature_importance.sort_values(by='KS statistic', ascending=False)