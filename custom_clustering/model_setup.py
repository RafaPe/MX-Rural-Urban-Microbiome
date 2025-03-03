import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score, silhouette_score, calinski_harabasz_score
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.base import BaseEstimator, ClusterMixin
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import make_scorer
from sklearn.cluster import HDBSCAN
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics.cluster import homogeneity_score, adjusted_mutual_info_score
import optuna

class CustomClustering(ClusterMixin, BaseEstimator):
    def __init__(self, algorithm='kmeans', filter_criteria='presence', filter_threshold=25, **kwargs):
        """
        Parameters:
        -----------
        k : int, default=2
            The number of clusters to generate.

        algorithm : str, default='kmeans'
            Clustering algorithm to use ('kmeans' or 'gaussian_mixture').

        filter_criteria : str, default='presence'
            Feature filtering criterion ('abundance' or 'presence').

        filter_threshold : int or float, default=25
            Threshold for feature filtering. Interpreted as:
            - Number of samples for 'presence'.
            - Mean abundance multiplier for 'abundance'.
        """
        self.algorithm = algorithm
        self.filter_criteria = filter_criteria
        self.filter_threshold = filter_threshold
        self.kwargs = kwargs

    def fit(self, X, y=None):
        """
        Fits the clustering algorithm to the data after applying the filtering criterion.
        """
        self.data_ = X.copy()
    
        # Feature filtering
        if self.filter_criteria == 'abundance':
            column_sums = X.sum()
            threshold = column_sums.mean() * self.filter_threshold
            self.filtered_data_ = X.loc[:, column_sums > threshold]
    
        elif self.filter_criteria == 'presence':
            binary_df = (X > 0).astype(int)
            presence = binary_df.sum(axis=0) > self.filter_threshold
            self.filtered_data_ = X.loc[:, presence]
    
        else:
            raise ValueError("Invalid filter_criteria. Choose 'abundance' or 'presence'.")

        self.selected_features_ = self.filtered_data_.columns
    
        # Scaling the data
        self.scaler_ = StandardScaler()
        self.data_scaled_ = self.scaler_.fit_transform(self.filtered_data_)
    
        # Clustering
        if self.algorithm == 'kmeans':
            self.model_ = KMeans(**self.kwargs)
            self.labels_ = self.model_.fit_predict(self.data_scaled_)
    
        elif self.algorithm == 'gaussian_mixture':
            self.model_ = GaussianMixture(**self.kwargs)
            self.labels_ = self.model_.fit_predict(self.data_scaled_)

        elif self.algorithm == 'hdbscan':
            self.model_ = HDBSCAN(**self.kwargs)
            self.labels_ = self.model_.fit_predict(self.data_scaled_)
    
        else:
            raise ValueError("Invalid algorithm. Choose 'kmeans' or 'gaussian_mixture'.")
    
        # Scoring metrics
        # self.silhouette_score_ = silhouette_score(self.data_scaled_, self.labels_)
        # self.davies_bouldin_score_ = davies_bouldin_score(self.data_scaled_, self.labels_)
        # self.calinski_harabasz_score_ = calinski_harabasz_score(self.data_scaled_, self.labels_)
    
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
        X_scaled = self.scaler_.transform(X_filtered)
        return self.model_.predict(X_scaled)


def optimization_function(trial, data, labels):
    algorithm = trial.suggest_categorical('algorithm', ['kmeans', 'gaussian_mixture'])
    filter_criteria = trial.suggest_categorical('filter_criteria', ['presence', 'abundance'])
    filter_threshold = trial.suggest_float('filter_threshold', 0, 120, step=5) if filter_criteria == 'presence' else trial.suggest_float('filter_threshold', 0, 2.7, step=0.1)

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