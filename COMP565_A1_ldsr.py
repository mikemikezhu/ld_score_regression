import pandas as pd
import numpy as np

# README: https://github.com/mikemikezhu/ld_score_regression

""" 
Constants 
"""

SAMPLE_SIZE = 1000
TOTAL_FEATURES = 4268

""" 
Data Loader
"""


class DataLoaderService:

    """ Load LD Matrix (Pearson Correlation) """

    def load_correlation(self, path: str,
                         compression: str) -> np.ndarray:

        result = pd.read_csv(path, compression=compression)
        print("LD matrix:")
        print(result.head())
        result = result.to_numpy()[:, 1:]
        print("LD matrix shape: {}".format(result.shape))
        return result

    """ Load marginal statistics """

    def load_marginal_statistics(self, path: str,
                                 compression: str) -> np.ndarray:

        result = pd.read_csv(path, compression=compression)
        print("Beta marginal statistics:")
        print(result.head())
        result = result.to_numpy()[:, 1]
        print("Beta marginal statistics shape: {}".format(result.shape))
        return result


""" 
LD Score
"""


class LdScoreCalculator:

    def calculate_ld_score(self, correlation: np.ndarray) -> np.ndarray:

        # LD score is the sum of the squared Pearson correlation between SNP j and SNP k
        result = np.sum(np.square(correlation), axis=1)
        print("LD Score:")
        print(result)
        print("LD Score shape: {}".format(result.shape))
        return result


""" 
LD Score Regression
"""


class LdScoreRegressionCalculator:

    def calculate_heritability(self, ld_score: np.ndarray,
                               marginal_statistics: np.ndarray,
                               sample_size: int,
                               total_features: int) -> float:

        numerator = np.sum(
            ld_score * (np.square(marginal_statistics) - 1 / sample_size))
        denominator = np.sum(np.square(ld_score)) / total_features
        return numerator / denominator


""" 
Main
"""


def main():

    # Load data
    data_loader = DataLoaderService()
    r = data_loader.load_correlation("data/LD.csv.gz", 'gzip')
    beta = data_loader.load_marginal_statistics(
        "data/beta_marginal.csv.gz", 'gzip')

    assert r.shape[0] == TOTAL_FEATURES
    assert r.shape[1] == TOTAL_FEATURES
    assert len(beta) == TOTAL_FEATURES

    # Calculate LD score
    ld_score_calculator = LdScoreCalculator()
    ld_score = ld_score_calculator.calculate_ld_score(r)

    assert len(ld_score) == TOTAL_FEATURES

    # Calculate heritability
    regression_calculator = LdScoreRegressionCalculator()
    heritability = regression_calculator.calculate_heritability(ld_score=ld_score,
                                                                marginal_statistics=beta,
                                                                sample_size=SAMPLE_SIZE,
                                                                total_features=TOTAL_FEATURES)
    print("Heritability: {}".format(heritability))


if __name__ == "__main__":
    main()
