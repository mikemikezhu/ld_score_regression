# Description

For a phenotype of interest, we have collected the marginal statistics $\tilde{\mathbf{\beta}}$ for $M = 4268$ SNPs and the $M × M$ LD matrix $\mathbf{R}$ (i.e., pairwise SNP-SNP Pearson correlation). The marginal statistics are based on $N = 1000$ individuals. Download the marginal statistics and LD matrix from here:

https://drive.google.com/file/d/119Wmw9ockQNssHel3CZ88L2GhWqvW8ZJ/view?usp=sharing

For this question, you may also assume there is no population stratiﬁcation in this dataset. Both phenotype and genotype were standardized.

Implement the very basic LD score regression algorithm with a programming language of your choice (preferably Python or R) to estimate the heritability of the phenotype.

# Solution

First, calculate LD score of SNP $j$
$$l_j = \sum_{k=1}^Mr_{jk}^2$$

Then, we define the SSE loss function between observed value $N \tilde{\beta}_j^2$ and predicted value $\mathbb{E}[\chi_j^2] = N \frac{h^2}{M} l_j + 1$:

$$L = \sum_{j=1}^M \left(N \tilde{\beta}_j^2 - N \frac{h^2}{M} l_j - 1\right)^2$$

Therefore, we can calculate the heritability by solving the following objective function is:

$$\hat{h}^2 = \argmin_{h^2} \sum_{j=1}^M \left(N \tilde{\beta}_j^2 - N \frac{h^2}{M} l_j - 1\right)^2$$

Let $\frac{\partial L}{\partial h^2} = 0$:

Then we can calculate heritability $\hat{h}^2$:

$$\hat{h}^2 = \frac{\sum_{j=1}^M l_j \left( \tilde{\beta}_j^2 - \frac{1}{N} \right)}{\sum_{j=1}^M l_j^2 / M}$$