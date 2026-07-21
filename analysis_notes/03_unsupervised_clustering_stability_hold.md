# Unsupervised clustering and cluster-stability analysis

Status: **analysis completed; retain for later manuscript integration after all planned additional analyses are complete. Do not generate a figure yet.**

## Dataset and preprocessing

Input: `SDY2583_integrated_clinical_immune_score_matrix_ALL10_with_CP23.csv`

- 829 complete cases across the 10 principal integrated scores, age, sex, disease group, immunotype, and cohort.
- Each principal score was residualized for age, sex, and disease group, then standardized.
- Primary clustering space: 10-dimensional disease-adjusted integrated-score space.
- Candidate solutions: k = 2–8.
- Methods: k-means and Ward hierarchical clustering; diagonal Gaussian-mixture information criteria as a supplementary check.
- Stability: repeated 80% subsampling and independent training/validation clustering.

## Cluster tendency

Hopkins statistic:

- mean = 0.757
- 95% empirical interval = 0.742–0.770

The score space is non-random and contains structure. However, Hopkins does not establish that the structure is composed of sharply discrete clusters; gradients and elongated manifolds can also produce high values.

## Internal validity

K-means results:

| k | Silhouette | Calinski–Harabasz | Davies–Bouldin | Minimum cluster n |
|---:|---:|---:|---:|---:|
| 2 | 0.207 | 187.99 | 1.971 | 286 |
| 3 | 0.175 | 143.31 | 2.031 | 185 |
| 4 | 0.110 | 122.77 | 2.070 | 81 |
| 5 | 0.106 | 111.72 | 2.000 | 70 |
| 6 | 0.104 | 102.90 | 1.944 | 41 |
| 7 | 0.110 | 95.30 | 1.835 | 34 |
| 8 | 0.100 | 88.93 | 1.910 | 31 |

The best silhouette and Calinski–Harabasz values occurred at k = 2, but the absolute silhouette remained modest, indicating broad within-cluster heterogeneity and overlap.

Diagonal Gaussian-mixture BIC decreased through k = 6 and then increased only slightly. The absence of a sharp isolated optimum is more consistent with continuous or multi-scale structure than with one uniquely supported discrete cluster count.

## Subsampling stability

| k | Mean ARI | Median ARI | 10th-percentile ARI | Mean cluster Jaccard |
|---:|---:|---:|---:|---:|
| 2 | 0.954 | 0.964 | 0.911 | 0.975 |
| 3 | 0.777 | 0.826 | 0.488 | 0.834 |
| 4 | 0.508 | 0.406 | 0.259 | 0.594 |
| 5 | 0.773 | 0.808 | 0.583 | 0.826 |

The k = 2 solution was highly stable under subsampling. The k = 3 and k = 5 solutions showed moderate average stability but substantial lower-tail instability. The k = 4 solution was poor.

## Training–validation reproducibility

Separate cluster solutions were fitted in the original training and validation subsets.

| k | Mean centroid cosine | Minimum centroid cosine | Validation assignment ARI |
|---:|---:|---:|---:|
| 2 | 0.691 | 0.616 | 0.666 |
| 3 | 0.354 | −0.352 | 0.405 |
| 4 | 0.449 | −0.514 | 0.180 |
| 5 | 0.416 | −0.470 | 0.255 |

Only k = 2 showed moderate cross-cohort reproducibility. Solutions with three or more clusters did not reproduce sufficiently to justify new biological subtypes.

## Relationship to the original five immunotypes

Agreement between new k-means clusters and the original five immunotypes was low:

- k = 2: ARI = 0.067; NMI = 0.071
- k = 3: ARI = 0.086; NMI = 0.101
- k = 4–8: ARI approximately 0.05–0.07

This is expected because the original immunotypes were derived from a much larger cellular-state feature space. The ten integrated scores summarize dominant biological axes but do not recreate the original clustering solution.

## k = 2 interpretation

The stable k = 2 solution consisted of:

- Cluster 1, n = 543: uniformly below-average remodeling across all ten integrated scores.
- Cluster 2, n = 286: uniformly above-average remodeling across all ten integrated scores.

PC1 alone classified the k = 2 solution with AUC = 0.9997. Mean PC1 values were −0.985 and +1.871, whereas mean PC2 values were approximately zero in both clusters.

Therefore, k = 2 is not a distinct new biological taxonomy. It is effectively a dichotomization of the continuous global multi-compartment remodeling gradient already captured by PC1.

The high-remodeling cluster was enriched for G4_Chronic and G5_Suppressive but remained heterogeneous:

- G1_Naive 9.8%
- G2_Primed 24.1%
- G3_Progressive 25.5%
- G4_Chronic 21.0%
- G5_Suppressive 19.6%

Disease composition was similar after residualization: 55.2% cancer in the high-remodeling cluster versus 50.5% in the low-remodeling cluster.

## k = 3 interpretation

The k = 3 solution produced biologically recognizable profiles:

1. Low-remodeling state, n = 451.
2. Myeloid/monocyte/NK-oriented state, n = 193.
3. CD8/T-NK/regulatory-checkpoint-oriented state, n = 185.

The two high-remodeling states aligned qualitatively with the PCA PC2 poles and with G5_Suppressive and G4_Chronic enrichment, respectively. However:

- silhouette was only 0.175;
- subsample ARI had a 10th percentile of 0.488;
- training–validation assignment ARI was 0.405;
- one validation centroid showed negative cosine similarity before label alignment.

Thus, k = 3 is useful as a descriptive decomposition of PC1/PC2 geometry but is not sufficiently reproducible to define new patient subtypes.

## Covariate sensitivity

Agreement between age/sex-adjusted and disease-adjusted solutions was only moderate:

- k = 2 ARI = 0.470
- k = 3 ARI = 0.423
- k = 4 ARI = 0.346
- k = 5 ARI = 0.325

This confirms that cluster membership is sensitive to whether the cancer-versus-healthy mean shift is retained. Disease-adjusted clustering is the appropriate basis for evaluating biology beyond case-control status, but the sensitivity further argues against presenting new clusters as robust natural classes.

## Primary conclusion

The ten-score space contains strong non-random structure, but that structure is better represented by continuous and partially branching remodeling dimensions than by a new set of discrete, reproducible clusters.

- k = 2 is highly stable but merely dichotomizes the continuous PC1 remodeling burden.
- k = 3 recovers low-remodeling, CD8/T-NK-oriented, and myeloid/NK-oriented configurations, but cross-cohort stability is inadequate for subtype claims.
- No solution with k >= 3 meets the combined requirements of internal separation, method agreement, and training–validation reproducibility.

Accordingly, the manuscript should **not introduce new immune-remodeling subtypes** from these ten scores. The PCA and original-immunotype mapping provide a more faithful and statistically defensible representation.

## Working Results text

Unsupervised clustering was used to test whether the ten integrated remodeling scores supported additional discrete immune phenotypes beyond the previously defined immunotypes. Although the disease-adjusted score space showed non-random structure, internal-validity metrics favored a two-cluster solution with only modest separation. This solution was highly stable under subsampling but almost completely represented a binary division of the continuous first principal component into low- and high-remodeling groups. A three-cluster solution separated a low-remodeling group from myeloid/monocyte/NK-oriented and CD8/T-NK-oriented configurations, but its silhouette coefficient was low and its cluster assignments reproduced only weakly across the original training and validation subsets. Solutions with four or more clusters showed poorer internal separation and limited cross-cohort reproducibility. These findings indicate that the integrated score space is structured but does not support a new set of sharply discrete and reproducible patient subtypes.

## Working Discussion text

The unsupervised analyses reinforced the interpretation of peripheral immune remodeling as a continuous multidimensional architecture rather than a collection of newly discoverable classes. The only highly stable clustering solution divided participants according to overall remodeling burden and was almost perfectly recoverable from PC1 alone. A biologically interpretable three-state solution reproduced the low-remodeling, CD8/T-NK-oriented, and myeloid/NK-oriented regions of the PCA space, but it lacked sufficient training–validation stability to justify formal subtype labels. Avoiding new cluster nomenclature is therefore important: the principal components and original immunotypes capture the underlying organization more faithfully, whereas additional clustering would impose artificial boundaries on overlapping immune configurations.

## Figure policy

Do not generate a clustering-specific figure. The clustering result is primarily a negative model-selection finding that protects the manuscript from overclaiming. It may be mentioned briefly in Results or Supplementary Methods/Results and used to justify the final continuous-network/PCA representation. If space is limited, retain the full analysis in Supplementary material or the reproducibility repository rather than adding a separate main-text figure.
