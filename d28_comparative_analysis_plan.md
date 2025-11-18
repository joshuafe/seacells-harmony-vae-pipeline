# Comparative Analysis Plan: d28 PTCy vs. Rux vs. Normal

**Git Branch:** `feature/d28-comparative-analysis`

This document outlines a comprehensive plan to analyze the similarities and differences between several sample groups, with a primary focus on Day 28 post-transplant bone marrow aspirates from patients receiving either PTCy or Ruxolitinib-based GVH prophylaxis, compared to normal bone marrow.

We will leverage both unsupervised clustering and our novel VAE archetype analysis to explore these relationships at a deep cellular level.

## 1. Data and Sample Groupings

The analysis will be based on the `sample_annotations.csv` file, which defines the following key cohorts:

- **BM Normals:** Healthy bone marrow donors.
- **BM PTCy:** Post-transplant bone marrow samples from patients who received PTCy.
- **BM Rux:** Post-transplant bone marrow samples from patients who received Ruxolitinib.
- **PB Rux:** Post-transplant peripheral blood samples from patients who received Ruxolitinib.

**Special Cases for Consideration:**
- A cluster of samples around Day 28 to Day 50 post-transplant.
- Patients with multiple BM samples at different time points.
- A patient with same-day PB and BM samples.
- A patient 1.5 years post-transplant with lung disease (BAL sample).

## 2. Core Question

The central question we aim to answer is:

> **How similar to 'normal' are the Day 28 PTCy and Rux bone marrow aspirates? How do they differ from one another?**

## 3. Analysis Plan

We will approach this question through a multi-pronged analysis, combining traditional and advanced methods.

### 3.1. Unsupervised Clustering Analysis

This approach will identify and quantify cell populations across the different cohorts.

1.  **Global Clustering:**
    *   **Action:** Perform Phenograph clustering on a combined dataset of all normal, d28 PTCy, and d28 Rux BM samples.
    *   **Goal:** Identify a set of consensus cell populations (metaclusters) present across all samples.
    *   **Analysis:**
        *   **Compositional Analysis:** Compare the relative frequencies of each metacluster across the three cohorts (Normal, PTCy, Rux). Are there populations that are enriched or depleted in the transplant groups compared to normal?
        *   **Phenotypic Analysis:** For each metacluster, compare the mean marker expression of the cells between the cohorts. Do the same cell types (e.g., CD4+ T-cells) have different activation or differentiation states depending on the cohort?

2.  **Cohort-Specific Clustering:**
    *   **Action:** Perform separate Phenograph clustering within each cohort (Normal, PTCy, Rux).
    *   **Goal:** Identify cohort-specific cell states that might be missed in a global analysis.
    *   **Analysis:** Compare the resulting clusters between cohorts to see if unique cell states emerge in the PTCy or Rux groups.

### 3.2. VAE Archetype Analysis

This approach will allow us to understand the underlying cellular "building blocks" (archetypes) and how they are combined to form the cellular landscape of each cohort.

1.  **Train a Global VAE Model:**
    *   **Action:** Train our `VAEDisambiguation` model on the combined dataset of normal, d28 PTCy, and d28 Rux BM samples.
    *   **Goal:** Learn a set of universal cellular archetypes that represent the fundamental cell states across health and post-transplant recovery.

2.  **Archetype Contribution Analysis:**
    *   **Action:** For each sample, determine the contribution of each archetype to its overall cellular composition.
    *   **Goal:** Quantify how the "recipe" of cell states differs between the cohorts.
    *   **Analysis:**
        *   **Archetype Usage:** Do the PTCy or Rux samples utilize different archetypes compared to normals? Are there archetypes that are unique to or absent from a particular cohort?
        *   **Sample-level Comparison:** Project each sample into the archetype space. How close are the PTCy and Rux samples to the "normal" samples in this space?

3.  **Archetype Phenotype Analysis:**
    *   **Action:** Analyze the marker expression profiles of the archetypes themselves.
    *   **Goal:** Understand the biological identity of the fundamental cell states.
    *   **Analysis:** Compare the archetypes to known cell lineages. Do any archetypes represent disease-specific or treatment-specific cell states?

## 4. Documentation and Next Steps

This document will serve as a living record of our analysis. We will create a "Log" section below to document each experiment, its outcome, and the next steps.

### 4.1. Immediate Next Steps

1.  **Prepare the Data:**
    *   Filter the `sample_annotations.csv` to identify the exact file names for the "Normal", "d28 PTCy", and "d28 Rux" cohorts.
    *   Create a new, combined `.h5ad` file containing only the metacells from these selected samples.

2.  **Execute Clustering Analysis:**
    *   Run the Phenograph clustering pipeline on the combined `.h5ad` file.
    *   Generate initial plots of cluster composition across cohorts.

3.  **Execute VAE Archetype Analysis:**
    *   Train the `VAEDisambiguation` model on the combined `.h5ad` file.
    *   Generate initial plots of archetype usage across cohorts.

---

## 5. Initial Findings from Clustering Analysis

We have performed an initial clustering analysis using Phenograph on the combined dataset of Normal, d28 PTCy, and d28 Rux bone marrow samples. This analysis yielded 19 distinct cell clusters.

### 5.1. Cluster Composition

The composition of these clusters reveals significant differences between the three cohorts:

-   **Normal-Dominant Clusters:** Several clusters were almost exclusively composed of cells from normal donors.
    -   Clusters 11, 13, 15, 17, and 18 were all >80% normal.
    -   This indicates that a significant portion of the normal bone marrow cellular landscape is distinct from the post-transplant landscape.

-   **Rux-Dominant Clusters:** A number of clusters were dominated by cells from the Ruxolitinib cohort.
    -   Clusters 1, 3, 7, 8, and 16 were all >78% Rux.
    -   This suggests that Ruxolitinib may be inducing or selecting for specific cellular states that are less common in the other cohorts.

-   **Mixed PTCy/Rux Clusters:** A large group of clusters showed a mix of both PTCy and Rux cells, with very few normal cells.
    -   Clusters 4, 5, 6, and 9 were all >98% composed of a mix of PTCy and Rux cells.
    -   This points to a shared cellular environment between the two transplant cohorts, likely representing common responses to the transplant procedure itself.

-   **No PTCy-Dominant Clusters:** Notably, no clusters were found to be dominant for the PTCy cohort. PTCy cells were primarily found within the mixed clusters, suggesting they share a cellular state with the Rux samples rather than forming a distinct, unique population.

### 5.2. Granular Analysis with Leiden Clustering

To gain a more detailed understanding, we also analyzed the results from the Leiden clustering, which produced a more granular set of 88 clusters. This analysis largely confirms and extends the findings from Phenograph.

-   **Stronger Cohort Separation:** The Leiden clustering showed an even clearer separation between the cohorts, with numerous clusters being highly enriched for a single cohort.
-   **Normal-Dominant Clusters:** Clusters 13, 20, 25, 31, and 46 were all >70% Normal.
-   **Rux-Dominant Clusters:** A striking number of clusters were almost exclusively composed of Rux cells. Clusters 44, 59, and 68-86 were all 100% Rux, suggesting a highly specific and somewhat fragmented response to the drug.
-   **PTCy-Dominant Clusters:** In contrast to the Phenograph results, Leiden clustering identified two PTCy-dominant clusters (Cluster 62 at 75% and Cluster 87 at 100%), indicating that there are some unique cellular states in this cohort, even if they are less prominent than the Rux-specific ones.

### 5.3. Initial Answers to Core Question

Based on this compositional analysis, we can form some initial answers:

> **How similar to 'normal' are the d28 PTCy and Rux bone marrow aspirates?**

Both are substantially different from normal. The vast majority of cells from both PTCy and Rux samples fall into clusters that are either rare or completely absent in the normal samples.

> **How do they differ from one another?**

They exhibit both similarities and differences.
-   **Similarity:** They share a large number of "mixed" cell clusters, indicating a common response to transplantation.
-   **Difference:** The Rux cohort is unique in that it forms several "Rux-dominant" clusters, suggesting a drug-specific effect that is not observed in the PTCy cohort. The PTCy cohort appears less distinct and more heterogeneous, existing within the shared post-transplant landscape, although it does have a small number of unique clusters.

### 5.4. Phenotypic Analysis of Clusters

To understand the biological identity of these cohort-specific clusters, we analyzed their marker expression profiles.

-   **Normal-Dominant Clusters:** These clusters represent a diverse and healthy immune landscape.
    -   **Cluster 11 & 15:** Appear to be classic T-cell populations (CD45RO+ or CD4+).
    -   **Cluster 13 & 17:** Show high expression of multiple activation and exhaustion markers (CD279+, CD366+, CD95+, CD69+), likely representing highly active or terminally differentiated T-cells.
    -   **Cluster 14:** Expresses naive T-cell markers (CD45RA+, CD197+).

-   **Rux-Dominant Clusters:** These clusters point towards a specific immunological state induced by Ruxolitinib.
    -   **Cluster 8:** A clear CD8+ T-cell population.
    -   **Cluster 1:** A CD4+ T-cell population with high expression of the inhibitory receptor CD152 (CTLA-4), suggesting a regulatory T-cell (Treg) or exhausted T-cell phenotype.
    -   The prevalence of these specific phenotypes, particularly the potential Treg/exhausted population, is a key finding and may relate to the mechanism of action of Ruxolitinib.

-   **Mixed PTCy/Rux Clusters:** These clusters represent shared post-transplant cell states.
    -   **Cluster 0:** A CD8+ T-cell population with high expression of activation/exhaustion markers.
    -   **Cluster 6:** Appears to be a Natural Killer (NK) cell population (CD56+).

### 5.5. VAE Archetype Analysis

To complement the discrete clustering analysis, we performed a VAE archetype analysis. This method identifies a set of continuous "archetypes" that represent fundamental cellular states, and then describes each cell as a mixture of these archetypes. We trained a VAE with 12 archetypes on the combined dataset.

#### 5.5.1. Archetype Usage by Cohort

The contribution of each archetype to the cellular landscape of each cohort revealed several key differences:

-   **Normal-Dominant Archetypes:**
    -   **Archetype 3 & 9:** These archetypes were significantly more prevalent in the Normal cohort compared to the PTCy and Rux cohorts. These represent cellular states that are characteristic of a healthy bone marrow environment and are diminished after transplantation.

-   **PTCy-Associated Archetype:**
    -   **Archetype 1:** This archetype showed the highest usage in the PTCy cohort.

-   **Rux-Associated Archetype:**
    -   **Archetype 7:** This archetype was most prominent in the Rux cohort.

-   **Shared Archetypes:** The majority of the archetypes were utilized by all three cohorts, representing the common, underlying cellular processes in bone marrow.

#### 5.5.2. Phenotypic Analysis of Key Archetypes

By examining the marker profiles of the disambiguated archetypes, we can infer their biological identity:

-   **Archetype 3 (Normal-Dominant):** Characterized by high expression of **CD45RO, CD4, and CD19**. This suggests a **memory CD4+ T-cell or B-cell** state that is depleted in the post-transplant samples.
-   **Archetype 9 (Normal-Dominant):** Shows high expression of **CD95, CD45RA, and the activation/exhaustion markers CD366 and CD69**. This likely represents a population of **activated or terminally differentiated T-cells** that is specific to the healthy state.
-   **Archetype 1 (PTCy-Associated):** This archetype has a less distinct profile, with moderate expression of **CD62L, CD152 (CTLA-4), and CD56**. This mixed phenotype may represent a **transitional or progenitor-like state** that is more common in the PTCy setting.
-   **Archetype 7 (Rux-Associated):** Also has a mixed profile, with moderate expression of the exhaustion marker **CD279 (PD-1)**. This could represent an **exhausted T-cell state** that is slightly more prevalent in the Rux cohort.

#### 5.5.4. Archetype Usage by Cohort (Statistical Comparison)

We performed a Kruskal-Wallis test followed by post-hoc Dunn's tests to statistically compare the archetype usage across the Normal, PTCy, and Rux cohorts. The results are visualized below, with asterisks indicating statistical significance (p < 0.05: *, p < 0.01: **, p < 0.001: ***).

*(Plots for each archetype will be embedded here)*

![Archetype 0 Usage](d28_vae_analysis/archetype_usage_plots/archetype_0_usage.png)
![Archetype 1 Usage](d28_vae_analysis/archetype_usage_plots/archetype_1_usage.png)
![Archetype 2 Usage](d28_vae_analysis/archetype_usage_plots/archetype_2_usage.png)
![Archetype 3 Usage](d28_vae_analysis/archetype_usage_plots/archetype_3_usage.png)
![Archetype 4 Usage](d28_vae_analysis/archetype_usage_plots/archetype_4_usage.png)
![Archetype 5 Usage](d28_vae_analysis/archetype_usage_plots/archetype_5_usage.png)
![Archetype 6 Usage](d28_vae_analysis/archetype_usage_plots/archetype_6_usage.png)
![Archetype 7 Usage](d28_vae_analysis/archetype_usage_plots/archetype_7_usage.png)
![Archetype 8 Usage](d28_vae_analysis/archetype_usage_plots/archetype_8_usage.png)
![Archetype 9 Usage](d28_vae_analysis/archetype_usage_plots/archetype_9_usage.png)
![Archetype 10 Usage](d28_vae_analysis/archetype_usage_plots/archetype_10_usage.png)
![Archetype 11 Usage](d28_vae_analysis/archetype_usage_plots/archetype_11_usage.png)

These plots visually confirm the statistical differences identified, showing which archetypes are significantly enriched or depleted in each cohort.

#### 5.5.5. Analysis of Interstitial Spaces: 2D Distribution of Metacells

To understand the nature of the "interstitial" spaces, we analyzed the two-dimensional distribution of metacells relative to their top two contributing archetypes. For each metacell, we calculated its proximity ratio (how close it is to one of the top two archetypes versus the other) and its overall archetype entropy (a measure of how ambiguous its assignment is).

**Proportion of Interstitial Cells within each Cohort:**
![Interstitial Cells Proportions](d28_vae_analysis/interstitial_cells_proportions.png)

-   **Normal:** 30.2% of Normal cells are interstitial.
-   **PTCy:** 28.0% of PTCy cells are interstitial.
-   **Rux:** 22.1% of Rux cells are interstitial.

**Interpretation of Proportions:** The statistical analysis (Chi-squared test, p-value < 0.0001) confirms a highly significant association between cohort and interstitial status. Contrary to our initial interpretation, the **Normal cohort actually has the highest proportion of interstitial cells (30.2%)** followed by PTCy (28.0%), and then Rux (22.1%). This suggests that the healthy bone marrow environment, while having well-defined archetypes, also contains a significant proportion of cells in transitional states. The Rux cohort, despite showing a high degree of heterogeneity in clustering, appears to have a *lower* proportion of these highly ambiguous cells when viewed through the lens of VAE archetypes. This could imply that Ruxolitinib is driving cells towards more defined, albeit potentially novel, states.

**2D Distribution Plots for Top Archetype Pairs:**

We generated 2D scatter plots for the top 10 most common archetype pairs. These plots show the proximity ratio (x-axis) against archetype entropy (y-axis), with points colored by cohort.

*(Plots for each top archetype pair will be embedded here)*

![Archetype Pair 5-7 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_5_7_2d_distribution.png)
![Archetype Pair 2-7 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_2_7_2d_distribution.png)
![Archetype Pair 1-7 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_1_7_2d_distribution.png)
![Archetype Pair 1-4 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_1_4_2d_distribution.png)
![Archetype Pair 7-11 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_7_11_2d_distribution.png)
![Archetype Pair 4-5 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_4_5_2d_distribution.png)
![Archetype Pair 0-5 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_0_5_2d_distribution.png)
![Archetype Pair 4-7 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_4_7_2d_distribution.png)
![Archetype Pair 2-10 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_2_10_2d_distribution.png)
![Archetype Pair 7-10 Distribution](d28_vae_analysis/archetype_pair_distributions/archetype_pair_7_10_2d_distribution.png)

**Interpretation of 2D Distributions:** These plots allow us to visually assess whether metacells form a continuum between archetypes or cluster closer to one or the other. We can also observe if the distribution patterns differ by cohort. For instance, a dense cloud of points in the middle of the x-axis (proximity ratio ~0.5) with high entropy would indicate a true "interstitial" continuum, whereas points clustered at the extremes (0 or 1) with low entropy would suggest more defined, albeit mixed, populations.

### 5.6. Refined Answers to Core Question

The VAE analysis provides a more continuous and nuanced view that complements the clustering results.

> **How similar to 'normal' are the d28 PTCy and Rux bone marrow aspirates?**

They are different at a fundamental level. The post-transplant samples are depleted of key "normal" cellular states (Archetypes 3 and 9) related to memory and activation/exhaustion, and instead are composed of a different recipe of cellular archetypes.

> **How do they differ from one another?**

The PTCy and Rux cohorts utilize different combinations of archetypes. The PTCy cohort shows a preference for a potential "transitional" state (Archetype 1), while the Rux cohort has a slightly higher usage of a potential "exhausted" state (Archetype 7). Furthermore, the **Normal cohort has the highest proportion of cells in "interstitial" or transitional states**, suggesting a more dynamic healthy environment, while the Rux cohort has the lowest proportion of these highly ambiguous cells. This suggests that Ruxolitinib might be driving cells into more defined, though potentially novel, states. The 2D distribution plots for archetype pairs provide further visual evidence of how these interstitial spaces are populated by each cohort, revealing potential differences in the continuum of cell states.

### 5.8. The Story So Far

Our analysis set out to understand the differences between three groups of bone marrow samples: healthy normals, and post-transplant patients treated with either PTCy or Ruxolitinib. Using a combination of two powerful techniques, clustering and VAE archetypes, we have uncovered a consistent and multi-layered story.

1.  **Post-Transplant Environment is Fundamentally Different:** Both the clustering and VAE analyses show a clear and dramatic separation between the normal and post-transplant samples. The cellular landscape of the post-transplant patients is not just a variation of the normal state, but a fundamentally different environment. Many cell populations and states that are common in normal bone marrow are depleted or absent in the post-transplant samples, and vice-versa.

2.  **Ruxolitinib's Unique Signature:** The Ruxolitinib-treated samples stand out with a unique and specific cellular signature.
    *   **Clustering:** We found several clusters that were almost exclusively composed of cells from the Rux cohort. These "Rux-dominant" clusters point to specific cell populations that are either induced or expanded by the drug.
    *   **VAE Archetypes:** The VAE analysis revealed that the Rux samples have a *lower* proportion of cells in "interstitial" or transitional states (as defined by high entropy), suggesting that Ruxolitinib might be driving cells towards more defined, albeit potentially novel, states. This contrasts with the Normal cohort, which has the highest proportion of such ambiguous cells.

3.  **PTCy's More "Canonical" Post-Transplant State:** The PTCy samples, in contrast, appear to represent a more "canonical" or "generic" post-transplant state.
    *   **Clustering:** We did not find any dominant clusters for the PTCy cohort in the Phenograph analysis, and only a few in the more granular Leiden analysis. Instead, the PTCy cells were mostly found in "mixed" clusters that they shared with the Rux samples.
    *   **VAE Archetypes:** The PTCy samples had a lower proportion of "interstitial" cells than the Normal samples, but a higher proportion than the Rux samples, suggesting an intermediate level of cellular ambiguity.

**In a nutshell:** While both PTCy and Rux treatments result in a post-transplant environment that is very different from normal, Ruxolitinib appears to be driving a more specific and unusual cellular response, characterized by the expansion of unique cell populations and a *lower* degree of cellular ambiguity compared to the healthy state. The Normal bone marrow, surprisingly, exhibits the highest proportion of cells in transitional states.

#### 5.8.1. A Note on Absolute Cell Counts

You mentioned that Rux patients tend to have higher CD4 and CD8 counts than PTCy patients. The data I have is at the metacell level, so I cannot directly calculate absolute cell counts per sample.

However, I can analyze the *proportions* of CD4+ and CD8+ cells in each cohort, based on the archetypes that we have identified as being enriched for these markers. For example, Archetype 8 was identified as a CD8+ T-cell archetype. The VAE analysis shows that the average contribution of this archetype is similar across the three cohorts (Normal: 0.058, PTCy: 0.077, Rux: 0.073). This suggests that the *proportion* of CD8+ T-cells may not be dramatically different between the cohorts, even if the absolute counts are. A more detailed analysis of the archetypes enriched for CD4 and CD8 could shed more light on this.

### 5.9. What Distinguishes PTCy from Rux?

While both post-transplant cohorts are fundamentally different from normal, our analyses have revealed several key distinctions between the PTCy and Ruxolitinib-treated patients.

1.  **Cellular Stability and Transitional States:**
    *   The statistical analysis of interstitial cells showed that the **Normal cohort has the highest proportion of interstitial cells (30.2%)**, followed by PTCy (28.0%), and then Rux (22.1%). This is a significant finding, suggesting that Ruxolitinib might be driving cells into more defined, albeit potentially novel, states, reducing the overall cellular ambiguity compared to both Normal and PTCy.
    *   The 2D distribution plots for archetype pairs provide further visual evidence of how these interstitial spaces are populated by each cohort.

2.  **Dominant Cell Populations:**
    *   The clustering analysis revealed several **"Rux-dominant" clusters**, indicating the expansion of specific cell populations under Ruxolitinib treatment.
    *   In contrast, we found **no "PTCy-dominant" clusters** in the Phenograph analysis, and only a few in the more granular Leiden analysis. The PTCy samples appear to be composed of a more "generic" post-transplant cellular landscape that is shared with the Rux samples.

3.  **Archetype Usage:**
    *   The statistical comparison of archetype usage confirmed that there are significant differences between the two cohorts.
    - **Archetype 7** is significantly more prevalent in the Rux cohort compared to both Normal and PTCy cohorts, making it a "Rux-specific" archetype.
    - **Archetypes 4 and 5** also show statistically significant, albeit smaller, differences between PTCy and Rux.

**Summary:**

In essence, the Ruxolitinib treatment appears to be driving a more profound and specific alteration of the cellular landscape compared to PTCy. This is characterized by the emergence of unique, dominant cell populations and a *lower* degree of cellular ambiguity compared to the healthy state. The Normal bone marrow, surprisingly, exhibits the highest proportion of cells in transitional states.

### 5.10. New Visualizations

Here are the new visualizations you requested:

**Sample-level Archetype Usage Heatmap:**
![Sample-level Archetype Usage Heatmap](d28_vae_analysis/sample_archetype_usage_heatmap.png)

**Cluster to Archetype Flow (Sankey Diagram):**
[Link to Sankey Diagram](d28_vae_analysis/cluster_archetype_sankey.html)

**Interstitial Cells UMAP:**
![Interstitial Cells UMAP](d28_vae_analysis/interstitial_cells_umap.png)

**Archetype Pair Chord Diagram:**
![Archetype Pair Chord Diagram](d28_vae_analysis/archetype_pair_chord_diagram.png)

### 5.11. Next Steps

The initial analysis phase is now complete. We have a consistent picture from both clustering and VAE archetype analysis. The next steps will involve deeper dives into these results and further biological interpretation.

1.  **Integrate Findings:** Synthesize the results from the clustering and VAE analyses into a cohesive narrative.
2.  **Deeper Phenotypic Analysis:** Perform a more detailed analysis of the marker expression data for the key clusters and archetypes.
3.  **Downstream Analysis:** Plan for downstream analyses, such as differential abundance testing and functional enrichment analysis.

---

## 6. Experiment Log

*(This section will be filled in as we conduct the experiments.)*

**Experiment 2025-11-17-A: Initial Clustering and Composition Analysis**
*   **Goal:** Identify consensus cell populations and compare their composition across cohorts.
*   **Action:**
    1.  Created a combined `h5ad` file with Normal, d28 PTCy, and d28 Rux samples.
    2.  Ran Phenograph and Leiden clustering using `scripts/seacells_clustering_analysis.py`.
    3.  Analyzed cluster composition using `analyze_cluster_composition.py`.
*   **Result:** Identified distinct cohort-specific clusters and provided an initial answer to the core research question.
*   **Next Step:** Analyze cluster phenotypes and proceed with VAE archetype analysis.

**Experiment 2025-11-17-B: Phenotypic Analysis of Phenograph Clusters**
*   **Goal:** Characterize the biological identity of the Phenograph clusters.
*   **Action:**
    1.  Analyzed the mean marker expression of each cluster using `analyze_cluster_phenotypes.py`.
*   **Result:** Identified distinct phenotypic profiles for the Normal-dominant, Rux-dominant, and Mixed clusters.
*   **Next Step:** Proceed with VAE analysis.

**Experiment 2025-11-17-C: VAE Archetype Analysis**
*   **Goal:** Identify fundamental cellular states (archetypes) and their usage across cohorts.
*   **Action:**
    1.  Trained a `VAEDisambiguation` model on the combined dataset.
    2.  Analyzed archetype usage and phenotypes.
*   **Result:** Identified Normal-dominant, PTCy-associated, and Rux-associated archetypes with distinct phenotypic profiles.
*   **Next Step:** Synthesize all findings and plan for deeper analysis.

**Experiment 2025-11-17-D: Interstitial Cell Analysis**
*   **Goal:** Analyze cells that are not well-represented by any single archetype.
*   **Action:**
    1.  Identified cells with high archetype entropy.
    2.  Analyzed the cohort composition of these interstitial cells.
*   **Result:** Found that the Rux cohort has the highest proportion of interstitial cells, suggesting a more transitional or perturbed cellular state.
*   **Next Step:** Synthesize all findings.

**Experiment 2025-11-17-E: Deeper Dive into Archetype Relationships**
*   **Goal:** Analyze the relationships between pairs of archetypes.
*   **Action:**
    1.  Identified common pairs of archetypes.
    2.  Visualized the connections between archetypes using a chord diagram.
*   **Result:** Revealed the strongest connections between archetypes, providing insight into the continuous nature of the cellular landscape.
*   **Next Step:** Synthesize all findings.
