![image](https://github.com/mango1213/Metabolic-heterogeneity-repository/assets/82602606/34208866-2159-4e55-99eb-b839f277351f)# Metabolic-heterogeneity-repository
Metabolic heterogeneity affects trastuzumab response and survival in HER2-positive advanced gastric cancer

# Abstract
Background: Trastuzumab is the only first-line treatment targeted against the human epidermal growth factor receptor 2 (HER2) approved for patients with HER2-positive advanced gastric cancer. The impact of metabolic heterogeneity on trastuzumab treatment efficacy remains unclear. /
Methods: Spatial metabolomics via high mass resolution imaging mass spectrometry was performed in pretherapeutic biopsies of patients with HER2-positive advanced gastric cancer in a prospective multicenter observational study. The mass spectra, representing the metabolic heterogeneity within tumor areas, were grouped by K–means clustering algorithm. Simpson’s diversity index was applied to compare the metabolic heterogeneity level of individual patients. /
Results: Clustering analysis revealed metabolic heterogeneity in HER2-positive gastric cancer patients and uncovered nine tumor subpopulations. High metabolic heterogeneity was shown as a factor indicating sensitivity to trastuzumab (p=0.008) and favorable prognosis at trend level. Two of the nine tumor subpopulations associated with favorable prognosis and trastuzumab sensitivity, and one subpopulation associated with poor prognosis and trastuzumab resistance. /
Conclusions: This work revealed that tumor metabolic heterogeneity associated with prognosis and trastuzumab response based on tissue metabolomics of HER2-positive gastric cancer. Tumor metabolic subpopulations may provide an association with trastuzumab therapy efficacy.

# code
## 1_df_clin_cluster_merged.xlsx
Excel of cluster presence of samples in different numbers of K from 2 to 10.
## main_for_BJC.py
Codes used for Cox proportional hazards regression model whose performace was evaluated by Akaike information criterion (AIC) and survival analyses.
