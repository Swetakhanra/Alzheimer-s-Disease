# Alzheimer's Disease Gene Expression Analysis
This project uses RNA-seq count data (GSE153873) to identify differentially expressed genes between Alzheimer's Disease (AD), Old, and Young samples. 

### 🔬 Methods
- STAR count matrix
- DESeq2 differential expression analysis
- Pathway enrichment using Enrichr (KEGG 2021 Human)
- Visualization and interpretation of top genes and pathways

### Key Findings
- CRH downregulated in AD (stress-related)
- TNF, RAGE signaling, LTD, and CAMs are altered in AD brains
  
### Visualizations 
Boxplots of normalized expression (VST) values were created for CRH and KLF15 across conditions (AD, Old, Young).

<img width="835" height="496" alt="image" src="https://github.com/user-attachments/assets/edb39e64-b05f-4cc7-a99f-d8cfad476f1b" />

In KLF15 gene-- AD is higher than old and young, in CRH, LOC101926975, in C2orf16 gene AD is lower than old and young, and in DUSP6 AD and young is overlapping and old is above.
 

<img width="859" height="511" alt="image" src="https://github.com/user-attachments/assets/b63d38f8-afe9-4cd1-b148-99c1712edeeb" />

In our differential gene expression analysis comparing Alzheimer's Disease (AD) with healthy aging (Old), two key genes stood out for further investigation, CRH and KLF15.

CRH (Left Plot)
- Lower in AD than in Old
- AD group expression is compressed around 4
- Very significant difference → This drop is not due to chance
- Suggests CRH is downregulated in Alzheimer’s
- This gene plays a central role in the hypothalamic-pituitary-adrenal (HPA) axis, which regulates stress response and has been previously implicated in neurodegenerative conditions. Its decreased expression may indicate a disrupted stress response system or neuroendocrine dysfunction in Alzheimer's.
  
KLF15 (Right Plot)
- Higher in AD than in Old
- AD samples show a big jump in expression
- Also significant → Supports that KLF15 is upregulated in AD
- May be activated in Alzheimer’s brain
- KLF15 is a transcription factor involved in metabolic regulation, circadian rhythm, and cellular stress response. Although its role in Alzheimer’s is less understood, its elevated expression could suggest a compensatory or pathological activation in diseased brain tissue.


<img width="936" height="557" alt="image" src="https://github.com/user-attachments/assets/d1bc5d39-6913-4ff9-ac2e-91a0c59b3655" />

Using the Wilcoxon rank-sum test (Mann–Whitney U test) we got the above result.
Both genes showed statistically significant differences in expression between AD and Old samples

Gene CRH
- Expression is lower in AD compared to Old
- p < 0.0001 (very significant)
- CRH is downregulated in Alzheimer’s disease.
- CRH controls stress hormone release → its reduction may reflect neuroendocrine dysfunction in AD.

Gene KLF15
- Expression is higher in AD than Old
- p < 0.0001
- KLF15 is upregulated in Alzheimer’s.
- It’s involved in metabolism and transcription regulation → possibly stress response or disease compensation.

### 🧾Result 
In this study, i analyzed gene expression differences between Alzheimer's Disease (AD), Old (non-AD aged), and Young brain tissue samples using DESeq2 on publicly available RNA-seq data.
Differential Expression Analysis--
Out of 27,135 genes, several showed significant differences between AD and Old samples. The top five genes ranked by p-value were
1.	CRH– Downregulated in AD
2.	KLF15– Upregulated in AD
3.	LOC101926975– Downregulated in AD
4.	C2orf16– Downregulated in AD
5.	DUSP6– Slightly downregulated in AD and Young compared to Old

CRH and KLF15 were selected for further investigation due to their consistent and biologically relevant expression patterns.
Differential expression between conditions was assessed using the Wilcoxon rank-sum test (Mann–Whitney U test), as implemented in the ggpubr package. Significance was based on p-values.
- CRH- Significantly lower in AD (p < 0.0001)
- KLF15- Significantly higher in AD (p < 0.0001)

Functional Interpretation--
- CRH (Corticotropin-Releasing Hormone) is involved in the stress response and is normally expressed in the brain. Its downregulation may reflect impaired neuroendocrine signaling or neuronal loss in AD.
- KLF15 (Kruppel-Like Factor 15) is a transcription factor linked to metabolism, inflammation, and circadian rhythm regulation. Its upregulation in AD suggests a potential compensatory metabolic response in diseased neurons.

Pathway Enrichment Analysis using Enrichr, the two genes were found to be associated with--
- Stress hormone signaling (CRH via NR3C1)
- Transcriptional regulation and metabolism (KLF15, co-regulated with KLF9, PPARA)
These findings suggest a possible shift in stress regulation and metabolic transcription in the Alzheimer's brain.

### Conclusion
CRH and KLF15 are significantly altered in Alzheimer’s Disease and may represent key molecular players in neurodegeneration. 
 

