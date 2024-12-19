# Gut Microbiome Population Identification  

This project focuses on distinguishing rural and urban populations in Mexico based on gut microbiome samples. By analyzing microbiome data at the taxonomic and functional levels, the aim is to identify unique characteristics of each population.  

## Project Overview  
- *Objective:* To characterize individuals into rural or urban populations based on their gut microbiome's taxonomic composition and functional pathways.  
- *Data Sources:* Microbiome samples processed through bioinformatics pipelines.  
- *Applications:* Understanding the relationship between microbiome composition and environmental or lifestyle factors, aiding in population-specific health research.  

---

## Data Overview  
### Dimensions of Processed Data  
This dataset is part of a project conducted under the **Latin Biota Initiative**, which aims to investigate the gut microbiome of diverse populations in Latin America. The analyzed data provides insights into both taxonomic diversity and functional capabilities.


- *Samples:* Total number of samples analyzed (208 samples; 123 from rural and 85 from urban populations).  
- *Features:*  
  - **Taxonomic Features:** The dataset includes **13,529 OTUs** (Operational Taxonomic Units) identified through the Kraken and Bracken tools. These OTUs represent different microbial species and are used to characterize the microbiome diversity of the populations. [Clic to see the data](datasets/bracken_taxonomy.csv)


  - **Functional Pathways:** Functional pathways were derived using the HUMAnN 3.0 pipeline, resulting in **532 functional pathways**. These pathways highlight the biochemical processes present in the gut microbiome of the analyzed samples. 
  The HUMAnN 3.0 pipeline calculates pathway abundances based on the abundances of component reactions, with each reaction's abundance determined by summing the abundances of genes catalyzing the reaction. Pathway abundance is computed both at the community level and for individual species (plus an "unclassified" stratum), leveraging gene abundances and pathway structures. [Clic to see the data](datasets/latinbiota_pathabundance_unstratified.tsv)




### Problem Task  
- *Task Type:* Clusters
- *Goal:* Form and characterize meaninful clusters    

---

## Data Preprocessing  
### Workflow  
1. *Data Cleaning:*  
   - Software: [Kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/)  
   - Steps: Removal of host contamination and low-quality reads.  

2. *Taxonomic Assignment:*  
   - Software: [Kraken2](https://ccb.jhu.edu/software/kraken2/) and [Bracken](https://ccb.jhu.edu/software/bracken/).  
   - Output: OTUs assigned to taxonomic categories.  

3. *Functional Pathway Analysis:*  
   - Software: [HUMAnN 3.0](https://huttenhower.sph.harvard.edu/humann/).  
   - Output: Abundances of functional pathways in the microbiome.  

---

## Analytical Objective  
The key focus is to utilize machine learning techniques to:  
1. *Identify:* Which taxonomic and/or functional features are most important for the clusters separation
2. *Characterize:* Unique microbiome traits associated with rural and urban populations.  

---


For Elvia and Isaac

Asignación taxonómica Kracken/Bracken [Aquí](datasets/bracken_taxonomy.csv)

Correspondencia de las especies a los OTUs [Aquí](datasets/taxonomy_table_otus.csv)

Metadatos de las muestras [Aquí](datasets/metadata_LATINBIOTA_MEXICO.xlsx)