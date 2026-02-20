# lncRNA_TF_pairs_analysis  

Bioinformatics analysis of regulatory relationships between long non-coding RNA genes and transcription factor coding genes.  
  
## Software Dependencies  

- R 4.5.2  
- Bioconductor packages: TissuEnrich, CAGEr  
- Additional dependencies are listed in `environment.yml`  
  
## Installation  
  
1. **Clone the repository:**  
    ```bash  
    git clone git@github.com:MartinLoza/lncRNA_TF_pairs_analysis.git
    ```  
2. **Create and activate the conda environment:**  
    ```bash  
    conda env create -f environment.yml  
    conda activate lncRNA_TF_pairs_analysis  
    ```  
  
## Usage Instructions  
  
The analysis is organized in folders and files to be run incrementally in numerical order:  
  
| Step | Folder/File                                 | Description                                    |  
|------|---------------------------------------------|------------------------------------------------|  
| 00   | 00_Setup_ENSEML_transcripts                 | Setup of ENSEMBL transcript data               |  
| 01   | 01_Get_PCG_distance_to_ncRNA.ipynb          | Calculate distance from PCGs to ncRNA          |  
| 02   | 02_Figure_1_comparison_species              | Generate Figure 1: Species comparison          |  
| 03   | 03_Annotate_TFs.ipynb                       | Annotate transcription factors                 |  
| 04   | 04_Figure_2_window_size_and_GO              | Generate Figure 2: Window size and GO analysis |  
| 05   | 05_Figure_3_Species_comparison_lncRNA_TF    | Generate Figure 3: lncRNA-TF species comparison|  
| 06   | 06_Setup_lncRNA_TF_pairs.ipynb              | Setup lncRNA-TF pairs                          |  
| 07   | 07_ENCODE_analyses                          | ENCODE dataset analyses                        |  
| 08   | 08_GTEx_analyses                            | GTEx dataset analyses                          |  
| 09   | 09_Figure_4_Correlation_ENCODE_GTEx         | Generate Figure 4: Correlation ENCODE/GTEx     |  
| 10   | 10_FANTOM5_time_course_analyses             | FANTOM5 time course analyses                   |  
| 11   | 11_Figure_5_FANTOM5                         | Generate Figure 5: FANTOM5                     |  


**Note:**    
Files inside each folder also follow the same numerical ordering.  
  
## Citation  
  
If you use this repository, please cite:    
[to be updated]  