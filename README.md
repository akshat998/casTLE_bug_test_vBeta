# casTLE: Systematic comparison of CRISPR/Cas9 and RNAi screens for essential genes

Reference: Morgens, David W., et al. "Systematic comparison of CRISPR/Cas9 and RNAi screens for essential genes." Nature biotechnology 34.6 (2016): 634-636.



## Installation: 
- Please make sure anaconda is installed on your machine. 
- Please run the following in order, in the command line: 
   ```
   conda create -name py2 python=2.7
   conda activate py2
   conda install -c bioconda bowtie=1.2.2
   conda install -c bioconda htseq
   conda install -c anaconda scipy=0.19.1
   ```  
- If all the installations went smoothly, you are ready to run the commands from 'Quick Start'. 

## File Descriptions:
- Scripts/makeIndices.py: Creation of a mapping for all the guides, used by future steps for association to a gene. Not requred if genome-wide alignment is prefered. 
- Scripts/makeCounts.py: Align FASTQ and make count file
- Scripts/analyzeCounts.py: Compares count files using casTLE
- Scripts/analyzeCombo.py: Combines data for two screens


## Quick Start Example: 
- Please activate your conda environment with: 
   ```
   conda activate py2
   ```  
- Please cd inside the downloaded casTLE directory.  
- Inside casTLE, please unzip the directory 'GenRef.zip'
- Creation of an index file: <br />
  Step 1: Please save the index file (in csv format) in the directory 'Align'. <br />
  Step 2: Please run the following command in the command line: <br />
  ```
  python ./Scripts/makeIndices.py -o ./Align/attardi1_align.csv attardi1_align_short attardi1_align_long -t
  ```  
- Run makeCounts.py: <br />
  ```
  python Scripts/makeCounts.py ./Sequencing_Data/220318_RVJS_NucKO_KS_colab/left_flank_S14_L00  Results/left_flank attardi1_align_short   
  python Scripts/makeCounts.py ./Sequencing_Data/220318_RVJS_NucKO_KS_colab/right_flank_S15_L00  Results/right_flank attardi1_align_short 
  python Scripts/makeCounts.py ./Sequencing_Data/220318_RVJS_NucKO_KS_colab/t0_S13_L00  Results/t0 attardi1_align_short 
  python Scripts/makeCounts.py ./Sequencing_Data/220318_RVJS_NucKO_KS_colab/lib_S16_L00  Results/lib attardi1_align_short 
  ``` 
- Run analyzeCounts.py: (Output will be saved within the directory 'Results') <br />
  ``` 
  python Scripts/analyzeCounts.py ./Data/Results/t0_counts.csv ./Data/Results/right_flank_counts.csv t0_v_right
  python Scripts/analyzeCounts.py ./Data/Results/t0_counts.csv ./Data/Results/left_flank_counts.csv t0_v_left
  python Scripts/analyzeCounts.py ./Data/Results/lib_counts.csv ./Data/Results/t0_counts.csv lib_v_t0
  ``` 
- Run analyzeCombo: (Output will be saved within the directory 'Results') <br />
  ```
  python Scripts/analyzeCombo.py ./Results/t0_v_right.csv  ./Results/t0_v_left.csv combo_left_right
  ```     
- Plot the guides for a particular gene: (Output will be saved within the directory 'Results') <br />
  ```
  python Scripts/plotGenes_RL.py /Results/t0_v_left.csv SPTAN1
  ```     
- Run analyzeCombo.py: 
  ```
  python Scripts/addCombo.py Results/combo_left_right.csv 1000000
  ```     


## Questions, problems?
Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (TODO: Email)


