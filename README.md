# **PROSAVA Tutorial**

## _Description_

PROSAVA is the abbreviation of **PRO**tein **S**equence **A**nalysis and **V**isualization **A**pplication. It is a standalone application for domain-centric protein analysis in prokaryotes, which can automatically dissect protein sequences into domains, and then used the domains for comparing and visualizing homologous sequences in the form of domain architectures. In addition, PROSAVA could also visualize prokaryotic proteomes in terms of life complexities and lifestyles from domain-centric view through identification of non-redundant domains and quantification of domain-domain interactions.

## _Feedback_

Please report any bugs or suggest enhancements at _healthscience@foxmail.com_

## _OS Requirements_

64-bit Linux or MacOS

## _License_

PROSAVA is a free software

## _System Configuration_

The development of PROSAVA mainly depends on the Anaconda environment. Anaconda integrates many software packages related to data analysis and scientific computing, which comes with the Conda package management system specifically designed to solve environmental dependency problems. Therefore, before using PROSAVA, please make sure that Anaconda environment is successfully installed on your computer. Anaconda installation program could be easily downloaded from the official website _https://www.anaconda.com_. It is worth noting that Anaconda with Python 3 version is required. After configuring the Anaconda environment, please use the pip command to install the **_ete3_** module on the terminal. Installation command is: **_pip install ete3_**. After completing the installation, the operating environment of PROSAVA has been successfully configured.

## _PROSAVA Usage_

1.Open the terminal, set the working path to where PROSAVA is located.

2.Type the following command in the terminal, set the corresponding parameter and the program will run automatically.

```
python prosava.py -E (1e-3 or 1e-10) -T (sequence or proteome) -I (path of input file) -M (path of Pfam model)
```

EXAMPLE

```
python prosava.py -E 1e-3 -T sequence -I sample.fasta -M Pfam-A.hmm-version33.1/Pfam-A.hmm
```

3.PROSAVA can classify homologous sequences according to domain compositions for homologous sequences, conduct corresponding phylogenetic analysis on specific categories, and compare the changes of domain compositions between different categories. For proteomes of prokaryotes, PROSAVA can quantify the complexity of life and analyze the associated network of co-occurrence of corresponding domains.

4.After the program runs, a **_domain_dissection_** folder is automatically generated under the PROSAVA directory, which stores domain dissection results after the automatic analysis of the input sequences. After all the analysis is completed, a **_result_** folder will be generated under the PROSAVA directory, and the corresponding result file will be stored. The specific content will be described in detail below. In particular, when the input data are homologous sequences, a folder named **_phylotree_** is additionally generated in the directory, and related intermediate files of phylogenetic analysis are stored therein.

5.When performing homologous sequence analysis, a total of 4 result files are generated under the result folder. The first is **_seqAnalysis.csv_**, which stores the domain analysis results of each sequence in the homologous sequence and gives the E-value corresponding to each domain. **_domainPic.pdf_** is the visualization result of the classification of the submitted homologous sequences organized by domains. **_domainDistribution.csv_** is the text information corresponding to the visualization result. **_phylogeneticTree.pdf_** is the visualization result of the phylogenetic analysis of the specific domain organization. The corresponding intermediate files during the execution of the phylogenetic analysis are stored in the phylotree folder under the PROSAVA directory. There are three subfolders and three files under the phylotree folder. The **_1_sequence_** folder contains the set of original sequences corresponding to each type of domain distribution. **_2_alignment_** folder contains the results of each kind of original sequence alignment. **_3_hmm_model_** folder contains corresponding hidden Markov model files generated according to sequence comparison. **_sequence.fasta_** is a sequence set that each type of sequence randomly generates according to the constructed HMM model. **_sequence.fasta.aln_** is the sequence alignment result of the sequence set, and **_sequence.fasta.phy_** is the phylogenetic result generated based on these alignment results.

6.When performing proteomic proteome analysis, a total of 6 result files are generated under the result folder. The first is **_Domain_Frequence.csv_**, which contains the frequency of occurrence of all domains in the proteome. The remaining five files are related data and visualization results in structural domain co-occurrence related network analysis. **_Domain_interaction.csv_** contains the co-occurrence related information of the structural domain and the corresponding frequency. **_Non-redundant_Domain.csv_** contains the frequency of occurrence of each node or structure domain in the network, and the corresponding maximum and minimum E-value. Domain_interaction.csv and Non-redundant_Domain.csv can be used as inputs into CytoScape as needed to visualize the cooccurrence network. The following **_domainInteract.jpg_** is a visualization of the highly connected part of the co-occurring network. **_Hub_Nodes.csv_** is the relevant information of the core structure domain in the visualization results, including the name of the core structure domain, Pfam ID, and the number of co-occurrence associations with other structure domains. The last **_networkAppendix.txt_** contains the total number of domains in the proteome, the total number of co-occurrences of all domains in the proteome, the number of non-redundant domains, the number of nodes and edges in the maximally connected network, and all Information related to the core structure domain.
