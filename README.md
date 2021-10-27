# CROPSR: An Automated Platform for Complex Genome-Wide CRISPR gRNA Design and Validation


## About CROPSR

**This repository is open sourced as specified in the LICENSE file. It is Apache License 2.0. For additional information, please check [LICENSE](LICENSE).**  
CROPSR is a python tool designed for genome-wide gRNA design and evaluation for CRISPR experiments, with special focus on complex genomes such as those found in energy-producing crops. CROPSR is a product of the DOE Center for Advanced Bioenergy and Bioproducts Innovation (CABBI).

### Citation
Please cite the following when utilizing CROPSR:

> H. Müller Paul, D. D. Istanto, J. Heldenbrand, and M. Hudson, "CROPSR: An automated platform for complex genome-wide CRISPR gRNA design and validation", BMC Bioinformatics, Oct 2021. [Online] https://doi.org/10.21203/rs.3.rs-927816/v1.


---
<details>
 <summary> Table of Contents </summary>

 [About](#about)
 - [Citation](#citation)

 [Prerequisites](#prerequisites)
 - [Dependencies](#dependencies)
 - [Installing dependencies](#installing-dependencies) 
 
 [Getting Started](#getting-started)
 - [Input Files](#input-files)
 - [First Steps](#first-steps)
 - [Output Files](#output-files) 
 
 [Example Data and Output](#example-data-and-output)
 - [Tutorial](#tutorial)
 
 <!-- [Contributing](#contributing) -->

 [Disclosures](#disclosures)
</details>

---
## Prerequisites

CROPSR does not require a separate Python environment for dependency management, as there are few dependencies and the code is updated to their current versions. Any required changes will be made to maintain compatibility with current versions of these dependencies. CROPSR is intended to be used on Python 3.7 or newer (newest stable recommended).

### **Dependencies:**

- Python version: *3.7 or newer*
- Python libraries: 
    - pandas
    - argparse

### Installing dependencies:

- If you have both Python 2.7 and Python 3 installed (e.g. Ubuntu 18.10 or older, MacOS Catalina or older), use pip3 to install the libraries to the correct path for Python 3:  
```bash
$ pip3 install pandas argparse
```
- If you have only Python 3 installed, or Python 3 is your default (e.g. Ubuntu 20.04 or newer, MacOS Big Sur or newer), the default instalation with :  
```bash
$ pip install pandas argparse 
```

#### [Back to top](#)
---
## Geting Started

### Input files

To perform a full genome analysis, the following two files are required:
- Fasta file containing whole genome sequence
- GFF file containing genome functional annotation
> **Imoprtant note for genomes downloaded from [Phytozome](https://phytozome-next.jgi.doe.gov)**:
> 
> If your genome comes from Phytozome, please make sure to also download the *annotation_info.txt* file.
> This is important because Phytozome GFF files contain a reference for a location within a separate file, so their genome browser can display the functional annotation as a separate layer. Without this file, CROPSR will output a database with no functional annotation when processing a Phytozome genome.

### First steps

CROPSR was developed as a CLI software, and requires a basic understanding of bash (or equivalent).

1. Download the project folder.

    ```bash
    $ git clone git@github.com:cabbi-bio/CROPSR.git
    ```

2. Navigate to the project folder.

    ```bash
    $ cd CROPSR
    ```

3. Place all input files on a folder (This **does not** need to be the CROPSR folder, as long as all input files are in the same folder. ***This is especially important for genomes downloaded from Phytozome***).

4. Running CROPSR (if Python 3 is set as your default `python` path, replace where it says `python3`):
    You can get the CROPSR help prompt by entering `python3 CROPSR.py -h` will return the list of arguments below:

    ```
    $ python3 CROPSR.py -h
    ```
    ```
    usage: CROPSR_6.py [-h] -f F -g G -p [-o O] [-l] [-L] [--cas9] [-v]

    optional arguments:
      -h, --help        show this help message and exit
      -f F, --fasta F   [required] path to input file in FASTA format
      -g G, --gff G     [required] path to input file in GFF format
      -p, --phytozome   path to input annotation info file in TXT format, default = None
      -o O, --output O  path to output file, default = data.csv
      -l , --length     length of the gRNA sequence, default = 20
      -L , --flanking   length of flanking region for verification, default = 200
      --cas9            specifies that design will be made for the Cas9 CRISPR system
      -v, --verbose     prints visual indicators for each iteration
    ```

    #### CROPSR arguments
    | Flag | Description |
    |---|---|
    | -h, --help | Quits the code and opens the help prompt |
    | -f, --fasta | Path to the FASTA file (`*.fasta`, `*.fa`) containing the genome sequence (always required) |
    | -g, --gff    | Path to the GFF file (`*.gff`, `*.gff3`) containing the functional annotation (always required) |
    | -p, --phytozome | Path to the `annotation_info.txt` file containing functional annotation (required for phytozome genomes) |
    | -o, --output | Path to save the output database file, including file name. The default is `data.csv` at the working directory |
    | -l, --length    | Desired length of the gRNA sequence. The default value is `20`, and this should not be changed unless required by a non-standard Cas protein. **Changing this value otherwise may cause the experiment to fail** |
    | -L, --flanking  | Desired length of flanking region for designing primers for PCR validation. The default value is `200` bases upstream and downstream of the cutsite |
    | --cas9           | Type of CRISPR system for the experiment. **At least one CRISPR system is required**. (Currently only `Cas9` is available, but other systems may be implemented in a future version) |
    | -v, --verbose    | Enables verbose mode and prints notes at several points of the process. Enable for debugging | 


### Output files

After completion, if `verbose` is enabled, a prompt will appear to inform the user that `The output file has been generated at example_data/cropsr_output.csv`. No temporary files are generated during the analysis.

CROPSR outputs a CSV (comma separated values) file by default. This file type was chosen due to the ease of handling, including importing it into the database manager of the user's preference. An option to output as a JSON following MongoDB formatting is also provided, requiring the `pymongo` library as an additional dependency.

<details>
<summary>Click here for MongoDB dependency instructions</summary>

  - If you have both Python 2.7 and Python 3 installed (e.g. Ubuntu 18.10 or older, MacOS Catalina or older), use pip3 to install the libraries to the correct path for Python 3:  
 ```bash
 $ pip3 install pymongo
 ```
 - If you have only Python 3 installed, or Python 3 is your default (e.g. Ubuntu 20.04 or newer, MacOS Big Sur or  newer), the default instalation with :  
 ```bash
 $ pip install pymongo
 ```

</details>

#### [Back to top](#)
---


## Example data and output

An example data set is provided to serve as a tutorial. This data set is comprised of the first chromosome of *Arabidopsis thaliana*, and includes both a `FASTA` and `GFF` files. The entire process will be described below.

The example data set is provided in a folder named [sample_data](sample_data) within the contents of this Git repository.
To follow along with this tutorial, no additional data should be required.

The data structure of the repository is represented below:
<!-- PLEASE CHECK THE FOLLOWING BEFORE UPLOADING -->
```
+-- README.md
+-- LICENCE.md
+-- CROPSR.py
+-- cropsr_functions.py
+-- prmrdsgn2.py
+-- sample_data/
    +-- sample_genome.fa
    +-- sample_genome.gff
+-- .DS_Store
+-- .gitignore
```

<details>
  <summary>Click here to preview FASTA file sample_genome.fa</summary>

  ```
  >NC_003070.9 Arabidopsis thaliana chromosome 1 sequence
  ccctaaaccctaaaccctaaaccctaaacctctGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATG
  AATCCCTAAATACCTAAttccctaaacccgaaaccggTTTCTCTGGTTGAAAATCATTGTGtatataatgataattttat
  CGTTTTTATGTAATTGCTTATTGTTGTGtgtagattttttaaaaatatcatttgagGTCAATACAAATCCTATTTCTTGT
  GGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATT
  TGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTATGTTT
  GGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGgaaaattatttagttg
  taGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAGCATTTAT
  TCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAAAAAAGCTA
  TCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACtcaaaaaagtatttttagatgtttgttttgc
  ttctttgaagTAGTTTCTCTTTGCAAAATTCCTCTTTTTTTAGAGTGATTTGGATGATTCAAGACTTCTCGGTACTGCAA
  AGTTCTTCCGCCTGATTAATTATCCATTTTACCTTTGTCGTAGATATTAGGTAATCTGTAAGTCAACTCATATACAActc
  ataatttaaaataaaattatgatcGACACACGTTTACACATAAAATCTGTAAATCAACTCATATACCCGTTATTCCCACA
  ATCATATGCTTTCTAAAAGCAAAAGTATATGTCAACAATTGGTTATAAATTATTAGAAGTTTTCCACTTATGACTTAAGA
  ACTTGTGAAGCAGAAAGTGGCAACACCCCCCAcctcccccccccccccccacccCCCAAATTGAGAAgtcaattttatat
  aatttaatcaaataaataagtttatgGTTAAGAGTTTTTtactctctttatttttctttttctttttgagacATACTGAA
  AAAAGTTGTAATTATTAATGATAGTTCTGTGATTCCTCCATGAATCACATCtgcttgatttttctttcataaatttataa
  gtaATACATTCTTATAAAATGGTCAGAGAAACACCAAAGATCCCGAGATTTCTTCTCacttactttttttctatctatCT
  AGATTATATAAATGAGATGTTGAATTAGAGGAACCTTTGATTCAATGATCATAGAAAAATTAGGTAAAGAGTCAGTGTCG
  TTATGTTATGGAAGATGTGAATGAAGTTTGACTTCTCATTGTATATGAGTAAAATCTTTTCTTACAAGGGAAGTCCCCAA
  ...
  TAGAAAGCGACCTCATAGATGGTAGGAATAGTGTGAAGAACGTTAGATTCTGTGACTAATTTCTTGTGCAGACTATAATA
  ACACTGGGTAGTAATACttatatcttcatcatcttaaACATGTAGTAATCCCATTACATCAGCTTAAATCGGTTAACCCT
  TCTACAAACACAGTCAATCCTGCAGAAAGGTACATCCAGATAATCTCAGTCGACGACCATGAGTTTTGGTTCATGTGTTT
  CTTAAACTACTGTGTGATATAATATTGGTGGATTCTTGATTAGAATAACATTATTTTGGTTGCaaactaaattatttgat
  tcttattcTTAATTAGTTACCATGTcttgattagggtttagggtttagggtttagggtttagggtttaggg
  ```

</details>

<details>
  <summary>Click here to preview GFF file sample_genome.gff</summary>

  ```
  ##gff-version 3
  #!gff-spec-version 1.21
  #!processor NCBI annotwriter
  #!genome-build TAIR10.1
  #!genome-build-accession NCBI_Assembly:GCF_000001735.4
  #!annotation-source TAIR and Araport 
  ##sequence-region NC_003070.9 1 30427671
  ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3702
  NC_003070.9     RefSeq  region  1       30427671        .       +       .       ID=NC_003070.9:1..30427671;Dbxref=taxon:3702;Name=1;chromosome=1;ecotype=Columbia;gbkey=Src;genome=chromosome;mol_type=genomic DNA
  NC_003070.9     RefSeq  gene    3631    5899    .       +       .       ID=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC domain containing protein 1,T25K16.1,T25K16_1;locus_tag=AT1G01010
  NC_003070.9     RefSeq  mRNA    3631    5899    .       +       .       ID=rna-NM_099983.2;Parent=gene-AT1G01010;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580,Genbank:NM_099983.2;Name=NM_099983.2;gbkey=mRNA;gene=NAC001;inference=similar to RNA sequence%2C mRNA:INSD:BT001115.1%2CINSD:AF439834.1%2CINSD:AK226863.1;locus_tag=AT1G01010;orig_protein_id=gnl|JCVI|AT1G01010.1;orig_transcript_id=gnl|JCVI|mRNA.AT1G01010.1;product=NAC domain containing protein 1;transcript_id=NM_099983.2
  NC_003070.9     RefSeq  exon    3631    3913    .       +       .       ID=exon-NM_099983.2-1;Parent=rna-NM_099983.2;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580,Genbank:NM_099983.2;gbkey=mRNA;gene=NAC001;inference=similar to RNA sequence%2C mRNA:INSD:BT001115.1%2CINSD:AF439834.1%2CINSD:AK226863.1;locus_tag=AT1G01010;orig_protein_id=gnl|JCVI|AT1G01010.1;orig_transcript_id=gnl|JCVI|mRNA.AT1G01010.1;product=NAC domain containing protein 1;transcript_id=NM_099983.2
  ...
  NC_003074.8     RefSeq  CDS     23459372        23459803        .       -       0       ID=cds-NP_001190172.1;Parent=rna-NM_001203243.2;Dbxref=TAIR:AT3G63540,GeneID:7922413,Genbank:NP_001190172.1,Araport:AT3G63540;Name=NP_001190172.1;Note=Mog1/PsbP/DUF1795-like photosystem II reaction center PsbP family protein%3B FUNCTIONS IN: calcium ion binding%3B INVOLVED IN: photosynthesis%3B LOCATED IN: oxygen evolving complex%2C extrinsic to membrane%3B CONTAINS InterPro DOMAIN/s: Photosystem II oxygen evolving complex protein PsbP (InterPro:IPR002683)%2C Mog1/PsbP/DUF1795%2C alpha/beta/alpha sandwich (InterPro:IPR016124)%2C Mog1/PsbP%2C alpha/beta/alpha sandwich (InterPro:IPR016123).;end_range=23459803,.;gbkey=CDS;inference=Similar to RNA sequence%2C EST:INSD:EH981858.1%2CINSD:AA728452.1%2CINSD:EG467302.1%2C INSD:BP618169.1%2CINSD:EL145784.1%2CINSD:EG443892.1%2C INSD:EG467305.1%2CINSD:EL079202.1%2CINSD:EG467300.1%2C INSD:EH839793.1%2CINSD:EG499296.1%2CINSD:BP607626.1%2C INSD:EG443886.1%2CINSD:EL285102.1%2CINSD:EG443894.1%2C INSD:BP789147.1%2CINSD:EL020616.1%2CINSD:EG467297.1%2C INSD:BP850338.1%2CINSD:BP828641.1%2CINSD:EL159236.1%2C INSD:EH950133.1%2CINSD:H76109.1%2CINSD:H76121.1%2C INSD:EL147984.1%2CINSD:BP561663.2%2CINSD:EG443885.1%2C INSD:EG467309.1%2CINSD:BP829784.1%2CINSD:EL123641.1%2C INSD:CK117773.1%2CINSD:H76108.1%2CINSD:EG467303.1%2C INSD:EH845308.1%2CINSD:H76445.1%2CINSD:EG499250.1%2C INSD:EH906185.1%2CINSD:EG499263.1%2CINSD:EG443896.1%2C INSD:DR368907.1%2CINSD:BP655570.1%2CINSD:EG443893.1%2C INSD:EH867494.1%2CINSD:AV796368.1%2CINSD:EG443888.1%2C INSD:EG443891.1%2CINSD:EH809220.1%2CINSD:ES174073.1%2C INSD:EL195858.1%2CINSD:EG443890.1%2CINSD:EG443960.1%2C INSD:EL170398.1%2CINSD:CB263188.1%2CINSD:EG443959.1%2C INSD:EL269452.1%2CINSD:BP831081.1%2CINSD:AA721800.1%2C INSD:EG499151.1%2CINSD:EL055832.1%2CINSD:BP625539.1%2C INSD:EH882173.1%2CINSD:EL268883.1%2CINSD:W43032.1%2C INSD:EL059043.1%2CINSD:BP651248.1%2CINSD:EL008352.1%2C INSD:EL063436.1%2CINSD:EG467290.1%2CINSD:ES094569.1%2C INSD:EG467311.1%2CINSD:BP651706.1%2CINSD:EG499307.1%2C INSD:EG443889.1%2CINSD:EG420141.1%2CINSD:BP623147.1%2C INSD:BP822126.1%2CINSD:EH805894.1%2CINSD:EH856030.1%2C INSD:BP808479.1%2CINSD:CK120153.1%2CINSD:H76442.1%2C INSD:EG467295.1%2CINSD:EL304565.1%2CINSD:ES181862.1%2C INSD:EG467281.1%2CINSD:EH930027.1%2CINSD:EL152252.1%2C INSD:EH994289.1%2CINSD:EL126605.1%2CINSD:EL249438.1%2C INSD:ES025201.1%2CINSD:AI099598.1%2CINSD:H76097.1%2C INSD:EH994460.1%2CINSD:CK120738.1%2CINSD:EG499128.1%2C INSD:H76446.1%2CINSD:EH840708.1%2CINSD:EG499274.1%2C INSD:H76100.1%2CINSD:EG467308.1%2CINSD:EL120114.1%2C INSD:H76094.1%2CINSD:EG443895.1%2CINSD:ES164938.1%2C INSD:EG467304.1%2CINSD:EG467299.1%2CINSD:EG499507.1%2C INSD:EL070078.1%2CINSD:EH883854.1%2CINSD:H76113.1%2C INSD:EH826245.1%2CINSD:EH956928.1%2CINSD:H76114.1%2C INSD:EL318399.1%2CINSD:H76101.1%2CINSD:EH903234.1%2C INSD:EG467282.1%2CINSD:EG420153.1%2CINSD:EL026196.1;locus_tag=AT3G63540;orig_transcript_id=gnl|JCVI|mRNA.AT3G63540.1;partial=true;product=thylakoid lumenal protein (Mog1/PsbP/DUF1795-like photosystem II reaction center PsbP family protein);protein_id=NP_001190172.1
  ```

</details>

### Tutorial

1. Navigate to the project folder.

    ```bash
    $ cd CROPSR
    ```

2. Run CROPSR with the sample data (code is available below).
    ```bash
    $ python3 CROPSR.py -f sample_data/sample_genome.fa -g sample_data/sample_genome.gff3 -o sample_data/sample_genome_output.csv --cas9 -v
    ```
    Note that the `verbose` flag was left on. This will cause the terminal to print notifications during the process, however, it means you will not be able to utilize the terminal window until it is finished. If you close the terminal window while the process is running, it will cause an interruption.
    
    <details>
      <summary>Click here to learn how to run this process in the background</summary>

      ```bash
      $ python3 CROPSR.py -f sample_data/sample_genome.fa -g sample_data/sample_genome.gff3 -o sample_data/sample_genome_output.csv --cas9 &
      ```
      In this variation, the `verbose` flag was removed, and a `&` was added at the end of the command. This will free your terminal window to perform other tasks or be closed. **This is the recommended approach when running real data, as the process may take more than a day to finish**. Make sure the computer remains powered on for the entirety of the process.
    </details>

3. While CROPSR is running (with the `verbose` flag active), you should see the following appear in your terminal:

    ```
    ################################################################################
    ##                                                                            ##
    ##                                                                            ##
    ##          .o88b.   d8888b.    .d88b.    d8888b.   .d8888.   d8888b.         ##
    ##         d8P  Y8   88  `8D   .8P  Y8.   88  `8D   88'  YP   88  `8D         ##
    ##         8P        88oobY'   88    88   88oodD'   `8bo.     88oobY'         ##
    ##         8b        88`8b     88    88   88ººº       `Y8b.   88`8b           ##
    ##         Y8b  d8   88 `88.   `8b  d8'   88        db   8D   88 `88.         ##
    ##          `Y88P'   88   YD    `Y88P'    88        `8888Y'   88   YD         ##
    ##                                                                            ##
    ##                                                                            ##
    ################################################################################
    U.S. Dept. of Energy's Center for Advanced Bioenergy and Bioproducts Innovation
    University of Illinois at Urbana-Champaign

            You are currently utilizing the following settings:

            CROPSR version:                                 1.11b
            Path to genome file in FASTA format:            /sample_data/sample_genome.fa
            Path to output file:                            /sample_data/sample_genome_output.csv
            Length of the gRNA sequence:                    20
            Length of flanking region for verification:     200
            Number of available CPUs:                       12
            Path to annotation file in GFF format:          /sample_data/sample_genome.gff3
            Path to annotation_info file in TXT format:     None
            Designing for CRISPR system:
                Streptococcus pyogenes Cas9                 True
    ```

4. After the process is complete, you should have access to the generated output file in `CSV` format.

    The folder structure should be similar to what is represented below:
    <!-- PLEASE CHECK THE FOLLOWING BEFORE UPLOADING -->
    ```
    +-- README.md
    +-- LICENCE.md
    +-- CROPSR.py
    +-- cropsr_functions.py
    +-- prmrdsgn2.py
    +-- sample_data/
        +-- sample_genome.fa
        +-- sample_genome.gff
        +-- sample_genome_output.csv
    +-- .DS_Store
    +-- .gitignore
    ```

    <!-- <details>
      <summary>Click here for a preview of sample_genome_output.csv</summary>
      
      CSV FILE OUTPUT
      
    </details> -->

#### [Back to top](#)
---

<!-- ## Contributing

#### [Back to top](#)
--- -->

## Disclosures

>This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research under Award Number DE-SC0018420). Any opinions, findings, and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy.

#### [Back to top](#)
---
