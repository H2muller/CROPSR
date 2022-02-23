# CROPSR: An Automated Platform for Complex Genome-Wide CRISPR gRNA Design and Validation


## About CROPSR

**This repository is open sourced as specified in the LICENSE file. It is Apache License 2.0. For additional information, please check [LICENSE](LICENSE).**  
CROPSR is a python tool designed for genome-wide gRNA design and evaluation for CRISPR experiments, with special focus on complex genomes such as those found in energy-producing crops. CROPSR is a product of the DOE Center for Advanced Bioenergy and Bioproducts Innovation (CABBI).

### Citation
Please cite the following when utilizing CROPSR:

> Müller Paul, H., Istanto, D.D., Heldenbrand, J. et al. CROPSR: an automated platform for complex genome-wide CRISPR gRNA design and validation. BMC Bioinformatics 23, 74 (2022). https://doi.org/10.1186/s12859-022-04593-2.


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
      -f , --fasta F   [required] path to input file in FASTA format
      -g , --gff G     path to input file in GFF format
      -p , --phytozome   path to input annotation info file in TXT format, default = None
      -o , --output O  path to output file, default = data.csv
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

An example data set is provided to serve as a tutorial. This data set is comprised of the first chromosome of *Saccharomyces cerevisiae*, and includes both a `FASTA` and `GFF` files. The entire process will be described below.

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
  >Chr01
  ccacaccacacccacacacccacacaccacaccacacaccacaccacacccacacacacacatCCTAACACTACCCTAAC
  ACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTTACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCAT
  TCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
  CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATACTGTTCTTCTACCCACCATAT
  TGAAACGCTAACAAATGATCGTAAATAACACACACGTGCTTACCCTACCACTTTATACCACCACCACATGCCATACTCAC
  CCTCACTTGTATACTGATTTTACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTC
  CACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGG
  TCTATACCCTGTGCCATTTACCCATAACGCCCATCATTATCCACATTTTGATATCTATATCTCATTCGGCGGTcccaaat
  attgtataaCTGCCCTTAATACATACGTTATACCACTTTTGCACCATATACTTACCACTCCATTTATATACACTTATGTC
  AATATTACAGAAAAATCCCCACAAAAATCacctaaacataaaaatattctacttttcaacaataataCATAAACATATTG
  GCTTGTGGTAGCAACACTATCATGGTATCACTAACGTAAAAGTTCCTCAATATTGCAATTTGCTTGAACGGATGCTATTT
  CAGAATATTTCGTACTTACACAGGCCATACATTAGAATAATATGTCACATCACTGTCGTAACACTCTTTATTCACCGAGC
  AATAATACGGTAGTGGCTCAAACTCATGCGGGTGCTATGATACAATTATATCTTATTTCCATTCCCATATGCTAACCGCA
  ATATCCTAAAAGCATAACTGATGCATCTTTAATCTTGTATGTGACACTACTCATACGAAGGGACTATATCTAGTCAAGAC
  GATACTGTGATAGGTACGTTATTTAATAGGATCTATAACGAAATgtcaaataattttacgGTAATATAACTTATCAGCGG
  CGTATACTAAAACGGACGTTACGATATTGTCTCACTTCATCTTACCACCCTCTATCTTATTGCTGATAGAACACTAACCC
  CTCAGCTTTATTTCTAGTTACAGTTACACAAAAAACTATGCCAACCCAGAAATCTTGATATTTTACGTGTCAAAAAATGA
  GGGTCTCTAAATGAGAGTTTGGTACCATGACTTGTAACTCGCACTGCCCTGATCTGCAATCTTGTTCTTAGAAGTGACGC
  ATATTCTATACGGCCCGACGCGACGCgccaaaaaatgaaaaacgAAGCAGCGactcatttttatttaagGACAAAGGTTG
  CGAAGCCGCACATTTCCAATTTCATTGTTGTTTATTGGACATACACTGTTAGCTTTATTACCGTCCACGTTTTTTCTACA
  ATAGTGTagaagtttctttcttatgTTCATCGTATTCATAAAATGCTTCACGAACACCGTCATTGATCAAATAGgtctat
  aatattaatatacatttatataaTCTACGGTATttatatcatcaaaaaaaagtagtttttttattttattttgttcgtta
  attttcaatttctatGGAAACCCGTTCGTAAAATTGGCGTTTGTCTCTAGTTTGCGATAGTGTAGATACCGTCCTTGGAT
  AGAGCACTGGAGATGGCTGGCTTTAATCTGCTGGAGTACCATGGAACACCGGTGATCATTCTGGTCACTTGGTCTGGAGC
  AATACCGGTCAACATGGTGGTGAAGTCACCGTAGTTGAAAACGGCTTCAGCAACTTCGACTGGGTAGGTTTCAGTTGGGT
  GGGCGGCTTGGAACATGTAGTATTGGGCTAAGTGAGCTCTGATATCAGAGACGTAGACACCCAATTCCACCAAGTTGACT
  CTTTCGTCAGATTGAGCTAGAGTGGTGGTTGCAGAAGCAGTAGCAGCGATGGCAGCGACACCAGCGGCGATTGAAGTTAA
  TTTGACCattgtatttgttttgtttgttaGTGCTGATATAAGCTTAACAGGAAAGGAAAGAATAAAGACATATTCTCAAA
  GGCATATAGTTGAAGCAGCTCTATTTATACCCATTCCCTCATGGGTTGTTGCTATTTAAACGATCGCTGACTGGCACCAG
  TTCCTCATCAAATATTCTCTATATCTCATCTTTCACACAATCTCATTATCTCTATGGAGATGCTCTTGTTTCTGAACGAA
  TCATAAATCTTTCATAGGTTTCGTATGTGGAGTACTGTTTTATGGCGCTTATGTGTATTCGTATGCGCAGAATGTGGGAA
  TGCCAATTATAGGGGTGCCGAGGTGCCTTATAAAACCCTTTTCTGTGCCTGTGacatttcctttttcggtcaaaaagaat
  atccGAATTTTAGATTTGGACCCTCGTACAGAAGCTTATTGTCTAAGCCTGAATTCAGTCTGCTTTAAACGGCTTCCGCG
  GAGGAAATATTTCCATCTCTTGAATTCGTACAACATTAAACGTGTGTTGGGAGTCGTATACTGTTAGGGTCTGTAAACTT
  GTGAACTCTCGGCAAATGCCTTGGTGCAATTACGTAATTTTAGCCGCTGAGAAGCGGATGGTAATGAGACAAGTTGATAT
  CAAACAGATACATATTTAAAAGAGGGTACCGCTAATTTAGCAGGGCAGTATTATTGTAGTTTGATATGTACGGCTAACTG
  AACCTAAGTAGGGATATGAGAGTAAGAACGTTCGGCTACTCTTCTTTCTAAGTGGGATTTTTCTTAATCCTTGGATTCTT
  AAAAGGTTATTAAAGTTCCGCACAAAGAACGCTTGGAAATCGCATTCATCAAAGAACAACTCTTCGTTTTCCAAACAATC
  TTCCCGAAAAAGTAGCCGTTCATTTCCCTTCCGATTTCATTCCTAGACtgccaaatttttcttgctcATTTATAATGATT
  GATAAGAATTGTATTTGTGTCCCATTCTCGTAGATAAAATTCTTGGAtgttaaaaaattattattttcttcataaagAAG
  CTTTCAAGATATAAGATACGAAATAGGGGTTGATAATTGCATGACAGTAGCTTTAgatcaaaaaggaaagcaTGGAGGGA
  AACAGTAAACAGTGAAAATTCTCTTGAGAACCAAAGTAAACCTTCATTGAAGAGCTTCcttaaaaaatttagaaTCTCCC
  ATGTCAACGGGTTTCCATACCTCCCCAGCATCATacatcttttttcaaagaaactTCAAATGCCTCTTTTATGCAAGGGG
  CAAAATCCTGAAATGACTTAAACTTAGCAGTttcgtcttttttcaaagagaatggttgaagaagaattgtttTGGACGCT
  TATTGACAATCTGTTGCATTGATAAAGTACCTACTATCCCAGACTATATTTGTATACAAGTACAAAATTAGGTTTGTTGA
  AACAACTTTCCGATCATTGGTGCCCGTATCTGATGTTTTTTTAGTAATTTCTTTGTAAATACAGGGAGTTGTTTCGAAAG
  CTTATGAGAAAAATACATGAATGACAGGTAAAAATATTGGCTCGAAAAAGAGGacaaaaagagaaatcaTAAATGAGTAA
  ACCCACTTGCTGGACATTATCCAGTAAAGGCTTGGTAGTAACCATAATATTACCCAGGTACGAAACGCTAAGAACTTGAA
  AGACTCATAAAACTTCCAGGTTAAgctatttttgaaaatattctgaGGTAAAAGCCATTAAGGTCCAGATAACCAAGGGA
  ...
  ```

</details>

<details>
  <summary>Click here to preview GFF file sample_genome.gff</summary>

  ```
  ##gff-version 3
  #!gff-spec-version 1.21
  #!processor NCBI annotwriter
  #!genome-build R64
  #!genome-build-accession NCBI_Assembly:GCF_000146045.2
  #!annotation-source SGD R64-3-1
  ##sequence-region Chr01 1 230218
  ##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=559292
  Chr01	RefSeq	region	1	230218	.	+	.	ID=Chr01:1..230218;Dbxref=taxon:559292;Name=I;chromosome=I;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=S288C
  Chr01	RefSeq	telomere	1	801	.	-	.	ID=id-Chr01:1..801;Dbxref=SGD:S000028862;Note=TEL01L%3B Telomeric region on the left arm of Chromosome I%3B composed of an X element core sequence%2C X element combinatorial repeats%2C and a short terminal stretch of telomeric repeats;gbkey=telomere
  Chr01	RefSeq	origin_of_replication	707	776	.	+	.	ID=id-Chr01:707..776;Dbxref=SGD:S000121252;Note=ARS102%3B Autonomously Replicating Sequence;gbkey=rep_origin
  Chr01	RefSeq	gene	1807	2169	.	-	.	ID=gene-YAL068C;Dbxref=GeneID:851229;Name=PAU8;end_range=2169,.;gbkey=Gene;gene=PAU8;gene_biotype=protein_coding;locus_tag=YAL068C;partial=true;start_range=.,1807
  Chr01	RefSeq	mRNA	1807	2169	.	-	.	ID=rna-NM_001180043.1;Parent=gene-YAL068C;Dbxref=GeneID:851229,Genbank:NM_001180043.1;Name=NM_001180043.1;end_range=2169,.;gbkey=mRNA;gene=PAU8;locus_tag=YAL068C;partial=true;product=seripauperin PAU8;start_range=.,1807;transcript_id=NM_001180043.1
  Chr01	RefSeq	exon	1807	2169	.	-	.	ID=exon-NM_001180043.1-1;Parent=rna-NM_001180043.1;Dbxref=GeneID:851229,Genbank:NM_001180043.1;end_range=2169,.;gbkey=mRNA;gene=PAU8;locus_tag=YAL068C;partial=true;product=seripauperin PAU8;start_range=.,1807;transcript_id=NM_001180043.1
  Chr01	RefSeq	CDS	1807	2169	.	-	0	ID=cds-NP_009332.1;Parent=rna-NM_001180043.1;Dbxref=SGD:S000002142,GeneID:851229,Genbank:NP_009332.1;Name=NP_009332.1;Note=hypothetical protein%3B member of the seripauperin multigene family encoded mainly in subtelomeric regions;experiment=EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695],EXISTENCE:mutant phenotype:GO:0045944 positive regulation of transcription by RNA polymerase II [PMID:12586695];gbkey=CDS;gene=PAU8;locus_tag=YAL068C;product=seripauperin PAU8;protein_id=NP_009332.1
  Chr01	RefSeq	gene	2480	2707	.	+	.	ID=gene-YAL067W-A;Dbxref=GeneID:1466426;Name=YAL067W-A;end_range=2707,.;gbkey=Gene;gene_biotype=protein_coding;locus_tag=YAL067W-A;partial=true;start_range=.,2480
  Chr01	RefSeq	mRNA	2480	2707	.	+	.	ID=rna-NM_001184582.1;Parent=gene-YAL067W-A;Dbxref=GeneID:1466426,Genbank:NM_001184582.1;Name=NM_001184582.1;end_range=2707,.;gbkey=mRNA;locus_tag=YAL067W-A;partial=true;product=uncharacterized protein;start_range=.,2480;transcript_id=NM_001184582.1
  Chr01	RefSeq	exon	2480	2707	.	+	.	ID=exon-NM_001184582.1-1;Parent=rna-NM_001184582.1;Dbxref=GeneID:1466426,Genbank:NM_001184582.1;end_range=2707,.;gbkey=mRNA;locus_tag=YAL067W-A;partial=true;product=uncharacterized protein;start_range=.,2480;transcript_id=NM_001184582.1
  Chr01	RefSeq	CDS	2480	2707	.	+	0	ID=cds-NP_878038.1;Parent=rna-NM_001184582.1;Dbxref=SGD:S000028593,GeneID:1466426,Genbank:NP_878038.1;Name=NP_878038.1;Note=hypothetical protein%3B identified by gene-trapping%2C microarray-based expression analysis%2C and genome-wide homology searching;gbkey=CDS;locus_tag=YAL067W-A;product=uncharacterized protein;protein_id=NP_878038.1
  Chr01	RefSeq	gene	7235	9016	.	-	.	ID=gene-YAL067C;Dbxref=GeneID:851230;Name=SEO1;end_range=9016,.;gbkey=Gene;gene=SEO1;gene_biotype=protein_coding;locus_tag=YAL067C;partial=true;start_range=.,7235
  ...
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
            Path to genome file in FASTA format:            sample_data/sample_genome.fa
            Path to output file:                            sample_data/output.csv
            Length of the gRNA sequence:                    20
            Length of flanking region for verification:     200
            Number of available CPUs:                       12
            Path to annotation file in GFF format:          /sample_data/sample_genome.gff3
            Path to annotation_info file in TXT format:     None
            Designing for CRISPR system:
                Streptococcus pyogenes Cas9                 True
            
    Genome file sample_data/sample_genome.fa successfully imported
    formatting genome
    Genome file sample_data/sample_genome.fa successfully formatted
    The genome was successfully converted to a dictionary
    Annotation file sample_data/sample_genome.gff successfully imported
    Annotation database successfully generated

            Initiating PAM site detection.
            
            Please wait, this may take a while...
            

                17314 Cas9 PAM sites were found on Chr01
                
    The output file has been generated at sample_data/output.csv

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

    <details>
      <summary>Click here for a preview of sample_genome_output.csv</summary>
      
      ```
      crispr_id,crispr_sys,sequence,long_sequence,chromosome,start_pos,end_pos,cutsite,strand,on_site_score,features
      A01NW7FGPN,cas9,GGUUAGAUUAGGGCUGUGUU,GCCAGGGUUAGAUUAGGGCUGUGUUAGGGU,Chr01,77,97,94,+,0.0388536288320188,,completed
      A01QLYYDXZ,cas9,GUGCGUACGUAAAAUCAGUA,UCCGUGUGCGUACGUAAAAUCAGUAUACAA,Chr01,411,431,428,+,0.5402860690510768,,completed
      A01RVDM36X,cas9,GGAGUGAAGUGGAAUCUGAG,GCCAUGGAGUGAAGUGGAAUCUGAGAGUAG,Chr01,471,491,488,+,0.6091206773441749,,completed
      A011BCL3O8,cas9,GCAUAAUGAUGUGAGUGCAU,GCCGUGCAUAAUGAUGUGAGUGCAUUUGGU,Chr01,521,541,538,+,0.05499563386100348,,completed
      A01J39N6WJ,cas9,UGAGGCAAGUGCCGUGCAUA,ACCGCUGAGGCAAGUGCCGUGCAUAAUGAU,Chr01,536,556,553,+,0.1199299810564199,,completed
      A01AG1OHUM,cas9,AUGAGAUAUAGAUAUCAAAA,GCCGAAUGAGAUAUAGAUAUCAAAAUGUGG,Chr01,605,625,622,+,0.7178898545540503,,completed
      A01HBGMKT6,cas9,CGAAUGAGAUAUAGAUAUCA,ACCGCCGAAUGAGAUAUAGAUAUCAAAAUG,Chr01,608,628,625,+,0.16003821718975292,,completed
      A01LPJZIOI,cas9,UAUGUUUAUGataataacaa,GCCAAUAUGUUUAUGataataacaactttt,Chr01,777,797,794,+,0.8846303027545066,,completed
      A019RQ1MNC,cas9,AAGCCAAUAUGUUUAUGata,ACCACAAGCCAAUAUGUUUAUGataataac,Chr01,784,804,801,+,0.2709374516526445,,completed
      ...
      ```
      
    </details>

#### [Back to top](#)
---

<!-- ## Contributing

#### [Back to top](#)
--- -->

## Disclosures

>This work was funded by the DOE Center for Advanced Bioenergy and Bioproducts Innovation (U.S. Department of Energy, Office of Science, Office of Biological and Environmental Research under Award Number DE-SC0018420). Any opinions, findings, and conclusions or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the U.S. Department of Energy.

#### [Back to top](#)
---
