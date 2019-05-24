# phyne (PHYlogenetic analysis pipeliNE)

## Workflow

![workflow_phyne](./image/workflow_phyne.png)

phyne는 phylogenetic analysis를 위한 파이프라인으로 다양한 인풋을 지원합니다.


## Development Environment
 - server : wolf
 - language : python 2.7

## Dependents

### Common dependents
 - [MAFFT](http://mafft.cbrc.jp/alignment/software) v7.123b
 - [MUSCLE](http://www.drive5.com/muscle) v3.8.31
 - [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html) 0.91b
 - [FastTree](http://www.microbesonline.org/fasttree/) 2.1.11
 - R package (RColorBrewer, ggplot2, phylogram)

### MLST dependents
 - [mlst](https://github.com/tseemann/mlst)
 - [PubMLST](https://pubmlst.org/general.shtml)

### ortholog dependents
 - [OrthoMCL](https://orthomcl.org/orthomcl/) v2.0.9.M
 - [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) 2.2.26

## Usage

### Help message

```
$python bin/phyne.py --help
usage: phyne.py [-h] [--mode {mlst,nssnp,ortholog}] [--config CONFIG]
                [--outdir OUTDIR] [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  --mode {mlst,nssnp,ortholog}
  --config CONFIG
  --outdir OUTDIR
  --prefix PREFIX
```


### MLST

#### Workflow

![workflow_mlst_mode](./image/workflow_mlst_mode.png)

#### General Command line
```
$python bin/phyne.py --mode mlst --config [mlst.conf] --outdir [result] --prefix [testset]
```

#### Input & Output

 - Input : Bacterial genome sequences (FASTA format)

 - Output
  1. {outdir}/{prefix}.mlst_profile.xls : MLST profile (a tab-separated line)


```
GCA_000439795.1 aphagocytophilum        64      pheS(42)        glyA(32)        fumC(28)        mdh(18) sucA(40)        dnaN(28)      atpA(26)
GCA_000013125.1 aphagocytophilum        161     pheS(42)        glyA(32)        fumC(28)        mdh(18) sucA(82)        dnaN(28)      atpA(26)
GCA_000439775.1 aphagocytophilum        64      pheS(42)        glyA(32)        fumC(28)        mdh(18) sucA(40)        dnaN(28)      atpA(26)
GCA_000689655.1 aphagocytophilum        215     pheS(103)       glyA(77)        fumC(28)        mdh(4)  sucA(94)        dnaN(28)      atpA(60)
GCA_000689635.2 aphagocytophilum        82      pheS(3) glyA(33)        fumC(29)        mdh(3)  sucA(2) dnaN(2) atpA(4)
GCA_000964685.1 aphagocytophilum        64      pheS(42)        glyA(32)        fumC(28)        mdh(18) sucA(40)        dnaN(28)      atpA(26)
GCA_000689615.1 aphagocytophilum        217     pheS(104)       glyA(78)        fumC(70)        mdh(53) sucA(95)        dnaN(77)      atpA(1)
GCA_000478425.1 aphagocytophilum        64      pheS(42)        glyA(32)        fumC(28)        mdh(18) sucA(40)        dnaN(28)      atpA(26)
GCA_000968455.1 aphagocytophilum        -       pheS(42)        glyA(32)        fumC(28)        mdh(18?)        sucA(40)        dnaN(28)       atpA(26)
```

 >  - col 1 : the genome name
 >  - col 2 : the matching PubMLST scheme name
 >  - col 3 : the ST (sequence type)
 >  - col 4 ~ : the allele IDs
 
 
 2. {outdir}/{prefix}.mlst_sequence.fa : Allele sequences are arranged in order
 3. {outdir}/{prefix}.mlst_sequence.fa.order : THE order
 4. {outdir}/{prefix}.mlst_align.fa : multiple sequence alignment ([fasta](https://en.wikipedia.org/wiki/FASTA_format))
 5. {outdir}/{prefix}.mlst_align.phylip : multiple sequence alignment ([phylip](http://rosalind.info/glossary/phylip-format/))
 6. {outdir}/{prefix}.mlst_align.newick : tree ([newick](https://en.wikipedia.org/wiki/Newick_format))
 7. {outdir}/{prefix}.mlst_align.dist : Sequence distance metrics
 8. {outdir}/{prefix}.PhylogeneticTree.png : very simple phylogenetic tree figure
 9. {outdir}/{prefix}.PCA.png : PCA plot


#### Config file

- Without target scheme

```
$ cat bin/phyne_mlst.conf
{
        "target_scheme" : "",
        "input_genome" : {
                "mygenome" : "denovogenome.fa",
                "relativegenome1" : "Salmonella_enterica_subsp_enterica_serovar_Typhimurium_DT104_v1.fa",
                "relativegenome2" : "Salmonella_enterica_subsp_enterica_serovar_Typhi_str_CT18_v1.fa",
                "relativegenome3" : "Salmonella_enterica_subsp_enterica_serovar_Weltevreden_str_10259_v0.2.fa"
        },
        "phyne_mlst_exe" : "/BiO/BioPeople/siyoo/phyne/bin/phyne_mlst.py"
}
```

- With target scheme [Scheme list](https://pubmlst.org/data/dbases.xml)
 
```
$ cat bin/phyne_mlst.conf
{
        "target_scheme" : "aphagocytophilum",
        "input_genome" : {
                "genome1" : "GCA_000013125.1_ASM1312v1_genomic.fna",
                "genome2" : "GCA_000439775.1_ASM43977v1_genomic.fna",
                "genome3" : "GCA_000439795.1_ASM43979v1_genomic.fna",
                "genome4" : "GCA_000478425.1_ASM47842v1_genomic.fna",
                "genome5" : "GCA_000689655.1_MRK1.0_genomic.fna"
        },
        "phyne_mlst_exe" : "/BiO/BioPeople/siyoo/phyne/bin/phyne_mlst.py"
}
```



### NSSNP (Non-Synoymous SNP)

```
python bin/phyne.py --mode nssnp --config --outdir --prefix
```

필터링옵션 case 별로 상세하게 적고.
input vcf 이다.

### Ortholog (Single copy genes)

Based on the evolutionary similarity of sequences, a pipeline that draws phylogenetic trees by selecting only single copy genes of two or more species/genomes

#### Workflow
![ortholog workflow](./image/ortholog_workflow.png)

 - Step1. Clustering ortholog Groups of Protein Sequences with orthoMCL
 - Step2. Single copy gene selection
 - Step3. Multiple sequence alignment
 - Step4. Remove ambiguously aligned regions 
 - Step5. Merge clean alignment and Make tree
 - Step6. Calculate the sequence distance

#### General Command line
```
python bin/phyne.py --mode ortholog --config [ortholog.conf] --outdir [result] --prefix [testset]
```

#### Input & Output
- Input : protein sequences (FASTA format)  

- Output : Intermediate directory
 1. {outdir}/1.orthoMCL : orthoMCL-Run Result files
 2. {outdir}/2.mafft : multiple sequence alignment of each single copy gene cluster (fasta)
 3. {outdir}/3.Gblock : The sequence from which the region with low quality is removed (fasta)
 4. {outdir}/4.fasttree : newick file

- Ooutput : Final files  
 1. {outdir}/Report/ortholog.fa : multiple sequnece alignmet merge
 2. {outdir}/Report/newick.txt : tree
 3. {outdir}/Report/SinglecopyGene.xls : single copy gene lists between species
 4. {outdir}/Report/distance.txt : Sequence distance metrics
 5. {outdir}/Report/PhylogeneticTree.png : very simple phylogenetic tree figure
 6. {outdir}/Report/PCA.png : PCA plot

#### Config file
- phyne_ortholog.conf :

```
#otrhoMCL
##########################################################################
orthoMCL        =       Y        
ortho_fasta_dir =       /BiO/BioPeople/boram/Projects/Phylogenetic/pepiteds
ortho_thread    =       40
parser  =       Y
ortholog_fa     =       Y
##########################################################################
Tree    =       fasttree
```

- orthoMCL을 실행 시키려면 Y, 다음 Step만 진행 하고자 한다면 N
- ortho_fasta_dir은 종들의 amino acid sequence (fasta) 파일이 있는 directory 경로 입력  
- ortho_thread는 orthoMCL의 thread를 설정해주는 것으로, 기본 40으로 설정  
- parser는 Single copy gene을 선별하는 script 실행, 이미 진행 했다면 N으로 설정  
- ortholog_fa는 종 별 Single copy gene의 서열을 모아놓은 fasta 생성

! 처음 실행시킬때에는 모두 Y로 실행 시키고, 원하는 step부터 실행시키고자 할 때 N으로 변경하여 사용하면 된다.


## Contributors
- seungil.yoo@theragenetex.com
- ingang.shin@theragenetex.com
- boram.choi@threagenetex.com

lastest update : 2019-05-24

