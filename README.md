# phyne (PHYlogenetic analysis pipeliNE)

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

#### Example Command line
```
$python bin/phyne.py -mode mlst --config [mlst.conf] --outdir [result] --prefix [testset]
```

#### Input & Output
- Input : Bacterial genome sequence (FASTA format)
- Output
 - {outdir}/{prefix}.mlst_profile.xls : MLST profile (a tab-separated line)
 - col 1 : the genome name
 - col 2 : the matching PubMLST scheme name
 - col 3 : the ST (sequence type)
 - col 4~ : the allele IDs

| example_1 | senterica | 365 | aroC(130) | dnaN(97) | hemD(25) | hisD(125) | purE(84) | sucA(9) | thrA(101) |
| example_2 | senterica | 2 | aroC(1) | dnaN(1) | hemD(2) | hisD(1) | purE(1) | sucA(1) | thrA(5) |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |


 - {outdir}/{prefix}.mlst_sequence.fa
 - {outdir}/{prefix}.mlst_sequence.fa.order
 - {outdir}/{prefix}.mlst_align.fa
 - {outdir}/{prefix}.mlst_align.phylip
 - {outdir}/{prefix}.mlst_align.newick
 - {outdir}/{prefix}.mlst_align.dist
 - {outdir}/{prefix}.PhylogeneticTree.png
 - {outdir}/{prefix}.PCA.png


#### Config file

- target scheme 이 없을 때, (de-novo mlst)

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

- target scheme 이 있을 때, (targeted mlst)
 - [Scheme list](https://pubmlst.org/data/dbases.xml)
 
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

### Orthologus gene (Single copy genes)

```
python bin/phyne.py --mode ortholog --config --outdir --prefix
```

configration file 정보
안하고 싶으면 N으로 바꿔라 



constact : seungil.yoo@theragenetex.com

lastest update : 2019-05-23