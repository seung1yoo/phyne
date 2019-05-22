

import os
import sys
import subprocess


class PHYNE_COMMON:
    phyne_mlst_exe = '/BiO/BioPeople/siyoo/phyne/bin/phyne_mlst.py'
    phyne_nssnp_exe = '/BiO/BioPeople/siyoo/phyne/bin/phyne_nssnp.py'
    phyne_ortholog_exe = '/BiO/BioPeople/siyoo/phyne/bin/phyne_ortholog.py'

    mafft_exe = '/usr/bin/mafft'
    muscle_exe = '/BiO/BioPeople/siyoo/phyne/tools/muscle3.8.31_i86linux64'
    gblocks_exe = '/BiO/BioPeople/siyoo/phyne/tools/Gblocks'
    fasttree_exe = '/BiO/BioPeople/siyoo/phyne/tools/FastTree'
    mlst_home = '/BiO/BioPeople/siyoo/phyne/tools/mlst'
    fa_to_phylip_exe = '/BiO/BioPeople/siyoo/phyne/lib/fasta-to-phylip.py'
    ninja_exe = '/BiO/BioPeople/siyoo/BioTools/ninja_1.2.2/ninja'
    newick_to_dist_exe = '/BiO/BioPeople/siyoo/phyne/lib/newick-to-distance.py'

    orthoMCL_exe = '/BiO/BioPeople/siyoo/phyne/lib/Orthomcl-Pipe.py'
    parser_exe = '/BiO/BioPeople/siyoo/phyne/lib/Orthomcl_Parsing.py'
    fommatting_exe = '/BiO/BioPeople/siyoo/phyne/lib/mclOutput_parser.py'

    _pca_template = '''\
library(RColorBrewer)
data <- read.table('[INPUT]', sep='\t', header=T, row.names=1)
pc <- princomp(as.matrix(data), scores=TRUE)
png(file='[OUTPUT]', height=500, width=500)
plot(pc$score[,1:2], type='n', xlab='PC 1', ylab='PC 2', main = 'PCA Plot')
text(pc$score[,1:2], labels=colnames(data), col=brewer.pal(length(colnames(data)),"Set1"))
dev.off()
'''
    pca_template = '''\
library(ggplot2)
library(RColorBrewer)
data <- read.table('[INPUT]', sep='\t', header=T, row.names=1)
df_pca <- prcomp(as.matrix(data))
df_out <- as.data.frame(df_pca$x)
df_out$samples <- row.names(df_out)
percentage <- round(df_pca$sdev / sum(df_pca$sdev) * 100, 2)
percentage <- paste( colnames(df_out), "(", paste( as.character(percentage), "%", ")", sep="") )
theme <- theme(panel.background = element_blank(), 
               panel.border=element_rect(fill=NA),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               strip.background=element_blank(),
               axis.text.x=element_text(colour="black"),
               axis.text.y=element_text(colour="black"),
               axis.ticks=element_line(colour="black"),
               plot.margin=unit(c(1,1,1,1),"line"))
p <- ggplot(df_out,aes(x=PC1,y=PC2,color=samples)) +
            geom_point(size=3) + 
            theme +
            xlab(percentage[1]) + ylab(percentage[2])
ggsave(filename = "[OUTPUT]", width = 8, height = 6)
'''

    drawtree_template = '''\
library(phylogram)
x <- read.dendrogram(file = '[INPUT]')
png(file='[OUTPUT]', height=500, width=500)
plot(x, yaxt = "n")
dev.off()
'''

    def set_dir(self, a_dir):
        if not os.path.isdir(a_dir):
            os.mkdir(a_dir)
        return a_dir

    def run_with_subprocess(self, args, stdout_fn, stderr_fn):
        print('[MASSAGE] RUN with subprocess : {0}'.format(' '.join(args)))
        p = subprocess.Popen(args, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        stdout, stderr = p.communicate()

        out = open(stdout_fn, 'w')
        out.write(stdout)
        out.close()

        out = open(stderr_fn, 'w')
        out.write(stderr)
        out.close()

        return stdout, stderr

    def run_with_ossystem(self, args):
        print('[MASSAGE] RUN with os.system : {0}'.format(' '.join(args)))
        cmd = ' '.join(args)
        os.system(cmd)

    def make_cmd_mafft_default(self, in_fa, out_align_fa):
        args = list()
        args.append(self.mafft_exe)
        args.append(in_fa)
        #args.append('>')
        #args.append(out_align_fa)
        return args

    def make_cmd_muscle_default(self, in_fa, out_align_fa):
        args = list()
        args.append(self.muscle_exe)
        args.append('-in')
        args.append(in_fa)
        args.append('-out')
        args.append(out_align_fa)
        return args

    def make_cmd_gblocks_dna(self, in_align_fa):
        args = list()
        args.append(self.gblocks_exe)
        args.append(in_align_fa)
        args.append('-t=d')
        return args

    def make_cmd_gblocks_protein(self, in_align_fa):
        args = list()
        args.append(self.gblocks_exe)
        args.append(in_align_fa)
        args.append('-t=p')
        return args

    def make_cmd_fasttree_dna_default(self, in_align_fa, out_newick):
        args = list()
        args.append(self.fasttree_exe)
        args.append('-nt')
        args.append(in_align_fa)
        #args.append('>')
        #args.append(out_newick)
        return args

    def make_cmd_fasttree_dna_default_for_ossystemrun(self, in_align_fa, out_newick):
        args = list()
        args.append(self.fasttree_exe)
        args.append('-nt')
        args.append(in_align_fa)
        args.append('>')
        args.append(out_newick)
        return args

    def make_cmd_fasttree_pep_default(self, in_align_fa, out_newick):
        args = list()
        args.append(self.fasttree_exe)
        args.append(in_align_fa)
        args.append('>')
        args.append(out_newick)
        return args

    def make_cmd_FaToPhylip(self, in_fa, out_phylip):
        args = list()
        args.append('python')
        args.append(self.fa_to_phylip_exe)
        args.extend(['--input-fasta', in_fa])
        args.extend(['--output-phy', out_phylip])
        return args

    def make_cmd_ninja_FaToDist(self, in_fa, out_dist, fa_type):
        args = list()
        args.append(self.ninja_exe)
        args.extend(['--in', in_fa])
        args.extend(['--out', out_dist])
        args.extend(['--method', 'default'])
        args.extend(['--in_type', 'a']) # an alignment in fasta format
        args.extend(['--out_type', 'd']) # d = distance matrix / t = tree (newick?)
        if fa_type in ['dna']:
            args.extend(['--alph_type', 'd']) # d = dna / a = amino acid
            args.extend(['--corr_type', 'k']) # k for dna / s for amino acid
        elif fa_type in ['protein']:
            args.extend(['--alph_type', 'a']) # d = dna / a = amino acid
            args.extend(['--corr_type', 's']) # k for dna / s for amino acid
        else:
            print('ERR : Check the fa_type (dna or protein)')
            sys.exit()
            
        # Correction for multiple same-site substitutions.
        #'n' no correction
        #'j' jukes-cantor correction  { dist = -3/4 * ln (1 - 4/3 * dist ) }
        #'k' kimura 2-parameter correction { dist = -1/2 * ln ( (1-2p-q)*sqrt(1-2q) ) }
        #'s' FastTree's scoredist-like correction { dist = -1.3 * ln (1.0 - dist) }
        return args

    def DistToTable(self, in_fn, out_fn):
        samples = list()
        itemss = list()
        for line in open(in_fn):
            items = line.rstrip().split()
            if len(items) in [1]:
                continue
            samples.append(items[0])
            itemss.append(items)
        out = open(out_fn, 'w')
        out.write('DISTANCE\t{0}\n'.format('\t'.join(samples)))
        for items in itemss:
            out.write('{0}\n'.format('\t'.join(items)))
        out.close()

    def make_cmd_distToPCA(self, in_dist, out_png, pca_rscript):
        out = open(pca_rscript, 'w')
        out.write(self.pca_template.replace("[INPUT]", in_dist).replace("[OUTPUT]", out_png))
        out.close()
        #
        args = list()
        args.append('/usr/bin/Rscript')
        args.append(pca_rscript)
        return args

    def make_cmd_newickToTree(self, in_newick, out_png, drawtree_rscript):
        out = open(drawtree_rscript, 'w')
        out.write(self.drawtree_template.replace("[INPUT]", in_newick).replace("[OUTPUT]", out_png))
        out.close()
        #
        args = list()
        args.append('/usr/bin/Rscript')
        args.append(drawtree_rscript)
        return args

    def make_cmd_newickToDist(self, in_newick, out_dist):
        args = list()
        args.append('python')
        args.append(self.newick_to_dist_exe)
        args.append(in_newick)
        args.append('>')
        args.append(out_dist)
        return args

            



