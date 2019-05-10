

import os
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

    orthoMCL_exe = '/BiO/BioPeople/siyoo/phyne/lib/Orthomcl-Pipe.py'
    parser_exe = '/BiO/BioPeople/siyoo/phyne/lib/Orthomcl_Parsing.py'
    fommatting_exe = '/BiO/BioPeople/siyoo/phyne/lib/mclOutput_seq.py'

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



