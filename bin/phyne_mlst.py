
import os
import sys
import subprocess
import json

from phyne_common import PHYNE_COMMON
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class PHYNE_MLST(PHYNE_COMMON):
    def __init__(self, outdir, prefix):
        print('[MASSAGE] Start phyne_mlst')

        self.mlst_home = PHYNE_COMMON.mlst_home
        self.outdir = self.set_dir(outdir)
        self.tmpdir = self.set_dir('{0}/tmp'.format(outdir))
        self.prefix = prefix

        self.mlst_json          = '{0}/{1}.json'.format(self.tmpdir, self.prefix)
        self.mlst_novel_fa      = '{0}/{1}.novel.fa'.format(self.tmpdir, self.prefix)
        self.mlst_table         = '{0}/{1}.table'.format(self.tmpdir, self.prefix)
        self.mlst_profile_fn    = '{0}/{1}.profile.xls'.format(self.outdir, self.prefix)
        self.mlst_multiseq_fn   = '{0}/{1}.multiseq.fa'.format(self.tmpdir, self.prefix)
        self.mlst_db_fa         = '{0}/db/blast/mlst.fa'.format(self.mlst_home)
        self.mlst_multialign_fn = '{0}/{1}.multialign.fa'.format(self.tmpdir, self.prefix)
        self.mlst_multialign_gb = '{0}/{1}.multialign.fa-gb'.format(self.tmpdir, self.prefix)
        self.mlst_newick_fn     = '{0}/{1}.newick'.format(self.outdir, self.prefix)


    def load_input(self, infns, labels):
        self.in_dic = dict()
        for idx, infn in enumerate(infns):
            if os.path.isfile(infn):
                print('[MASSAGE] load input : {0} ==> {1}'.format(infn, labels[idx]))
                self.in_dic.setdefault(infn, labels[idx])
            else:
                print('[ERROR] Can not found {0}'.format(infn))
                sys.exit()
        return self.in_dic

    def make_cmd_mlst(self, target_scheme):
        _cmd_s = list()
        _cmd_s.append('{0}/bin/mlst'.format(self.mlst_home))
        if target_scheme:
            _cmd_s.append('--scheme')
            _cmd_s.append(target_scheme)
        _cmd_s.append('--json')
        _cmd_s.append(self.mlst_json)
        _cmd_s.append('--novel')
        _cmd_s.append(self.mlst_novel_fa)
        _cmd_s.extend(self.in_dic.keys())
        return _cmd_s

    def write_mlst_profiling(self):
        out_fh = open(self.mlst_profile_fn, 'w')
        for line in open(self.mlst_table):
            items = line.rstrip('\n').split('\t')
            new_items = [self.in_dic[items[0]]]
            new_items.extend(items[1:])
            out_fh.write('{0}\n'.format('\t'.join(new_items)))
        out_fh.close()

    def write_mlst_multifa(self):
        self.profile_dic = dict()
        for fn, info_dic in json.load(open(self.mlst_json)).items():
            sample_id = self.in_dic[fn]
            self.profile_dic.setdefault(sample_id, info_dic)

        self.allele_seq_dic = dict()
        for sample_id, info_dic in self.profile_dic.items():
            scheme = info_dic['scheme']
            for gene, allele_num in info_dic['alleles'].items():
                allele_seq_id = self.get_allele_seq_id(scheme, gene, allele_num)
                if allele_seq_id:
                    self.allele_seq_dic.setdefault(allele_seq_id, '')
                else:
                    print('[WARNING] Can not get allele_seq_id in {0} {1} {2}'.format(sample_id, gene, allele_num))
        self.get_allele_seq()

        out_fh = open(self.mlst_multiseq_fn, 'w')
        out_log_fh = open('{0}.order'.format(self.mlst_multiseq_fn), 'w')
        for sample_id, info_dic in self.profile_dic.items():
            scheme = info_dic['scheme']
            _id_s = []
            _seq_s = []
            _len_s = []
            for gene, allele_num in sorted(info_dic['alleles'].items()):
                allele_seq_id = self.get_allele_seq_id(scheme, gene, allele_num)
                if allele_seq_id in self.allele_seq_dic:
                    allele_seq = self.allele_seq_dic[allele_seq_id]
                    _id_s.append(allele_seq_id)
                    _seq_s.append(allele_seq)
                    _len_s.append(str(len(allele_seq)))
                else:
                    _id_s.append('{0}.{1}_-'.format(scheme, gene))
                    _seq_s.append('')
                    _len_s.append('0')
                #
            if len(''.join(_seq_s)):
                multi_seq = Seq(''.join(_seq_s))
                record = SeqRecord(multi_seq, id=sample_id, description='')
                SeqIO.write(record, out_fh, 'fasta')
            else:
                pass

            log_items = [sample_id]
            for idx, _id in enumerate(_id_s):
                log_items.append('{0}(len{1})'.format(_id, _len_s[idx]))
            out_log_fh.write('{0}\n'.format('\t'.join(log_items)))

        out_fh.close()
        out_log_fh.close()


    def get_allele_seq_id(self, scheme, gene, allele_num):
        allele_seq_id = None
        try:
            int(allele_num)
        except ValueError:
            if allele_num in ['-']:
                pass
            elif allele_num.startswith('~'):
                print('[WARNING] There are Missing data. type "~"')
                allele_seq_id = '{0}.{1}{2}'.format(scheme, gene, allele_num)
            elif allele_num.endswith('?'):
                print('[ERROR] There are Missing data. type "?"')
                print('Call to Seung-il Yoo')
                sys.exit()
            elif ',' in allele_num:
                print('[ERROR] There are Missing data. type ","')
                print('Call to Seung-il Yoo')
                sys.exit()
            else:
                print('[ERROR] There are Missing data.')
                print('Call to Seung-il Yoo')
                sys.exit()
        else:
            allele_seq_id = '{0}.{1}_{2}'.format(scheme, gene, allele_num)
        return allele_seq_id

    def get_allele_seq(self):

        if os.path.isfile(self.mlst_db_fa):
            for record in SeqIO.parse(open(self.mlst_db_fa), 'fasta'):
                if record.id in self.allele_seq_dic:
                    self.allele_seq_dic[record.id] = str(record.seq)

        if os.path.isfile(self.mlst_novel_fa):
            for record in SeqIO.parse(open(self.mlst_novel_fa), 'fasta'):
                if record.id in self.allele_seq_dic:
                    self.allele_seq_dic[record.id] = str(record.seq)

        # below lines are made for check
        for allele_seq_id, allele_seq in self.allele_seq_dic.items():
            if not allele_seq:
                print('[ERROR] allele_seq is missing. {0}'.format(allele_seq_id))
                sys.exit()


def main(args):
    phyne_mlst = PHYNE_MLST(args.outdir, args.prefix)
    phyne_mlst.load_input(args.infn_s, args.label_s)

    cmd = phyne_mlst.make_cmd_mlst(args.target_scheme)
    phyne_mlst.run_with_subprocess(cmd, phyne_mlst.mlst_table, '{0}/mlst.run.log'.format(phyne_mlst.tmpdir))

    phyne_mlst.write_mlst_profiling()

    phyne_mlst.write_mlst_multifa()

    cmd = phyne_mlst.make_cmd_mafft_default(phyne_mlst.mlst_multiseq_fn, phyne_mlst.mlst_multialign_fn)
    phyne_mlst.run_with_subprocess(cmd, phyne_mlst.mlst_multialign_fn, '{0}/mafft.run.log'.format(phyne_mlst.tmpdir))

    cmd = phyne_mlst.make_cmd_gblocks_dna(phyne_mlst.mlst_multialign_fn)
    phyne_mlst.run_with_subprocess(cmd, '{0}/gblocks.run.out'.format(phyne_mlst.tmpdir), '{0}/gblocks.run.log'.format(phyne_mlst.tmpdir))

    cmd = phyne_mlst.make_cmd_fasttree_dna_default(phyne_mlst.mlst_multialign_gb, phyne_mlst.mlst_newick_fn)
    phyne_mlst.run_with_subprocess(cmd, phyne_mlst.mlst_newick_fn, '{0}/fasttree.run.log'.format(phyne_mlst.tmpdir))



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--target-scheme', help='recommanded, if you do not specify, will be rised error.')
    parser.add_argument('--infn-s', nargs='+',
            help='Simply just give it a genome files in FASTA or GenBank format, optionally compressed with gzip.')
    parser.add_argument('--label-s', nargs='+',
            help='Must be same order with --infn-s')
    parser.add_argument('--outdir', default='./Results')
    parser.add_argument('--prefix', default='phyne_mlst')
    args = parser.parse_args()
    main(args)

