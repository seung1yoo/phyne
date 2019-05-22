
import os
import sys
from pandas import Series
from pandas import DataFrame
from phyne_common import PHYNE_COMMON

class NSSNP(PHYNE_COMMON):
	def __init__(self):
		pass

	def read_config(self, config_fn):
		self.conf_dic = dict()
		for line in open(config_fn):
			if not line.strip():
				continue
			if line.startswith('#'):
				continue
			splitted = line.rstrip('\n').split('=')
			conf_key = splitted[0]
			conf_value = splitted[1]
			if conf_key in ['SAMPLE_IDS']:
				samples = conf_value.split(',')
				self.conf_dic.setdefault(conf_key, samples)
			elif conf_key in ['INPUT']:
				self.conf_dic.setdefault(conf_key, conf_value)
			elif conf_key in ['GQ', 'GT_DP']:
				self.conf_dic.setdefault(conf_key, int(conf_value))
			elif conf_key in ['REF_ratio', 'ALT_ratio']:
				self.conf_dic.setdefault(conf_key, float(conf_value))
			elif conf_key in ['HETERO_ratio']:
				hetero_ratio_s = conf_value.split(',')
				self.conf_dic.setdefault(conf_key, [float(x) for x in hetero_ratio_s])
			else:
				print('un-expected conf_key')
				sys.exit()

	def read_afile(self, infn):
		ofh_dic = dict()
		for line in open(infn) :
				word = line.rstrip()
				splitted = word.split('\t')
				if word.startswith('##') :
					continue
				elif word.startswith('#CHROM'):
					S_INFO = '{0}'.format('\t'.join(splitted[0:9]))
					S_ID = '{0}'.format('\t'.join(splitted[9:]))
					ofh_dic.setdefault(S_INFO,S_ID)
					continue
				if 'missense_variant' in splitted[7] :
					S_INFO = '{0}'.format('\t'.join(splitted[0:9]))
					S_ID = '{0}'.format('\t'.join(splitted[9:]))
					if len(splitted[3]) != 1 or len(splitted[4]) != 1 :
						continue
					ofh_dic.setdefault(S_INFO,S_ID)
		return ofh_dic


	def filter_ofh_dic(self, ofh_dic):
		ID_lst = self.conf_dic['SAMPLE_IDS']
		lst_idx = []
		sample_idx = []
		filter_gt_dic = dict()
		for key,value in ofh_dic.iteritems():
			key_splitted = key.split('\t')
			v_split = value.split('\t')
			if key.startswith('#CHROM'):
				headder = v_split
				for i in v_split :
					if i in ID_lst :
						i_idx = headder.index(i)
						lst_idx.append(i_idx)
						sample_idx.append(v_split[i_idx])
				filter_gt_dic.setdefault(key,'\t'.join(sample_idx))
				continue
			Len = len(ID_lst)
			gt_lst = []
			cnt = 0
			for i_idx in lst_idx:
				n = v_split[i_idx]
				if n.startswith('./.'):
					continue
				ref_gt = key_splitted[3]
				alt_gt = key_splitted[4]
				gt_info = n.split(':')
				gt = gt_info[0]
				ref = gt_info[1].split(',')[0]
				alt = gt_info[1].split(',')[1]
				DP = gt_info[2]
				GQ = gt_info[3]
				if int(GQ) >= self.conf_dic['GQ'] and int(DP) >= self.conf_dic['GT_DP']:
					if gt == '0/0' :
						if (int(ref) * 100.0 / int(DP)) < self.conf_dic['REF_ratio']:
							continue
						gt = '{0}{1}'.format(ref_gt, ref_gt)
					elif gt == '0/1' :
						if (int(ref) * 100.0 / int(DP)) < self.conf_dic['HETERO_ratio'][0] or \
                        (int(ref) * 100.0 / int(DP)) > self.conf_dic['HETERO_ratio'][1]:
							continue
						gt = '{0}{1}'.format(ref_gt, alt_gt)
					elif gt == '1/1' :
						if (int(alt) * 100.0 / int(DP)) < self.conf_dic['ALT_ratio']:
							continue
						gt = '{0}{1}'.format(alt_gt, alt_gt)
				else:
					continue
				gt_lst.append(gt)
			gt_len = len(gt_lst)
			if Len != gt_len :
				continue
			filter_gt_dic.setdefault(key, '\t'.join(gt_lst))

		return filter_gt_dic

	def output_make_file(self, filter_gt_dic):
		make_dic = dict()
		ID_lst = self.conf_dic['SAMPLE_IDS']
		lst_idx = []
		for key,value in filter_gt_dic.iteritems():
			key_splitted = key.split('\t')
			v_split = value.split('\t')
			if key.startswith('#CHROM'):
				id_split = value.split('\t')
				headder = id_split
				for i in id_split :
					if i in ID_lst :
						i_idx = headder.index(i)
						lst_idx.append(i_idx)
				continue
			else:
				for i_idx in lst_idx :
					make_dic.setdefault(id_split[i_idx], []).append(v_split[i_idx])
		return make_dic

	def make_outfile_handle_dic(self, make_dic, outdir, prefix):
		self.set_dir(outdir)
		out_fn = '{0}/{1}.nssnp.fasta'.format(outdir, prefix)
		out_fh = open(out_fn, 'w')
		for key, value in make_dic.iteritems():
			out_fh.write('>{0}\n{1}\n'.format(key, ''.join(value)))
                out_fh.close()
                return out_fn


def main(args):
	nssnp = NSSNP()
	nssnp.read_config(args.config)
	nssnp.read_afile(nssnp.conf_dic['INPUT'])
	ofh_dic = nssnp.read_afile(nssnp.conf_dic['INPUT'])
	filter_gt_dic = nssnp.filter_ofh_dic(ofh_dic)
	make_dic = nssnp.output_make_file(filter_gt_dic)

	nssnp_fa_fn = nssnp.make_outfile_handle_dic(make_dic, args.outdir, args.prefix)
	nssnp_newick_fn = '{0}/{1}.nssnp.newick'.format(args.outdir, args.prefix)
	fasttree_log = '{0}/fasttree.log'.format(args.outdir)

	cmd = nssnp.make_cmd_fasttree_dna_default(nssnp_fa_fn, nssnp_newick_fn)
	nssnp.run_with_subprocess(cmd, nssnp_newick_fn, fasttree_log)

	drawtree_rscript = os.path.join(args.outdir, 'DrawTree.R')
	tree_png = os.path.join(args.outdir, 'PhylogeneticTree.png')
	cmd = nssnp.make_cmd_newickToTree(nssnp_newick_fn, tree_png, drawtree_rscript)
	nssnp.run_with_ossystem(cmd)

	dist_fn = os.path.join(args.outdir, 'distance.txt')
	cmd = nssnp.make_cmd_newickToDist(nssnp_newick_fn, dist_fn)
	nssnp.run_with_ossystem(cmd)

	pca_rscript = os.path.join(args.outdir, 'PCA.R')
	pca_png = os.path.join(args.outdir, 'PCA.png')
	cmd = nssnp.make_cmd_distToPCA(dist_fn, pca_png, pca_rscript)
	nssnp.run_with_ossystem(cmd)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser(description='NonSyn SNP CALL')
	parser.add_argument('--config')
	parser.add_argument('--outdir')
	parser.add_argument('--prefix')
	args = parser.parse_args()
	main(args)
