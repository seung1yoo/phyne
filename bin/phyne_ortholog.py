import os
import pickle
from phyne_common import PHYNE_COMMON


class PHYNE_ORTHOLOGOUS(PHYNE_COMMON):
	def __init__(self):
		self.set_dir(args.outdir)
		self.set_dir('{0}/{1}'.format(args.outdir, 'tmp'))
		self.set_dir('{0}/{1}'.format(args.outdir, 'Report'))


class READ_CONF:
	def __init__(self):
		self.toolDic = self.make_ToolList(args.conf)

	def make_ToolList(self, conf):
		toolDic = dict()
		toolList = list()
		for line in open(conf):
			if line.startswith('#'):
				pass
			else :
				lines = line.split('\t')
				key = lines[0]
				value = lines[2].strip()
				toolDic.setdefault(key, value)

		return toolDic


class ORTHOMCL(READ_CONF, PHYNE_COMMON, PHYNE_ORTHOLOGOUS):
	def __init__(self):
		READ_CONF.__init__(self)

	def make_cmd_orthoMCL(self, toolDic):
		cmd = list()

		cmd.append('python')
		cmd.append(self.orthoMCL_exe)
		cmd.append(toolDic['ortho_fasta_dir'])
		cmd.append('{0}/{1}'.format(args.outdir, '1.orthoMCL'))
		cmd.append('40')

		return cmd

	def make_orthoMCL_xls(self):
		cmd = list()
		cmd.append('python')
		cmd.append(self.parser_exe)
		cmd.append('{0}/{1}/{2}'.format(args.outdir, '1.orthoMCL', 'mclOutput_I2_0.group'))
		cmd.append('{0}/{1}/{2}'.format(args.outdir, '1.orthoMCL', 'mclOutput_I2_0.group'))
		return cmd

	def orthoMCL_parser(self):
		cmd = list()
		cmd.append('python')
		cmd.append(self.fommatting_exe)
		cmd.append('-f')
		cmd.append('{0}/{1}'.format(args.outdir, '1.orthoMCL/2.All_fasta/goodProteins.fasta'))
		cmd.append('-p')
		cmd.append('{0}/{1}'.format(args.outdir, '1.orthoMCL'))
		cmd.append('-o')
		cmd.append('{0}/{1}'.format(args.outdir, 'tmp'))
		return cmd


class RESULT_FASTA(PHYNE_COMMON):
	def __init__(self):
		self.c_file = '{0}/{1}'.format(args.outdir, '1.orthoMCL/mclOutput_I2_0.group.count.xls')
		self.g_file = '{0}/{1}'.format(args.outdir, '1.orthoMCL/mclOutput_I2_0.group.gene.xls')


	def make_speciesDic(self):
		speciesDic = dict()
		for line in open(self.c_file, 'r'):
			if line.startswith('#'):
				items = line.rstrip('\n').split('\t')[1:]
			else :
				break
		for num, item in enumerate(items) :
			speciesDic.setdefault(num, item)

		return speciesDic

	def make_seqDic(self):
		with open('{0}/{1}/{2}'.format(args.outdir, 'tmp', 'seqDic.txt'), 'rb') as handle:
			seqDic = pickle.load(handle)
		return seqDic

	def make_tmpFile(self):
		self.seqDic = self.make_seqDic()
		self.speciesDic = self.make_speciesDic()

		singlecopy = '{0}/{1}/{2}'.format(args.outdir, 'Report', 'SinglecopyGene.xls')
		singleout = open(singlecopy, 'w')
		headers = ['ClusterID', 'Species', 'Gene', 'Sequence']
		singleout.write('{0}\n'.format('\t'.join(headers)))

		n = 1
		for groupId, geneDic in self.seqDic.iteritems():
			tmpfile = '{0}/{1}/{2}{3}.fa'.format(args.outdir, 'tmp', 'tmp', str(n))
			tmp = open(tmpfile, 'w')

			for num in range(len(geneDic.keys())):
				tmp.write(geneDic[num])
				cluster = '{0}{1}'.format('Cluster', n)
				species = self.speciesDic[num]
				gene_seq = geneDic[num].lstrip('>').split('\n')
				gene = gene_seq[0]
				seq = gene_seq[1]

				preunits = [cluster, species, gene, seq]
				units = [x.strip() for x in preunits]
				singleout.write('{0}\n'.format('\t'.join(units)))
			n += 1
		singleout.close()

		return n

	def run_mafft(self):
		self.seqDic = self.make_seqDic()
		self.mafft_path = '{0}/{1}'.format(args.outdir, '2.mafft')
		self.set_dir(self.mafft_path)

		for num in range(len(self.seqDic)):
			tmpfile = '{0}/{1}/{2}{3}.fa'.format(args.outdir, 'tmp', 'tmp', str(num+1))
			mafft_tmp = '{0}/{1}{2}.fa'.format(self.mafft_path, 'mafft_tmp', str(num+1) )
			cmd = self.make_cmd_mafft_default(tmpfile, mafft_tmp)
			self.run_with_subprocess(cmd, mafft_tmp, '{0}/{1}'.format(self.mafft_path, 'mafft.log'))

	def run_Gblock(self):
		self.seqDic = self.make_seqDic()
		self.Gblock_path = '{0}/{1}'.format(args.outdir, '3.Gblock')
		self.set_dir(self.Gblock_path)
		for num in range(len(self.seqDic)):
			mafft_tmp = '{0}/{1}{2}.fa'.format(self.mafft_path, 'mafft_tmp', str(num+1))
			cmd = self.make_cmd_gblocks_protein(mafft_tmp)
			self.run_with_ossystem(cmd)
		
		gbfiles = 'gb'
		cmd = list()
		cmd.append('mv')
		cmd.append('{0}/*gb*'.format(self.mafft_path))
		cmd.append(self.Gblock_path)
		self.run_with_ossystem(cmd)

	def make_orthologDic(self): 
		self.speciesDic = self.make_speciesDic()
		self.seqDic = self.make_seqDic()
		self.gblock_path = '{0}/{1}'.format(args.outdir, '3.Gblock')

		orthologDic = dict()
		for num in range(len(self.seqDic)):
			gblock_file = '{0}/{1}{2}.fa-gb'.format(self.gblock_path, 'mafft_tmp', str(num+1))
			spec_list = list()
			for n, line in enumerate(open(gblock_file)):
				if line.startswith('>'):
					spec_list.append(n)
					spec_key = len(spec_list) -1
				else :
					seq = line
					if  spec_key in orthologDic.keys():
						orthologDic[spec_key] += seq
					else :
						orthologDic.setdefault(spec_key, seq)	

		return orthologDic

	def make_result_fa(self):
		self.make_tmpFile()
		self.run_mafft()
		self.run_Gblock()
		self.make_orthologDic()

		orthologDic = self.make_orthologDic()

		ortholog = '{0}/{1}/{2}'.format(args.outdir, 'Report', 'ortholog.fa')
		outfile = open(ortholog, 'w')
		for key, seq in orthologDic.iteritems():
			species = '{0}{1}'.format('>', self.speciesDic[key])
			outfile.write(species)
			seq = '{0}\n'.format(seq.replace(' ', '').replace('\n', ''))
			for n, string in enumerate(seq):
				if n % 60 != 0 :
					outfile.write(string)
				else :
					outfile.write('\n')
					outfile.write(string)
		return outfile


class RUN_PIPELINE(PHYNE_ORTHOLOGOUS, ORTHOMCL, RESULT_FASTA):
	def __init__(self):
		PHYNE_ORTHOLOGOUS.__init__(self)
		ORTHOMCL.__init__(self)
		RESULT_FASTA.__init__(self)

	def run_Tools(self, toolDic):

		if 'Y' in toolDic['orthoMCL']:
			cmd = self.make_cmd_orthoMCL(toolDic)
			self.run_with_subprocess(cmd, 'mclOutput_I2_0.group', 'orthoMCL.log')

		if 'Y' in toolDic['parser']:
			self.run_with_ossystem(self.make_orthoMCL_xls())
			self.run_with_ossystem(self.orthoMCL_parser())

		if 'Y' in toolDic['ortholog_fa']:
			self.make_result_fa()

		if 'fasttree' in toolDic['Tree']:
			path = '{0}/{1}'.format(args.outdir, '4.Fasttree')
			self.set_dir(path)
			ortholog_fn = os.path.join(args.outdir, 'Report', 'ortholog.fa')
			cmd = self.make_cmd_fasttree_pep_default(ortholog_fn, '{0}/{1}'.format(path, 'out_newick'))
			self.run_with_ossystem(cmd)

		newick_fn = os.path.join(args.outdir, '4.Fasttree', 'out_newick')
		report_fn = os.path.join(args.outdir, 'Report', 'newick.txt')
		if not os.path.isfile(report_fn):
			os.link(newick_fn, report_fn)

		drawtree_rscript = os.path.join(args.outdir, 'DrawTree.R')
		tree_png = os.path.join(args.outdir, 'Report', 'PhylogeneticTree.png')
		cmd = self.make_cmd_newickToTree(newick_fn, tree_png, drawtree_rscript)
		self.run_with_ossystem(cmd)


		dist_fn = os.path.join(args.outdir, 'Report', 'distance.txt')
		cmd = self.make_cmd_newickToDist(newick_fn, dist_fn)
		self.run_with_ossystem(cmd)

		pca_rscript = os.path.join(args.outdir, 'PCA.R')
		pca_png = os.path.join(args.outdir, 'Report', 'PCA.png')
		cmd = self.make_cmd_distToPCA(dist_fn, pca_png, pca_rscript)
		self.run_with_ossystem(cmd)

		os.system('cp {0} {1}'.format(self.c_file, os.path.join(args.outdir, 'Report', 'orthologous_group.count.xls')))
		os.system('cp {0} {1}'.format(self.g_file, os.path.join(args.outdir, 'Report', 'orthologous_group.gene.xls')))


def main(args):
	phyne_ortho = PHYNE_ORTHOLOGOUS()
	toolDic = READ_CONF().make_ToolList(args.conf)
	run = RUN_PIPELINE()
	run.run_Tools(toolDic)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', '--conf', help = 'input orthology config file')
	parser.add_argument('-o', '--outdir', help = 'Result directory')
	args = parser.parse_args()
	main(args)

