import os
import pickle
from phyne_common import PHYNE_COMMON


class PHYNE_ORTHOLOGOUS(PHYNE_COMMON):
	def __init__(self):

		self.set_dir(args.outdir)

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
		self.outdir = args.outdir


	def make_cmd_orthoMCL(self, toolDic, outdir):
		cmd = list()

		cmd.append('python')
		cmd.append(self.orthoMCL_exe)
		cmd.append(toolDic['ortho_fasta_dir'])
		cmd.append(outdir + '/1.orthoMCL')
		cmd.append('40')

		return cmd

	def make_orthoMCL_xls(self, outdir):
		cmd = list()
		cmd.append('python')
		cmd.append(self.parser_exe)
		cmd.append('mclOutput_I2_0.group')
		cmd.append(outdir + '/1.orthoMCL/mclOutput_I2_0.group')
		return cmd

	def orthoMCL_parser(self, outdir):
		cmd = list()
		cmd.append('python')
		cmd.append(self.fommatting_exe)
		cmd.append('-f')
		cmd.append(outdir + '/1.orthoMCL/2.All_fasta/goodProteins.fasta')
		cmd.append('-p')
		cmd.append(outdir + '/1.orthoMCL')
		return cmd


class RESULT_FASTA(PHYNE_COMMON):
	def __init__(self):
		self.c_file = self.outdir + '/1.orthoMCL/mclOutput_I2_0.group.count.xls'

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
		with open('seqDic.txt', 'rb') as handle:
			seqDic = pickle.load(handle)

		return seqDic

	def make_tmpFile(self):
		self.seqDic = self.make_seqDic()
		self.set_dir('tmp')

		n = 0
		for groupId, geneDic in self.seqDic.iteritems():
			tmpfile = 'tmp/tmp' + str(n) + '.fa'
			tmp = open(tmpfile, 'w')
			for num in range(len(geneDic.keys())):
				tmp.write(geneDic[num])
			n += 1

		return n

	def run_mafft(self):
		self.seqDic = self.make_seqDic()
		self.mafft_path = args.outdir + '/2.mafft'
		self.set_dir(self.mafft_path)

		for num in range(len(self.seqDic)):
			tmpfile = 'tmp/tmp' + str(num) + '.fa'
			mafft_tmp = self.mafft_path + '/mafft_tmp' + str(num) 
			cmd = self.make_cmd_mafft_default(tmpfile, mafft_tmp)
			self.run_with_subprocess(cmd, mafft_tmp, 'mafft.log')


	def run_Gblock(self):
		self.seqDic = self.make_seqDic()
		self.Gblock_path = args.outdir + '/3.Gblock'
		self.set_dir(self.Gblock_path)
		for num in range(len(self.seqDic)):
			mafft_tmp = self.mafft_path + '/mafft_tmp' + str(num)
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
		self.Gblock_path = args.outdir + '/3.Gblock'

		orthologDic = dict()
		for num in range(len(self.seqDic)):
			gblock_file = self.Gblock_path + '/mafft_tmp' + str(num) + '-gb'
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

		outfile = open('ortholog.fa', 'w')
		for key, seq in orthologDic.iteritems():
			species = '{0}{1}\n'.format('>', self.speciesDic[key])
			outfile.write(species)
			outfile.write(seq)
		return outfile


class RUN_PIPELINE(PHYNE_ORTHOLOGOUS, ORTHOMCL, RESULT_FASTA):
	def __init__(self):
		PHYNE_ORTHOLOGOUS.__init__(self)
		ORTHOMCL.__init__(self)
		RESULT_FASTA.__init__(self)

	def run_Tools(self, toolDic, outdir):

		if 'Y' in toolDic['orthoMCL']:
			cmd = self.make_cmd_orthoMCL(toolDic, outdir)
			self.run_with_subprocess(cmd, 'mclOutput_I2_0.group', 'orthoMCL.log')

		if 'Y' in toolDic['parser']:
			self.run_with_ossystem(self.make_orthoMCL_xls(self.outdir))
			self.run_with_ossystem(self.orthoMCL_parser(self.outdir))

		if 'Y' in toolDic['ortholog_fa']:
			self.make_result_fa()

		if 'fasttree' in toolDic['Tree']:
			path = outdir + '/4.Fasttree'
			self.set_dir(path)
			cmd = self.make_cmd_fasttree_pep_default('ortholog.fa', path +'/out_newick')
			self.run_with_ossystem(cmd)

		rm_cmd = ['rm', 'seqDic.txt']
		self.run_with_ossystem(rm_cmd)

def main(args):
	phyne_ortho = PHYNE_ORTHOLOGOUS()
	toolDic = READ_CONF().make_ToolList(args.conf)
	run = RUN_PIPELINE()
	run.run_Tools(toolDic, args.outdir)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', '--conf', help = 'input orthology config file')
	parser.add_argument('-o', '--outdir', help = 'Result directory')
	args = parser.parse_args()
	main(args)

