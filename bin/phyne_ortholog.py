
import os
from phyne_common import PHYNE_COMMON



class PHYNE_ORTHOLOGOUS(PHYNE_COMMON):
	def __init__(self):
		self.outdir = self.set_dir(args.outdir)


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


class ORTHOMCL(READ_CONF, PHYNE_COMMON):
	def __init__(self):
		READ_CONF.__init__(self)

	def make_cmd_orthoMCL(self, toolDic, outdir):
		cmd = list()

		cmd.append('python')
		cmd.append(self.orthoMCL_exe)
		cmd.append(toolDic['ortho_fasta_dir'])
		cmd.append(outdir + '/orthoMCL')
		cmd.append('40')

		return cmd

	def make_orthoMCL_xls(self, outdir):
		cmd = list()
		cmd.append('python')
		cmd.append(self.parser_exe)
		cmd.append('mclOutput_I2_0.group')
		cmd.append(outdir + '/orthoMCL/mclOutput_I2_0.group')
		return cmd

	def orthoMCL_parser(self, outdir):
		cmd = list()
		cmd.append('python')
		cmd.append(self.fommatting_exe)
		cmd.append('-f')
		cmd.append(outdir + '/orthoMCL/2.All_fasta/goodProteins.fasta')
		cmd.append('-p')
		cmd.append(outdir + '/orthoMCL')
		cmd.append('-o')
		cmd.append(outdir + '/orthoMCL/orthoMCL_group.fasta')

		return cmd


class RUN_PIPELINE(PHYNE_ORTHOLOGOUS, ORTHOMCL):
	def __init__(self):
		PHYNE_ORTHOLOGOUS.__init__(self)
		ORTHOMCL.__init__(self)

	def run_Tools(self, toolDic, outdir):

		if 'Y' in toolDic['orthoMCL']:
			cmd = self.make_cmd_orthoMCL(toolDic, outdir)
			self.run_with_subprocess(cmd, 'mclOutput_I2_0.group', 'orthoMCL.log')

		if 'Y' in toolDic['orthoMCL_parsing']:
			self.run_with_ossystem(self.make_orthoMCL_xls(outdir))
			self.run_with_ossystem(self.orthoMCL_parser(outdir))

		if 'mafft' in toolDic['aligner']:
			path = outdir + '/mafft'
			self.set_dir(path)
			cmd = self.make_cmd_mafft_default(outdir + '/orthoMCL/orthoMCL_group.fasta', path + '/mafft_align.fa')
			self.run_with_subprocess(cmd, path + '/mafft_align.fa', path + '/mafft.log')
			align_fa = path + '/mafft_align.fa'

		if 'muscle' in toolDic['aligner']:
			path = outdir + '/muscle'
			self.set_dir(path)
			cmd = self.make_cmd_muscle_default(outdir + '/orthoMCL/orthoMCL_group.fasta', path + '/muscle_align.fa')
			self.run_with_ossystem(cmd)
			align_fa = path + '/muscle_align.fa'

		if 'Y' in toolDic['Gblock']:
			path = outdir + '/Gblock'
			self.set_dir(path)
			os.system('cp {0} {1}'.format( align_fa, path))
			align_fa = align_fa.split('/')[-1]
			align_fa = path + '/' + align_fa
			cmd = self.make_cmd_gblocks_protein(align_fa)
			self.run_with_ossystem(cmd)
			align_fa = align_fa = align_fa + '-gb'

		if 'fasttree' in toolDic['Tree']:
			path = outdir + '/Fasttree'
			self.set_dir(path)
			cmd = self.make_cmd_fasttree_pep_default(align_fa, path +'/out_newick')
			self.run_with_ossystem(cmd)


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



