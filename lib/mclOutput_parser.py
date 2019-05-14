import glob 


class INPUTFILES:
	def __init__(self):
		self.path = args.path
		self.f_file = args.fasta
		self.c_file = '{0}/{1}'.format(self.path, '/mclOutput_I2_0.group.count.xls')
		self.g_file = '{0}/{1}'.format(self.path, '/mclOutput_I2_0.group.gene.xls')

class ORTHOMCL(INPUTFILES):
	def __init__(self):
		#PHYNE_ORTHOLOGOUS.__init__(self)
		INPUTFILES.__init__(self)

		self.speciesDic, self.groups = self.extract_groupId()
#		self.groupDic = self.makeGroupDic(self.groups)
		
	def extract_groupId(self):
		groups = list()
		speciesDic = dict()
		for line in open(self.c_file):
			if line.startswith("#"):
				species = line.rstrip('\n').split('\t')[1:]
				speciesNum = len(line.split('\t'))-1
				find_singleCopy = '\t1'*speciesNum
				speciesDic = self.makeValueDic(species)
			else:
				if find_singleCopy in line:
					groupId = line.split('\t')[0]
					groups.append(groupId)

		return speciesDic, groups

	def makeValueDic(self, values):
		valueDic = dict()
		for n, value in enumerate(values):
			valueDic.setdefault(n, value)

		return valueDic

	def makeGroupDic(self, groups):
		groupDic = dict()
		for line in open(self.g_file):
			for group in groups:
				if group in line:
					items = line.rstrip('\n').split('\t')[1:]
					for n, item in enumerate(items):
						groupDic.setdefault(group, {}).setdefault(n, item)

		return groupDic

	def parsingSeqFile(self, f_file):
		p_file = '{0}/{1}'.format(args.output, 'parsing.txt')
		parsing = open(p_file, 'w')
		for n, line in enumerate(open(f_file)):
			if line.startswith('>'):
				parsing.write('\n'+line)
			else:
				parsing.write(line.rstrip('\n').rstrip('*'))
		return p_file

	def makeSeqDic(self):
		import pickle
		groupDic = self.makeGroupDic(self.groups)
		p_file = self.parsingSeqFile(self.f_file)
		seqDic = dict()
		geneNum = -2
		for groupId, geneDic in groupDic.iteritems():
			for num, gene in geneDic.iteritems():
				for n, line in enumerate(open(p_file)):
					if '>' in line:
						if gene in line:
							geneNum = n
							seq = '> {0}{1}'.format(gene, '\n')
						else:
							geneNum = -2
					else :
						if n == geneNum + 1:
							seq += line
							seqDic.setdefault(groupId, {}).setdefault(num, seq)
						else:
							pass

		with open('{0}/{1}'.format(args.output, 'seqDic.txt'), 'wb') as handle:
			pickle.dump(seqDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

		return seqDic


def main(args):
#	phyne_ortho = PHYNE_ORTHOLOGOUS():
	run = ORTHOMCL()	

	seqDic = run.makeSeqDic()

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-p", "--path", help = " OrhoMCL PATH")
	parser.add_argument("-f", "--fasta", help = "goodProteins.fasta")
	parser.add_argument("-o", "--output", help = "output directory")
	args = parser.parse_args()
	main(args)
