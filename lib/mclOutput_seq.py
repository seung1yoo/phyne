import glob 

def files(args):
	path = args.path
	f_file = args.fasta
	c_file = path + '/mclOutput_I2_0.group.count.xls'
	g_file = path + '/mclOutput_I2_0.group.gene.xls'
	out_file = args.outfile
	return f_file, c_file, g_file, out_file

def extract_groupId(c_file):
	groups = list()
	speciesDic = dict()
	for line in open(c_file):
		if line.startswith("#"):
			species = line.rstrip('\n').split('\t')[1:]
			speciesNum = len(line.split('\t'))-1
			find_singleCopy = '\t1'*speciesNum
			speciesDic = makeValueDic(species)
		else:
			if find_singleCopy in line:
				groupId = line.split('\t')[0]
				groups.append(groupId)
	
	return speciesDic, groups

def makeValueDic(values):
	valueDic = dict()
	for n, value in enumerate(values):
		valueDic.setdefault(n, value)

	return valueDic

def makeGeneDic(groups, g_file):
	geneDic = dict()
	for line in open(g_file):
		for group in groups:
			if group in line:
				items = line.rstrip('\n').split('\t')[1:]
				for n, item in enumerate(items):
					if n in geneDic.keys():
						geneDic[n].append(item)
					else :
						geneDic.setdefault(n, [item])
	return geneDic

def parsingSeqFile(f_file):
	p_file = 'parsing.txt'
	parsing = open(p_file, 'w')
	for n, line in enumerate(open(f_file)):
		if line.startswith('>'):
			parsing.write('\n'+line)
		else:
			parsing.write(line.rstrip('\n').rstrip('*'))
	return p_file

def makeSeqDic(geneDic, f_file):
	p_file = parsingSeqFile(f_file)
	seqDic = dict()
	geneNum = -2
	for num, genes in geneDic.iteritems():
		seq = ''
		seqDic.setdefault(num, seq)
		for gene in genes:
			for n, line in enumerate(open(p_file)):
				if '>' in line:
					if gene in line:
						geneNum = n
					else:
						geneNum = -2
				else :
					if n == geneNum + 1:
						seqDic[num] = seqDic[num] + line.rstrip('\n')
					else:
						pass
	return seqDic

def makeOutfile(speciesDic, seqDic, o_file):
	outfile = open(o_file, 'w')
	
	for num, species in speciesDic.iteritems():
		outfile.write('>'+speciesDic[num]+'\n')
		outfile.write(seqDic[num]+'\n')

	outfile.close()

def main(args):
	f_file, c_file, g_file, out_file = files(args)
	speciesDic, groups = extract_groupId(c_file)
	geneDic = makeGeneDic(groups, g_file)
	seqDic = makeSeqDic(geneDic, f_file)
	makeOutfile(speciesDic, seqDic, out_file)

if __name__=='__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-p", "--path", help = " OrhoMCL PATH")
	parser.add_argument("-f", "--fasta", help = "goodProteins.fasta")
	parser.add_argument("-o", "--outfile", help = "out fasta file") 
	args = parser.parse_args()
	main(args)
