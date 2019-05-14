import sys, glob, os

if len(sys.argv) != 4 :
	print "Usage : python %s <in.fasta.dir> <out.dir> <thread>" % sys.argv[0]
	sys.exit()

def makeDirectory(file):
	if not os.path.exists(file):
		os.system("mkdir -p %s" % file)

INDIR = sys.argv[1]
OUTDIR = sys.argv[2]
THREAD = sys.argv[3]

ORTHOMCLPATH = '/BiO/sgpark/DenovoAssemblyPipeline/OrthologousGensCluster/orthomclSoftware-v2.0.9.M/bin'
MCL = "/BiO/sgpark/DenovoAssemblyPipeline/OrthologousGensCluster/mcl-14-137/bin"
BLASTPATH = "/BiO/BioTools/ncbi-blast/blast-2.2.26/bin"

CONFIGUREFILE = '/BiO/sgpark/DenovoAssemblyPipeline/OrthologousGensCluster/orthomclSoftware-v2.0.9.M/config.example.txt'


INFASTDIR = "%s/1.compliant_fasta" % OUTDIR
makeDirectory(INFASTDIR)


fastafiles = glob.glob("%s/*.fasta" % INDIR)

for seqfile in fastafiles :
	sampleid = seqfile.split("/")[-1].split(".fasta")[0]
	print "%s/orthomclAdjustFasta %s %s 1 %s" % (ORTHOMCLPATH, sampleid, seqfile, INFASTDIR)
	os.system("%s/orthomclAdjustFasta %s %s 1 %s" % (ORTHOMCLPATH, sampleid, seqfile, INFASTDIR))

print "%s/orthomclFilterFasta %s 10 20" % (ORTHOMCLPATH, INFASTDIR)
os.system("%s/orthomclFilterFasta %s 10 20" % (ORTHOMCLPATH, INFASTDIR))


ALLFASTA = "%s/2.All_fasta" % OUTDIR
makeDirectory(ALLFASTA)

print "mv goodProteins.fasta %s/goodProteins.fasta" % ALLFASTA
os.system("mv goodProteins.fasta %s/goodProteins.fasta" % ALLFASTA)

print "mv poorProteins.fasta %s/poorProteins.fasta" % ALLFASTA
os.system("mv poorProteins.fasta %s/poorProteins.fasta" % ALLFASTA)

print "%s/formatdb -i %s/goodProteins.fasta -p T" % (BLASTPATH, ALLFASTA)
os.system("%s/formatdb -i %s/goodProteins.fasta -p T" % (BLASTPATH, ALLFASTA))


BLASTRESULTDIR = "%s/3.blast_result" % OUTDIR
makeDirectory(BLASTRESULTDIR)

print "%s/blastall -d %s/goodProteins.fasta -i %s/goodProteins.fasta -p blastp -e 1e-5 -m 8 -o %s/goodProteins.fasta.blastp -a %s" % (BLASTPATH, ALLFASTA, ALLFASTA, BLASTRESULTDIR, THREAD)
os.system("%s/blastall -d %s/goodProteins.fasta -i %s/goodProteins.fasta -p blastp -e 1e-5 -m 8 -o %s/goodProteins.fasta.blastp -a %s" % (BLASTPATH, ALLFASTA, ALLFASTA, BLASTRESULTDIR, THREAD))

print "cp %s %s/config.txt" % (CONFIGUREFILE, OUTDIR)
os.system("cp %s %s/config.txt" % (CONFIGUREFILE, OUTDIR))

print "%s/orthomclInstallSchema %s/config.txt %s/config.log" % (ORTHOMCLPATH, OUTDIR, OUTDIR)
os.system("%s/orthomclInstallSchema %s/config.txt %s/config.log" % (ORTHOMCLPATH, OUTDIR, OUTDIR))

print "%s/orthomclBlastParser %s/goodProteins.fasta.blastp %s > %s/similarSequences.txt" % (ORTHOMCLPATH, BLASTRESULTDIR, INFASTDIR, BLASTRESULTDIR)
os.system("%s/orthomclBlastParser %s/goodProteins.fasta.blastp %s > %s/similarSequences.txt" % (ORTHOMCLPATH, BLASTRESULTDIR, INFASTDIR, BLASTRESULTDIR))

print "%s/orthomclLoadBlast %s/config.txt %s/similarSequences.txt" % (ORTHOMCLPATH, OUTDIR, BLASTRESULTDIR)
os.system("%s/orthomclLoadBlast %s/config.txt %s/similarSequences.txt" % (ORTHOMCLPATH, OUTDIR, BLASTRESULTDIR))

print "cp orthomcl %s/orthomcl.load" % (BLASTRESULTDIR)
os.system("cp orthomcl %s/orthomcl.load" % (BLASTRESULTDIR))

print "%s/orthomclPairs %s/config.txt %s/orthomcl.pair.log cleanup=no" % (ORTHOMCLPATH, OUTDIR, BLASTRESULTDIR)
os.system("%s/orthomclPairs %s/config.txt %s/orthomcl.pair.log cleanup=no" % (ORTHOMCLPATH, OUTDIR, BLASTRESULTDIR))

print "cp orthomcl %s/orthomcl.pair" % (BLASTRESULTDIR)
os.system("cp orthomcl %s/orthomcl.pair" % (BLASTRESULTDIR))

print "%s/orthomclDumpPairsFiles %s/config.txt" % (ORTHOMCLPATH, OUTDIR)
os.system("%s/orthomclDumpPairsFiles %s/config.txt" % (ORTHOMCLPATH, OUTDIR))

print "%s/mcl mclInput --abc -I 1.5 -o %s/mclOutput_I1_5" % (MCL, OUTDIR)
os.system("%s/mcl mclInput --abc -I 1.5 -o %s/mclOutput_I1_5" % (MCL, OUTDIR))
print "%s/orthomclMclToGroups H 000000 < %s/mclOutput_I1_5 > %s/mclOutput_I1_5.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR)
os.system("%s/orthomclMclToGroups H 000000 < %s/mclOutput_I1_5 > %s/mclOutput_I1_5.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR))

print "%s/mcl mclInput --abc -I 2.0 -o %s/mclOutput_I2_0" % (MCL, OUTDIR)
os.system("%s/mcl mclInput --abc -I 2.0 -o %s/mclOutput_I2_0" % (MCL, OUTDIR))
print "%s/orthomclMclToGroups H 000000 < %s/mclOutput_I2_0 > %s/mclOutput_I2_0.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR)
os.system("%s/orthomclMclToGroups H 000000 < %s/mclOutput_I2_0 > %s/mclOutput_I2_0.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR))

print "%s/mcl mclInput --abc -I 4.0 -o %s/mclOutput_I4_0" % (MCL, OUTDIR)
os.system("%s/mcl mclInput --abc -I 4.0 -o %s/mclOutput_I4_0" % (MCL, OUTDIR))
print "%s/orthomclMclToGroups H 000000 < %s/mclOutput_I4_0 > %s/mclOutput_I4_0.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR)
os.system("%s/orthomclMclToGroups H 000000 < %s/mclOutput_I4_0 > %s/mclOutput_I4_0.group" % (ORTHOMCLPATH, OUTDIR, OUTDIR))
