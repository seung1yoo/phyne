import sys

if len(sys.argv) != 3 :
	print "Usage : python %s <OrthoMclResult.group> <OrthoMcl.outprefix>" % sys.argv[0]
	sys.exit()

f = open(sys.argv[1], 'r')

LISTSPE = {}
GROUPCNT = {}
GROUPGENE = {}

for lines in f.xreadlines() :
	words = lines.rstrip().split()
	GROUPCNT[words[0][:-1]] = {}
	GROUPGENE[words[0][:-1]] = {}
	for num in range(1,len(words)) :
		sp = words[num].split("|")
		if not LISTSPE.has_key(sp[0]) :
			LISTSPE[sp[0]] = 1
		if not GROUPCNT[words[0][:-1]].has_key(sp[0]) :
			GROUPCNT[words[0][:-1]][sp[0]] = 0
			GROUPGENE[words[0][:-1]][sp[0]] = []
		GROUPCNT[words[0][:-1]][sp[0]] += 1
		GROUPGENE[words[0][:-1]][sp[0]].append(sp[1])

f.close()

w1 = open(sys.argv[2] + ".count.xls", 'w')
w2 = open(sys.argv[2] + ".gene.xls", 'w')

LIST = sorted(LISTSPE.keys())

w1.write("#GroupName\t%s\n" % "\t".join(LIST))
w2.write("#GroupName\t%s\n" % "\t".join(LIST))

for group in sorted(GROUPCNT.keys()) :
	w1.write(group)
	w2.write(group)
	for spname in LIST :
		if GROUPCNT[group].has_key(spname) :
			w1.write("\t%d" % GROUPCNT[group][spname])
			w2.write("\t%s" % (" ".join(GROUPGENE[group][spname])))
		else :
			w1.write("\t0")
			w2.write("\t.")
	w1.write("\n")
	w2.write("\n")

w1.close()
w2.close()
