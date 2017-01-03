from collections import defaultdict

def readfile_myb(filename):
    d = defaultdict(list)
    ## read in files and return list of genes and logFC and pval
    for line in open(filename):
        args = line.strip().split(" ")
        gene_name = args[2].split(".")[0]
        logfc = float(args[0])
        pval = float(args[1])
        d[gene_name].append((logfc,pval))
    return d


def readfile_mybsig(filename):
    d = defaultdict(list)
    ## read in files and return list of genes and logFC and pval
    for line in open(filename):
        args = line.strip().split(",")
        gene_name = args[1]
        logfc = float(1000.0)
        pval = float(10.0)
        d[gene_name.upper()].append((logfc,pval))
    return d

def readfile_vnd(filename):
    d = defaultdict(list)
    ## read in files and return list of genes and logFC and pval
    genes = []
    for line in open(filename):
        if line[0] == ",": continue
        args = line.strip().split(",")
        if len(args) < 6: continue
        for g in args[2].split(";"):
            gene_name = g.upper()
            if gene_name == "NA": continue
            #genes.append(gene_name)
            try:
                logfc = float(args[-6])
                pval = float(args[-2])
                if pval >= 0.05: continue
                #if logfc < 0 : continue
                d[gene_name].append((logfc,pval))
            except ValueError: continue
            genes.append(gene_name)
            #print len(set(genes))
            #print len(genes)
    return d

def readfile_vnda(filename):
    d = defaultdict(list)
    ## read in files and return list of genes and logFC and pval
    genes = []
    for line in open(filename):
        if line[0] == ",": continue
        args = line.strip().split(",")
        if len(args) < 6: continue
        for g in args[5].split(";"):
            gene_name = g.upper()
            if gene_name == "NA": continue
            #genes.append(gene_name)
            try:
                logfc = float(args[-6])
                pval = float(args[-2])
                #print logfc
                #print pval
                if pval >= 0.05: continue
                d[gene_name].append((logfc,pval))
            except ValueError: continue
            genes.append(gene_name)
            #print gene_name
            #print len(set(genes))
            #print len(genes)
            ## split these
            ###also some lines are sep by commas
    return d


def read_vas(vas_file):
    genes = []
    for line in open(vas_file):
        genes.append(line.strip())
    print line
    return  set(genes)


def find_shared_genes(vnd7,myb46,vas):
    vnd7_genes = set(vnd7.keys())
    myb46_genes = set(myb46.keys())
    shared = vnd7_genes.intersection(myb46_genes)
    vnd7_vas = vnd7_genes.intersection(vas)
    myb46_vas = vas.intersection(myb46_genes)
    all_shared = shared.intersection(vas)
    
    print "AT4G18770" in vnd7_genes
    print len(vas), len(vnd7_genes), len(myb46_genes)
    print len(shared), len(myb46_vas), len(vnd7_vas), len(all_shared)
    print all_shared

def main(vnd7_microarray, myb46_microarray, vas_genes):
    vnd7 = readfile_vnd(vnd7_microarray)
    #myb46 = readfile_myb(myb46_microarray)
    myb46 = readfile_mybsig(myb46_microarray)
    vas = read_vas(vas_genes)
    find_shared_genes(vnd7, myb46, vas)
    #print len(shared_genes), vnd7_n, myb46_n
    #print shared_genes


#main("GSE23555.csv", "myb46_ramires.txt")
#main("GSE23555.csv", "GSE16143.csv")
#main("GSE23555.csv", "GSE25838.csv")
main("../data/GSE23555.csv", "../data/TPJ_3989_sm_Table S1.csv","../data/At_Vas_genes.txt")


### the vnd7 dataset has half pos half neg
## but the myb46 dataset does not

### why only 4 genes in common... maybe think about using a diff dataset...
### double check code is working properly

