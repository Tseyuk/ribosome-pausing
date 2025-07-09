# written by Wang,Zeyu for the ribopausing analyze

from numpy import *
import argparse
parser=argparse.ArgumentParser(description="to find the tricodon window from riboseq")
parser.add_argument("-i","--input",type=str,metavar="input",required=True,help="please input the A site bedgraph file's fullname")
parser.add_argument("-r","--reference",type=str,metavar="reference",required=True,help="the longest CDS reference  file's  fullname")
args=parser.parse_args()
input=args.input
reference=args.reference

fgene = open(reference,'r')
fbed = open(input,'r')
fo = open(input[:-3]+"_tricodon_pasuing.tab", 'w')
fo.write("Codons"+"\t"+"AA"+"\t"+"Whole_window"+"\t"+"Signals"+"\n")

#translate codon to aa
def protein(dna):#
    dna = dna.upper()
    rna = dna.replace('T','U')
    table = {'UGA':'*','UAA':'*','UAG':'*','GCA': 'A', 'GCG': 'A', 'GCC': 'A', 'GAC': 'D', 'GAG': 'E', 'GAA': 'E', 'GGG': 'G', 'GGC': 'G', 'GGA': 'G', 'CAC': 'H', 'AAA': 'K', 'AAG': 'K', 'AAC': 'N', 'CCG': 'P', 'CCA': 'P', 'CCC': 'P', 'CAA': 'Q', 'CAG': 'Q', 'CGG': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'AGA': 'R', 'AGC': 'S', 'ACG': 'T', 'ACC': 'T', 'ACA': 'T', 'GCU': 'A', 'UGU': 'C', 'UGC': 'C', 'GAU': 'D', 'UUU': 'F', 'UUC': 'F', 'GGU': 'G', 'CAU': 'H', 'AUC': 'I', 'AUU': 'I', 'AUA': 'I', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'AUG': 'M', 'AAU': 'N', 'CCU': 'P', 'CGU': 'R', 'AGU': 'S', 'UCG': 'S', 'UCC': 'S', 'UCA': 'S', 'UCU': 'S', 'ACU': 'T', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GUA': 'V', 'UGG': 'W', 'UAU': 'Y', 'UAC': 'Y'}
    global aa
    aa = ''
    for i in range(0, int(len(rna)/3)):
        codon = rna[i * 3:i * 3 + 3]
        if codon in table:
            aa+=table.get(codon)
        else:
            continue
    return aa

#get the median value of a list
def get_median(data):
   data = sorted(data)
   size = len(data)
   if size % 2 == 0: # Determine that the list length is even number
    median = (data[size//2]+data[size//2-1])/2
    data[0] = median
   if size % 2 == 1: # Determine that the list length is an odd number
    median = data[(size-1)//2]
    data[0] = median
   return data[0]

#sliding
def sliding(i):
    window = 96#choose the window : 15codon + dicodon + 15codon
    ori = i[1]
    name = i[0]
    del_10codon_cds = ori[30:-30] #delete the 10 codons of the each end
    if len(del_10codon_cds) < window:
        print('not avalable')
    elif len(del_10codon_cds) > window:
        slidecounts = int(round((len(del_10codon_cds) - 96)/3))
        for count in range(slidecounts):
            start = count*3
            stop = count*3+window
            seq = del_10codon_cds[start:stop]
            dicodon= seq[45:54]
            zero = ['0']
            windownum = 96
            kongls = zero*windownum
            l = [seq[i:i + 3] for i in range(0, len(seq), 3)]
            if  name in bedposdic:
                refls = bedposdic.get(name)
                for bed3 in refls:
                    normvalue = average_dic.get(name)
                    #normvalue = median_dic.get(name)
                    position = bed3[0]
                    reads = bed3[1]
                    relapos =((int(position)-130) - start)
                    if 0 < relapos <= 96  and normvalue!= None:  # only the signal from the 0 frame of each codon could be recorded
                        lspos = relapos-1 # 0 +1 -1 frame all TAKING INTO ACCOUNT
                        # lspos = int(relapos / 3)
                        normed_val= round(reads / normvalue ,3)#normalize with average is better
                        #normed_val= reads * normvalue
                        kongls[lspos] = str(normed_val)
                    else:
                        continue
                a = 0
                for read in kongls:
                    if read != '0':
                        a = a + 1
                if a >= 15:
                    fo.write(dicodon + '\t' + protein(dicodon) + '\t' + ','.join(l) + '\t' + ','.join(kongls) + '\n')
            else:break

#get the dic of gene name and its' length and the specific sequence
genedic = {}
genelist = []
for key in fgene:
    key = (key.replace('\n', ''))[1:]#the name of each gene
    value = ((next(fgene, None)).strip())[100:-100]#the specific sequence
    genels = []
    if value is None:
        # oops, we got to the end of the file early
        raise ValueError('Invalid file format, missing value')
    genedic[key] = [value,len(value)]
    genels.append(key)
    genels.append(value)
    genels.append(len(value))
    genelist.append(genels)

#read the bg file and read signal above 0 of each gene and release the reads of continuouse site(like 351 353 1)
ori = ()
bedposls = []
bedposdic = {}
for sig in fbed:
    sig = sig.replace('\n', '')
    # sig = sig.split('   ')
    sig = sig.split('\t')
    name = sig[0]
    read = int(sig[-1])
    pos = sig[1].strip()
    posend = sig[2].strip()
    kongtotal = list(ori)
    kposls= [name,pos,read]#the signal is in the first nt of a codon
    if kposls[-1] == 0:
        continue
    else:
        if int(posend) - int(pos) == 1:
            bedposls.append(kposls)
            newls = [pos,read]
            kongtotal.append(newls)
        else:
            count = int(posend) - int(pos)
            for copo in range(count):
                kposls[1] = int(pos) + int(copo)
                add = str(kposls[1])
                lsn = [name, add, read]
                bedposls.append(lsn)
                newls = [add,read]
                kongtotal.append(newls)
        if name in bedposdic:
            olddicls = bedposdic.get(name)
            olddicls=olddicls+kongtotal
            bedposdic[name] = olddicls
        if name not in bedposdic:
            bedposdic[name] = kongtotal


# total signal above 0 of each gene
totalsigdic = {}
for s in bedposls:
    name = s[0]
    read = int(s[-1])
    if name in totalsigdic:
        totalsigdic.get(name).append(read)
    if name not in totalsigdic :  # the avarage signal
        totalsigdic[name] = [read]
#get the average normalizing index
average_dic = {}
#median_dic = {}
for key, value in totalsigdic.items():
    if key in genedic and len(value) > 30:
        average_dic[key] = float(mean(value))
        #median_dic[key] = (get_median(value) / len(value))
    else:
        continue
for k in genelist:
    sliding(k)
print('finished!')
fo.close()
