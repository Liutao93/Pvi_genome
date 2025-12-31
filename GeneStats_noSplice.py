#!/usr/bin/env python3
# -*- coding: utf-8 -*- 
from Bio import SeqIO
import sys, os, argparse, glob, os.path,time,glob,re,gzip
#################################################
#北京组学生物科技有限公司
#python 语言入门与基础绘图：
# https://zzw.xet.tech/s/2bdd89
#生信基础课:
#https://zzw.xet.tech/s/1bWTkv
####################################################

parser = argparse.ArgumentParser(description='This script was used to stat gene attributes. Note one gene must have one transcripts')
parser.add_argument('-g', '--gff', dest='gff', required=True, help='input gff  file')
#parser.add_argument('-f', '--fa', dest='fa', required=True, help='input genome fasta file')
parser.add_argument('-s', '--species', dest='species', required=True, help='Input species name and Key of output file')
parser.add_argument('-S', '--soft', dest='soft', required=False, default=None ,help='predict software name')
parser.add_argument('-o','--outDir',dest='outDir',required=False,default=os.getcwd(),help='specify the output file dir,default is current dir')

args = parser.parse_args()
###################################################################
#if not os.path.exists(args.outDir): os.mkdir(args.outDir)
#dout=os.path.abspath(args.outDir)
#cp_template(dout,args.temp)

args.gff=os.path.abspath(args.gff)
#args.fa=os.path.abspath(args.fa)
if not os.path.exists(args.outDir): os.mkdir(args.outDir)
args.outDir=os.path.abspath(args.outDir)


class GenomeFeature(object):
    def __init__(self):
        self.scaffold = None
        self.start = None
        self.end = None
        self.strand = None
        self.attribute =None
        self.id=None

    def __len__(self):
        return self.end-self.start+1

    def load_gff_record(self,gff_record):
        self.scaffold = gff_record.scaffold
        self.start = gff_record.start
        self.end = gff_record.end
        self.strand = gff_record.strand
        self.attribute = gff_record.attribute
        self.id=gff_record.id


class Gff_Record():
    def __init__(self,record_string):
        record = record_string.strip().split("\t")
        id=None
        pid=None
        g=re.search(r"ID=([^;]+)",record[8])
        if g :
            id=g.group(1)
        g1=re.search(r"Parent=([^;]+)",record[8])
        if g1 :
            pid=g1.group(1)
        self.scaffold = record[0]
        self.source = record[1]
        self.feature = record[2]
        self.start = int(record[3])
        self.end = int(record[4])
        self.score = record[5]
        self.strand = record[6]
        self.frame = record[7]
        self.attribute = record[8]
        self.id=id
        self.parent_id=pid
    def __len__(self):
        return self.end-self.start+1       


class Transcript(GenomeFeature):
    def __init__(self,genome):
        super(Transcript,self).__init__()
        self.genome=genome
        self.exons=[]
        self.cds=[]
        self.introns = []
        

    def add_exon(self,exon):
        self.exons.append(exon)
    def add_cds(self,cds):
        self.cds.append(cds)
    def retrive_introns(self):
        if len(self.exons) == 0:
            print("fail to find an exon in this Transcript: %s"%(self.attribute))
            self.add_intron(self.start, self.end)
        elif len(self.exons) > 1:
            self.exons.sort(key=lambda x:x.start)
            for i in range(1,len(self.exons)):
                if self.exons[i-1].end +1 < self.exons[i].start -1:
                    self.add_intron(self.exons[i-1].end + 1, self.exons[i].start - 1)

    def add_intron(self,start,end):
        new_intron = GeneFeature(self, "intron")
        new_intron.scaffold = self.scaffold
        new_intron.start = start
        new_intron.end = end
        new_intron.strand = self.strand
        new_intron.attribute = "intron at %s[%d:%d] %s" % (new_intron.scaffold, new_intron.start, new_intron.end, new_intron.strand)
        self.introns.append(new_intron)


class Gene(GenomeFeature):
    def __init__(self,genome):
        super(Gene,self).__init__()
        self.genome=genome
        self.exons=[]
        self.cds=[]
        self.introns = []
        self.transcripts = []
        self.id=None
        

    def add_exons(self):
        for transcript in self.transcripts:
            self.exons += transcript.exons

    def add_introns(self):
        for transcript in self.transcripts:
            self.introns += transcript.introns
    def add_cds(self):
        for transcript in self.transcripts:
            self.cds += transcript.cds

class GeneFeature(GenomeFeature):

    def __init__(self,gene,feature_type):
        super(GeneFeature,self).__init__()
        self.gene=gene
        self.feature_type=feature_type

#     def get_sequence(self):
# 
#         # would get 1 more base at 3' end
#         genome=self.gene.genome
#         scaffold = genome.genomeid_dict[self.scaffold]
# 
#         if self.strand == "+":
#             return scaffold[self.start-1:self.end+1]
#         elif self.strand == "-":
#             return scaffold[self.start-2:self.end].reverse_complement()


class Node():
    def __init__(self,number,ntype):
        self.number=number
        self.type=ntype

    def return_number(self):
        return self.number


class IntergenicRegions():
    def __init__(self,genome,scaffold):
        self.genome=genome
        self.scaffold = scaffold
        self.nodes=[Node(1,"scaffold_start"),Node(len(genome.genomeid_dict[scaffold]),"scaffold_end")]

    def update(self,gene_inteval):
        gene_start=Node(gene_inteval[0],"gene_start")
        gene_end=Node(gene_inteval[1],"gene_end")
        self.nodes.insert(0,gene_start)
        self.nodes.append(gene_end)
        self.nodes.sort(key=Node.return_number)
        idx_gene_start=self.nodes.index(gene_start)
        idx_gene_end=self.nodes.index(gene_end)

        if idx_gene_start%2 == 0:
            nodes_lst1 = self.nodes[0:idx_gene_start] or []
        else:
            new_node=Node(gene_inteval[0]-1,"inter_end")
            nodes_lst1 = self.nodes[0:idx_gene_start]+ [new_node]

        if idx_gene_end%2 == 1:
                try:
                        nodes_lst2 = self.nodes[idx_gene_end+1:]
                except IndexError:
                        nodes_lst2 = []

        else:
            new_node=Node(gene_inteval[1]+1,"inter_start")
            nodes_lst2 = [new_node] + self.nodes[idx_gene_end+1:]

        self.nodes=nodes_lst1+nodes_lst2

    def call_length_inter(self):
        inter_len_lst = []
        for i in range(0,len(self.nodes)//2):
            if self.nodes[2*i+1].type=="scaffold_end" or self.nodes[2*i].type=="scaffold_start":
                continue
            elif self.nodes[2*i+1].type!="inter_end" or self.nodes[2*i].type!="inter_start":
                print ("got errors in call intergenetic length")
            length = self.nodes[2*i+1].number - self.nodes[2*i].number + 1
            inter_len_lst.append(length)
        return inter_len_lst


class Genome():
    def __init__(self):
        #self.genomeid_dict = SeqIO.to_dict(SeqIO.parse(file,"fasta"))
        self.genes=[]
        self.transcripts=[]
        #self.intergenicregions_all={}
        self.cnt_donors={}
        self.cnt_acceptors={}
        self.genes_with_intron=[]
#         for scaffold in self.genomeid_dict:
#             self.intergenicregions_all[scaffold]=IntergenicRegions(self,scaffold)

    def annotate_genes(self,anno):
#         with open(args.outDir+"/"+"introns.fa","w"):
#             pass
        if(anno.endswith(".gz")):
            annotation=gzip.open(anno,"r")
        else:
            annotation=open(anno,"r")
        for line in annotation:
            if not line.strip():
                continue
            elif line.strip()[0] == "#":
                continue
            record = Gff_Record(str(line))
            if record.feature == "gene":
                new_gene = Gene(self)
                new_gene.load_gff_record(record)
                self.genes.append(new_gene)
            elif record.feature in ["transcript","mRNA"]:
                new_transcript=Transcript(self)
                new_transcript.load_gff_record(record)
                self.transcripts.append(new_transcript)
                if new_transcript.start<new_gene.start or new_transcript.end>new_gene.end:
                    print ("unrecognised transcript\n")
                new_gene.transcripts.append(new_transcript)
            elif record.feature == "exon":
                new_exon=GeneFeature(new_transcript,"exon")
                new_exon.load_gff_record(record)
                if new_exon.start<new_transcript.start or new_exon.end>new_transcript.end:
                    print ("unrecognised exon\n")
                else:
                    new_transcript.add_exon(new_exon)
            elif record.feature == "CDS":
                new_cds=GeneFeature(new_transcript,"cds")
                new_cds.load_gff_record(record)
                if new_cds.start<new_transcript.start or new_cds.end>new_transcript.end:
                    print ("unrecognised cds\n")
                else:
                    new_transcript.add_cds(new_cds)
        annotation.close()
                    
        for transcript in self.transcripts:
            transcript.retrive_introns()

            #for intron in transcript.introns:
#                sequence = intron.get_sequence()
#                 with open(args.outDir+"/"+"introns.fa", "a") as intr:
#                     intr.write(">"+intron.attribute + "\n")
#                     intr.write(str(sequence.seq[:-1])+"\n")
#                 donor = sequence[:2]
#                 self.cnt_donors[str(donor.seq).upper()] = self.cnt_donors.get(str(donor.seq).upper(), 0) + 1
#                 acceptor = sequence[-3:]
#                 self.cnt_acceptors[str(acceptor.seq).upper()] = self.cnt_acceptors.get(str(acceptor.seq).upper(), 0) + 1
        #print(dir(self.genes[0].cds))   
        for gene in self.genes:
            gene.add_exons()
            gene.add_introns()
            gene.add_cds()
            gene_inteval = (gene.start, gene.end)
            #self.intergenicregions_all[gene.scaffold].update(gene_inteval)
        #print(dir(self.genes[0].cds))   

    def get_number_of_genes(self):
        return len(self.genes)

    def get_average_gene_length(self):
        with open(args.outDir+"/"+"gene_length.tsv", "w"):
            pass
        with open(args.outDir+"/"+"cds_length.tsv", "w"):
             pass
        length_distri = list(map(len,self.genes))
        with open(args.outDir+"/"+"gene_length.tsv","a") as genlen:
            genlen.write(("\t"+args.species+"\n").join(["gene_length\t"+str(x) for x in length_distri])+("\t"+args.species+"\n"))
            F.write(("\t"+args.species+"\n").join(["gene_length\t"+str(x) for x in length_distri])+("\t"+args.species+"\n"))
        total_length = sum(length_distri)
        length_cds_distri = list(map(lambda x: sum(list(map(len,x.cds))), self.genes))
        #print(list(map(lambda x: sum(list(map(len,x.cds))), self.genes[0])))
        
        with open(args.outDir+"/"+"cds_length.tsv","a") as genlen:
            genlen.write(("\t"+args.species+"\n").join(["cds_length\t"+str(x) for x in length_cds_distri])+("\t"+args.species+"\n"))
            F.write(("\t"+args.species+"\n").join(["cds_length\t"+str(x) for x in length_cds_distri])+("\t"+args.species+"\n"))
        cds_length = sum(length_cds_distri)
        number_of_genes=self.get_number_of_genes()
        return total_length/(number_of_genes*1.0),cds_length/(number_of_genes*1.0)

    def get_number_of_scaffolds(self):
        return len(self.genomeid_dict)

#     def get_length_and_number_of_exons(self):
#         total = 0
#         number = 0
#         for gene in self.genes:
#             len_exons=sum(map(len,gene.exons))
#             total += len_exons
#             number += len(gene.exons)
#         return total,number

    def get_length_and_number_of_exons(self):
        total = 0
        number = 0
        with open(args.outDir+"/"+"exons_length.tsv", "w"):
            pass
        with open(args.outDir+"/"+"exon_number.tsv", "w"):
            pass   
        for gene in self.genes:
            len_exons_distribution = list(map(len,gene.exons))
            with open(args.outDir+"/"+"exons_length.tsv","a") as inlen:
                inlen.write(("\t"+args.species+"\n").join(["exons_length\t"+str(x) for x in len_exons_distribution])+("\t"+args.species+"\n"))
                F.write(("\t"+args.species+"\n").join(["exons_length\t"+str(x) for x in len_exons_distribution])+("\t"+args.species+"\n"))
            with open(args.outDir+"/"+"exon_number.tsv","a") as inlenn:
                inlenn.write("exon_number\t"+str(len(gene.exons))+"\t"+args.species+"\n")
                F.write("exon_number\t"+str(len(gene.exons))+"\t"+args.species+"\n")
            len_exons=sum(map(len,gene.exons))
            total += len_exons
            number += len(gene.exons)
        return total,number
    def get_length_and_number_of_cds(self):
        total = 0
        number = 0
        with open(args.outDir+"/"+"cds_number.tsv", "w"):
            pass   
        for gene in self.genes:
            len_cds_distribution = list(map(len,gene.cds))
            with open(args.outDir+"/"+"cds_number.tsv","a") as inlenn:
                inlenn.write("cds_number\t"+str(len(gene.cds))+"\t"+args.species+"\n")
                F.write("cds_number\t"+str(len(gene.cds))+"\t"+args.species+"\n")
                
            len_cds=sum(map(len,gene.cds))
            total += len_cds
            number += len(gene.cds)
        return total,number
    def get_number_of_genes_with_introns(self):
        number = 0
        for gene in self.genes:
            if len(gene.introns) > 0 :
                number += 1
                self.genes_with_intron.append(gene)
        return number

    def get_length_and_number_of_introns(self):
        total = 0
        number = 0
     #  empty the file if exist
        with open(args.outDir+"/"+"introns_length.tsv", "w"):
            pass
        for gene in self.genes_with_intron:
            len_introns_distribution = list(map(len,gene.introns))
            with open(args.outDir+"/"+"introns_length.tsv","a") as inlen:
                inlen.write(("\t"+args.species+"\n").join(["introns_length\t"+str(x) for x in len_introns_distribution])+("\t"+args.species+"\n"))
                F.write(("\t"+args.species+"\n").join(["introns_length\t"+str(x) for x in len_introns_distribution])+("\t"+args.species+"\n"))

            len_introns=sum(len_introns_distribution)
            total += len_introns
            number += len(gene.introns)
        return total,number

#     def get_inter_length_and_number(self):
#         total = 0
#         number = 0
#         for inter_per_scaffold in self.intergenicregions_all.values():
#             length_lst = inter_per_scaffold.call_length_inter()
#             total += sum(length_lst)
#             number += len(length_lst)
#         return total,number

    def print_stats(self):
        number_genes= self.get_number_of_genes()
        average_gene_length,average_cds_length = self.get_average_gene_length()
        #number_of_scaffolds = self.get_number_of_scaffolds()
        length_exons,number_exons = self.get_length_and_number_of_exons()
        length_cds,number_cds = self.get_length_and_number_of_cds()
        number_genes_with_introns = self.get_number_of_genes_with_introns()
        length_introns,number_introns = self.get_length_and_number_of_introns()
#        length_inters,number_inters = self.get_inter_length_and_number()
#         donors_stats="/".join(map(lambda x:"%s:%f"%(x,self.cnt_donors[x]/(number_introns*1.0)*100),self.cnt_donors))
#         acceptors_stats="/".join(map(lambda x:"%s:%f"%(x[:2]+"|"+x[2:],self.cnt_acceptors[x]/(number_introns*1.0)*100),self.cnt_acceptors))
        
        with open(args.outDir+"/"+"gene_stats.tsv","w") as f:
            #if(args.soft):f.write("%s \t %s\n"%("Category",args.soft))
            # f.write("Number of genes: \t %d\n"%(number_genes))
            # f.write ("Average gene length: \t %f\n"%(average_gene_length))
            # f.write ("Average cds length: \t %f\n"%(average_cds_length))
            # #f.write ("Number of sequences in genome: \t %d\n"%(number_of_scaffolds))
            # 
            # f.write ("Total number of exons: \t %d\n"%(number_exons))
            # f.write ("Average number of exons per gene: \t %f\n"%(number_exons/(number_genes*1.0)))
            # f.write ("Average exon length: \t %f\n"%(length_exons/(number_exons*1.0)))
            # 
            # f.write ("Number of genes with introns: \t %d\n"%(number_genes_with_introns))
            # f.write ("Total number of introns: \t %d\n"%(number_introns))
            # f.write ("Total intron length: \t %d\n"%(length_introns))
            # f.write ("Average intron length: \t %f\n"%(length_introns/(number_introns*1.0)))
            # f.write ("Average number of introns per gene: \t %f\n"%(number_introns/(number_genes*1.0)))




#             f.write ("Splice donors (%%): \t %s\n"%(donors_stats))
#             f.write ("Splice acceptors (%%): \t %s\n"%(acceptors_stats))

#             f.write ("Number of intergenic regions: \t %d\n"%(number_inters))
#             f.write ("Total intergenic region length:  \t %d\n"%(length_inters))
#             f.write ("Average intergenic region length: \t %f\n"%(length_inters/(number_inters*1.0)))
            if(args.soft):f.write("%s \t %s\n"%("Category",args.soft))
            f.write("Number of genes: \t %s\n"%("{:,d}".format(number_genes)))
            f.write ("Average gene length: \t %s\n"%("{:,.2f}".format(average_gene_length)))
            #f.write ("Average cds length: \t %s\n"%("{:,.2f}".format(average_cds_length)))
            #f.write ("Number of sequences in genome: \t %s\n"%("{:,d}".format(number_of_scaffolds)))
            
            f.write ("Total number of cds: \t %s\n"%("{:,d}".format(number_cds)))
            f.write ("Total cds length: \t %s\n"%("{:,d}".format(length_cds)))
            f.write ("Average cds length: \t %s\n"%("{:,.2f}".format(length_cds/(number_cds*1.0))))
            f.write ("Average number of cds per gene: \t %s\n"%("{:,.2f}".format(number_cds/(number_genes*1.0))))

            f.write ("Total number of exons: \t %s\n"%("{:,d}".format(number_exons)))
            f.write ("Total exon length: \t %s\n"%("{:,d}".format(length_exons)))
            f.write ("Average exon length: \t %s\n"%("{:,.2f}".format(length_exons/(number_exons*1.0))))
            f.write ("Average number of exons per gene: \t %s\n"%("{:,.2f}".format(number_exons/(number_genes*1.0))))

            f.write ("Number of genes with introns: \t %s\n"%("{:,d}".format(number_genes_with_introns)))
            f.write ("Total number of introns: \t %s\n"%("{:,d}".format(number_introns)))
            f.write ("Total intron length: \t %s\n"%("{:,d}".format(length_introns)))
            f.write ("Average intron length: \t %s\n"%("{:,.2f}".format(length_introns/(number_introns*1.0))))
            f.write ("Average number of introns per gene: \t %s\n"%("{:,.2f}".format(number_introns/(number_genes*1.0))))

            #f.write ("Splice donors (%%): \t %s\n"%("{:,.2f}".format(donors_stats)))
            #f.write ("Splice acceptors (%%): \t %s\n"%("{:,.2f}".format(acceptors_stats)))
# 
#             f.write ("Number of intergenic regions: \t %s\n"%("{:,d}".format(number_inters)))
#             f.write ("Total intergenic region length:  \t %s\n"%("{:,d}".format(length_inters)))
#             f.write ("Average intergenic region length: \t %s\n"%("{:,.2f}".format(length_inters/(number_inters*1.0))))
            f.close()


#genome1 = Genome(args.fa)

F=open(args.outDir+"/"+"all_gene_structure.tsv","w")
genome1 = Genome()
genome1.annotate_genes(args.gff)
genome1.print_stats()
F.close()
