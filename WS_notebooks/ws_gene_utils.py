# ws_gene_utils.py
# WESmith 11/10/22, 12/02/22
# tailored utility functions for bioinformatics


from collections import defaultdict
import gffutils
from Bio import Seq, SeqIO
import matplotlib.pyplot as plt


# WS convenience function 04/24/23
def print_dict(d, leader='', keywidth=8, indent=4, sep=3):
    '''
    Convenience function to print nested dicts.
    d:        dictionary to print (dict)
    leader:   optional initial identifying char or string for each line of output (str)
    keywidth: presumed max text width of keys, so that they print right-aligned (int>0);
    indent:   the indentation value for the next nested dict (int>0)
    sep:      the separation from the key to the value (int>0)
    '''
    keywidthstr = str(int(keywidth))
    txt = '{:>' + keywidthstr + '}:' # right-justify keys
    sepstr    = ' ' * int(sep)
    nextident = ' ' * int(indent)
    for k, v in d.items():
        if isinstance(v, dict):
            print(leader + txt.format(k))
            print_dict(v, leader=leader + nextident, keywidth=keywidth, indent=indent, sep=sep)
        else:
            print(leader + txt.format(str(k)) + sepstr + '{}'.format(v))


def attrs(obj, skip=True, token='__'):
    '''
    Convenience function to print an objects attributes.
    '''
    attr = ['OBJECT TYPE: {}'.format(type(obj))]
    for k in dir(obj):
        if skip and k.__contains__(token): continue
        attr.append(k)
    return attr


def phred_to_percent_accurate(phred):
    return 100 * (1. - 10**(-phred/10))


def amino_histo(seq, order=True):
    num_amino = defaultdict(int)
    for k in seq:
        num_amino[k] += 1
    if order:
        amino_list = list(num_amino.keys())
        amino_list.sort()
        return {k: num_amino[k] for k in amino_list}
    return num_amino


def sequence_from_CDSs(chrom_seq, CDSs, strand, frame=False):
    # CDSs can bea generator or a list
    # frame: True means use frame info from each CDS
    #        False means ignore frame info (default): apparently 
    #        frame info isn't required: in WS tests, only 1 stop
    #        signal produced in protein with False, which is correct; 
    #        get many stop stop with True, which is indicating an 
    #        incorrect framing
    seq = Seq.Seq('')
    for CDS in CDSs:
        off = int(CDS.frame) if frame else 0
        my_cds = Seq.Seq(str(chrom_seq[CDS.start - 1 + off: CDS.end + off]))
        seq += my_cds
        
    return seq if strand == '+' else seq.reverse_complement()


def seq_protein_from_gene(fasta_file, gff_db, gene_id):
    recs = SeqIO.parse(open(fasta_file, 'rt', encoding='utf-8'), 'fasta')
    db   = gffutils.FeatureDB(gff_db)
    gene = db[gene_id]
    
    for rec in recs:
        not_found = True
        if rec.id == gene.seqid:
            not_found = False
            print('\n', rec.description, '//', rec.id, ' strand {}'.\
                  format(gene.strand))
            gene_seq = rec.seq
            break
            
    if not_found:
        raise AssertionError('gene_id not found in fasta file')
        
    mRNAs = db.children(gene, featuretype='mRNA')  # a generator
    
    seqs = defaultdict()
    prot = defaultdict()
    for mRNA in mRNAs:
        CDSs = db.children(mRNA, featuretype='CDS', order_by='start')
        seqs[mRNA.id] = sequence_from_CDSs(gene_seq, CDSs, 
                                           gene.strand, frame=False)
        prot[mRNA.id] = seqs[mRNA.id].translate()
        
    for key, val in prot.items():
        hist = amino_histo(val)
        print('for mRNA {}, stop signal occurs {} time(s) (it should be 1)'.\
              format(key, hist['*']))
        
    return {'seq': seqs, 'prot': prot}


def get_child_of_parent(db, parent, feature):
    return list(db.children(parent, featuretype=feature))


def get_object_coords(obj, y_offset, name='', scale=1): # WS
    out = []
    for k in obj:
        #print(k.id, k.start, k.stop, y_offset)  # WS diagnostic
        out.append([(k.start/scale, k.stop/scale), 
                    (y_offset, y_offset), name, k.frame])
    return out


def plot_mRNA_exon_cds(db, obj, y_offset, scale=1, delta=0.02): # WS
    out = [[(obj.start, obj.stop), (y_offset, y_offset), obj.id, obj.frame]]
    #print('mRNA: ', obj.id, obj.start, obj.stop, y_offset)  # WS diagnostic
    y_offset -= delta

    exon      = get_child_of_parent(db, obj, 'exon')
    cds       = get_child_of_parent(db, obj, 'CDS')
    utr_five  = get_child_of_parent(db, obj, 'five_prime_UTR')
    utr_three = get_child_of_parent(db, obj, 'three_prime_UTR')

    out.extend(get_object_coords(exon, y_offset, scale=scale, 
                                 name='exon'))
    y_offset -= delta
    out.extend(get_object_coords(cds, y_offset, scale=scale, 
                                 name='CDS'))
    y_offset -= delta
    out.extend(get_object_coords(utr_five, y_offset, scale=scale, 
                                 name="5'-UTR"))
    y_offset -= delta
    out.extend(get_object_coords(utr_three, y_offset, scale=scale, 
                                 name="3'-UTR"))
    return out


class PlotGFF():
    def __init__(self, db, gene_id, w=7, y_offset=1, lmarg=0.15, rmarg=0.05, 
                 tmarg=0.1, bmarg=0.1, xwin=16, ywin=None, txtoff=0.01, 
                 scale=1, sep=0.02, mRNA_sep=0.12, gene_sep=0.1):
        '''
        Plot sequence bounds for a gene, its mRNA transcriptions, and their 
        associated exons, CDSs, and UTRs, from a gff version-3 file for the 
        organism of interest.
        
        db  is derived from a gff version 3 file using gffutils:
            eg: db = gffutils.create_db(gff3_file, saved_db) to create 
                     'saved_db' on disk and load into memory as 'db'
                db = gffutils.FeatureDB(saved_db) to load existing 
                     'saved_db' into memory as 'db'
        '''  
        self.db       = db       # gff3 database produced using gffutils
        self.gene_id  = gene_id  # gene ID as per gff convention
        self.w        = w        # width of seq line in vertical direction
        self.y_offset = y_offset # vertical location of primary gene
        self.lmarg    = lmarg    # left   margin of figure
        self.rmarg    = rmarg    # right  margin of figure
        self.tmarg    = tmarg    # top    margin of figure
        self.bmarg    = bmarg    # bottom margin of figure
        self.xwin     = xwin     # width  of figure
        self.ywin     = ywin     # height of figure: if None, auto calculated
        self.txtoff   = txtoff   # offset for inline text before seq graphic
        self.scale    = scale    # defaults to 1, not needed at present
        self.sep      = sep      # sep between exon, CDS, UTRs in mRNA group
        self.mRNA_sep = mRNA_sep # sep between different mRNA groups
        self.gene_sep = gene_sep # separation from chromosome and first mRNA


    def plot_all_sequences(self):
        gene = self.db[self.gene_id]
        self.desc = gene.attributes['description'][0]
        txt = '{}\nstrand: {}'.format(self.gene_id, gene.strand)
        self.out = [[(gene.start/self.scale, gene.stop/self.scale), 
                     (self.y_offset, self.y_offset), txt, gene.frame]]
        #print('gene: ', gene.id, gene.start, gene.stop, y_offset)
        self.y_offset -= self.gene_sep
        mRNAs = get_child_of_parent(self.db, gene, 'mRNA')
        for mRNA in mRNAs:
            self.out.extend(plot_mRNA_exon_cds(self.db, mRNA, 
                                               self.y_offset, 
                                               scale=self.scale, 
                                               delta=self.sep))
            self.y_offset -= self.mRNA_sep   


    def plot(self):
        self.plot_all_sequences()
        xmin, xmax = self.out[0][0]
        ymin = self.out[-1][1][0] - self.bmarg
        ymax = self.out[0][1][0]  + self.tmarg
        xdel = xmax - xmin
        ydel = ymax - ymin
        if self.ywin is None:
            ywin = int(10 * ydel) + 1 # plot window height
        else:
            ywin = self.ywin
        xmin = xmin - self.lmarg * xdel
        xmax = xmax + self.rmarg * xdel

        fig, axs = plt.subplots(1, 1, figsize=(self.xwin, ywin))
        for k in self.out:
            axs.plot(k[0], k[1], linewidth=self.w)
            axs.text(xmin + self.txtoff * xdel, k[1][0], k[2], 
                     fontsize='medium')
            if k[3] != '.':
                axs.text(xmin + 5 * self.txtoff * xdel, k[1][0], 
                         'with frame #', fontsize='medium')
                axs.text(k[0][1], k[1][0] - self.sep/2, k[3], 
                         fontsize='small')

        axs.grid()
        axs.set_xlim(xmin, xmax)
        axs.set_ylim(ymin, ymax)
        txt = 'GENE {}: {}'.format(self.gene_id, self.desc)
        axs.set_title(txt)
        plt.show()