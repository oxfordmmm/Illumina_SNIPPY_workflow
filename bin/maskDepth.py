#!/usr/bin/env python3 
import sys
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import pandas as pd

class depthMasker:
    def __init__(self, fasta, output, pysamstats, depth, mw,mask='N', renamer=None):
        self.fasta = fasta
        self.output = output
        self.pysamstats = pysamstats
        self.depth = depth
        self.mask = mask
        self.maskWeak = mw
        if renamer != None:
            self.renamer={i.split(':')[0]:i.split(':')[1] for i in renamer}
        else:
            self.renamer={}

    def run(self):
        self.loadDepth()
        self.loadFasta()
        self.maskDepth()
        self.makeOutput()
        self.depthStats()

    def ps(self,row):
        if row['seq']=='N':
            x=None
        else:
            x = float(row[row['seq']+'_PS'])
        return x

    def depthStats(self):
        avDepth=self.depthDF.groupby('chrom')['reads_all'].mean()
        medianDepth=self.depthDF.groupby('chrom')['reads_all'].median()
        # calculate the number of position above 10x
        tenX=self.depthDF[self.depthDF['reads_all']>=10].groupby('chrom')['reads_all'].count()

        avDepth.rename("Average_Depth", inplace=True)
        medianDepth.rename("Median_Depth", inplace=True)
        tenX.rename("10x coverage breadth", inplace=True)
        print(avDepth)
        print(medianDepth)

        df=pd.concat([avDepth,medianDepth,tenX],axis=1)
        df.to_csv('depthStats.csv')


    def loadDepth(self):
        df=pd.read_csv(self.pysamstats,sep='\t')
        df=df[['chrom','pos','reads_all','A','T','C','G']]
        df['total']=df['A']+df['C']+df['G']+df['T']
        bases=['A','T','C','G']
        for b in bases:
            df['{0}_PS'.format(b)]=(df[b]/df['total'])
        self.depthDF=df

    def loadFasta(self):
        dfs=[]
        for seq in SeqIO.parse(open(self.fasta,'rt'),'fasta'):
            l=[seq.seq[i] for i in range(0,len(seq.seq))]
            df=pd.DataFrame( {'pos' : list( range( 1, len( seq.seq ) + 1 )), 'seq' : l } )
            df['chrom']=seq.id
            dfs.append(df)
        self.seqsDF=pd.concat(dfs)

    def maskDepth(self):
        df=self.seqsDF.merge(self.depthDF,left_on=['chrom','pos'],right_on=['chrom','pos'],how='left')
        df['ps']=df.apply(self.ps, axis=1)
        df['seq']=np.where(df['reads_all']>=int(self.depth), df['seq'], self.mask)
        if self.maskWeak != False:
            df['seq']=np.where(df['ps']>=float(self.maskWeak), df['seq'], self.mask)
        print(df[df['ps']<0.80])
        self.df=df

    def _getSeqs(self):
        chroms=list(self.df['chrom'].unique())
        for chrom in chroms:
            df=self.df[self.df['chrom']==chrom]
            seq=''.join(list(df['seq']))
            if chrom in self.renamer:
                chrom=chrom.replace(chrom,self.renamer[chrom])
            yield SeqRecord(Seq(seq), 
                id=chrom, 
                description='mininum_depth={0} mask_character={1}'.format(self.depth, self.mask))

    def makeOutput(self):
        _seqs=self._getSeqs()
        SeqIO.write(_seqs, self.output, 'fasta')

def run(opts):
    dm=depthMasker(opts.fasta,opts.output, opts.pysamstats, opts.depth,opts.maskWeak,renamer=opts.renamer)
    dm.run()

if __name__ == "__main__":
    parser = ArgumentParser(description='Mask positions below depth (d) with Ns')
    parser.add_argument('-p', '--pysamstats', required=True,
            help='pysam stats file')
    parser.add_argument('-f', '--fasta', required=True,
            help='input fasta')
    parser.add_argument('-o', '--output', required=True,
            help='output fasta')
    parser.add_argument('-d', '--depth', required=True,
            help='depth below which to mask with Ns')
    parser.add_argument('-r', '--renamer', required=False,nargs='+',
            help='sequence id to rename:name to give')
    parser.add_argument('-mw','--maskWeak', required=False, default=False,
            help='mask weakly supported positions from pysam (float 0-1) default=False')
    opts, unknown_args = parser.parse_known_args()
    run(opts)
