# Import libraries
from Bio import AlignIO
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import os
import pandas as pd
import re

# Import de novo assembly data
sqanti_table = pd.read_csv('../flair.collapse.isoforms_classification.filtered_lite_classification.txt', sep='\t', header=0)

# Preprocessing

# Transcript - gene correspondence
manifest = pd.DataFrame(sqanti_table, columns=['isoform', 'associated_gene'])
manifest['isoform'] = manifest['isoform'].str.split('.', expand=True)[0]
manifest['associated_gene'] = manifest['associated_gene'].str.split('.', expand=True)[0]

os.mkdir('SA')
os.mkdir('SB')
# Sequences gathering (needed for next step)
for index, row in manifest.iterrows():
    # Get SA & SB sequences
    # Transcripts sequences
    SB = (r for r in SeqIO.parse('flair.collapse.isoforms.fa', 'fasta') if row['isoform'] in r.id)
    SeqIO.write(SB, 'SB/' + row['isoform'] + '.fa', "fasta")
    # Canonical sequences
    SA = (r for r in SeqIO.parse('Ensembl_Canonical_CDS.fasta', 'fasta') if row['associated_gene'] in r.id)
    SeqIO.write(SA, 'SA/' + row['associated_gene'] + '.fa', "fasta")

# Pairwise nucleotide alignment
os.mkdir('needle')
for index, row in manifest.iterrows():
    afile = 'SA/' + row['associated_gene'] + '.fa'
    bfile = 'SB/' + row['isoform'] + '.fa'
    out_file = 'needle/' + row['associated_gene'] + '_' + row['isoform'] + '.needleall'
    if os.path.exists(afile) and os.path.exists(bfile):
        # NOTE: NeedleCommandline does not support new needle syntax, using os.system()
        cmd = "needle -outfile " + out_file + " -asequence " + afile + " -bsequence " + bfile + " -gapopen 10.0 -gapextend 0.5"
        os.system(cmd)

# Translation using canonical M codon

os.mkdir('SA/Prot')
os.mkdir('SB/Prot')
for nucleotide_alignment_file in os.listdir('needle'):
    canon_name = re.sub(r'_.*$', '', nucleotide_alignment_file)
    isoform_name = re.sub(r'^.*?_', '', nucleotide_alignment_file.rstrip('.needleall'))
    try:
        # Read needle alignments
        alignment = AlignIO.read(str('needle/' + nucleotide_alignment_file), 'emboss')
    except:
        pass
    canon = alignment[0,:].seq
    try:
        # Translate canonical form
        canon_prot = canon.ungap('-').translate(to_stop=True)
        canon_prot = SeqRecord.SeqRecord(canon_prot, id=canon_name)
        SeqIO.write(canon_prot, str('SA/Prot/' + canon_name) + '.fa', 'fasta')
    except:
            pass
    # Use canonical M codon as M codon for all transcripts
    AUG_index = canon.find('ATG')
    isoform = alignment[1,AUG_index:].seq
    # Remove intermediate gaps
    isoform = re.sub('(?<=[A-Z])-+(?=[A-Z])', '', str(isoform))
    # Only translate forms preserving canonical ATG
    if not isoform.startswith('-'):
        try:
            # Translate transcript
            # Some forms do not produce protein
            isoform_prot = SeqRecord.SeqRecord(Seq.Seq(Seq.translate(isoform, to_stop=True)), id=isoform_name)
            SeqIO.write(isoform_prot, str('SB/Prot/' + isoform_name + '.fa'), 'fasta')
        except:
            pass