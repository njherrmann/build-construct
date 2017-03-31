
import log
import csv


class OutputFormatter(object):
  """Writes potential construct data to a csv file."""

  def __init__(self):
    self.logger = log.getLogger(__name__)

    self.gene_block_constants = {}
    self._read_gene_block_constants('gene_block_constants.dat')



  def _read_gene_block_constants(self, constsfilepath):
    """Reads constant elements of the gene block construct."""
    
    with open(constsfilepath, 'r') as constsfile:
      for line in constsfile:
        tokens = line.strip().split()
        self.gene_block_constants[tokens[0]] = tokens[1]



  def _assemble_construct(self, guidepair):
    """Assembles the gene block constants and the sequences into a full construct string."""

    return self.gene_block_constants['u6prom__cln_hom_arm'] + \
              guidepair.seq1.sequence + \
              self.gene_block_constants['crispr_scaff_2'] + \
              self.gene_block_constants['csy4_clvg'] + \
              guidepair.seq2.sequence + \
              self.gene_block_constants['crispr_scaff_1__cln_hom_arm']



  def write(self, guidepairs, outfilepath):

    with open(outfilepath, 'w') as constsfile:
      writer = csv.writer(constsfile, dialect='excel')

      writer.writerow(['gRNA 1', 'genomic loc 1', 
                        'gene loc frac 1', 'exon 1',
                        'strand 1', 'GC% 1', 
                        'off-targets 1', 
                        '',
                       'gRNA 2', 'genomic loc 2',
                        'gene loc frac 2', 'exon 2',
                        'strand 2', 'GC% 2', 
                        'off-targets 2', 
                        '',
                       'gRNA separation (bp)', 
                       'del count (bp)', 'del%', 
                       '',
                       'full construct'])

      for pair in guidepairs:
        writer.writerow([pair.seq1.sequence, pair.seq1.gnm_loc,
                          pair.seq1.gene_loc_frac, pair.seq1.exon_num,
                          pair.seq1.strand, pair.seq1.gc_content, 
                          pair.seq1.offtargets,
                         '',
                         pair.seq2.sequence, pair.seq2.gnm_loc,
                          pair.seq2.gene_loc_frac, pair.seq2.exon_num,
                          pair.seq2.strand, pair.seq2.gc_content, 
                          pair.seq2.offtargets,
                         '',
                         pair.genomic_separation,
                         pair.deletion_count,
                         pair.deletion_pct,
                         '',
                         self._assemble_construct(pair)])

    self.logger.info('Successfully wrote constructs to %s' % outfilepath)