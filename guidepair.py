
#import targetsequence as ts


class GuidePair(object):
  """Container for a pair of guide RNAs (each a ts.TargetSequence object)."""

  def __init__(self, seq1, seq2, exon_edges=None):
    self.seq1 = seq1
    self.seq2 = seq2

    # reorders seq1 and seq2 so that seq1 is upstream
    if seq1.cut_site > seq2.cut_site:
      self.seq1, self.seq2 = seq2, seq1


    self.genomic_separation = abs(self.seq2.gen_loc - self.seq1.gen_loc)


    self.deletion_content = None
    self.deletion_fraction = None
    if exon_edges is not None:
      self.compute_deletion_stats(exon_edges)



  def compute_deletion_stats(self, exon_edges):
    """Computes the number of deleted bps and the fraction of 
    the full gene deleted then sets those attributes."""

    # The size of the gene in bps
    gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0], 
                                  exon_edges, 0)

    self.deletion_content = 0
    
    for edge in exon_edges:

      if edge[1] <= self.seq1.cut_site:
        # XXXXX----|---|--
        continue

      elif edge[0] <= self.seq1.cut_site < edge[1] <= self.seq2.cut_site:
        # XXX|XXXXX----|--
        self.deletion_content += (edge[1] - self.seq1.cut_site)

      elif edge[0] <= self.seq1.cut_site and self.seq2.cut_site < edge[1]:
        # XXX|XXXXXX|XXX
        self.deletion_content += (self.seq2.cut_site - self.seq1.cut_site)

      elif self.seq1.cut_site < edge[0] and edge[1] <= self.seq2.cut_site:
        # -|--XXXXXX---|--
        self.deletion_content += (edge[1] - edge[0] + 1)

      elif self.seq1.cut_site < edge[0] <= self.seq2.cut_site < edge[1]:
        # --|---XXXXX|XXX
        self.deletion_content += (self.seq2.cut_site - edge[0] + 1)

      else:   # self.seq2.cut_site < edge[0]
        break
        

    self.deletion_fraction = float(self.deletion_content) / gene_size





