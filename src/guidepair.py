
#import targetsequence as ts


class GuidePair(object):
  """Container for a pair of guide RNAs (each a ts.TargetSequence object)."""

  def __init__(self, seq1, seq2, exon_edges=None):
    self.seq1 = seq1
    self.seq2 = seq2

    self.genomic_separation = abs(self.seq2.cut_site - self.seq1.cut_site)


    self.deletion_count = None
    self.deletion_fraction = None
    self.deletion_pct = None
    if exon_edges is not None:
      self.compute_deletion_stats(exon_edges)



  def compute_deletion_stats(self, exon_edges):
    """Computes the number of deleted bps and the fraction of 
    the full gene deleted then sets those attributes."""

    # The size of the gene in bps
    gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0], 
                                  exon_edges, 0)

    self.deletion_count = 0

    left_cut = min(self.seq1.cut_site, self.seq2.cut_site)
    right_cut = max(self.seq1.cut_site, self.seq2.cut_site)
    
    for edge in exon_edges:

      if edge[1] <= left_cut:
        # XXXXX----|---|--
        continue

      elif edge[0] <= left_cut < edge[1] <= right_cut:
        # XXX|XXXXX----|--
        self.deletion_count += (edge[1] - left_cut)

      elif edge[0] <= left_cut and right_cut < edge[1]:
        # XXX|XXXXXX|XXX
        self.deletion_count += (right_cut - left_cut)

      elif left_cut < edge[0] and edge[1] <= right_cut:
        # -|--XXXXXX---|--
        self.deletion_count += (edge[1] - edge[0] + 1)

      elif left_cut < edge[0] <= right_cut < edge[1]:
        # --|---XXXXX|XXX
        self.deletion_count += (right_cut - edge[0] + 1)

      else:   # self.seq2.cut_site < edge[0]
        break


    self.deletion_fraction = float(self.deletion_count) / gene_size
    self.deletion_pct = int(100 * self.deletion_fraction)





