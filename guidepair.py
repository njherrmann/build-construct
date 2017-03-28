
#import targetsequence as ts


class GuidePair(object):
  """Container for a pair of guide RNAs (each a ts.TargetSequence object)."""

  def __init__(self, seq1, seq2):
    self.seq1 = seq1
    self.seq2 = seq2

