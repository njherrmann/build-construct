#!/usr/bin/env python

import log


class TargetSequence(object):
  """A data container for a target sequence and its properties."""

  def __init__(self, **kwargs):
    self.sequence = str(kwargs['sequence'])

    # the genomic location of the first base of the sequence - [int]
    self.gen_loc = self._convert_gen_loc(kwargs['gen_loc'])

    self.exon_num = int(kwargs['exon_num'])

    # strand identified by '+' (reads left-right) or '-' (reads right-left)
    self.strand = kwargs['strand']

    self.gc_content = int(kwargs['gc_content'])

    self.offsites = tuple((int(i) for i in kwargs['offsites']))

    self.target = self._compute_target()



  def _convert_gen_loc(self, gen_loc):
    """The gen_loc attribute should be an int but it will probably
     be served as a string in the form 'chr##:#######' with the gen_loc 
     given by the number following the colon."""

    if type(gen_loc) == int:
      return gen_loc

    elif type(gen_loc) == str:
      return int(gen_loc.split(':')[1].strip())


  def _compute_target(self):
    """Computes the genomic location of the cut site (target). The location of the
    cut site is the index of the base preceding the cut on the + strand."""

    if not(self.strand == '+' or self.strand == '-'):
      return
    return self.gen_loc + (16 if self.strand == '+' else 5)


  def bare_sequence(self):
    """Returns the target sequence without the PAM site."""
    return self.sequence[:-3]






class GeneParser(object):
  """An object that parses the ChopChop results.txt file"""

  def __init__(self, filepath=None):
    self.file = None
    self.sequences = []

    if filepath is not None:
      self.read(filepath)

  def read(self, filepath):

    self.sequences = []
    self.file = open(filepath, 'r')

    for i, line in enumerate(self.file):

      if i == 0:  # The first line is a text header
        continue

      # All ChopChop results tables have these columns in order:
      #   Rank, Target_sequence, Genomic_location, Exon, Strand, GC_content,
      #   Self_complementarity, MM0, MM1, MM2, MM3, Efficiency
      # We discard the Rank, Self_complementarity, and Efficiency fields.
      tokens = line.split()
      self.sequences.append(
                      TargetSequence(sequence=tokens[1], gen_loc=tokens[2],
                                     exon_num=tokens[3], strand=tokens[4],
                                     gc_content=tokens[5], offsites=tokens[7:11]))




if __name__ == '__main__':
  gp = GeneParser('NM_005308_results.txt')
  for i in range(10):
    print '{} {}'.format(gp.sequences[i].sequence, gp.sequences[i].gen_loc)