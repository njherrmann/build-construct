#!/usr/bin/env python

import log
import targetsequence as ts
import guidepair as gp



class GuideBuilder(object):
  """An object that parses the ChopChop results.txt file"""
  
  # A container for default settings
  DEFAULTS = { 'max_offtargets'   : (0, 0, 0, 0) ,
               'gRNA2_start_G'    : True ,  # True if seq 2 must start with G
               'separation_limit' : 10000 , # Max separation between cut sites in bp
               'latest_gRNA1'     : 0.5     # Latest location for gRNA1 as a
                                            #   as a fraction of gene sequence
               }


  def __init__(self, settings=DEFAULTS):
    
    self.logger = log.getLogger(__name__)

    self.file = None
    self.sequences = []
    self.guidepairs = []
    self.exon_edges = None
    self.gene_size = None

    self.settings = settings

    # Copies unspecified fields to settings from DEFAULTS
    for key, value in DEFAULTS.items():
      if key not self.settings.keys():
        self.settings[key] = value

    if 'input_path' in self.settings.keys():
      self.read()



  def read(self, filepath=None):
    """Reads and processes the ChopChop results file in filepath.
    Uses self.settings['filepath'] if no filepath argument is given."""

    if filepath is not None:
      self.settings['input_path'] = filepath
    elif 'input_path' not in self.settings.keys():
      self.logger.warning('Cannot read file: specify a filepath')
      return


    self.sequences = []
    self.guidepairs = []
    self.exon_edges = exon_edges


    self.file = open(self.settings['input_path'], 'r')

    for i, line in enumerate(self.file):

      if i == 0:  # The first line is a text header
        continue

      # All ChopChop results tables have these columns in order:
      #   Rank, Target_sequence, Genomic_location, Exon, Strand, GC_content,
      #   Self_complementarity, MM0, MM1, MM2, MM3, Efficiency
      # We discard the Rank, Self_complementarity, GC_content and Efficiency fields.
      # We also discard the PAM sites from the target sequences.
      tokens = line.split()
      self.sequences.append(
                   ts.TargetSequence(sequence=tokens[1][:-3], gen_loc=tokens[2],
                                     exon_num=tokens[3], strand=tokens[4],
                                     offtargets=tokens[7:11]))
      
      # The sequence table is sorted by ascending genomic location
      self.sort_sequences(self.sequences, keystr="gen_loc")

      self._filter_offtargets()

      if self.exon_edges is not None:
        self.gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0],
                          self.exon_edges, 0)
        self._filter_targets_in_exons()
        for seq in self.sequences:
          seq.set_gene_loc_frac(self.exon_edges)

      




  def sort_sequences(self, sequencelist, keystr=None):
    """Sorts the list of sequences in place by a given field. 
    Sequences are sorted by ascending genomic location by default (no keystr)"""

    if keystr is None or keystr == "gen_loc":
      sequencelist.sort(key=lambda seq: seq.gen_loc)
    else:
      self.logger.warning("Cannot sort pair list: invalid key string")



  def sort_pairs(self, pairlist, keystr=None):
    """Sorts the list of GuidePairs in place by a given field. 
    Pairs are sorted by ascending genomic location of seq1 by default (no keystr)"""

    if keystr is None or keystr == "gen_loc" or keystr == "genomic_location":
      pairlist.sort(key=lambda pair: pair.seq1.gen_loc)
    elif keystr == "gen_sep" or keystr == "genomic_separation":
      pairlist.sort(key=lambda pair: pair.genomic_separation)
    elif keystr == "del_cont" or keystr == "deletion_content":
      pairlist.sort(key=lambda pair: pair.deletion_content)
    elif keystr == "del_frac" or keystr == "deletion_fraction":
      pairlist.sort(key=lambda pair: pair.deletion_fraction)
    else:
      self.logger.warning("Cannot sort pair list: invalid key string")






  def set_exon_edges(self, exon_edges):
    self.exon_edges = exon_edges

    self.gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0],
                          self.exon_edges, 0)
    
    self._filter_targets_in_exons()
    
    for seq in self.sequences:
      seq.set_gene_loc_frac(self.exon_edges)



  def _filter_targets_in_exons(self):
    """Removes the sequences for which the target site falls outside an exon."""

    # The exon_edges input must be a list of tuples which mark the exon boundaries.
    # The boundaries for an exon are the indices of the first and last base pair
    #   of that exon (INCLUSIVE)

    if self.exon_edges is None:
      self.logger.warning("Cannot filter targets: assign exon edges first")
      return

    # A target site is indexed by the preceding base on the + strand.
    #   i.e. for a target site i, the cut falls between bases i and i+1.
    # These computations are done by each TargetSequence instance
    self.sequences = filter(lambda seq: seq.cut_in_range(self.exon_edges),
                            self.sequences)



  def _filter_offtargets(self, max_offtargets=None):
    """Removes sequences for which any offsite count exceeds the given 
    max value. Filters to sequences with no offtargets by default."""

    if max_offtargets is None:
      max_offtargets = self.settings['max_offtargets']

    if len(max_offtargets) != 4:
      self.logger.warning("Cannot filter targets: specify all four max offsite values")
      return

    self.sequences = filter(lambda seq: all([seq.offtargets[i] <= max_offtargets[i] for i in range(4)]),
                            self.sequences)


  

  def build_pairs(self):
    """Builds GuidePair objects for all valid sequence pairs."""

    if self.exon_edges is None:
      self.logger.warning("Cannot filter targets: assign exon edges first")
      return

    self.guidepairs = []


    # Sorts sequences by ascending genomic location
    self.sort_sequences(self.sequences, keystr="gen_loc")

    # A queue of gRNAs awaiting second sequences to pair
    start_seqs = []

    for seq in self.sequences:

      # Iterates across a *copy* of start_seqs so that some sequences
      # can be removed from the original during iteration
      for start in list(start_seqs):
        # Pops start sequences that are too far upstream
        if abs(seq.cut_site - start.cut_site) > self.settings['separation_limit']:
          start_seqs.pop(0)
        
        else:
          if start.cut_site == seq.cut_site:
            continue
          if self.settings['gRNA2_start_G']:
            self.guidepairs.append(gp.GuidePair(start, seq))
          else:
            self.guidepairs += [gp.GuidePair(start, g_seq) 
                                  for g_seq in seq.find_G_starts()]

        start_seq_queue += seq.find_G_starts()


      if seq.gene_loc_frac > self.settings['latest_gRNA1']:
        break








if __name__ == '__main__':
  gp = BlockBuilder()
  gp.read('NM_005308_results.txt')
  print 'first 10 unfiltered sequences:'
  for i in range(10):
    print '{} ({}) {}'.format(gp.sequences[i].gen_loc%10000, gp.sequences[i].strand, gp.sequences[i].sequence)

  