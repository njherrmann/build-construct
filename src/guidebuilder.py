#!/usr/bin/env python

import log
import targetsequence as ts
import guidepair as gp



class GuideBuilder(object):
  """An object that parses the ChopChop results.txt file.
  Requires a settings dict with a filepath in 'input_file':
  Usage:
    builder = GuideBuilder()
    builder.read(settings)
    builder.set_exon_edges(edges)
    builder.build_pairs()
    pairs = builder.get_pairs()
  """

  # A GuideBuilder can be initialized with settings which are used for
  # naked read() command, or a GuideBuilder can be initalized bare
  # then settings are as an argument: read(settings).
  
  # A container for default settings
  _DEFAULTS = { 'max_offtargets'    : (0, 0, 0, 0) ,
                'gRNA2_start_G'     : True ,  # True if seq 2 must start with G
                'separation_limit'  : 10000 , # Max separation between cut sites in bp
                'latest_gRNA2'      : 0.5 ,   # Latest location for gRNA2 as a
                                              #   as a fraction of gene sequence
                'min_exon_deletion' : 0       # min count of exon bps to delete
               }


  def __init__(self, settings=_DEFAULTS):
    
    self.logger = log.getLogger(__name__)

    self.sequences = []
    self.guidepairs = []
    self.exon_edges = None
    self.gene_size = None

    self.settings = settings

    if not('input_file' in self.settings.keys() and
           'CCDS_ID' in self.settings.keys()):
      self.logger.error('Cannot instantiate GuideBuilder without ' +
                        'input_file and CCDS_ID.')

    # Copies unspecified fields to settings from _DEFAULTS
    for key, value in GuideBuilder._DEFAULTS.items():
      if key not in self.settings.keys():
        self.settings[key] = value

    if 'input_file' in self.settings.keys():
      self.read()



  def clear(self):
    """Clears instance attributes."""
    self.sequences = []
    self.guidepairs = []
    self.exon_edges = None
    self.gene_size = None
    self.settings = GuideBuilder._DEFAULTS



  def read(self, filepath=None):
    """Reads and processes the ChopChop results file in filepath.
    Uses self.settings['input_file'] if no filepath argument is given."""

    if filepath is not None:
      self.settings['input_file'] = filepath
    elif 'input_file' not in self.settings.keys():
      self.logger.warning('Cannot read file: specify a filepath')
      return


    self.sequences = []
    self.guidepairs = []

    file = open(self.settings['input_file'], 'r')

    for i, line in enumerate(file):

      if i == 0:  # The first line is a text header
        continue

      # All ChopChop results tables have these columns in order:
      #   Rank, Target_sequence, Genomic_location, Exon, Strand, GC_content,
      #   Self_complementarity, MM0, MM1, MM2, MM3, Efficiency
      # We discard the Rank, Self_complementarity, GC_content and Efficiency fields.
      # We also discard the PAM sites from the target sequences.
      tokens = line.split()
      self.sequences.append(
                   ts.TargetSequence(sequence=tokens[1][:-3], gnm_loc=tokens[2],
                                     exon_num=tokens[3], strand=tokens[4],
                                     offtargets=tokens[7:11]))
      
      # The sequence table is sorted by ascending genomic location
      self.sort_sequences(keystr="cut_site")

      self._filter_offtargets()

      # Filters and sets locations for sequences if exon_edges already set
      if self.exon_edges is not None:
        self.gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0],
                          self.exon_edges, 0)
        self._filter_targets_in_exons()
        for seq in self.sequences:
          seq.set_gene_loc_frac(self.exon_edges, self.settings['strand'])

    self.logger.info('Successfully read input file %s' % self.settings['input_file'].split('/')[-1])

      


  def sort_sequences(self, keystr=None):
    """Sorts the list of sequences in place by a given field. 
    By default (no keystr), sequences are sorted by genomic location.
    If the strand setting is '-', the sequences are sorted in descending order."""

    if keystr is None or keystr == "gnm_loc":
      self.sequences.sort(key=lambda seq: seq.gnm_loc * (-1 if self.settings['strand'] == '-' else 1))
    if keystr == "cut_site":
      self.sequences.sort(key=lambda seq: seq.cut_site * (-1 if self.settings['strand'] == '-' else 1))
    else:
      self.logger.warning("Cannot sort pair list: invalid key string")



  def sort_pairs(self, keystr=None):
    """Sorts the list of GuidePairs in place by a given field. 
    Pairs are sorted by genomic location of seq1 by default (no keystr)"""

    if keystr is None or keystr == "gnm_loc" or keystr == "genomic_location":
      self.guidepairs.sort(key=lambda pair: pair.seq1.gnm_loc * (-1 if self.settings['strand'] == '-' else 1))
    elif keystr == "gen_sep" or keystr == "genomic_separation":
      self.guidepairs.sort(key=lambda pair: pair.genomic_separation)
    elif keystr == "del_count" or keystr == "deletion_count":
      self.guidepairs.sort(key=lambda pair: -1 * pair.deletion_count)
    elif keystr == "del_frac" or keystr == "deletion_fraction":
      self.guidepairs.sort(key=lambda pair: -1 * pair.deletion_fraction)
    else:
      self.logger.warning("Cannot sort pair list: invalid key string")






  def set_exon_edges(self, exon_edges):
    self.exon_edges = exon_edges

    self.gene_size = reduce(lambda tot, edge: tot + edge[1]+1 - edge[0],
                          self.exon_edges, 0)
    
    self._filter_targets_in_exons()
    
    for seq in self.sequences:
      seq.set_gene_loc_frac(self.exon_edges, self.settings['strand'])



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


    # Sorts sequences by location of cut_site
    # Ascending if gene falls on +strand, descending on -strand
    self.sort_sequences(keystr="cut_site")

    # A queue of gRNAs awaiting second sequences to pair
    start_seqs = []

    for seq2 in self.sequences:

      # Iterates across a *copy* of start_seqs so that some sequences
      # can be removed from the original during iteration
      for seq1 in list(start_seqs):
        # Pops start sequences that are too far upstream
        if abs(seq2.cut_site - seq1.cut_site) > self.settings['separation_limit']:
          start_seqs.pop(0)
          continue
        
        # Skips overlapping sequences
        elif seq1.overlap_Q(seq2):
          continue

        else:
          if self.settings['gRNA2_start_G']:
            self.guidepairs += filter(lambda pair: pair.deletion_count >= self.settings['min_exon_deletion'],
                        [gp.GuidePair(seq1, g_seq2, self.exon_edges) for g_seq2 in seq2.find_G_starts()])
          else:
            new_guidepair = gp.GuidePair(seq1, seq2, self.exon_edges)
            if new_guidepair.deletion_count >= self.settings['min_exon_deletion']:
              self.guidepairs.append(new_guidepair)
            else:
              continue

      start_seqs += seq2.find_G_starts()


      if seq2.gene_loc_frac > self.settings['latest_gRNA2']:
        break

    self.sort_pairs(keystr="del_count")
    self.logger.info('Guide pairs compiled')


  def get_pairs(self):
    return self.guidepairs

  def get_sequences(self):
    return self.sequences

  def get_gene_size(self):
    return self.gene_size  
