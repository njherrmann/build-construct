import log


class TargetSequence(object):
  """A data container for a target sequence and its properties."""

  def __init__(self, **kwargs):
    self.logger = log.getLogger(__name__)

    if 'sequence' not in kwargs.keys():
      self.logger 

    self.sequence = str(kwargs['sequence'])

    # the genomic location of the first base of the sequence - [int]
    self.gen_loc = self._convert_gen_loc(kwargs['gen_loc'])

    self.exon_num = int(kwargs['exon_num'])

    # strand identified by '+' (reads left-right) or '-' (reads right-left)
    self.strand = kwargs['strand']

    self.offtargets = tuple((int(i) for i in kwargs['offtargets']))

    self.gc_content = self._compute_gc_content()

    self.target = self._compute_target()

    self.contains_seed_T = self._check_seed_T()



  def _convert_gen_loc(self, gen_loc):
    """The gen_loc attribute should be an int but it will probably
     be served as a string in the form 'chr##:#######' with the gen_loc 
     given by the number following the colon."""

    if type(gen_loc) == int:
      return gen_loc

    elif type(gen_loc) == str:
      return int(gen_loc.split(':')[1].strip())


  def _compute_gc_content(self):
    """Sets the gc_content member to the integer-rounded percentage
    GC content of the sequence."""

    # This calculation includes the PAM site. Is this correct?
    return int(round(
            100*(self.sequence.count('G') + self.sequence.count('C'))/len(self.sequence)))


  def _compute_target(self):
    """Computes the genomic location of the cut site (target). The location of the
    cut site is the index of the base preceding the cut on the + strand."""

    if not(self.strand == '+' or self.strand == '-'):
      return
    return self.gen_loc + (len(self.sequence)-4 if self.strand == '+' else 2)


  def _check_seed_T(self):
    """Returns True if there is at least one T in the seed region
    (4bps before the PAM site)"""

    return self.sequence[-7 : -3].count('T') > 0


  def get_bare_sequence(self):
    """Returns the target sequence without the PAM site."""
    return self.sequence[:-3]


  def truncate_front(self, n_trunc):
    """Returns a TargetSequence of this sequence with 
    the first n_trunc characters removed."""

    if n_trunc == 0:
      return self

    if len(self.sequence) - n_trunc < 18:
      self.logger.warning('Cannot truncate: resulting guide would have fewer than 18 bases.')
      return

    return TargetSequence(sequence=self.sequence[n_trunc : ], 
                          gen_loc=self.gen_loc + (n_trunc if self.strand == '+' else 0),
                          exon_num=self.exon_num,
                          strand=self.strand,
                          offtargets=self.offtargets)


  def find_start_guides(self):
    """Returns a list containing all possible starting guides
    derived from this sequence."""

    guides = []

    for i in xrange(2):
      if self.sequence[i] == 'G':
        guides.append(self.truncate_front(i))
        # NB: truncate_front(0) returns self

    guides = filter(lambda seq: 40 <= seq.gc_content <= 70, guides)

    return guides


  def cut_in_range(self, boundarylist):
    """Returns True if the target site of this sequence falls
    between any of the given boundaries (each boundary a two-element list)."""


    for bnd in boundarylist:
      # Each boundary is a 2-element list or tuple of ints
      if len(bnd) < 2:
        self.logger.warning('Cannot validify cut site: improper boundary shape')
        continue
      if bnd[0] >= bnd[1]:
        self.logger.error('Cannot validify cut site: invalid boundary specification')
        continue

      # Returns true if the target site is to the right of the first boundary
      #   or to the right of the penultimate boundary.
      if (bnd[0] <= self.target <= bnd[1]-1):
        return True

    return False