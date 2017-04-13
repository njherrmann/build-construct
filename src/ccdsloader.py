#!/usr/bin/env python

import log
import requests
from bs4 import BeautifulSoup


class CcdsLoader(object):
  """An object that loads CCDS pages for proteins and reports
  data such as exon edges, gene strand, and reference genome."""

  # The base URL for the CCDS lookup page.
  # The CCDS code for the desired protein is appended to the URL
  #   WITHOUT the 'CCDS' prefix (already built into the URL)
  URL_BASE = "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA=CCDS"

  def __init__(self, settings=None):

    self.logger = log.getLogger(__name__)

    self.settings = settings

    self.ccds_id = None
    self.valid_id = None
    self.soup = None
    self.exon_edges = None
    self.strand = None

    # init with a protein id is optional
    if settings is not None:
      self.load(self.settings)


  def load(self, settings):
    self.soup = None
    self.exon_edges = None

    self.settings = settings

    self.ccds_id = self.settings['CCDS_ID']
    self._clean_id()
    if self.valid_id:
      self._get_soup()
      self.logger.info("CCDS%s entry loaded" % self.ccds_id)


  def _clean_id(self):
    """Cleans the ccds_id attr to just the id digits and removes the version suffix."""

    self.valid_id = False

    if self.ccds_id is None:
      self.logger.warning("failed to clean CCDS id: value not set")
      return

    new_id = self.ccds_id

    new_id = str(new_id).strip().lower()

    # Strips everything after the '.' character
    if new_id.find('.') > 0:
      new_id = new_id[ : new_id.find('.')]

    if new_id[:4] == "ccds":
      new_id = new_id[4:]

    if not new_id.isdigit():
      self.logger.warning("Failed to clean CCDS id: invalid format")
      return

    self.ccds_id = new_id
    self.valid_id = True



  def _get_soup(self):
    """Retrieves and sets the soup attr to a BeautifulSoup object for CCDS page."""

    if not self.valid_id:
      self.logger.warning("Cannot get soup: invalid CCDS id")

    self.soup = BeautifulSoup(requests.get(CcdsLoader.URL_BASE + self.ccds_id).content, 'lxml')

    self._add_strand_to_settings()


  def _add_strand_to_settings(self):
    """Gets strand from soup and assigns settings['strand'] to the fetched value."""

    # This is not strictly necessary because get_strand() assigns the local strand attr.
    self.strand = self.get_strand()

    if 'strand' in self.settings.keys() and self.settings['strand'] != self.strand:
      use_ccds_strand = raw_input("\n> Strand specified in settings does not match CCDS info.\n" + \
                                   "> Use '%s' strand from CCDS? [Y/n]\n> " % self.strand)
      print
      if use_ccds_strand.lower() == 'y' or use_ccds_strand.lower() == 'yes':
        self.settings['strand'] = self.strand
        self.logger.info("Overwriting strand to '%s'." % self.strand)
    else:
      self.settings['strand'] = self.strand





  def get_exon_edges(self):
    """Scrapes soup for and reports exon boundaries in a list of tuples."""

    if self.exon_edges is not None:
      return self.exon_edges

    if self.soup is None or not(self.valid_id):
      self.logger.warning("Cannot get exon edges: no soup")
      return
    
    extable = self.soup.find('small', string="Chromosome").parent.parent.parent

    # Caches exon_edges until a new load command is executed.
    self.exon_edges = [tuple([int(cell.get_text()) for cell in row.find_all('td')[1:3]])
                            for row in extable.find_all('tr')[1:]]
    return self.exon_edges



  def get_strand(self):
    """Scrapes soup for and reports the gene strand (either '+' or '-')."""

    # This method is called when the soup is loaded and the return value
    # is assigned to settings['strand']

    if self.strand is not None:
      return self.strand

    if self.soup is None or not(self.valid_id):
      self.logger.warning("Cannot get strand: no soup")
      return

    strandline = filter(lambda s: 'strand of chromosome' in s.lower(), 
                            map(lambda tag: str(tag.get_text()), self.soup.findAll('b')))[0]
    self.strand = strandline.split("'")[1]
    return self.strand





if __name__ == "__main__":
  
  bcor = 'CCDS48093.1'
  grk5 = 'CCDS7612.1'

  loader = CcdsLoader()
  loader.load(bcor)
  
  loader.get_strand()
  print bcor
  print "strand: {}".format(loader.get_strand())
  print "exon edges"
  print loader.get_exon_edges()
  