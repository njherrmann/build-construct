#!/usr/bin/env python

import log
import requests
from bs4 import BeautifulSoup


class CcdsLoader(object):
  """An object that loads CCDS pages for proteins and reports
  data such as exon edges and reference genome."""

  # The base URL for the CCDS lookup page.
  # The CCDS code for the desired protein is appended to the URL
  #   WITHOUT the 'CCDS' prefix (already built into the URL)
  URL_BASE = "https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&GO=MainBrowse&DATA=CCDS"

  def __init__(self, ccds_id=None):

    self.logger = log.getLogger(__name__)

    self.ccds_id = None
    self.valid_id = None
    self.soup = None
    self.exon_edges = None

    # init with a protein id is optional
    if ccds_id is not None:
      self.load(ccds_id)


  def load(self, ccds_id):
    self.soup = None
    self.exon_edges = None

    self.ccds_id = ccds_id
    self._clean_id()
    if self.valid_id:
      self._get_soup()
      self.logger.info("CCDS entry loaded")


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
      self.logger.warning("failed to clean CCDS id: invalid format")
      return

    self.ccds_id = new_id
    self.valid_id = True



  def _get_soup(self):
    """Retrieves and sets the soup attr to a BeautifulSoup object for CCDS page."""

    if not self.valid_id:
      self.logger.warning("cannot get soup: invalid CCDS id")

    self.soup = BeautifulSoup(requests.get(CcdsLoader.URL_BASE + self.ccds_id).content, 'lxml')


  def get_exon_edges(self):
    if self.exon_edges is not None:
      return self.exon_edges

    if self.soup is None or not(self.valid_id):
      self.logger.warning("cannot get exon edges: no soup")
      return
    
    extable = self.soup.find('small', text="Chromosome").parent.parent.parent

    # Caches exon_edges until a new load command is executed.
    self.exon_edges = [tuple([int(cell.get_text()) for cell in row.find_all('td')[1:3]])
                            for row in extable.find_all('tr')[1:]]
    return self.exon_edges


if __name__ == "__main__":
  loader = CcdsLoader()
  loader.load('CCDS7612.1')
  print "CCDS7612"
  loader.get_exon_edges()
  print
  print "CCDS2"
  loader.load(2)
  loader.get_exon_edges()