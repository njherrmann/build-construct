
import log
import guidebuilder as gb
import ccdsloader as ccds



CCDS_ID = 'CCDS7612'
CHOPCHOP_RESULTS_FILE = 'NM_005308_results.txt'



if __name__ == '__main__':
  loader = ccds.CcdsLoader()
  loader.load(CCDS_ID)
  edges = loader.get_exon_edges()
  print 'exon edges'
  for edge in edges:
    print '{} - {}'.format(edge[0], edge[1])
  print

  parser = gb.GuideBuilder()
  parser.read(CHOPCHOP_RESULTS_FILE)

  parser.filter_targets_in_exons(edges)
  print 'first 20 unfiltered sequences:'
  for i in range(20):
    print '[{}] {} {:>20} {} - ({}, {}, {}, {})'.format(parser.sequences[i].strand,
                                 parser.sequences[i].gen_loc,
                                 parser.sequences[i].sequence,
                                 parser.sequences[i].gen_loc + len(parser.sequences[i].sequence),
                                 parser.sequences[i].offtargets[0],
                                 parser.sequences[i].offtargets[1],
                                 parser.sequences[i].offtargets[2],
                                 parser.sequences[i].offtargets[3]
                             )

  parser.filter_offtargets()
  print 'all {} filtered sequences by offtargets:'.format(len(parser.sequences))
  for i in range(len(parser.sequences)):
    print '[{}] {} {:>23} {} - ({}, {}, {}, {})'.format(parser.sequences[i].strand,
                                 parser.sequences[i].gen_loc,
                                 parser.sequences[i].sequence,
                                 parser.sequences[i].gen_loc + len(parser.sequences[i].sequence),
                                 parser.sequences[i].offtargets[0],
                                 parser.sequences[i].offtargets[1],
                                 parser.sequences[i].offtargets[2],
                                 parser.sequences[i].offtargets[3]
                             )
  
  
  


