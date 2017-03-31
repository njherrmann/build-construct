#!/usr/bin/env python


import log
import guidebuilder as gb
import ccdsloader as ccds
import settingsreader as sr
import outputformatter as of


SETTINGS = {'CCDS_ID'     : 'CCDS7612' ,
            'input_file'  : 'NM_005308_results.txt' }


if __name__ == '__main__':
  
  settings = sr.SettingsReader('settings.inp')

  loader = ccds.CcdsLoader()
  loader.load(settings.settings['CCDS_ID'])
  edges = loader.get_exon_edges()
  
  #print 'exon edges'
  #for edge in edges:
  #  print '{} - {}'.format(edge[0], edge[1])
  #print

  builder = gb.GuideBuilder(settings.settings)
  builder.set_exon_edges(edges)
  builder.build_pairs()
  builder.sort_pairs(keystr="deletion_count")
  pairs = builder.get_pairs()

  outputter = of.OutputFormatter()
  outputter.write(pairs, settings.settings['output_file'])
    
#  print 'number of pairs: {}'.format(len(pairs))
#  print ' n  |  genm_loc_1  loc_frac_1  |  genm_loc_2  loc_frac_2  |  separation  |  del_count   del_pct    |    seq_1                 seq_2'
#  for i in range(len(pairs)):
#    print pairs[i].seq2
#    print '{:>2}  |  {:>10}  {:>8.3}    |  {:>10}  {:>8.3}    |  {:>8}    |  {:>7}   {:>8}     |  {:>20}  {:>20}'.format(
#                    i+1,
#                    pairs[i].seq1.gnm_loc,
#                    pairs[i].seq1.gene_loc_frac,
#                    pairs[i].seq2.gnm_loc,
#                    pairs[i].seq2.gene_loc_frac,
#                    pairs[i].seq2.cut_site - pairs[i].seq1.cut_site,
#                    pairs[i].deletion_count,
#                    pairs[i].deletion_pct,
#                    pairs[i].seq1,
#                    pairs[i].seq2 )
