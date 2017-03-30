
import log
import guidebuilder as gb
import ccdsloader as ccds



SETTINGS = {'CCDS_ID'     : 'CCDS7612' ,
            'input_file'  : 'NM_005308_results.txt' }


if __name__ == '__main__':
  loader = ccds.CcdsLoader()
  loader.load(SETTINGS['CCDS_ID'])
  edges = loader.get_exon_edges()
  
  #print 'exon edges'
  #for edge in edges:
  #  print '{} - {}'.format(edge[0], edge[1])
  #print

  builder = gb.GuideBuilder(SETTINGS)
  builder.read()
  builder.set_exon_edges(edges)
  builder.build_pairs()
  builder.sort_pairs(keystr="deletion_count")
  pairs = builder.get_pairs()
  
  print 'number of pairs: {}'.format(len(pairs))
  print 'genm_loc_1  loc_frac_1  |  genm_loc_2  loc_frac_2  |  separation  |  del_count   del_pct    |    seq_1                 seq_2'
  for i in range(len(pairs)):
 #   print pairs[i].seq2
    print '{:>10}  {:>8.3}    |  {:>10}  {:>8.3}    |  {:>8}    |  {:>7}   {:>8}     |  {:>20}  {:>20}'.format(
                    pairs[i].seq1.gen_loc,
                    pairs[i].seq1.gene_loc_frac,
                    pairs[i].seq2.gen_loc,
                    pairs[i].seq2.gene_loc_frac,
                    pairs[i].seq2.cut_site - pairs[i].seq1.cut_site,
                    pairs[i].deletion_count,
                    int(pairs[i].deletion_fraction * 100),
                    pairs[i].seq1,
                    pairs[i].seq2 )
