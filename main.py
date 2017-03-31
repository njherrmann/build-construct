#!/usr/bin/env python


import log
import sys
import os
import guidebuilder as gb
import ccdsloader as ccds
import settingsreader as sr
import outputformatter as of


SETTINGS = {'CCDS_ID'     : 'CCDS7612' ,
            'input_file'  : 'NM_005308_results.txt' }


if __name__ == '__main__':

  logger = log.getLogger("main")
  logger.info('Building constructs in %s' % os.getcwd())


  settings = sr.SettingsReader(sys.argv[1])
  logger.info('Reading settings from %s' % sys.argv[1])

  loader = ccds.CcdsLoader()
  loader.load(settings.settings['CCDS_ID'])
  edges = loader.get_exon_edges()
  

  builder = gb.GuideBuilder(settings.settings)
  builder.set_exon_edges(edges)
  builder.build_pairs()
  builder.sort_pairs(keystr="deletion_count")
  pairs = builder.get_pairs()

  outputter = of.OutputFormatter()
  outputter.write(pairs, settings.settings['output_file'])