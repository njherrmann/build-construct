#!/usr/bin/env python


import log
import sys
import os
import guidebuilder as gb
import ccdsloader as ccds
import settingsreader as sr
import outputformatter as of



if __name__ == '__main__':

  logger = log.getLogger("main")
  logger.info('Pairing gRNAs in %s' % os.getcwd())

  # Binds project_dir the base project directory
  #   where the gene_block_constants.const file should be stored.
  project_dir = os.path.dirname( os.path.realpath(__file__) )
  if os.path.basename(project_dir) == 'src':
    project_dir = os.path.normpath( os.path.join(project_dir, "..") )


  # The settings file can be specified as a cmd-line input
  # otherwise, the script will search the working dir for
  # a file named settings.inp or gene_block_settings.inp
  settings_file = None
  if len(sys.argv) >= 2:
    settings_file = sys.argv[1]
  elif os.path.exists('gene_block_settings.inp'):
    settings_file = 'gene_block_settings.inp'
  elif os.path.exists('settings.inp'):
    settings_file = 'settings.inp'
  else:
    logger.error('Cannot find a settings file. Specify a file as input or ' +
                  'run pair-guides in a directory with gene_block_settings.inp')
    exit()

  settings = sr.SettingsReader(settings_file)
  logger.info('Reading settings from %s' % settings_file)

  loader = ccds.CcdsLoader()
  loader.load(settings.settings['CCDS_ID'])
  edges = loader.get_exon_edges()
  

  builder = gb.GuideBuilder(settings.settings)
  builder.set_exon_edges(edges)
  builder.build_pairs()
  builder.sort_pairs(keystr="deletion_count")
  pairs = builder.get_pairs()

  outputter = of.OutputFormatter(project_dir + '/gene_block_constants.const')
  outputter.write(pairs, settings.settings['output_file'])