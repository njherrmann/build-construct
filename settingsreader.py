import log


class SettingsReader(object):
  """Reads the settings file for the grna-block-builder."""

  def __init__(self, filepath):

    self.logger = log.getLogger(__name__)

    self.settings = {}
    self.filepath = filepath

    self._read(self.filepath)


  def _read(self, filepath):

    with open(filepath, 'r') as settingsfile:
      
      for line in settingsfile:
        
        # Comments are prefixed by hash marks
        if line[0] == '#':
          continue
        
        tokens = line.strip().split()

        if len(tokens) == 0:
          continue
        if len(tokens) < 2:
          self.logger.warning('No value for %s' % str(tokens[0]))


        if tokens[0].lower() == 'input_file':
          self.settings['input_file'] = tokens[1]

        elif tokens[0].lower() == 'ccds_id':
          self.settings['CCDS_ID'] = tokens[1]

        elif tokens[0].lower() == 'output_file':
          self.settings['output_file'] = tokens[1]

        elif tokens[0].lower() == 'grna2_start_g':
          if (str(tokens[1]).lower() == 't' or str(tokens[1]).lower() == 'true' or
              str(tokens[1]).lower() == '1' or str(tokens[1]).lower() == 'y' or
              str(tokens[1]).lower() == 'yes'):
            self.settings['gRNA2_start_G'] = True
          elif (str(tokens[1]).lower() == 'f' or str(tokens[1]).lower() == 'false' or
                str(tokens[1]).lower() == '0' or str(tokens[1]).lower() == 'n' or
                str(tokens[1]).lower() == 'no'):
            self.settings['gRNA2_start_G'] = False
          else:
            self.logger.warning('Value for %s must be True or False.' % tokens[0])

        elif tokens[0].lower() == 'separation_limit':
          self.settings['separation_limit'] = int(tokens[1]) * 1000 # converts to kbp

        elif tokens[0].lower() == 'latest_grna2':
          self.settings['latest_gRNA2'] = float(tokens[1])

        elif tokens[0].lower() == 'min_exon_deletion':
          self.settings['min_exon_deletion'] = int(tokens[1])

        elif tokens[0].lower() == 'max_offtargets':
          self.settings['max_offtargets'] = tuple(map(int, tokens[1:5]))

        else:
          self.logger.warning('Invalid key: %s' % str(tokens[0]))


    if 'output_file' not in self.settings.keys():
      self.settings['output_file'] = ".".join(self.settings['input_file'].split('.')[:-1]) \
                                             + '_constructs.csv'