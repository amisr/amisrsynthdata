# create_synth_data.py
# Create an AMISR data file with synthetic data

from .SyntheticData import SyntheticData
from argparse import ArgumentParser, RawDescriptionHelpFormatter


def main():
    config_file_help = 'Some help string'

    # Build the argument parser tree
    parser = ArgumentParser(description=config_file_help,
                            formatter_class=RawDescriptionHelpFormatter)
    arg = parser.add_argument('synth_config_file',help='Configuration file for synthetic data set.')
    args = vars(parser.parse_args())

    sd = SyntheticData(args['synth_config_file'])

if __name__=='__main__':
    main()
