import argparse
import src.tool_functions as dlbclass
import warnings
warnings.filterwarnings('ignore')


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--id',
                    help='Cohort set name.',
                    required=True, type=str)
parser.add_argument('-g', '--gsm_file',
                    help='The gsm file to classify samples from. ',
                    required=True, type=str)
parser.add_argument('-o', '--output_dir',
                    help='Output directory.',
                    required=False, type=str, default='./')

args = parser.parse_args()

dlbclass.classify_samples(args.gsm_file, args.id, args.output_dir)
