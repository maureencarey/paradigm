import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import logging
import argparse
from datetime import datetime

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='check_logs{}.log'.format(day), level=logging.INFO) # ADDING TO LOG FILE
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description='Read in the log')
parser.add_argument('log_file')
args = parser.parse_args()
log_fname = vars(args)['log_file']

# parse arguments for global variables
INFO_name = log_fname

logging.info('these steps were infeasible:')

with open(log_fname) as f:
    prev_line = ' '
    for i, line in enumerate(f):
        if 'Infeasible' in line or 'infeasible' in line:
            logging.info('infeasible step in {}'.format(log_fname))
            logging.info(prev_line)
            logging.info(line)
            logging.info('__________________________________')
        if 'error' in line or 'Error' in line:
            logging.info('error message in {}'.format(log_fname))
            logging.info(prev_line)
            logging.info(line)
            logging.info('__________________________________')
        prev_line = line

