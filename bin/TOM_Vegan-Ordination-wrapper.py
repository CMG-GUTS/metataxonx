#!/usr/bin/python

"""
TOM_Vegan-Ordination-wrapper.py
---------------------------------

:module:	TOM_Vegan-Ordination-wrapper.py ;

:syntax:	TOM_Vegan-Ordination-wrapper.py --inomics=taxa.xlsx --osheet=genera --inpheno=mapping.xlsx --psheet=mapping --contrast=MyContrast [options] ;

:synopsis:	Performs automated ordination [..] ;

:example:	TOM_Vegan-Ordination-wrapper.py -v --inomics=input.xlsx --osheet=genera --inpheno=input.xlsx --psheet=pheno --log=1 --contrast=\#BBC_BaCo --output=ResultsDir ;

"""

import os
import re
import subprocess
import sys

from optparse import OptionParser


def parse_args():
	"""
    Parse CLI options.

    Returns
    -------
    dict
        User CLI options

    """

	usage = "usage: %prog --inomics=taxa.xlsx --osheet=genera --inpheno=mapping.xlsx --psheet=mapping --contrast=MyContrast [options]"
	required=["inomics", "osheet", "inpheno", "psheet", "contrast"]
	parser = OptionParser(usage)
	parser.add_option('-i', '--inomics', 
		              dest="inomics", 
					  help="Excel with omics data, i.e. Filtered relative abundance feature table(s) of a specific taxonomic level, e.g. 'taxa.xlsx'. [required]",
					  default="taxa.xlsx"
					  )
	parser.add_option('-s', '--osheet', 
		              dest="osheet", 
					  help="Sheet name in the Excel with omics data, e.g. 'genera'. [required]", 
					  default="genera"
					  )
	parser.add_option('-p', '--inpheno', 
		              dest="inpheno", 
					  help="Excel with phenotypic data, i.e. Study data file with sample characteristics and study conditions/contrasts, e.g. 'mapping.xlsx'. [required]",
					  default="mapping.xlsx"
					  )
	parser.add_option('-m', '--psheet', 
		              dest="psheet", 
					  help="Sheet name in the Excel with phenotypic data, e.g. 'mapping'. [required]", 
					  default="mapping"
					  )
	parser.add_option('-c', '--contrast',
		              dest="contrast",
		              help="Header name in the Mapping file of the study contrast to be used for the ordination plot. [required]",
					  default=None
		              )
	parser.add_option('-l', '--log',
		              dest="logtransf",
		              help="This sets the log-transformation (int) of the omics data based on formula [ Y' = log ( A * Y + 1 ) ], \
					  where A is the 'strength' of the log-transformation : 1, 10, 100, 1000, etc., default = 1. [optional]",
					  default=1
		              )
	parser.add_option('-o', '--output',
		              dest="output", 
					  help="Directory to send output to, default is (./) current working directory. [optional]",
					  default="./"
					  )
	parser.add_option('-v', '--verbose',
		              dest="verbose",
		              default=False,
		              action="store_true"
		              )
	options, args = parser.parse_args()
	for r in required:
		if options.__dict__[r] is None:
			parser.print_help()
			parser.error("parameter '%s' is required !"%r)
			sys.exit(1)
			
	if options.verbose:
		print ('VERBOSE :', options.verbose)
		print ('INPUT OMICS FILE :', options.inomics, '( sheet :',options.osheet,')')
		print ('INPUT PHENO FILE :', options.inpheno, '( sheet :',options.psheet,')')
		print ('CONTRAST HEADER :', options.contrast)
		print ('OUTPUT :', options.output)
		verbose = options.verbose
	else:
		verbose = False
	return options, verbose


def run_Vegan( inomics, osheet, inpheno, psheet, contrast, logtransf, output ):

	Rcommand = ( 'Rscript TOM_Vegan-Ordination.R --slave', \
	'--inomics=', inomics, '--osheet=', osheet, '--inpheno=', inpheno, '--psheet=', psheet, '--contrast=', contrast, '--log=', str(logtransf), '--output=', output )
	
	Rcommand = ' '.join(Rcommand)
	Rcommand = Rcommand.replace('= ','=')
	
	print('# R command will be executed :', Rcommand)
	
	subprocess.call(Rcommand, shell=True)



def main():
	options, verbose = parse_args()
	
	run_Vegan( options.inomics, options.osheet, options.inpheno, options.psheet, options.contrast, options.logtransf, options.output )
	
	
main()
