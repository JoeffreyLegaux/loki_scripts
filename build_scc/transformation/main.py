# (C) Copyright 2023- ECMWF.
# (C) Copyright 2023- Meteo-France.

from openacc_setup import call_openacc_trans

import click

@click.command()
#@click.option('--pathr', help='path of the file to open')
#@click.option('--pathw', help='path of the file to write to')
@click.option('--pathpack', help='absolute path to the pack')
@click.option('--pathview', help='path to src/local/... or src/main/...')
@click.option('--pathfile', help='path to the file, with the file name in the path')
@click.option('--pathacc', help='path to the place where acc files are stored')

@click.option('--horizontal_opt', default=None, help='additionnal possible horizontal idx')
@click.option('--inlined', '-in', default=None, multiple=True, help='names of the routine to inline')

def main_function(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined):
   call_openacc_trans(pathpack, pathview, pathfile, pathacc, horizontal_opt, inlined)

main_function()

