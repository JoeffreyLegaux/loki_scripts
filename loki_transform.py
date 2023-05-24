#!/usr/bin/env python

# (C) Copyright 2018- ECMWF.
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

"""
Loki head script for source-to-source transformations concerning ECMWF
physics, including "Single Column" (SCA) and CLAW transformations.
"""

import sys
from pathlib import Path
import click

from loki import (
    Sourcefile, Transformation, Scheduler, SchedulerConfig,
    Frontend, as_tuple, auto_post_mortem_debugger, flatten, info
)

# Get generalized transformations provided by Loki
from loki.transform import (
    DependencyTransformation, FortranCTransformation, FileWriteTransformation
)

# pylint: disable=wrong-import-order
from transformations.argument_shape import (
    ArgumentArrayShapeAnalysis, ExplicitArgumentArrayShapeTransformation
)
from transformations.derived_types import DerivedTypeArgumentsTransformation

from single_column_coalesced_seq import SingleColumnCoalescedTransformationSeq
# from array_syntax import ExplicitArraySyntaxes
# from logicals_preproc import
from cray_stack import InsertCrayPointers
from mf_transforms import DerivedTypeDimensions, UniformizeLoops, LogicalsPreproc, ExplicitArraySyntaxes, CleanSpec

from loki import Subroutine



@click.group()
@click.option('--debug/--no-debug', default=False, show_default=True,
              help=('Enable / disable debug mode. This automatically attaches '
                    'a debugger when exceptions occur'))
def cli(debug):
    if debug:
        sys.excepthook = auto_post_mortem_debugger


@cli.command()
@click.option('--out-path', '-out', type=click.Path(),
              help='Path for generated souce files.')
@click.option('--path', '-p', type=click.Path(),
              help='Path to search during source exploration.')
@click.option('--header', '-h', type=click.Path(), multiple=True,
              help='Path for additional header file(s).')
@click.option('--cpp/--no-cpp', default=False,
              help='Trigger C-preprocessing of source files.')
@click.option('--include', '-I', type=click.Path(), multiple=True,
              help='Path for additional header file(s)')
@click.option('--define', '-D', multiple=True,
              help='Additional symbol definitions for the C-preprocessor')
@click.option('--omni-include', type=click.Path(), multiple=True,
              help='Additional path for header files, specifically for OMNI')
@click.option('--xmod', '-M', type=click.Path(), multiple=True,
              help='Path for additional .xmod file(s) for OMNI')
@click.option('--data-offload', is_flag=True, default=False,
              help='Run transformation to insert custom data offload regions.')
@click.option('--remove-openmp', is_flag=True, default=False,
              help='Removes existing OpenMP pragmas in "!$loki data" regions.')
@click.option('--mode', '-m', default='sca',
              type=click.Choice(['idem', 'sca', 'claw', 'scc', 'scc-hoist']),
              help='Transformation mode, selecting which code transformations to apply.')
@click.option('--frontend', default='fp', type=click.Choice(['fp', 'ofp', 'omni']),
              help='Frontend parser to use (default FP)')
@click.option('--config', default=None, type=click.Path(),
              help='Path to custom scheduler configuration file')
def convert(out_path, path, header, cpp, include, define, omni_include, xmod,
            data_offload, remove_openmp, mode, frontend, config):
    """
    Single Column Abstraction (SCA): Convert kernel into single-column
    format and adjust driver to apply it over in a horizontal loop.

    Optionally, this can also insert CLAW directives that may be use
    for further downstream transformations.
    """
    if config is None:
        config = SchedulerConfig.from_dict(cloudsc_config)
    else:
        config = SchedulerConfig.from_file(config)

    build_args = {
        'preprocess': cpp,
        'includes': include,
        'defines': define,
        'xmods': xmod,
        'omni_includes': omni_include,
    }

    frontend = Frontend[frontend.upper()]
    frontend_type = Frontend.FP if frontend == Frontend.OMNI else frontend

    # Note, in order to get function inlinig correct, we need full knowledge
    # of any imported symbols and functions. Since we cannot yet retro-fit that
    # after creation, we need to make sure that the order of definitions can
    # be used to create a coherent stack of type definitions.
    definitions = []
    for h in header:
        sfile = Sourcefile.from_file(filename=h, frontend=frontend_type, **build_args)
        definitions = definitions + list(sfile.modules)

    # Create a scheduler to bulk-apply source transformations
    paths = [Path(p).resolve() for p in as_tuple(path)]
    paths += [Path(h).resolve().parent for h in as_tuple(header)]
    scheduler = Scheduler(paths=paths, config=config, frontend=frontend,
                          definitions=definitions, **build_args)


    # First, remove all derived-type arguments; caller first!
    scheduler.process(transformation=DerivedTypeArgumentsTransformation())


    # Now we instantiate our transformation pipeline and apply the main changes
    horizontal = scheduler.config.dimensions['horizontal']
    vertical = scheduler.config.dimensions['vertical']
    block_dim = scheduler.config.dimensions['block_dim']

    scheduler.process(transformation=LogicalsPreproc([], ['LHOOK', 'LMUSCLFA']))

    scheduler.process(ExplicitArraySyntaxes(horizontal=horizontal))
    
    # Temporary tranformation for dimensions that are part of a derived type
    # The member variables are added to the list of routine variables
    # This ensures proper excution of subsequent scripts
    derived_type_dimensions_transformation = DerivedTypeDimensions(horizontal=horizontal, vertical=vertical, block_dim=block_dim)

    scheduler.process(transformation=derived_type_dimensions_transformation)

    scheduler.process(transformation=UniformizeLoops(horizontal=horizontal))

    transformation = SingleColumnCoalescedTransformationSeq(
        horizontal=horizontal, vertical=vertical, block_dim=block_dim,
        directive='openacc', hoist_column_arrays='hoist'
    )

    if transformation:
        scheduler.process(transformation=transformation)
    else:
        raise RuntimeError('[Loki] Convert could not find specified Transformation!')

    # Activate flag for second pass that will remove exactly what has been added
    derived_type_dimensions_transformation.first_pass = False
    scheduler.process(transformation=derived_type_dimensions_transformation)

    scheduler.process(transformation=InsertCrayPointers())    

    # Write out all modified source files into the build directory
    
    scheduler.process(transformation=FileWriteTransformation(builddir=out_path, mode='', suffix='F90'))

    scheduler.process(transformation=CleanSpec())

    scheduler.process(DependencyTransformation(suffix='', mode='strict', include_path=out_path))

if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
