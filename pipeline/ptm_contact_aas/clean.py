#!/usr/bin/env python3
"""
Clean: Clean up temporary, input, and output files as desired
"""

# Initialize

import os
from blang import Run

# Start

# Safety prompt (input)

inclean = 'disabled'
# inclean = 'enabled'
# inclean = 'always'

if os.path.isdir('./input/') and (inclean != 'disabled'):

    print("\nAlso remove input files? (y/n)")

    i = input()
    if i.lower() == 'y':
        inclean = 'all'

    if inclean in ('all', 'always'):
        Run("Clean up input files", "rm -f input-*")
        
        Run("Clean up input directory", "mv ./input ./deleting_input")
        Run("Clean up input directory", "rm -rf ./deleting_input &")
        Run("Recreate input directory", "mkdir ./input")

# Safety prompt (output)

# outclean = 'disabled'
outclean = 'some'
# outclean = 'all'

if os.path.isdir('./output/'):
    print("\nAlso remove output files? (y/n)")
    
    i = input()
    if i.lower() == 'y':
        outclean = 'all'

# Clean

Run("Clean up temporary files", "rm -f tmp-*")

if os.path.isdir('./tmp/'):
    Run("Clean up temporary directory", "mv ./tmp ./deleting_tmp")
    Run("Clean up temporary directory", "rm -rf ./deleting_tmp &")
    Run("Recreate temporary directory", "mkdir ./tmp")


if os.path.isdir('./output/'):
    if outclean == 'all':
        Run("Clean up output files", "rm -f output-*")
        
        Run("Clean up output directory", "mv ./output ./deleting_output")
        Run("Clean up output directory", "rm -rf ./deleting_output &")
        Run("Recreate output directory", "mkdir ./output")
else:
    Run("Clean up output files", "rm -f output-*")

# print("Done!")
