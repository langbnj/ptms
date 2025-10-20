#!/usr/bin/env python3 -u
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import Run, Waitforjobs

# Start

# Disorder calls (type='AlphaFold')
# # Run("Disorder without smoothing", "~/scripts/qsub.sh disorder.py dis -humanonly")
# # Run("Disorder based on ±10 aa smoothing", "~/scripts/qsub.sh disorder.py dis10 -humanonly")
# # Run("Disorder based on ±10 aa smoothing (with pLDDT < 70 as blank)", "~/scripts/qsub.sh disorder.py dis10 -humanonly -plddt70")
# # Run("Disorder based on ±10 aa smoothing (with pLDDT < 90 as blank)", "~/scripts/qsub.sh disorder.py dis10 -humanonly -plddt90")
# Run("Disorder based on ±10 aa smoothing (with pLDDT < 70 or PAE > 2 as blank)", "~/scripts/qsub.sh disorder.py dis10 -humanonly -plddt70 -pae2")
# Run("Disorder based on ±10 aa smoothing (with pLDDT < 90 or PAE > 1 as blank)", "~/scripts/qsub.sh disorder.py dis10 -humanonly -plddt90 -pae1")
Run("Disorder based on ±10 aa smoothing (with pLDDT < 70 or PAE > 2 as blank)", "disorder.py dis10 -humanonly -plddt70 -pae2")
Run("Disorder based on ±10 aa smoothing (with pLDDT < 90 or PAE > 1 as blank)", "disorder.py dis10 -humanonly -plddt90 -pae1")

# Note: The subsequent scripts (e.g. pipeline/evolutionary_rate_analysis and snps_contingency_coresurf) are robust to having values other than 'C' or 'S' in the 'CoreSurf' sequence from uniseq. Can safely use 'M' rather than ' '.
# Core/Surface calls (type='CoreSurf')
# # Run("Core/Surface without smoothing", "~/scripts/qsub.sh coresurf.py surf -humanonly")
# # Run("Core/Surface without smoothing (with membrane as 'M')", "~/scripts/qsub.sh coresurf.py surf -humanonly")
# # Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 70 as blank)", "~/scripts/qsub.sh coresurf.py surf -humanonly -plddt70")
# # Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 90 as blank)", "~/scripts/qsub.sh coresurf.py surf -humanonly -plddt90")
# Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 70 or PAE > 2 as blank)", "~/scripts/qsub.sh coresurf.py surf -humanonly -plddt70 -pae2")
# Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 90 or PAE > 1 as blank)", "~/scripts/qsub.sh coresurf.py surf -humanonly -plddt90 -pae1")
Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 70 or PAE > 2 as blank)", "coresurf.py surf -humanonly -plddt70 -pae2")
Run("Core/Surface without smoothing (with membrane as 'M' and with pLDDT < 90 or PAE > 1 as blank)", "coresurf.py surf -humanonly -plddt90 -pae1")
# Run("Core/Surface based on ±10 aa smoothing", "~/scripts/qsub.sh coresurf.py surf10 -humanonly")

# Waitforjobs()

print("Done!")
