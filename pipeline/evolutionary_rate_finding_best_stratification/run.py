#!/usr/bin/env python3
"""
Run: Run entire pipeline
"""

# Initialize

# from blang_mysql import *
from blang import *

# Start

combinations = []

# # All (2*3*5*6*7*5*7*4 = 176,400 combinations, though, so we need to scan them individually instead)
# nosurfs = ("-coresurf",)
# tests = ("wilcox_test", "oneway_test",)
# alternatives = ("two.sided", "less", "greater",)
# sources = ("all", "Ochoa", "UniProt", "PhosphoSitePlus", "dbPTM",)
# evorates = ("rate4site", "lichtarge", "capra0", "norm_rate4site", "norm_lichtarge", "norm_capra0",)
# mafftmodes = ("einsi_tree_1para", "einsi_tree", "einsi", "linsi_tree", "linsi", "ginsi_tree", "ginsi",)
# preds = ("AlphaFold", "AlphaFold_pLDDT70", "AlphaFold_pLDDT70_PAE2", "AlphaFold_pLDDT90", "AlphaFold_pLDDT90_PAE1",)
# minsizes = (2, 3, 5, 10, 20, 50, 100,)
# disfilts = ("all", "strcore", "strsurf", "dissurf",)
# for nosurf in nosurfs:
#     for test in tests:
#         for alternative in alternatives:
#             for source in sources:
#                 for evorate in evorates:
#                     for mafftmode in mafftmodes:
#                         for pred in preds:
#                             for minsize in minsizes:
#                                 for disfilt in disfilts:
#                                     combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best test & alternative:
nosurfs = ("-coresurf",)
tests = ("wilcox_test", "oneway_test",)
alternatives = ("two.sided", "less", "greater",)
sources = ("all",)
evorates = ("rate4site",)
mafftmodes = ("einsi_tree_1para",)
preds = ("AlphaFold_pLDDT70_PAE2",)
minsizes = (2,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best source:
nosurfs = ("-coresurf",)
tests = ("wilcox_test",)
alternatives = ("less",)
sources = ("all", "Ochoa", "UniProt", "PhosphoSitePlus", "dbPTM",)
evorates = ("rate4site",)
mafftmodes = ("einsi_tree_1para",)
preds = ("AlphaFold_pLDDT70_PAE2",)
minsizes = (2,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best evorate:
nosurfs = ("-coresurf",)
tests = ("wilcox_test",)
alternatives = ("less",)
sources = ("all",)
evorates = ("rate4site", "lichtarge", "capra0", "norm_rate4site", "norm_lichtarge", "norm_capra0",)
mafftmodes = ("einsi_tree_1para",)
preds = ("AlphaFold_pLDDT70_PAE2",)
minsizes = (2,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best mafftmode:
nosurfs = ("-coresurf",)
tests = ("wilcox_test",)
alternatives = ("less",)
sources = ("all",)
evorates = ("rate4site",)
mafftmodes = ("einsi_tree_1para", "einsi_tree", "einsi", "linsi_tree", "linsi", "ginsi_tree", "ginsi",)
preds = ("AlphaFold_pLDDT70_PAE2",)
minsizes = (2,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best pred:
nosurfs = ("-coresurf",)
tests = ("wilcox_test",)
alternatives = ("less",)
sources = ("all",)
evorates = ("rate4site",)
mafftmodes = ("einsi_tree_1para",)
preds = ("AlphaFold", "AlphaFold_pLDDT70", "AlphaFold_pLDDT70_PAE2", "AlphaFold_pLDDT90", "AlphaFold_pLDDT90_PAE1",)
minsizes = (2,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])

# Finding out best minsize:
nosurfs = ("-coresurf",)
tests = ("wilcox_test",)
alternatives = ("less",)
sources = ("all",)
evorates = ("rate4site",)
mafftmodes = ("einsi_tree_1para",)
preds = ("AlphaFold_pLDDT70_PAE2",)
minsizes = (2, 3, 5, 10, 20, 50, 100,)
disfilts = ("all", "strcore", "strsurf", "dissurf",)
for nosurf in nosurfs:
    for test in tests:
        for alternative in alternatives:
            for source in sources:
                for evorate in evorates:
                    for mafftmode in mafftmodes:
                        for pred in preds:
                            for minsize in minsizes:
                                for disfilt in disfilts:
                                    combinations.append([nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt])


# Unique & sort
combinations = sorted(set(tuple(i) for i in combinations))
print(f"\nNumber of unique combinations: {len(combinations)}\n")

Cd("tmp")
# d()
for (nosurf, test, alternative, source, evorate, mafftmode, pred, minsize, disfilt) in combinations:
    Run("Submit job", f"~/scripts/qsub.sh ../main.R {nosurf} {test} {alternative} {source} {evorate} {mafftmode} {pred} {minsize} {disfilt}")
Cd("..")

Waitforjobs()

# Compile output TSV files into one
Run("Write header to output table", "echo 'nosurf	test	alternative	source	evorate	mafftmode	pred	minsize	disfilt	strat	pval	sig	min' > output-table.tsv")
Run("Compile output TSV files into one", "cat output/output-*.tsv | grep -vP '^nosurf	test	alternative	source	evorate	mafftmode	pred	minsize	disfilt	strat	pval	sig	min' >> output-table.tsv")

# Import into mysql
Run("Import into table 'evorate_strat'", "~/scripts/import.pl output-table.tsv evorate_strat -overwrite -allindices")

print("Done!")
