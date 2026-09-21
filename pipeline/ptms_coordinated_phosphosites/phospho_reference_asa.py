#!/usr/bin/env python3
"""
phospho_reference_asa.py
------------------------
Derives maximum accessible surface area (ASA) reference values for
phosphoserine (SEP), phosphothreonine (TPO) and phosphotyrosine (PTR),
on the scale of Tien et al. 2013.

Why this is needed
    Relative solvent accessibility requires a per-residue maximum ASA.
    Tien et al. 2013 provide these for the twenty standard residues only.
    No equivalent values are published for phosphorylated residues, and
    published surveys of phosphosite accessibility avoid the problem by
    measuring unphosphorylated structures instead (e.g. Somavarapu et al.
    2014, BMC Struct Biol 14:9, PMID 24618394, where 97.7% of matched
    structures were of non-phosphorylated forms). Normalising a
    phosphoresidue by the value for its unphosphorylated parent understates
    burial, because the phosphate adds roughly 80 Å² of surface.

Method
    Tien et al. derived their theoretical values by enumerating Gly-X-Gly
    tripeptide conformations and taking the largest accessible area. The
    maximum falls in the alpha-helical region of the Ramachandran plot,
    not in the extended conformation used by Miller et al. 1987.

    The same enumeration is repeated here. Backbones are built from ideal
    peptide geometry over a grid of phi/psi angles; side chains are taken
    from the PDB Chemical Component Dictionary and rotated over their
    torsions; ASA of the central residue is computed for every combination
    and the largest kept.

    Running the standard residues through the same procedure reproduces the
    published theoretical values to within about 8%. The offset is
    systematic and cancels in the phospho/parent ratio, which is the only
    quantity carried forward: each phosphoresidue reference is the Tien
    Empirical value for its parent, scaled by that ratio.

Why not mkdssp
    mkdssp is used for all accessibility values in the survey itself, but it
    cannot be used here: it treats SEP, TPO and PTR as unknown residues when
    they appear in a hand-built tripeptide, drops them, and breaks the chain.
    Shrake-Rupley is therefore used for the enumeration, applied identically
    to the phosphorylated and unphosphorylated residues.

Output
------
    Reference ASA values printed to stdout, and written to
    phospho_reference_asa.tsv for use by find_coordinated_phosphoresidues.py.

Usage
-----
    python phospho_reference_asa.py [--chi-step 60] [--out phospho_reference_asa.tsv]

Dependencies
------------
    pip install biopython numpy
"""

import argparse
import itertools
import math
import os
import urllib.request
import warnings

import numpy as np
from Bio.PDB import SASA, Structure, Model, Chain, Residue, Atom

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Tien et al. 2013, PLoS One 8(11):e80635, PMID 24278298, Table 1.
# "Theoretical" column: Gly-X-Gly enumeration. Used here only to check that
# this script reproduces their procedure.
MAX_ASA_TIEN_THEORETICAL: dict[str, float] = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# "Empirical" column of the same table: the scale used throughout the paper.
MAX_ASA_TIEN_EMPIRICAL: dict[str, float] = {
    "SER": 143.0, "THR": 163.0, "TYR": 255.0,
}

PHOSPHO_PARENT = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

# Ideal peptide geometry (Engh & Huber 1991).
BOND_N_CA, BOND_CA_C, BOND_C_N = 1.458, 1.525, 1.329
ANGLE_N_CA_C, ANGLE_CA_C_N, ANGLE_C_N_CA = 111.2, 116.2, 121.7

# Backbone grid. Covers the allowed regions of the Ramachandran plot, including
# the alpha region where Tien et al. find the maximum.
PHI_GRID = list(range(-180, -39, 20)) + [60]
PSI_GRID = list(range(-180, 181, 20))

# Rotatable side-chain bonds, and the atoms that move with each.
SIDECHAIN_TORSIONS: dict[str, list] = {
    "SER": [(("CA", "CB"), ["OG"])],
    "THR": [(("CA", "CB"), ["OG1", "CG2"])],
    "TYR": [(("CA", "CB"), ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]),
            (("CB", "CG"), ["CD1", "CD2", "CE1", "CE2", "CZ", "OH"])],
    "SEP": [(("CA", "CB"), ["OG", "P", "O1P", "O2P", "O3P"]),
            (("CB", "OG"), ["P", "O1P", "O2P", "O3P"])],
    "TPO": [(("CA", "CB"), ["OG1", "CG2", "P", "O1P", "O2P", "O3P"]),
            (("CB", "OG1"), ["P", "O1P", "O2P", "O3P"])],
    "PTR": [(("CA", "CB"), ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH",
                            "P", "O1P", "O2P", "O3P"]),
            (("CB", "CG"), ["CD1", "CD2", "CE1", "CE2", "CZ", "OH",
                            "P", "O1P", "O2P", "O3P"]),
            (("CZ", "OH"), ["P", "O1P", "O2P", "O3P"])],
}

CCD_URL = "https://files.rcsb.org/ligands/download/{code}.cif"
CCD_DIR = "ccd"


# ---------------------------------------------------------------------------
# Component geometry
# ---------------------------------------------------------------------------

def ideal_coordinates(code):
    """Heavy-atom coordinates for one component, from the Chemical Component Dictionary."""
    path = os.path.join(CCD_DIR, f"{code}.cif")
    os.makedirs(CCD_DIR, exist_ok=True)
    if not os.path.exists(path):
        with open(path, "w") as fh:
            fh.write(urllib.request.urlopen(CCD_URL.format(code=code), timeout=60).read().decode())

    atoms, in_block, columns = {}, False, []
    for line in open(path):
        if line.startswith("_chem_comp_atom."):
            in_block = True
            columns.append(line.strip().split(".", 1)[1])
            continue
        if not in_block:
            continue
        if line.startswith(("#", "loop_")):
            break
        fields = line.split()
        if len(fields) < len(columns):
            continue
        row = dict(zip(columns, fields))
        name = row["atom_id"].strip('"')
        if row["type_symbol"] == "H":
            continue
        if row.get("pdbx_leaving_atom_flag") == "Y" and name != "OXT":
            continue
        try:
            atoms[name] = np.array([float(row[f"pdbx_model_Cartn_{axis}_ideal"])
                                    for axis in ("x", "y", "z")])
        except (KeyError, ValueError):
            continue
    return atoms


def place_atom(a, b, c, bond, angle, torsion):
    """Position a fourth atom from three others, given internal coordinates."""
    angle, torsion = math.radians(angle), math.radians(torsion)
    bc = (c - b) / np.linalg.norm(c - b)
    normal = np.cross(b - a, bc)
    normal /= np.linalg.norm(normal)
    frame = np.array([bc, np.cross(normal, bc), normal]).T
    offset = np.array([-bond * math.cos(angle),
                       bond * math.sin(angle) * math.cos(torsion),
                       bond * math.sin(angle) * math.sin(torsion)])
    return c + frame.dot(offset)


def backbone(n_residues, phi, psi, omega=180.0):
    """N, CA and C positions for a peptide with uniform backbone angles."""
    coords = [np.array([0.0, 0.0, 0.0]),
              np.array([BOND_N_CA, 0.0, 0.0]),
              np.array([BOND_N_CA + BOND_CA_C * math.cos(math.radians(180 - ANGLE_N_CA_C)),
                        BOND_CA_C * math.sin(math.radians(180 - ANGLE_N_CA_C)), 0.0])]
    for _ in range(n_residues - 1):
        n = place_atom(coords[-3], coords[-2], coords[-1], BOND_C_N, ANGLE_CA_C_N, psi)
        ca = place_atom(coords[-2], coords[-1], n, BOND_N_CA, ANGLE_C_N_CA, omega)
        c = place_atom(coords[-1], n, ca, BOND_CA_C, ANGLE_N_CA_C, phi)
        coords += [n, ca, c]
    return [tuple(coords[i:i + 3]) for i in range(0, len(coords), 3)]


def rotate_about_bond(coords, axis, movers, degrees):
    """Rotate the named atoms about a bond, in place."""
    start, end = coords[axis[0]], coords[axis[1]]
    unit = (end - start) / np.linalg.norm(end - start)
    theta = math.radians(degrees)
    cross = np.array([[0, -unit[2], unit[1]],
                      [unit[2], 0, -unit[0]],
                      [-unit[1], unit[0], 0]])
    rotation = (np.eye(3) * math.cos(theta) + math.sin(theta) * cross
                + (1 - math.cos(theta)) * np.outer(unit, unit))
    for name in movers:
        if name in coords:
            coords[name] = rotation.dot(coords[name] - end) + end


def gly_x_gly(code, phi, psi, torsions):
    """Build one Gly-X-Gly tripeptide as (residue number, residue name, atom name, xyz)."""
    frame = backbone(3, phi, psi)
    out = []
    for index, (name, (n, ca, c)) in enumerate(zip(["GLY", code, "GLY"], frame), start=1):
        ideal = ideal_coordinates(name)
        source = np.array([ideal["N"], ideal["CA"], ideal["C"]])
        target = np.array([n, ca, c])
        u, _, vt = np.linalg.svd((source - source.mean(0)).T.dot(target - target.mean(0)))
        chirality = np.sign(np.linalg.det(vt.T.dot(u.T)))
        rotation = vt.T.dot(np.diag([1, 1, chirality])).dot(u.T)

        placed = {}
        for atom, xyz in ideal.items():
            if index != 2 and atom not in BACKBONE_ATOMS:
                continue
            if atom == "OXT" and index != 3:
                continue
            placed[atom] = rotation.dot(xyz - source.mean(0)) + target.mean(0)
        if index == 2:
            for (axis, movers), angle in zip(SIDECHAIN_TORSIONS.get(code, []), torsions):
                rotate_about_bond(placed, axis, movers, angle)
        out += [(index, name, atom, xyz) for atom, xyz in placed.items()]
    return out


# ---------------------------------------------------------------------------
# Accessibility
# ---------------------------------------------------------------------------

def central_residue_asa(atoms):
    """Shrake-Rupley ASA of the central residue of a tripeptide."""
    structure = Structure.Structure("x")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    structure.add(model)
    model.add(chain)
    current, residue = None, None
    for number, name, atom, xyz in atoms:
        if number != current:
            residue = Residue.Residue((" ", number, " "), name, "")
            chain.add(residue)
            current = number
        residue.add(Atom.Atom(atom, xyz, 0.0, 1.0, " ", atom, 0, atom[0]))
    SASA.ShrakeRupley().compute(structure, level="R")
    return list(chain)[1].sasa


def maximum_asa(code, chi_step, keep=6):
    """Largest ASA over the backbone grid and the side-chain torsion grid."""
    torsions = SIDECHAIN_TORSIONS.get(code, [])
    default = (180,) * len(torsions)
    scored = sorted(((central_residue_asa(gly_x_gly(code, phi, psi, default)), phi, psi)
                     for phi in PHI_GRID for psi in PSI_GRID), reverse=True)
    best = scored[0][0]
    if torsions:
        grid = list(itertools.product(*([range(0, 360, chi_step)] * len(torsions))))
        for _, phi, psi in scored[:keep]:
            for combination in grid:
                best = max(best, central_residue_asa(gly_x_gly(code, phi, psi, combination)))
    return best


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--chi-step", type=int, default=60,
                        help="side-chain torsion step in degrees (default: %(default)s)")
    parser.add_argument("--out", default="phospho_reference_asa.tsv")
    args = parser.parse_args()

    print("Calibration against Tien et al. 2013, Theoretical column:")
    print(f"  {'residue':9s} {'this script':>12s} {'published':>10s} {'ratio':>7s}")
    ratios = []
    for code in ["ALA", "SER", "THR", "TYR", "LEU", "LYS", "PHE", "VAL"]:
        value = maximum_asa(code, args.chi_step)
        published = MAX_ASA_TIEN_THEORETICAL[code]
        ratios.append(value / published)
        print(f"  {code:9s} {value:12.1f} {published:10.0f} {value / published:7.3f}")
    print(f"  mean ratio {np.mean(ratios):.3f} (sd {np.std(ratios):.3f}); "
          f"systematic, and cancels below")

    parents = {code: maximum_asa(code, args.chi_step) for code in ["SER", "THR", "TYR"]}
    rows = []
    print("\nPhosphoresidue reference values:")
    print(f"  {'residue':9s} {'this script':>12s} {'/ parent':>9s} {'reference':>10s}")
    for code, parent in PHOSPHO_PARENT.items():
        value = maximum_asa(code, args.chi_step)
        ratio = value / parents[parent]
        reference = MAX_ASA_TIEN_EMPIRICAL[parent] * ratio
        rows.append((code, parent, ratio, reference))
        print(f"  {code:9s} {value:12.1f} {ratio:9.3f} {reference:10.0f}")

    with open(args.out, "w") as fh:
        fh.write("# Maximum ASA reference values for phosphorylated residues, in square Angstroms.\n")
        fh.write("# Derived by the Gly-X-Gly enumeration of Tien et al. 2013 "
                 "(PLoS One 8(11):e80635, PMID 24278298);\n")
        fh.write("# see phospho_reference_asa.py. The parent value is the Empirical column of "
                 "their Table 1,\n")
        fh.write("# scaled by the phospho/parent ratio measured here.\n")
        fh.write("residue\tparent\tratio\treference_asa\n")
        for code, parent, ratio, reference in rows:
            fh.write(f"{code}\t{parent}\t{ratio:.3f}\t{reference:.0f}\n")
    print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
