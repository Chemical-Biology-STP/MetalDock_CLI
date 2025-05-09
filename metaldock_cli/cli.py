#!/usr/bin/env python3

import os
import re
import subprocess
import sys

import numpy as np
import questionary
from ase.calculators.emt import EMT
from ase.io import read, write
from ase.optimize import BFGS
from Bio.PDB import PDBParser
from mendeleev import element
from prompt_toolkit.completion import PathCompleter as PromptPathCompleter


# === Helper ===
def safe_ask(q):
    try:
        answer = q.ask()
        if answer is None:
            raise KeyboardInterrupt
        return answer
    except KeyboardInterrupt:
        print("\n👋 Exiting by user request.")
        sys.exit(0)


def sanitize_filename(name):
    return re.sub(r"[^a-zA-Z0-9]", "_", name).lower()


def detect_metal_symbols(xyz_path):
    metals = {
        "Sc",
        "Ti",
        "V",
        "Cr",
        "Mn",
        "Fe",
        "Co",
        "Ni",
        "Cu",
        "Zn",
        "Y",
        "Zr",
        "Nb",
        "Mo",
        "Tc",
        "Ru",
        "Rh",
        "Pd",
        "Ag",
        "Cd",
        "Hf",
        "Ta",
        "W",
        "Re",
        "Os",
        "Ir",
        "Pt",
        "Au",
        "Hg",
    }
    found = set()
    with open(xyz_path, "r") as f:
        lines = f.readlines()[2:]
        for line in lines:
            parts = line.split()
            if parts and parts[0] in metals:
                found.add(parts[0])
    return list(found)


def parse_xyz(file_path):
    atoms = []
    with open(file_path, "r") as f:
        lines = f.readlines()
        natoms = int(lines[0].strip())
        for line in lines[2 : 2 + natoms]:
            parts = line.split()
            if parts:
                atoms.append(parts[0])
    return atoms


def calculate_valence(atoms, charge):
    total = sum(element(atom).electrons for atom in atoms)
    return total - charge


def ase_validate_and_optimize_structure(path):
    atoms = read(path)
    atoms.center()
    z_coords = [atom.position[2] for atom in atoms]
    if all(abs(z) < 1e-3 for z in z_coords):
        atoms.set_calculator(EMT())
        dyn = BFGS(atoms, logfile=None)
        dyn.run(fmax=0.05)
    write(path, atoms)


def find_ligands_in_pdb(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    ligands = {}
    for model in structure:
        for chain in model:
            for res in chain:
                hetflag, resseq, _ = res.id
                if hetflag.startswith("H") and res.resname not in ("HOH",):
                    coords = [
                        atom.coord for atom in res if atom.element != "H"
                    ]
                    if coords:
                        center = np.mean(coords, axis=0)
                        label = f"{res.resname}_{chain.id}_{resseq}"
                        ligands[label] = center
    return ligands


def check_orca_and_metaldock():
    errors = []
    orca_path = subprocess.run(
        ["which", "orca"], capture_output=True, text=True
    )
    metaldock_path = subprocess.run(
        ["which", "metaldock"], capture_output=True, text=True
    )

    if orca_path.returncode != 0:
        errors.append(
            "❌ ORCA not found in PATH. Contact Yew Mun (yewmun.yip@crick.ac.uk)"
        )
    if metaldock_path.returncode != 0:
        errors.append(
            "❌ MetalDock not found in PATH. Contact Yew Mun (yewmun.yip@crick.ac.uk)"
        )

    return errors, orca_path.stdout.strip(), metaldock_path.stdout.strip()


def get_docking_center(protein_path):
    is_holo = safe_ask(
        questionary.select(
            "Is this a holo or apo structure?", choices=["holo", "apo"]
        )
    )
    if is_holo == "holo":
        ligands = find_ligands_in_pdb(protein_path)
        if not ligands:
            print("⚠️ No ligands found. Switching to apo mode.")
            is_holo = "apo"
        elif len(ligands) == 1:
            return tuple(ligands[list(ligands.keys())[0]])
        else:
            selected = safe_ask(
                questionary.select(
                    "Select ligand for docking center:",
                    choices=list(ligands.keys()),
                )
            )
            return tuple(ligands[selected])

    # apo mode
    method = safe_ask(
        questionary.select(
            "How to define docking center?",
            choices=["Manual", "Reference ligand (XYZ)"],
        )
    )
    if method == "Manual":
        x = float(safe_ask(questionary.text("Dock X", default="10.0")))
        y = float(safe_ask(questionary.text("Dock Y", default="10.0")))
        z = float(safe_ask(questionary.text("Dock Z", default="10.0")))
        return x, y, z
    else:
        path = safe_ask(questionary.path("Reference ligand path:"))
        mol = read(path)
        mol.center()
        return tuple(mol.get_center_of_mass())


def run_single_ligand_flow():
    path_completer = PromptPathCompleter()

    protein = os.path.normpath(
        safe_ask(questionary.path("Protein PDB:", completer=path_completer))
    )
    ligand = os.path.normpath(
        safe_ask(questionary.path("Ligand XYZ:", completer=path_completer))
    )

    base = sanitize_filename(os.path.splitext(os.path.basename(ligand))[0])
    job_dir = os.path.abspath(base)
    os.makedirs(job_dir, exist_ok=True)

    shutil_protein = os.path.join(job_dir, os.path.basename(protein))
    shutil_ligand = os.path.join(job_dir, os.path.basename(ligand))
    with open(protein, "r") as src, open(shutil_protein, "w") as dst:
        dst.write(src.read())
    with open(ligand, "r") as src, open(shutil_ligand, "w") as dst:
        dst.write(src.read())

    charge = int(safe_ask(questionary.text("Ligand charge:", default="0")))
    ncpu = int(safe_ask(questionary.text("CPUs:", default="60")))
    ph = float(safe_ask(questionary.text("Protein pH:", default="7.5")))
    box_size = int(safe_ask(questionary.text("Box size:", default="20")))
    num_poses = int(
        safe_ask(questionary.text("Number of poses:", default="10"))
    )
    vacant = safe_ask(
        questionary.confirm("Vacant site on metal?", default=False)
    )
    geom = safe_ask(
        questionary.confirm("Perform geometry optimization?", default=True)
    )
    orca_input = safe_ask(
        questionary.text("ORCA input line:", default="B3LYP def2-TZVP")
    )
    orca_block = safe_ask(
        questionary.text(
            "ORCA block:", default='%basis newECP {{metal}} "def2-SD" end end'
        )
    )
    submit = safe_ask(questionary.confirm("Submit job?", default=False))

    errors, orca_cmd, metaldock_cmd = check_orca_and_metaldock()
    if errors:
        for e in errors:
            print(e)
        return

    try:
        ase_validate_and_optimize_structure(shutil_ligand)
    except Exception as e:
        print(f"❌ ASE optimization failed: {e}")
        return

    metals = detect_metal_symbols(shutil_ligand)
    if not metals or len(metals) > 1:
        print(f"❌ Metal detection failed: found {metals}")
        return
    metal = metals[0]

    atoms = parse_xyz(shutil_ligand)
    valence = calculate_valence(atoms, charge)
    spin = 1 if valence % 2 else 0

    dock_x, dock_y, dock_z = get_docking_center(protein)

    ini_path = os.path.join(job_dir, "input.ini")
    with open(ini_path, "w") as f:
        f.write(
            f"""
[DEFAULT]
metal_symbol = {metal}
method = dock
ncpu = {ncpu}

[PROTEIN]
pdb_file = {os.path.basename(protein)}
pH = {ph}
clean_pdb = True

[QM]
engine = ORCA
orcasimpleinput = {orca_input}
orcablocks = {orca_block.replace('{{metal}}', metal)}

[METAL_COMPLEX]
geom_opt = {geom}
xyz_file = {os.path.basename(ligand)}
charge = {charge}
spin = {spin}
vacant_site = {vacant}

[DOCKING]
ga_dock = True
dock_x = {dock_x:.2f}
dock_y = {dock_y:.2f}
dock_z = {dock_z:.2f}
box_size = {box_size}
random_pos = True
num_poses = {num_poses}
"""
        )

    sbatch_path = os.path.join(job_dir, f"metaldock_{base}.sbatch")
    with open(sbatch_path, "w") as f:
        f.write(
            f"""#!/bin/bash
#SBATCH --job-name=metaldock_{base}
#SBATCH --partition=ncpu
#SBATCH --ntasks={ncpu}
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=480G
#SBATCH --output=metaldock_{base}_%j.out
#SBATCH --error=metaldock_{base}_%j.err

module load OpenMPI/4.1.6-GCC-13.2.0
export ASE_ORCA_COMMAND='{orca_cmd} PREFIX.inp > PREFIX.out'
{metaldock_cmd} -i input.ini -m dock
"""
        )

    if submit:
        subprocess.run(["sbatch", sbatch_path], cwd=job_dir)


def run_batch_ligand_flow():
    path_completer = PromptPathCompleter()

    protein = os.path.normpath(
        safe_ask(questionary.path("Protein PDB:", completer=path_completer))
    )
    ligand_dir = os.path.normpath(
        safe_ask(
            questionary.path(
                "Directory of ligand XYZs:", completer=path_completer
            )
        )
    )
    ligand_files = [f for f in os.listdir(ligand_dir) if f.endswith(".xyz")]

    if not ligand_files:
        print("❌ No .xyz ligands found.")
        return

    charge = int(safe_ask(questionary.text("Ligand charge:", default="0")))
    ncpu = int(safe_ask(questionary.text("CPUs:", default="60")))
    ph = float(safe_ask(questionary.text("Protein pH:", default="7.5")))
    box_size = int(safe_ask(questionary.text("Box size:", default="20")))
    num_poses = int(
        safe_ask(questionary.text("Number of poses:", default="10"))
    )
    vacant = safe_ask(
        questionary.confirm("Vacant site on metal?", default=False)
    )
    geom = safe_ask(
        questionary.confirm("Perform geometry optimization?", default=True)
    )
    orca_input = safe_ask(
        questionary.text("ORCA input line:", default="B3LYP def2-TZVP")
    )
    orca_block = safe_ask(
        questionary.text(
            "ORCA block:", default='%basis newECP {{metal}} "def2-SD" end end'
        )
    )
    submit = safe_ask(questionary.confirm("Submit jobs?", default=False))

    errors, orca_cmd, metaldock_cmd = check_orca_and_metaldock()
    if errors:
        for e in errors:
            print(e)
        return

    dock_x, dock_y, dock_z = get_docking_center(protein)
    failed = []

    for file in ligand_files:
        full_path = os.path.join(ligand_dir, file)
        name = sanitize_filename(os.path.splitext(file)[0])
        job_dir = os.path.abspath(name)
        os.makedirs(job_dir, exist_ok=True)

        dst = os.path.join(job_dir, file)
        with open(full_path, "r") as src, open(dst, "w") as out:
            out.write(src.read())

        try:
            ase_validate_and_optimize_structure(dst)
        except Exception as e:
            failed.append((file, f"ASE failed: {e}"))
            continue

        metals = detect_metal_symbols(dst)
        if len(metals) != 1:
            failed.append(
                (file, f"{'No metal' if not metals else 'Multiple metals'}")
            )
            continue
        metal = metals[0]
        atoms = parse_xyz(dst)
        valence = calculate_valence(atoms, charge)
        spin = 1 if valence % 2 else 0

        ini = os.path.join(job_dir, "input.ini")
        with open(ini, "w") as f:
            f.write(
                f"""
[DEFAULT]
metal_symbol = {metal}
method = dock
ncpu = {ncpu}

[PROTEIN]
pdb_file = {os.path.basename(protein)}
pH = {ph}
clean_pdb = True

[QM]
engine = ORCA
orcasimpleinput = {orca_input}
orcablocks = {orca_block.replace('{{metal}}', metal)}

[METAL_COMPLEX]
geom_opt = {geom}
xyz_file = {file}
charge = {charge}
spin = {spin}
vacant_site = {vacant}

[DOCKING]
ga_dock = True
dock_x = {dock_x:.2f}
dock_y = {dock_y:.2f}
dock_z = {dock_z:.2f}
box_size = {box_size}
random_pos = True
num_poses = {num_poses}
"""
            )
        sbatch = os.path.join(job_dir, f"metaldock_{name}.sbatch")
        with open(sbatch, "w") as f:
            f.write(
                f"""#!/bin/bash
#SBATCH --job-name=metaldock_{name}
#SBATCH --partition=ncpu
#SBATCH --ntasks={ncpu}
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem=480G
#SBATCH --output=metaldock_{name}_%j.out
#SBATCH --error=metaldock_{name}_%j.err

module load OpenMPI/4.1.6-GCC-13.2.0
export ASE_ORCA_COMMAND='{orca_cmd} PREFIX.inp > PREFIX.out'
{metaldock_cmd} -i input.ini -m dock
"""
            )
        if submit:
            subprocess.run(["sbatch", sbatch], cwd=job_dir)

    if failed:
        print("❌ Failed ligands:")
        for f, reason in failed:
            print(f"  - {f}: {reason}")


# --- Entry ---
def main():
    print("💡 You can quit this script at any time by pressing Ctrl + C.\n")
    mode = questionary.select(
        "Run MetalDock for:", choices=["Single ligand", "Batch directory"]
    ).ask()
    if mode == "Single ligand":
        run_single_ligand_flow()
    else:
        run_batch_ligand_flow()


if __name__ == "__main__":
    print("💡 You can quit this script at any time by pressing Ctrl + C.\n")
    main()
