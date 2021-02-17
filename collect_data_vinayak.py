import os
import re
import gc
import json
from rdkit.Chem import MolToSmiles
from pymatgen.core import Molecule
from pymatgen.core import structure
from pymatgen.io.gaussian import GaussianOutput
from ocelot.routines.conformerparser import pmgmol_to_rdmol

def MakeJSON(mol_dir,omega_file,json_file):
    """
    For a molecule rotation directory, mol_dir, this finds the energy, xyz, and
    degree of rotation and places those into  json_file.
    """
    mol_name = str(mol_dir.split('/')[-1])
    energy = 0
    mol_xyz = 0
    try:
    # if 0 == 0:
        log_fn = [x for x in os.listdir(mol_dir) if x.endswith('opt_0.log')][0]
        log_path = os.path.join(mol_dir,log_fn)
        mol = GaussianOutput(log_path)
        if mol.properly_terminated:
            num_electrons = mol.electrons[0]
            eigens = list(mol.eigenvalues.values())[0]
            homo = eigens[num_electrons - 1] * 27.2114
            lumo = eigens[num_electrons] * 27.2114
            homo_lumo_gap = lumo - homo
            pymat_mol = Molecule.from_file(log_path)
            smi_mol = pmgmol_to_rdmol(pymat_mol)[1]
            with open(omega_file,'r') as fn:
                omega_data = fn.readlines()[-2].split()[1]
                omega = "0{}".format(omega_data.split('.')[1])
            normal = 1
        else:
            normal = 0
            print("Error. wtuning for did not properly terminate for {}".format(mol_name))
    except Exception as e:
        normal = 0
        print("Error. Data NOT collected for {}. Energy calculations for dihedral rotations may not have finished!".format(mol_name))
        print(e)
    if normal == 1:
        json_data = {
        "molecule_name" : mol_name,
        "smiles" : smi_mol,
        "tuned_omega" : omega,
        "homo" : homo,
        "lumo" : lumo,
        "homo_lumo_gap" : homo_lumo_gap
        }
        with open(json_file,'w+') as fn:
            json.dump(json_data, fn)
            print("Data collected for {}".format(mol_name))
    else:
        pass
    try:
        json_file.close()
        omega_file.close()
    except:
        pass

cwd = os.getcwd()
gather_home = os.path.join(cwd,'rotated_dihed')
out_home = os.path.join(cwd,'vjsons')
for m in os.listdir(gather_home):
    mpath = os.path.join(gather_home,m)
    if os.path.isdir(mpath):
        json_file = "{}/{}.json".format(out_home,m)
        omega_file = os.path.join(mpath,'output.log')
        if os.path.isfile(json_file) == False:
            MakeJSON(mpath,omega_file,json_file)
