import os
import re
import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt


class MoleculeRot:
    """
    MoleculeRot Object

    Creates an object for each molecule with a variety of properties including a dictionary (and a normalized
    dictionary) for it's energies at each dihedral rotation. Each object also has descriptive properties such as name,
    smile ring number, unit number, polymer number, etc. There are only a couple class functions, but these draw the
    molecule structure and plot the PE curve.

    """

    def __init__(
            self,
            name,
            master_dir,
            unified_unconst = False,
            energy_fn="master_energy.json",
            homo_fn="master_homos.json",
            lumo_fn="master_lumos.json",
            unconst_fn="master_geomopt.json",
            omega_fn="master_omega.json",
            smiles_fn="master_smiles.json"
    ):
        # Note: when entering this analysis, all energies should be in eV
        self.name = name
        self.master_dir = master_dir
        self.unified_unconst = unified_unconst

        self.ring_num = int(self.name.split('_')[1])
        self.unit_num = int(self.name.split('_')[2])
        self.polymer_num = int(self.name.split('_')[3])
        self.substituents = self.name.split('_')[4]

        self.energy_dict = self.make_dict_floats(self.get_data_from_master(os.path.join(master_dir, energy_fn)))
        self.homo_dict = self.make_dict_floats(self.get_data_from_master(os.path.join(master_dir, homo_fn)))
        self.lumo_dict = self.make_dict_floats(self.get_data_from_master(os.path.join(master_dir, lumo_fn)))
        self.tuned_omega = self.get_data_from_master(os.path.join(master_dir, omega_fn))
        self.smiles = self.get_data_from_master(os.path.join(master_dir, smiles_fn))
        self.unconst_path = os.path.join(master_dir, unconst_fn)
        self.unconst_data = self.get_data_from_master(self.unconst_path)


    def __str__(self):
        return f'name: {self.name}\n{self.ring_num} ring type, {self.unit_num} monomer units, {self.substituents} ' \
               f'substituents\nenergy dictionary: {self.energy_dict}'

    def get_data_from_master(self, master_file):
        with open(master_file, 'r') as fn:
            _dict = json.load(fn)
        try:
            _data = _dict[self.name]
            return _data
        except KeyError:
            try:
                _data = [v for k, v in _dict.items() if k.startswith(self.name)][0]
                return _data
            except IndexError:
                return None

    @staticmethod
    def make_dict_floats(_dict):
        if _dict is not None:
            return dict(sorted({float(key): float(value) for key, value in _dict.items()}.items()))

    @property
    def unified_unconst_data(self):
        """
        This returns the unconstrained energy data from the minimum energy length out of all lengths of this polymer,
        instead of the unconstrained energy data from this polymer at this particular length
        """
        units = [1,3,5,7]
        angle_dict, energy_dict = {}, {}
        with open(self.unconst_path, 'r') as fn:
            _dict = json.load(fn)
        for unit in units:
            mol_name = 'mols_{}_{}_{:02d}_{}'.format(self.ring_num, unit, self.polymer_num, self.substituents)
            try:
                angle_dict[unit] = _dict[mol_name][0]
                energy_dict[unit] = _dict[mol_name][1]
            except KeyError:
                try:
                    angle_dict[unit], energy_dict[unit] = [v for k, v in _dict.items() if k.startswith(mol_name)][0]
                except IndexError:
                    pass
        try:
            min_energy_unit = min(energy_dict, key=energy_dict.get)
            return [angle_dict[min_energy_unit], energy_dict[min_energy_unit]]
        except ValueError:
            return None

    @property
    def unconst_energy(self):
        if self.unified_unconst:
            _dict = self.unified_unconst_data
        else:
            _dict = self.unconst_data
        if _dict is None:
            return None
        else:
            return _dict[1]

    @property
    def unconst_angle(self):
        if self.unified_unconst:
            _dict = self.unified_unconst_data
        else:
            _dict = self.unconst_data
        _dict = self.unconst_data
        if _dict is None:
            return None
        else:
            return abs(_dict[0])

    @property
    def min_e(self):
        return min(self.energy_dict.values())

    @property
    def max_e_norm(self):
        return max(self.norm_energy_dict.values())

    @property
    def norm_energy_dict(self):
        # This dictionary is the only one with kcal/mol energies
        # if self.unconst_energy is None:
        _norm_edict = {float(deg): 23.06 * (float(eng) - self.min_e) for deg, eng in self.energy_dict.items()}
        _sorted_norm_edict = dict(sorted(_norm_edict.items()))
        return _sorted_norm_edict
        # else:
        #     _norm_edict = {float(deg): 23.06 * (float(eng) - self.unconst_energy) for deg, eng in
        #                    self.energy_dict.items()}
        #     _sorted_norm_edict = dict(sorted(_norm_edict.items()))
        #     return _sorted_norm_edict

    def draw_structure(self, out_dir=None):
        molecule = Chem.MolFromSmiles(self.smiles)
        AllChem.Compute2DCoords(molecule)
        if out_dir is not None:
            Draw.MolToFile(molecule, out_dir + 'img_{}.png'.format(self.name))
        return Draw.MolToImage(molecule)

    def plot_torsion_energy(self, out_dir=None):
        """
        Plots the energy vs the torsion angle for a molecule given a torsion energy
        file. Saves the plot to a png file.
        """
        fig, ax = plt.subplots()
        phi, energy = self.norm_energy_dict.keys(), self.norm_energy_dict.values()
        plt.scatter(phi, energy)

        # ax.set_xlim(0, 3)
        # ax.set_ylim(0, 3)

        plt.xlim(min(phi) - 3, max(phi) + 3)
        # plt.xticks(np.linspace(start=0, stop=180, num=7))
        plt.ylim(top=max(energy)+5, bottom=min(energy)-5)
        # plt.yticks(np.linspace(start=-10, stop=14, num=5))
        plt.xlabel("dihedral angle (degrees)")
        plt.ylabel("energy (kcal/mol)")
        plt.title("Energy for " + self.name)
        fig.set_facecolor('w')
        if out_dir is not None:
            plt.savefig(out_dir + 'torsionE_plt_{}.png'.format(self.name), dpi=300, bbox_inches='tight')
            plt.close('all')

    def plot_homo_lumo(self, out_dir=None):
        """
        Plots the homo and lumo vs the torsion angle for a molecule given a torsion energy
        file. Saves the plot to a png file.
        """
        fig, ax = plt.subplots()
        phi_h, homo = self.homo_dict.keys(), self.homo_dict.values()
        phi_l, lumo = self.lumo_dict.keys(), list(self.lumo_dict.values())

        # eng_max, eng_min = max(max([homo, lumo])), min(min([homo, lumo]))

        ax.scatter(phi_h, homo, label='HOMO')
        ax.plot(phi_h, homo)
        ax.scatter(phi_l, lumo, label='LUMO')
        ax.plot(phi_l, lumo)

        ax.set_xlim(-3, 183)
        ax.set_xticks(np.linspace(start=0, stop=180, num=7))
        ax.set_ylim(top=5,bottom=-20)
        ax.set_yticks(np.linspace(start=-20, stop=5, num=6))
        ax.set_xlabel("dihedral angle (degrees)")
        ax.set_ylabel("energy (kcal/mol)")
        ax.set_title("HOMO/LUMO for " + self.name)
        fig.set_facecolor('w')
        fig.legend(loc='lower left', bbox_to_anchor=(0.9, 0.2))
        if out_dir is not None:
            plt.savefig(out_dir + 'tortion_homo_lumo_{}.png'.format(self.name), dpi=300, bbox_inches='tight')
            plt.close('all')

    def get_angle_from_energy(self, energy):
        for angle, eg in self.energy_dict.items():
            if energy == eg:
                return angle
        return "Dihedral angle doesn't exist for {}.".format(energy)

    @staticmethod
    def write_json(data, filename):
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)

    def write_nnff_json(self, json_outdir):
        for angle in self.norm_energy_dict.keys():
            json_data = {
                "molecule_name": self.name,
                "degree": angle,
                "smile": self.smiles,
                "energy": self.norm_energy_dict[angle]
            }
            self.write_json(json_data, "{}/{}_{}deg.json".format(json_outdir, self.name, int(angle)))
