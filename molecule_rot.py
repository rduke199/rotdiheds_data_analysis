import os


class MoleculeRot:
    def __init__(self, origin_fn):
        self.origin_fn = os.path.abspath(origin_fn)
        self.name = str(origin_fn.split('energies_')[-1]).strip('.cvs')
        self.ring_num = int(self.name.split('_')[1])
        self.unit_num = int(self.name.split('_')[2])
        self.polymer_num = int(self.name.split('_')[3])
        self.substituents = self.name.split('_')[4]
        self.zero_energy = self.energy_dict[0]

    def __str__(self):
        return  f'name: {self.name}\n{self.ring_num} ring type, {self.unit_num} monomer units, {self.substituents} substituents\nenergy dictionary: {self.energy_dict}'

    @property
    def energy_dict(self):
        data = np.genfromtxt(fname=self.origin_fn,delimiter=',', dtype='float')
        return dict(data)

    @property
    def norm_energy_dict(self):
        _norm_energy_dict = {deg: eng - self.zero_energy for deg, eng in self.energy_dict.items()}
        return _norm_energy_dict

    def Plot(self, title, out_dir, label=None):
        """
        Plots the energy vs the torsion angle for a molecule given a torsion energy
        file. Saves the plot to a png file.
        """
        phi, energy = self.norm_energy_dict.keys(), self.norm_energy_dict.values()
        min_energy, max_energy = min(energy), max(energy)
        energy_range = max_energy - min_energy

        plt.scatter(phi, energy, label = label)
        plt.xlim(min(phi)-3, max(phi)+3)
        plt.xticks(np.linspace(start=0, stop=180, num=7))
        plt.ylim(top = max_energy + energy_range*0.15,bottom = min_energy - energy_range*0.15)
        plt.yticks(np.linspace(start=min_energy, stop=max_energy, num=5))
        plt.xlabel("dihedral angle (degrees)")
        plt.ylabel("energy (Hartrees)")
        plt.title(title)
        plt.legend()
        plt.savefig(out_dir+'torsionE_plt_{}.png'.format(title),dpi=300)
        plt.close('all')


def AllMolList(energy_dir):
    all_molecules = []
    for f in os.listdir(energy_dir):
        fpath = os.path.join(energy_dir,f)
        if f.startswith("energies"):
            try:
                mol = MoleculeRot(fpath)
                all_molecules.append(mol)
            except:
                print("Error. Did not convert {} to Molecule Object.".format(f))
    return all_molecules

def MakeMolList(all_molecules,ring_num=None,unit_num=None,polymer_num=None,substituents=None):
    mol_list = []
    for mol in all_molecules:
        if ring_num is not None:
            if mol.ring_num != ring_num:
                continue
        if unit_num is not None:
            if mol.unit_num != unit_num:
                continue
        if polymer_num is not None:
            if mol.polymer_num != polymer_num:
                continue
        if substituents is not None:
            if bool(re.search(substituents, mol.substituents)) == False:
                continue
        mol_list.append(mol)
    return mol_list


def OverlayPlts(mol_list,title,out_dir,varying_atrib='unit_num'):
    """
    Makes overlay plot of allenergies files in cwd.
    """
    min_energy, max_energy = 1000, -1000
    for mol in mol_list:
        phi, energy = mol.norm_energy_dict.keys(), mol.norm_energy_dict.values()
        plt.scatter(phi, energy, label = eval('mol.'+varying_atrib))
        min_e, max_e = min(energy), max(energy)
        if min_e < min_energy:
            min_energy = min_e
        if max_e > max_energy:
            max_energy = max_e
    energy_range = max_energy - min_energy

    plt.xlim(min(phi)-3, max(phi)+3)
    plt.xticks(np.linspace(start=0, stop=180, num=7))
    plt.ylim(top = max_energy + energy_range*0.15,bottom = min_energy - energy_range*0.15)
    plt.yticks(np.linspace(start=min_energy, stop=max_energy, num=5))
    plt.xlabel("dihedral angle (degrees)")
    plt.ylabel("energy (Hartrees)")
    plt.title(title)
    plt.legend(loc=(1.04,0))
    plt.savefig(out_dir+'torsionE_plt_{}.png'.format(title),dpi=300)

def AveragesPlts(mol_list,title,out_dir):
    """
    Makes overlay plot of polymer unit energy averages (with error bars) in cwd.
    """
    data = np.zeros([19,3])
    for n,d in enumerate(np.linspace(0,180,19)):
        data[n][0] = d
        energies = []
        for mol in mol_list:
            phi, energy = mol.norm_energy_dict.keys(), mol.norm_energy_dict.values()
            try:
                energies.append(mol.norm_energy_dict[d])
                np_energies = np.array(energies)
                data[n][1] = np.average(np_energies)
                data[n][2] = np.std(np_energies)/np.sqrt(len(np_energies))
            except:
                print('Error. Did not find energy for {} at {} degree rotation.'.format(mol.name,d))
    min_energy = np.amin(data[:,1]) - np.amax(data[:,2])
    max_energy = np.amax(data[:,1]) + np.amax(data[:,2])
    energy_range = max_energy - min_energy

    print('\n\nNumber of molecules averaged: {}'.format(len(mol_list)))
    plt.scatter(data[:,0], data[:,1])
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,2], linestyle="None")
    plt.xlim(min(phi)-3, max(phi)+3)
    plt.xticks(np.linspace(start=0, stop=180, num=7))
    plt.ylim(top = max_energy + energy_range*0.15,bottom = min_energy - energy_range*0.15)
    plt.yticks(np.linspace(start=min_energy, stop=max_energy, num=5))
    plt.xlabel("dihedral angle (degrees)")
    plt.ylabel("energy (Hartrees)")
    plt.title(title)
    plt.savefig(out_dir+'torsionE_plt_{}.png'.format(title),dpi=300)
