import os
import re
import json
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from molecule_rot import MoleculeRot


def make_all_mol_list(master_dir):
    """
    Make molecule object list.
    This function makes a list of MoleculeRot objects for all molecules in the energy_dir.

    :param master_dir: str path to master directory
    :return: list of all molecules as MoleculeRot objects.
    """
    all_molecules = []
    master_files = [x for x in os.listdir(master_dir) if x.startswith('master')]
    master_path = os.path.join(master_dir, master_files[0])
    with open(master_path, 'r') as fn:
        data = json.load(fn)
    mol_names = data.keys()
    for name in mol_names:
        try:
            mol = MoleculeRot(name, master_dir)
            all_molecules.append(mol)
        except:
            print("Error. Did not convert {} to Molecule Object.".format(name))
    return all_molecules


def make_small_mol_list(all_molecules, ring_num=None, unit_num=None, polymer_num=None, substituents=None):
    """
    Parse small list
    This function parses a list of molecules with specified attributes from a list of MoleculeRot objects.

    :param all_molecules: list of all molecules as MoleculeRot objects
    :param ring_num: str specified ring number
    :param unit_num: str specified unit number
    :param polymer_num: str specified polymer number
    :param substituents: str specified substituents
    :return: list of molecules with specified attributes as MoleculeRot objects
    """
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
            if not bool(re.search(substituents, mol.substituents)):
                continue
        mol_list.append(mol)
    return mol_list


def overlay_plot(mol_list, title, out_dir=None, varying_attribute='unit_num', legend_outside=True, draw_1unit=False):
    """
    Plot overlay plot.
    This function plots an overlay plot of all molecule PE curves in the mol_list

    :param mol_list: list of molecules as MoleculeRot objects
    :param title: str title
    :param out_dir: path to directory in which to save figure
    :param varying_attribute: str attribute that will be specified in legend
    :param legend_outside: boolean. If Ture, the legend will be outside the plot
    :param draw_1unit: boolean. If Ture, the monomer structure will be drawn outside the plot
    :return: None
    """
    fig, ax = plt.subplots()
    for mol in mol_list:
        if mol.unit_num == 1:
            mol_image = mpimg.pil_to_array(mol.draw_structure())
        try:
            phi, energy = mol.norm_energy_dict.keys(), mol.norm_energy_dict.values()
            ax.scatter(phi, energy, label=eval('mol.' + varying_attribute))
            ax.plot(phi, energy)
        except:
            print('Error. Did not plot {} {} for {}.'.format(varying_attribute, eval('mol.'+varying_attribute), title))

    ax.set_xlim(-3, 183)
    ax.set_xticks(np.linspace(start=0, stop=180, num=7))
    ax.set_ylim(top=15, bottom=-1)
    ax.set_yticks(np.linspace(start=0, stop=14, num=8))
    ax.set_xlabel("dihedral angle (degrees)")
    ax.set_ylabel("energy (kcal/mol)")
    ax.set_title(title)
    fig.patch.set_facecolor("w")
    if legend_outside:
        fig.legend(loc='lower left', bbox_to_anchor=(1, 0))
    if draw_1unit:
        new_ax = fig.add_axes([0.7, .2, 0.4, 0.4], anchor='NE')
        img = new_ax.add_artist(AnnotationBbox(OffsetImage(mol_image, zoom=0.5), (1, 1)))
        new_ax.axis('off')
        fig.patch.set_facecolor("w")
        if out_dir is not None:
            fig.savefig(out_dir + 'torsionE_OverlayPlt_{}.png'.format(title), dpi=300, bbox_inches='tight',
                        bbox_extra_artists=(img,))
            plt.close('all')
    elif out_dir is not None:
        fig.savefig(out_dir + 'torsionE_OverlayPlt_{}.png'.format(title), dpi=300)
        plt.close('all')


def average_plot(mol_list, title, out_dir):
    """
    Plot average plot.
    This function plots the average PE curve (with error bars) of all molecules objects in mol_list

    :param mol_list: list of molecules as MoleculeRot objects
    :param title: str title for plot
    :param out_dir: path to directory in which to save figure
    :return: None
    """
    data = np.zeros([19, 3])
    for n, d in enumerate(np.linspace(0, 180, 19)):
        data[n][0] = d
        energies = []
        for mol in mol_list:
            phi, energy = mol.norm_energy_dict.keys(), mol.norm_energy_dict.values()
            try:
                energies.append(mol.norm_energy_dict[d])
                np_energies = np.array(energies)
                data[n][1] = np.average(np_energies)
                data[n][2] = np.std(np_energies) / np.sqrt(len(np_energies))
            except KeyError:
                print('Error. Did not find energy for {} at {} degree rotation.'.format(mol.name, d))

    print('\n\nNumber of molecules averaged: {}'.format(len(mol_list)))
    fig = plt.figure()
    plt.scatter(data[:, 0], data[:, 1], label=title)
    plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], linestyle="None")
    plt.xlim(min(phi) - 3, max(phi) + 3)
    plt.xticks(np.linspace(start=0, stop=180, num=7))
    plt.ylim(top=15, bottom=-1)
    plt.yticks(np.linspace(start=0, stop=14, num=8))
    plt.xlabel("dihedral angle (degrees)")
    plt.ylabel("energy (kcal/mol)")
    plt.title(title)
    fig.patch.set_facecolor("w")
    plt.savefig(out_dir + 'torsionE_AvgPlt_{}.png'.format(title), dpi=300, bbox_inches='tight')
    plt.close('all')
