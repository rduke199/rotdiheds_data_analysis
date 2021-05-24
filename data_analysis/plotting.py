import os
import re
import json
import numpy as np
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from molecule_rot import MoleculeRot


def make_all_mol_list(master_dir, name_file='master_smiles.json'):
    """
    Make molecule object list.
    This function makes a list of MoleculeRot objects for all molecules in the energy_dir.

    :param master_dir: str path to master directory
    :param name_file: name of file from which to draw molecule names
    :return: list of all molecules as MoleculeRot objects.
    """
    all_molecules = []
    master_files = [x for x in os.listdir(master_dir) if x.startswith('master')]
    master_path = os.path.join(master_dir, name_file)
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

def make_small_mol_list(all_molecules, ring_num=None, unit_num=None, polymer_num=None, substituents=None,
                        chromophore=None, side_chain=None):
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
        if chromophore is not None:
            if mol.chromophore != chromophore:
                continue
        if side_chain is not None:
            if mol.side_chains != side_chain:
                continue
        if substituents is not None:
            if not bool(re.search(substituents, mol.substituents)):
                continue
        mol_list.append(mol)
    return mol_list

def overlay_plot(mol_list,
                 title,
                 x_values,
                 y_values,
                 x_labels,
                 y_labels,
                 plot_type,
                 x_min=None,
                 x_max=None,
                 y_min=None,
                 y_max=None,
                 out_dir=None,
                 varying_attribute='unit_num',
                 draw_1unit=False
    ):
    """
    Plot overlay plot.
    This function plots an overlay plot of all molecule PE curves in the mol_list

    :param mol_list: list of molecules as MoleculeRot objects
    :param title: str title
    :param
    :param out_dir: path to directory in which to save figure
    :param varying_attribute: str attribute that will be specified in legend
    :param draw_1unit: boolean. If Ture, the monomer structure will be drawn outside the plot
    :return: None
    """
    color_dict = {1: 'k', 3: 'c', 5: 'm', 7: 'y'}
    fig, ax = plt.subplots()
    for mol in mol_list:
        if mol.unit_num == 1:
            mol_image = mpimg.pil_to_array(mol.draw_structure())
        for x_val, y_val in zip(x_values, y_values):
            try:
                x = eval('mol.' + x_val)
                y = eval('mol.' + y_val)
                ax.scatter(x, y, label=eval('mol.' + varying_attribute), color=color_dict[mol.unit_num])
                ax.plot(x, y, color=color_dict[mol.unit_num])
            except:
                print('Error. Did not plot {} and {} for {} {} in {}.'.format(x_val, y_val, varying_attribute, eval('mol.'+varying_attribute), title))

    if x_min and x_max:
        ax.set_xlim(x_min, x_max)
    if y_min and y_max:
        ax.set_ylim(top=y_max, bottom=y_min)
    ax.set_xlabel(x_labels)
    ax.set_ylabel(y_labels)
    ax.set_title(title)
    fig.patch.set_facecolor("w")
    # sort both labels and handles by labels
    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles, labels)
    if draw_1unit:
        new_ax = fig.add_axes([0.7, .2, 0.4, 0.4], anchor='NE')
        img = new_ax.add_artist(AnnotationBbox(OffsetImage(mol_image, zoom=0.5), (1, 1)))
        new_ax.axis('off')
        fig.patch.set_facecolor("w")
        if out_dir is not None:
            fig.savefig(out_dir + '{}_OverlayPlt_{}.png'.format(plot_type, title), dpi=300, bbox_inches='tight',
                        bbox_extra_artists=(img,))
            plt.close('all')
    elif out_dir is not None:
        fig.savefig(out_dir + '{}_OverlayPlt_{}.png'.format(plot_type, title), dpi=300, bbox_inches='tight')
        plt.close('all')


def overlay_energy_plot(mol_list, title, out_dir=None, varying_attribute='unit_num', draw_1unit=False):
    x_values = ["norm_energy_dict.keys()"]
    y_values = ["norm_energy_dict.values()"]
    x_labels = "dihedral angle (degrees)"
    y_labels = "energy (kcal/mol)"
    y_min = -1
    y_max = 15
    plot_type = "torsionE"
    overlay_plot(mol_list=mol_list,
                 title=title,
                 x_values=x_values,
                 y_values=y_values,
                 x_labels=x_labels,
                 y_labels=y_labels,
                 plot_type=plot_type,
                 y_min=y_min,
                 y_max=y_max,
                 out_dir=out_dir,
                 varying_attribute=varying_attribute,
                 draw_1unit=draw_1unit
                 )

def overlay_homo_lumo_plot(mol_list, title, out_dir=None, varying_attribute='unit_num', draw_1unit=False):
    x_values = ['homo_dict.keys()','lumo_dict.keys()']
    y_values = ['homo_dict.values()','lumo_dict.values()']
    x_labels = "dihedral angle (degrees)"
    y_labels = "energy (kcal/mol)"
    plot_type = "HomoLumo"
    y_min = -20
    y_max = 5
    overlay_plot(mol_list=mol_list,
                 title=title,
                 x_values=x_values,
                 y_values=y_values,
                 x_labels=x_labels,
                 y_labels=y_labels,
                 plot_type=plot_type,
                 y_min=y_min,
                 y_max=y_max,
                 out_dir=out_dir,
                 varying_attribute=varying_attribute,
                 draw_1unit=draw_1unit
                 )


def plot_omega(mol_list, title, out_dir=None, varying_attribute='name', draw_1unit=False):
    x_values = ["unit_num"]
    y_values = ["tuned_omega"]
    x_labels = "Number of Monomer Units"
    y_labels = "Omega Value"
    plot_type = "omega"

    fig, ax = plt.subplots()
    plot_data = {}
    for mol in mol_list:
        if mol.unit_num == 1:
            mol_image = mpimg.pil_to_array(mol.draw_structure())
        for x_val, y_val in zip(x_values, y_values):
            try:
                x = eval('mol.' + x_val)
                y = eval('mol.' + y_val)
                plot_data[x] = y
            except:
                print('Error. Did not plot {} and {} for {} {} in {}.'.format(x_val, y_val, varying_attribute, eval('mol.'+varying_attribute), title))

    ax.scatter(plot_data.keys(), plot_data.values())
    ax.plot(plot_data.keys(), plot_data.values())

    ax.set_ylim(top=0.4, bottom=0.0)
    ax.set_xlabel(x_labels)
    ax.set_ylabel(y_labels)
    ax.set_title(title)
    fig.patch.set_facecolor("w")
    # sort both labels and handles by labels
    if draw_1unit:
        new_ax = fig.add_axes([0.7, .2, 0.4, 0.4], anchor='NE')
        img = new_ax.add_artist(AnnotationBbox(OffsetImage(mol_image, zoom=0.5), (1, 1)))
        new_ax.axis('off')
        fig.patch.set_facecolor("w")
        if out_dir is not None:
            fig.savefig(out_dir + '{}_OverlayPlt_{}.png'.format(plot_type, title), dpi=300, bbox_inches='tight',
                        bbox_extra_artists=(img,))
            plt.close('all')
    elif out_dir is not None:
        fig.savefig(out_dir + '{}_OverlayPlt_{}.png'.format(plot_type, title), dpi=300, bbox_inches='tight')
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
            try:
                phi, energy = mol.norm_energy_dict.keys(), mol.norm_energy_dict.values()
                energies.append(mol.norm_energy_dict[d])
                np_energies = np.array(energies)
                data[n][1] = np.average(np_energies)
                data[n][2] = np.std(np_energies) / np.sqrt(len(np_energies))
            except:
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
