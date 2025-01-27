import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patheffects import withStroke
import pandas as pd
import pubchempy as pcp
# from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import sys


# Branched:
    # Carboxylic Acid: (C(=O)O)

# Multiple Functionalization:
    # Fluorination (F3): FC(F)(F)



# ========================================== User Inputs =========================================
# Input 1: Target Molecule
inSMILES = ['C#CC(=O)O',
            'C=CC(=O)O']
inAddHydrogens = False

# Input 2: 3D Model
inPlotMolecule3D = False
inLabelAtoms = False
inPlotGrid = False
inFigureColor = '#FFFFFF'
inPlotSpheres = True
inResolution = 1000

# Input 3: 2D Image
inPlotMolecule2D = False
inImageSizeX = 1000
inImageSizeY = inImageSizeX
inSave2DFigure = False



# ======================================== Set Parameters ========================================
# Colors: Console
white = '\033[38;2;255;255;255m'
silver = '\033[38;2;204;204;204m'
purple = '\033[38;2;189;22;255m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'


# Define: Atom diameters
scale = 5
diameters = {
    'H': 53 * scale,
    'C': 77 * scale,
    'N': 75 * scale,
    'O': 66 * scale,
    'F': 64 * scale,
    'P': 110 * scale,
    'S': 104 * scale,
    'Cl': 99 * scale,
    'Br': 114 * scale
}

# Define: Atomic properties
radius = 0.4
atoms = {
    'H': {'radius': radius, 'color': '#CCCCCC'},
    'C': {'radius': radius, 'color': '#353535'},
    'N': {'radius': radius, 'color': '#52B7FF'},
    'O': {'radius': radius, 'color': '#CE1D1D'},
    'F': {'radius': radius, 'color': '#FF8000'},
    'P': {'radius': radius, 'color': '#7E1DA0'},
    'S': {'radius': radius, 'color': '#E1B831'},
    'Cl': {'radius': radius, 'color': '#FF0080'},
    'Br': {'radius': radius, 'color': '#03428E'}
}
print(f'Avalible Atoms:')
for atom, properties in atoms.items():
    print(f'     {red}{atom}{resetColor}:{silver} {properties}')
print(f'{resetColor}\n')

# Print options
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 2000)
pd.set_option('display.float_format', '{:,.6f}'.format)



# ======================================= Define Functions =======================================
def pressKey(event):
    if event.key == 'escape':
        plt.close()



def plot3DBond(ax, coord1, coord2, radius, color, resolution=20):
    # Parameters:
    # - ax: Matplotlib 3D axis.
    # - coord1: Coordinates of the first atom as a tuple (x, y, z).
    # - coord2: Coordinates of the second atom as a tuple (x, y, z).
    # - radius: Radius of the cylinder (bond thickness).
    # - color: Color of the bond.
    # - resolution: Number of faces for the cylinder (default: 20).

    from matplotlib.colors import to_rgba
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection



    # Vector from coord1 to coord2
    v = np.array(coord2) - np.array(coord1)
    length = np.linalg.norm(v)

    # Normalize vector
    v /= length

    # Generate orthogonal vectors for the circular cross-section
    not_v = np.array([1, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1, 0])
    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)

    # Create the circle
    theta = np.linspace(0, 2 * np.pi, resolution)
    circle = radius * np.cos(theta)[:, None] * n1 + radius * np.sin(theta)[:, None] * n2

    # Create the cylinder
    bottom = np.array(coord1) + circle
    top = np.array(coord2) + circle
    for i in range(resolution):
        vertices = [
            bottom[i], bottom[(i + 1) % resolution],
            top[(i + 1) % resolution], top[i]
        ]
        ax.add_collection3d(Poly3DCollection([vertices], color=to_rgba(color, alpha=0.8)))



def stickFigure(compound):
    # Load molecule from SMILES
    mol = Chem.MolFromSmiles(compound)

    # Generate 2D coordinates
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

    # Draw the molecule
    img = Draw.MolToImage(mol, size=(inImageSizeX, inImageSizeY))
    img.show()
    if inSave2DFigure:
        img.save('molecule.png')



def plotMolecule(atomicCoords, properties, compound, labelMolecule):
    # Process data
    x = atomicCoords['X']
    y = atomicCoords['Y']
    z = atomicCoords['Z']


    # Make a figure
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Molecular Structure:\n{labelMolecule}', color='limegreen', fontsize=16)
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    fig.patch.set_facecolor(inFigureColor)
    ax.set_facecolor(inFigureColor)


    if inPlotSpheres:
        elements = list(atomicCoords['Elements'])
        # Draw spheres
        for index, (cx, cy, cz) in enumerate(zip(x, y, z)):
            element = elements[index]
            radius = atoms[element]['radius']  # Get the atomic radius
            color = atoms[element]['color']  # Get the atomic color
            # print(f'Plot atom:{white} {element}{resetColor}')

            # Create a 3D sphere
            u = np.linspace(0, 2 * np.pi, inResolution)
            v = np.linspace(0, 2 * np.pi, inResolution)
            sphereX = radius * np.outer(np.cos(u), np.sin(v)) + cx
            sphereY = radius * np.outer(np.sin(u), np.sin(v)) + cy
            sphereZ = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + cz

            # Plot the sphere
            ax.plot_surface(
                sphereX, sphereY, sphereZ,
                color=color, edgecolor='k', linewidth=0, alpha=1
            )


        # Plot bonds as cylinders
        for bond in compound.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()

            # Get coordinates for the two atoms involved in the bond
            coord1 = compound.GetConformer().GetAtomPosition(atom1)
            coord2 = compound.GetConformer().GetAtomPosition(atom2)

            # Determine bond type and radius
            bondType = bond.GetBondType()
            if bondType == Chem.rdchem.BondType.SINGLE:
                radius = 0.1
            elif bondType == Chem.rdchem.BondType.DOUBLE:
                radius = 0.15
            elif bondType == Chem.rdchem.BondType.TRIPLE:
                radius = 0.2
            else:
                radius = 0.1  # Default radius

            # Plot the cylinder bond
            plot3DBond(ax,
                       (coord1.x, coord1.y, coord1.z),
                       (coord2.x, coord2.y, coord2.z),
                       radius, color='black', resolution=20)
    else:
        # Draw circles
        ax.scatter(x, y, z, color=elementColors, s=elementSizes, label='atoms')


        # Plot the bonds as lines
        for bond in compound.GetBonds():
            atom1 = bond.GetBeginAtomIdx()
            atom2 = bond.GetEndAtomIdx()

            # Get coordinates for the two atoms involved in the bond
            coord1 = compound.GetConformer().GetAtomPosition(atom1)
            coord2 = compound.GetConformer().GetAtomPosition(atom2)

            # Determine bond type
            bondType = bond.GetBondType()
            if bondType == Chem.rdchem.BondType.SINGLE:
                linewidth = 3  # Standard line for single bonds
            elif bondType == Chem.rdchem.BondType.DOUBLE:
                linewidth = 4
            elif bondType == Chem.rdchem.BondType.TRIPLE:
                linewidth = 5

            # Plot bond
            ax.plot([coord1.x, coord2.x], [coord1.y, coord2.y], [coord1.z, coord2.z],
                    color='black', linewidth=linewidth)

    if inLabelAtoms:
        # Add atom labels
        outline = withStroke(linewidth=2, foreground='black')
        for i, atom in enumerate(compound.GetAtoms()):
            ax.text(x[i], y[i], z[i],  atom.GetSymbol(), color=properties[atom]['color'][i],
                    fontsize=18, fontweight='bold', path_effects=[outline])


    # Disable the grid
    if not inPlotGrid:
        ax.grid(False)
        ax.set_axis_off()

    # Get figure coordinates
    manager = plt.get_current_fig_manager()
    window = manager.window
    windowPosX = 605
    windowPosY = 65
    window.geometry(f"+{int(windowPosX)}+{int(windowPosY)}")
    window.title('Simulation')  # Change title
    window.configure(bg=inFigureColor)


    # # Set axis limits
    # Calculate: Molecular lenghts
    lenX = max(x) - min(x)
    lenY = max(y) - min(y)
    lenZ = max(z) - min(z)
    range = max(lenX, lenY, lenZ)

    # Find the center of the plot
    centerX, centerY, centerZ = ((max(x) + min(x)) / 2,
                                 (max(y) + min(y)) / 2,
                                 (max(z) + min(z)) / 2)

    # Set: Axis limits
    spacer = 0
    axisLimitX = (centerX + range / 2) + spacer
    axisLimitY = (centerY + range / 2) + spacer
    axisLimitZ = (centerZ + range / 2) + spacer

    # Set axis limits to be equal
    ax.set_xlim([-axisLimitX, axisLimitX])
    ax.set_ylim([-axisLimitY, axisLimitY])
    ax.set_zlim([-axisLimitZ, axisLimitZ])

    fig.canvas.mpl_connect('key_press_event', pressKey)
    fig.tight_layout()
    plt.show()



def evaluateMolecule(compound):
    print('================================= Molecular Coordinates '
          '=================================')
    # Load a molecule from a SMILES string
    if inAddHydrogens:
        mol = Chem.MolFromSmiles(compound)
        mol = Chem.AddHs(mol)
    else:
        from rdkit import RDLogger

        mol = Chem.MolFromSmiles(compound)

        # Suppress RDKit warnings
        logger = RDLogger.logger()
        logger.setLevel(RDLogger.ERROR)
    if mol is None:
        raise ValueError(f"{orange}ERROR: Invalid SMILES string{cyan} {compound}")

    # Get compound name
    resetCompoundName = False
    compoundPubChemPy = pcp.get_compounds(compound, 'smiles')[0]
    compoundName = compoundPubChemPy.iupac_name
    if compoundName is None:
        resetCompoundName = True
        compoundName = compound


    # Generate 3D coordinates for the molecule
    AllChem.EmbedMolecule(mol)

    # Get atomic coordinates
    coords = mol.GetConformer().GetPositions()
    coords = pd.DataFrame(coords,
                          columns=['X', 'Y', 'Z'])
    coords['Elements'] = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = coords[['Elements', 'X', 'Y', 'Z']] # Reorder columns
    if resetCompoundName:
        print(f'Coordinates:{pink} {compoundName}\n'
              f'{silver}{coords}{resetColor}\n\n')
    else:
        print(f'Coordinates:{pink} {compoundName}{resetColor} ({pink}{compound}{resetColor})\n'
              f'{silver}{coords}{resetColor}\n\n')


    # Plot the data
    if inPlotMolecule3D:
        plotMolecule(atomicCoords=coords, properties=atoms,
                     compound=mol, labelMolecule=compoundName)

    if inPlotMolecule2D:
        stickFigure(compound=compound)

    return coords, mol, compoundName



def alignMolecules(mol1, mol2):
    print('==================================== Align Sequences '
          '====================================')
    from Bio.PDB import PDBParser, Superimposer

    parser = PDBParser()
    structure1 = parser.get_structure(inSMILES[0], "molecule1.pdb")
    structure2 = parser.get_structure(inSMILES[1], "molecule2.pdb")

    atoms1 = [atom for atom in structure1.get_atoms()]
    atoms2 = [atom for atom in structure2.get_atoms()]

    # Align and calculate RMSD
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())

    print(f"RMSD: {super_imposer.rms}")


# ========================================= Run The Code =========================================
# Evaluate the molecule
coordinates1, molecule1, name1 = evaluateMolecule(compound=inSMILES[0])
coordinates2, molecule2, name2 = evaluateMolecule(compound=inSMILES[1])

# Align structures
alignMolecules(mol1=molecule1, mol2=molecule2)
