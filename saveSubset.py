import pandas as pd



# ===================================== User Inputs ======================================
loadSubset = False # True: load a smaller set, False: load the full set
saveSubset = True # Save a subset of the full dataset
saveCount = 100 # Save this many molecules
pathMolecules = '/path/structure.tsv' # Path: full dataset
pathSubset = 'path/structureSubset.tsv' # Path: partial dataset



# =================================== Define Functions ===================================
def loadMolecules(path):
    print(f'Loading Molecules:\n'
          f'     {path}\n\n')

    # Load the TSV file
    try:
        structures = pd.read_csv(path, sep='\t')  # Use '\t' for tab-separated values
        print(f'Loaded Molecules:\n'
             f'{structures.head(10)}')
        return structures
    except FileNotFoundError:
        print(f'File not found.')
    except Exception as e:
        print(f'ERROR:\n{e}')


def saveMolSubset(saveMol, saveCount):
    # Save a subset of molecules to a new TSV file
    subset = saveMol.head(saveCount)
    subset.to_csv(pathSubset, sep='\t', index=False)  # Save as TSV
    print(f'First {saveCount} rows saved to:\n'
          f'     {pathSubset}')



# ===================================== Run The Code =====================================
pd.set_option('display.max_columns', None)

# Load the data
if loadSubset:
    molecules = loadMolecules(path=pathSubset)
else:
    molecules = loadMolecules(path=pathMolecules)
    if saveSubset:
        saveMolSubset(saveMol=molecules, saveCount=100)
