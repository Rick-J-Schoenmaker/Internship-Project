import os
import shutil
import rdkit
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors

sdfFiles = []
smilesList = []
listRicksPDB = []
pdbPocketFile = []
pdbFiles = []

# folder with all subfolders
folders = []  # list with all the folders within the PDBBind 2018version folder.
files = []  # list with all the file paths


def openfolder():
    for entry in os.scandir("/home/rick/Documenten/Rick'spdb"): # Give pathway to PDBBind dataset.
        if entry.is_dir():
            folders.append(entry.path)  # fill list with folderpaths

    for value in folders:  # loop through the folders
        for entry in os.scandir(value):
            if entry.is_file():
                files.append(entry.path)  # fill list with filepaths

    for value in files:  # loop through filepaths
        if ".sdf" in value:
            sdfFiles.append(value)
        elif ".pdb" in value:
            pdbFiles.append(value)
    for value in pdbFiles:  # loop through pdbfiles and seperate on pocket or full file
        if "pocket.pdb" in value:
            pdbPocketFile.append(value)

    return;


def rdkit():

    HAcceptor = 12
    HDonor = 7
    LogP = 6
    Mol_Weight_Max = 650
    Mol_Weight_Min = 180

    lip_HA_count = 0
    lip_HD_count = 0
    lip_LogP_count = 0
    lip_MolW_count = 0
    count_element = 0

    family = ""
    pos = ""
    smiles = ""

    for x in sdfFiles:
        suppl = Chem.SDMolSupplier(x)   # Get supplement from SDF File.
        for mol in suppl:
            if mol is None: continue
            lip_HA_count = Lipinski.NumHAcceptors(mol) #set Hydrogen Acceptor count
            lip_HD_count = Lipinski.NumHDonors(mol)     # Get Hydrogen Donor count
            lip_LogP_count = Chem.Crippen.MolLogP(mol)  # Get LogP count.
            lip_MolW_count = Chem.Descriptors.ExactMolWt(mol)   # Get Molecule weight count.
            feats = []
            if lip_HA_count < HAcceptor and lip_HD_count < HDonor and lip_LogP_count < LogP and \
                Mol_Weight_Max > lip_MolW_count > Mol_Weight_Min:   # select ligand on custom cutoff.
                fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
                factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
                feats = factory.GetFeaturesForMol(mol)
                smiles = Chem.MolToSmiles(mol)  # Retrieve Smiles from mol.
                smilesList.append(Chem.MolToSmiles(mol))
                listRicksPDB.append(x)
            for f in feats:
                atom = mol.GetAtomWithIdx(f.GetAtomIds()[0]).GetSymbol()   # Get atoms from molecule.
                pos = f.GetPos()  # Get X,Y,Z positions from atoms in molecule.
                family = str(f.GetFamily())  # Get type of atom
    print(listRicksPDB)
    for value in listRicksPDB:
        shutil.move(str(value)[:-15], "/home/rick/Documenten/Rick'spdb")   # Give pathway to directory
    return;


def main():
    openfolder()
    rdkit()


main()