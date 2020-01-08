### Author: Rick Schoenmaker
### Version 1.9
### Last updated: 26-11-19
import math

import numpy
import psycopg2
import shutil

from Bio.PDB import *

import rdkit
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolfiles

import os

# List for file paths of specific files
mol2Files = []
sdfFiles = []
pdbFiles = []

listLigandName = []
listNumAtoms = []
listMolType = []
listChargeType = []
listBondSort = []
listSubstructures = []
mol2bondfiles = []
mol2substructurefiles = []
mol2Atomfiles = []

mol2threeletter =[]

listChain = []
listPocketfile = []
listPath = []
pdbProteinFile = []
pdbProteinName = []

pdbPocketFile = []
mol2File = []
smilesList = []

Names = []
atom_list = []
result_Three_Letter = []

listRicksPDB = []
group_type = []
dict_groupid = {}
dict_protein = {}
position_protein = []
group_protein = []
dict_ligand = {}
position_ligand = []
group_ligand = []
groups_protein_ligand = []
interaction_protein = []
interaction_ligand = []

Hydrophobic_list = ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"]

def openfolder():
    # folder with all subfolders
    folders = []  # list with all the folders within the PDBBind 2018version folder.
    files = []  # list with all the file paths
    for entry in os.scandir("/home/rick/Documenten/Rick'spdb"):
        if entry.is_dir():
            folders.append(entry.path)  # fill list with folderpaths

    for value in folders:  # loop through the folders
        for entry in os.scandir(value):
            if entry.is_file():
                files.append(entry.path)  # fill list with filepaths

    for value in files:  # loop through filepaths
        if ".mol2" in value:
            mol2Files.append(value)
        elif ".sdf" in value:
            sdfFiles.append(value)
        elif ".pdb" in value:
            pdbFiles.append(value)

    for value in pdbFiles:  # loop through pdbfiles and seperate on pocket or full file
        if "protein.pdb" in value:
            pdbProteinFile.append(value)
        elif "pocket.pdb" in value:
            pdbPocketFile.append(value)

    return;

# This function implements RDkit to retrieve atom binding information from sdf files.
def Protein_LigandData():
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
    count = 0
    family = ""
    pos = ""
    smiles = ""
    connection = ()

    for value in mol2Files:  # loop through all mol2 files in the mol2Files list.
        with open(value) as infile:
            copy = False
            for line in infile:
                if line.strip() == "@<TRIPOS>SUBSTRUCTURE":
                    copy = True
                    continue
                elif line.strip() == "@<TRIPOS>MOLECULE":
                    copy = False
                    continue
                elif copy:
                    mol2threeletter.append(
                        str(line[7:10]).strip("\n"))  # append all lines between substructure and molecule.
    while '' in mol2threeletter:
        mol2threeletter.remove('')
    while '   ' in mol2threeletter:
        mol2threeletter.remove('   ')
    try:
        connection = psycopg2.connect(user="postgres",
                                      password="Gerbils1",
                                      host="localhost",
                                      port="5433",
                                      database="postgres")
        cursor = connection.cursor()
        cursor.execute("DELETE FROM end_concept.atoms *")
        cursor.execute("DELETE FROM end_concept.group_type *")
        cursor.execute("DELETE FROM end_concept.lig_conform *") # SQL statement to delete rows from ligand_atoms table.
        cursor.execute("DELETE FROM end_concept.group_atoms *")
        cursor.execute("DELETE FROM end_concept.pro_conform *") # SQL statement to delete rows from protein_atoms table.
        sqlla = "INSERT INTO end_concept.lig_conform(pdb_id, lig_id, atom_id, x, y, z)" \
                " VALUES (%s,%s,%s,%s,%s,%s)"  # SQL statement to fill the ligand_atoms.
        sqlg = "INSERT INTO end_concept.group_type(group_type_id, type)  VALUES(%s,%s) "
        sqlatom = "INSERT INTO end_concept.atoms(atom_id, atom) VALUES (%s,%s)"
        sqlga = "INSERT INTO end_concept.group_atoms(group_type_id_fk, atom_id_fk,group_id,pseudo_x,pseudo_y,pseudo_z) VALUES (%s,%s,%s,%s,%s,%s)"
        sqlpro = "INSERT INTO end_concept.pro_conform(pdb_id,chain_id,x,y,z,atom_id) VALUES (%s,%s,%s,%s,%s,%s)"
        sqlinter = "INSERT INTO end_concept.interactions(group_id_li, group_id_pro) VALUES (%s,%s)"
        count_ligand_atoms = 0
        endfamily = []
        countgroup = -1
        countligandid = -1
        listPDBSub= []
        countdonor = 0
        countacceptor = 1
        countaromatic = 2
        counthydrophobe = 3
        countlumpedhydrophobe = 4
        countposionizable = 5
        countznbinder = 6
        countnegionizable =7
        countamino = 0
        for i in pdbPocketFile:
            count += 1
            print(count)
            pdbmol = rdkit.Chem.rdmolfiles.MolFromPDBFile(i)
            feats_protein = []
            fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
            factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
            feats_protein = factory.GetFeaturesForMol(pdbmol)
            for fp in feats_protein:
                atom_protein = pdbmol.GetAtomWithIdx(fp.GetAtomIds()[0]).GetSymbol()  # Get atoms from molecule.
                pos_protein = fp.GetPos()  # Get X,Y,Z positions from atoms in molecule.
                family_protein = str(fp.GetFamily())  # Get type of atom.
                conformer_protein = pdbmol.GetConformer()
                for atomid in fp.GetAtomIds():
                    count_ligand_atoms += 1
                    valatompro = (count_ligand_atoms, atom_protein)
                    cursor.execute(sqlatom, valatompro)
                    coords_protein = conformer_protein.GetAtomPosition(atomid)
                    at = pdbmol.GetAtomWithIdx(0)
                    rwm = Chem.RWMol(pdbmol)
                    idx = rwm.AddAtom(at)
                    newat = rwm.GetAtomWithIdx(idx)
                    aa = rwm.GetAtomWithIdx(atomid).GetPDBResidueInfo().GetResidueName()
                    countamino += 1
                    valpro = (i[37:41].upper(), str(newat.GetPDBResidueInfo().GetChainId()), coords_protein.x,
                              coords_protein.y, coords_protein.z, count_ligand_atoms)  # Values for SQL statement.
                    cursor.execute(sqlpro, valpro)  # Execute SQL statement.
                    if "Donor" in family_protein:
                        valga = (1, count_ligand_atoms, countdonor, float(pos_protein.x), float(pos_protein.y), float(pos_protein.z))  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countdonor)
                        groups_protein_ligand.append(countdonor)
                        interaction_protein.append(i[35:39])
                        group_type.append(1)
                        countdonor += 8
                    if "Acceptor" in family_protein:
                        valga = (2, count_ligand_atoms, countacceptor, float(pos_protein.x), float(pos_protein.y), float(pos_protein.z))  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countacceptor)
                        groups_protein_ligand.append(countacceptor)
                        interaction_protein.append(i[35:39])
                        group_type.append(2)
                        countacceptor += 8
                    if "Aromatic" in family_protein:
                        valga = (3, count_ligand_atoms, countaromatic, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countaromatic)
                        groups_protein_ligand.append(countaromatic)
                        interaction_protein.append(i[35:39])
                        group_type.append(3)
                    if "Hydrophobe" in family_protein and aa in Hydrophobic_list:
                        valga = (4, count_ligand_atoms, counthydrophobe, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(counthydrophobe)
                        groups_protein_ligand.append(counthydrophobe)
                        interaction_protein.append(i[35:39])
                        group_type.append(4)
                        counthydrophobe += 8
                    if "LumpedHydrophobe" in family_protein and aa in Hydrophobic_list:
                        valga = (5, count_ligand_atoms, countlumpedhydrophobe, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countlumpedhydrophobe)
                        groups_protein_ligand.append(countlumpedhydrophobe)
                        interaction_protein.append(i[35:39])
                        group_type.append(5)
                    if "PosIonizable" in family_protein:
                        valga = (6, count_ligand_atoms, countposionizable, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countposionizable)
                        groups_protein_ligand.append(countposionizable)
                        interaction_protein.append(i[35:39])
                        group_type.append(6)
                    if "ZnBinder" in family_protein:
                        valga = (7, count_ligand_atoms, countznbinder, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countznbinder)
                        groups_protein_ligand.append(countznbinder)
                        interaction_protein.append(i[35:39])
                        group_type.append(7)
                    if "NegIonizable" in family_protein:
                        valga = (8, count_ligand_atoms, countnegionizable, pos_protein.x, pos_protein.y, pos_protein.z)  # Values for SQL statement.
                        cursor.execute(sqlga, valga)  # Execute SQL statement.
                        position_protein.append(pos_protein)
                        group_protein.append(countnegionizable)
                        groups_protein_ligand.append(countnegionizable)
                        interaction_protein.append(i[35:39])
                        group_type.append(8)
                countaromatic += 8
                countlumpedhydrophobe += 8
                countposionizable += 8
                countznbinder += 8
                countnegionizable += 8
        for x in sdfFiles:
            countligandid += 1
            suppl = Chem.SDMolSupplier(x)   # Get supplement from SDF File.
            for mol in suppl:
                if mol is None: continue
                lip_HA_count = Lipinski.NumHAcceptors(mol) #set Hydrogen Acceptor count
                lip_HD_count = Lipinski.NumHDonors(mol)     # Get Hydrogen Donor count
                lip_LogP_count = Chem.Crippen.MolLogP(mol)  # Get LogP count.
                lip_MolW_count = Chem.Descriptors.ExactMolWt(mol)   # Get Molecule weight count.
                feats_ligand = []

                if lip_HA_count < HAcceptor and lip_HD_count < HDonor and lip_LogP_count < LogP and \
                    Mol_Weight_Max > lip_MolW_count > Mol_Weight_Min:   # select ligand on custom cutoff.
                    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
                    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
                    feats_ligand = factory.GetFeaturesForMol(mol)
                    smiles = Chem.MolToSmiles(mol)  # Retrieve Smiles from mol.
                    smilesList.append(Chem.MolToSmiles(mol))
                    listPocketfile.append(str(x).replace("ligand.sdf", "pocket.pdb"))
                    mol2File.append(str(x).replace(".sdf", ".mol2"))
                    listRicksPDB.append(x)
                for fl in feats_ligand:
                    atom = mol.GetAtomWithIdx(fl.GetAtomIds()[0]).GetSymbol()   # Get atoms from molecule.
                    pos_ligand = fl.GetPos()  # Get X,Y,Z positions from atoms in molecule.
                    family_ligand = str(fl.GetFamily())  # Get type of atom.
                    conformer_ligand = mol.GetConformer()
                    if family_ligand not in endfamily:
                        countgroup += 1
                        endfamily.append(family_ligand)
                        valg = (countgroup+1, endfamily[countgroup])
                        cursor.execute(sqlg, valg)
                    for atomid in fl.GetAtomIds():
                        count_ligand_atoms += 1
                        valatom = (count_ligand_atoms, atom)
                        cursor.execute(sqlatom, valatom)
                        coords_ligand = conformer_ligand.GetAtomPosition(atomid)
                        vall = (x[37:41].upper(), mol2threeletter[countligandid], count_ligand_atoms, coords_ligand.x, coords_ligand.y, coords_ligand.z)   # Values for SQL statement.
                        cursor.execute(sqlla, vall)  # Execute SQL statement.
                        if "Donor" in family_ligand:
                            valga = (1, count_ligand_atoms, countdonor, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countdonor)
                            groups_protein_ligand.append(countdonor)
                            interaction_ligand.append(x[35:39])
                            group_type.append(1)
                            countdonor += 8
                        if "Acceptor" in family_ligand:
                            valga = (2, count_ligand_atoms, countacceptor, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countacceptor)
                            groups_protein_ligand.append(countacceptor)
                            interaction_ligand.append(x[35:39])
                            group_type.append(2)
                            countacceptor += 8
                        if "Aromatic" in family_ligand:
                            valga = (3, count_ligand_atoms, countaromatic, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countaromatic)
                            groups_protein_ligand.append(countaromatic)
                            interaction_ligand.append(x[35:39])
                            group_type.append(3)
                        if "Hydrophobe" in family_ligand and aa in Hydrophobic_list:
                            valga = (4, count_ligand_atoms, counthydrophobe, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(counthydrophobe)
                            groups_protein_ligand.append(counthydrophobe)
                            interaction_ligand.append(x[35:39])
                            group_type.append(4)
                            counthydrophobe += 8
                        if "LumpedHydrophobe" in family_ligand and aa in Hydrophobic_list:
                            valga = (5, count_ligand_atoms, countlumpedhydrophobe, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countlumpedhydrophobe)
                            groups_protein_ligand.append(countlumpedhydrophobe)
                            interaction_ligand.append(x[35:39])
                            group_type.append(5)
                        if "PosIonizable" in family_ligand:
                            valga = (6, count_ligand_atoms, countposionizable, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countposionizable)
                            groups_protein_ligand.append(countposionizable)
                            interaction_ligand.append(x[35:39])
                            group_type.append(6)
                        if "ZnBinder" in family_ligand:
                            valga = (7, count_ligand_atoms, countznbinder, pos_ligand.x, pos_ligand.y, pos_ligand.z)  # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countznbinder)
                            groups_protein_ligand.append(countznbinder)
                            interaction_ligand.append(x[35:39])
                            group_type.append(7)
                        if "NegIonizable" in family_ligand:
                            valga = (8, count_ligand_atoms, countnegionizable, pos_ligand.x, pos_ligand.y, pos_ligand.z) # Values for SQL statement.
                            cursor.execute(sqlga, valga)  # Execute SQL statement.
                            position_ligand.append(pos_ligand)
                            group_ligand.append(countnegionizable)
                            groups_protein_ligand.append(countnegionizable)
                            interaction_ligand.append(x[35:39])
                            group_type.append(8)
                    countaromatic += 8
                    countlumpedhydrophobe += 8
                    countposionizable += 8
                    countznbinder += 8
                    countnegionizable += 8
        connection.commit()
    except (Exception, psycopg2.Error) as error:
        print("Error while connecting to PostgreSQL", error)


    finally:
        # closing database connection.
        if (connection):
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")


def Mol2Splitter():
    mol2moleculefiles = []
    filteredmol2moleculefiles = []
    list_length = 0

    for value in mol2File:  # loop through all mol2 files in the mol2Files list.
        with open(value) as infile:
            copy = False
            for line in infile:
                if line.strip() == "@<TRIPOS>MOLECULE":
                    copy = True
                    continue
                elif line.strip() == "@<TRIPOS>ATOM":
                    copy = False
                    continue
                elif copy:
                    mol2moleculefiles.append(line.strip("\n"))  # append all lines between molecule and atom.

        with open(value) as infile:
            copy = False
            for line in infile:
                if line.strip() == "@<TRIPOS>ATOM":
                    mol2Atomfiles.append("###")
                    copy = True
                    continue
                elif line.strip() == "@<TRIPOS>BOND":
                    copy = False
                    continue
                elif copy:
                    mol2Atomfiles.append(line.strip("\n"))  # append all lines between atom and bond.

        with open(value) as infile:
            mol2bondfiles.append("###")
            copy = False
            for line in infile:
                if line.strip() == "@<TRIPOS>BOND":
                    copy = True
                    continue
                elif line.strip() == "@<TRIPOS>SUBSTRUCTURE":
                    copy = False
                    continue
                elif copy:
                    mol2bondfiles.append(line.strip("\n"))  # append all lines between bond and substructure.
        with open(value) as infile:
            copy = False
            for line in infile:
                if line.strip() == "@<TRIPOS>SUBSTRUCTURE":
                    copy = True
                    continue
                elif line.strip() == "@<TRIPOS>MOLECULE":
                    copy = False
                    continue
                elif copy:
                    mol2substructurefiles.append(
                        str(line[7:10]).strip("\n"))  # append all lines between substructure and molecule.
    while '' in mol2substructurefiles:
        mol2substructurefiles.remove('')


    # Filter out empty objects from the list.
    filteredmol2moleculefiles = [x for x in mol2moleculefiles if len(x.strip('')) > 0]
    # Get length of filtered molecule files
    list_length = len(filteredmol2moleculefiles)
    for protein in filteredmol2moleculefiles[0:list_length:4]:
        listLigandName.append(protein)  # fill list with ligandnames.

    for atoms in filteredmol2moleculefiles[1:list_length:4]:
        listNumAtoms.append(atoms)  # fill list with atom numbers
    for mol_type in filteredmol2moleculefiles[2:list_length:4]:
        listMolType.append(mol_type)  # fill list with molecule types

    for charge_type in filteredmol2moleculefiles[3:list_length:4]:
        listChargeType.append(charge_type)  # fill list with charge type

    return;

# Function that parses PDB protein files.
def PDBProteinParser():
    parser = PDBParser()  # Call the PDBParser function
    pdbPocketFile.sort()
    for value in pdbPocketFile:
        with open(value) as infile:
            for line in infile:
                if line.__contains__("HEADER"):
                    listPath.insert(0, line.strip("HEADER").strip("\n"))
                    listPath.insert(1, value)
                    continue


    return;

def Interactions():
    try:
        connection = psycopg2.connect(user="postgres",
                                      password="Gerbils1",
                                      host="localhost",
                                      port="5433",
                                      database="postgres")
        cursor = connection.cursor()
        cursor.execute("DELETE FROM end_concept.interactions *")
        sqlinter = "INSERT INTO end_concept.interactions(group_id_li, group_id_pro) VALUES (%s,%s)"
        dict_protein = {kp: vp for kp, vp in zip(position_protein, group_protein)}
        dict_ligand = {kl: vl for kl, vl in zip(position_ligand, group_ligand)}
        dict_protein_name = {kp: vp for kp, vp in zip(position_protein, interaction_protein)}
        dict_ligand_name = {kl: vl for kl, vl in zip(position_ligand, interaction_ligand)}
        dict_groupid = {kg: vg for kg, vg in zip(groups_protein_ligand, group_type)}
        for kp, vp in dict_protein_name.items():
            for kl, vl in dict_ligand_name.items():
                if vl == vp:
                    a = numpy.array(kp)
                    b = numpy.array(kl)
                    dist = numpy.linalg.norm(a - b)
                    # Selects aromatic interactions with a 5 angstrom cut-off
                    if dict_groupid.get(dict_ligand.get(kl)) == 3 and dict_groupid.get(dict_protein.get(kp)) == 3:
                        if dist <= 5:
                            valinter = (dict_ligand.get(kl), dict_protein.get(kp))
                            cursor.execute(sqlinter, valinter)
                    #  Selects hydrophobic interactions with a 4.5 angstrom cut-off
                    if dict_groupid.get(dict_ligand.get(kl)) == 4 and dict_groupid.get(dict_protein.get(kp)) == 4 or \
                    dict_groupid.get(dict_ligand.get(kl)) == 5 and dict_groupid.get(dict_protein.get(kp)) == 5 or \
                    dict_groupid.get(dict_ligand.get(kl)) == 4 and dict_groupid.get(dict_protein.get(kp)) == 5 or \
                    dict_groupid.get(dict_ligand.get(kl)) == 5 and dict_groupid.get(dict_protein.get(kp)) == 4:
                        if dist <= 4.5:
                            valinter = (dict_ligand.get(kl), dict_protein.get(kp))
                            cursor.execute(sqlinter, valinter)
                    # Selects donor-acceptor interactions with a 4.4 angstrom cut-off
                    if dict_groupid.get(dict_ligand.get(kl)) == 1 and dict_groupid.get(dict_protein.get(kp)) == 2 or\
                       dict_groupid.get(dict_ligand.get(kl)) == 2 and dict_groupid.get(dict_protein.get(kp)) == 1:
                        if dist <= 4.4:
                            valinter = (dict_ligand.get(kl), dict_protein.get(kp))
                            cursor.execute(sqlinter, valinter)
                    # Selects cation - anion interactions and cation - aromatic interactions:
                    if dict_groupid.get(dict_ligand.get(kl)) == 8 and dict_groupid.get(dict_protein.get(kp)) == 6 or\
                       dict_groupid.get(dict_ligand.get(kl)) == 3 and dict_groupid.get(dict_protein.get(kp)) == 6 or\
                       dict_groupid.get(dict_ligand.get(kl)) == 6 and dict_groupid.get(dict_protein.get(kp)) == 8 or\
                       dict_groupid.get(dict_ligand.get(kl)) == 6 and dict_groupid.get(dict_protein.get(kp)) == 3:
                        if dist <= 4.5:
                            valinter = (dict_ligand.get(kl), dict_protein.get(kp))
                            cursor.execute(sqlinter, valinter)
        connection.commit()
    except (Exception, psycopg2.Error) as error:
        print("Error while connecting to PostgreSQL", error)


    finally:
        # closing database connection.
        if (connection):
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")

# This function fills a database with retrieved and filtered data from the PDB-BIND data files.
def ProteinData():
    connection = ()
    # Try to make connection and fill a database.
    try:
        connection = psycopg2.connect(user="postgres",
                                      password="Gerbils1",
                                      host="localhost",
                                      port="5433",
                                      database="postgres")

        cursor = connection.cursor()
        cursor.execute("DELETE FROM end_concept.protein *")  # SQL statement to delete rows from protein table.
        cursor.execute("DELETE FROM end_concept.interacting_ligand_protein *")  # SQL statement to delete rows from ligand table.
        cursor.execute("DELETE FROM end_concept.pro_conform *")  # SQL statement to delete rows from protein_atoms table.

        sqll = "INSERT INTO end_concept.interacting_ligand_protein(ligand_id,protein_id) VALUES (%s,%s)"  # SQL statement for ligand data.
        sqlP = "INSERT INTO end_concept.protein(protein_id, seq) VALUES (%s,%s)"  # SQL statement for protein data.
        sqla = "INSERT INTO end_concept.pro_conform(pdb_id,chain_id,x,y,z,atom_id) VALUES (%s,%s,%s,%s,%s,%s)"
        sqls = "INSERT INTO end_concept.ligand(ligand_id,smiles) VALUES (%s,%s)"  # SQL statement to fill smiles columns
        sqlatom = "INSERT INTO end_concept.atoms(atom_id, atom_type, atom) VALUES (%s,%s,%s)"

        parser = PDBParser()
        countheader = 0
        countpath = 1
        lenlistpath = len(listPath)
        lenlistpath = lenlistpath / 2
        length = len(listPath)
        count = 0
        count_atom_id = 215607

        ###Ligand_variables
        dict_three_letter = {}
        get_id = []
        get_str = ""
        id_list_filler = []
        id_list = []
        countli = 5
        counter = 0
        chain_count = 0
        with open("/home/rick/Documenten/Internship-project/Internship-project/cc-to-pdb.tdd") as infile:
            rows = (line.strip("\n").split('\t') for line in infile)
            dict_three_letter = {row[1]: row[0:] for row in rows}
        while count < lenlistpath:
            # countligandname += 1
            length -= 1
            data = parser.get_structure(id=listPath[countheader], file=listPath[countpath])
            atom_list = []
            count = count + 1
            count_atom = 0

            model = data.get_models()
            models = list(model)
            chains = list(models[0].get_chains())
            residue = list(chains[0].get_residues())
            atoms = list(residue[0].get_atoms())

            for model in data:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_list.append(atom)
            ppb = PPBuilder()
            peptides = []
            for pp in ppb.build_peptides(data):
                peptides.append(str(pp.get_sequence()))
            str_peptides = ''.join(str(e) for e in peptides)
            peptide_list = str_peptides.split("\n")
            #PROTEINS
            valp = (
                str(data.get_full_id()).replace("_POCKET", "").strip("(").strip(")").replace("'", "").replace(",", ""),
                peptide_list[0:count].__str__().strip("{").strip("}").strip("[").strip("]")
                    .strip("'"))
            cursor.execute(sqlP, valp)
            str(atom_list[count_atom]).strip("<>").strip("Atom").strip("'").strip(",")
            #PROTEIN-ATOMS
            with open(listPath[countpath]) as infile:
                for chainpdb in infile:
                    if chainpdb.__contains__("ATOM"):
                        listChain.append(chainpdb[21:22])
                    if chainpdb.__contains__("HETATM"):
                        listChain.append(chainpdb[21:22])

            countheader = countheader + 2
            countpath = countpath + 2
            for f in atom_list:
                coordinates = list(atom_list[count_atom].get_coord())
                if count_atom != len(atom_list):
                    vala = (str(data.get_full_id()).replace(" ", "").replace("_POCKET", "").strip("(").strip(")")
                            .replace("'", "").replace(",", ""), listChain[chain_count], float(coordinates[0]),
                            float(coordinates[1]), float(coordinates[2]), count_atom_id)
                    cursor.execute(sqla, vala)
                    chain_count +=1
                    valatom = (count_atom_id, str(data.get_full_id()).replace(" ", "").replace("_POCKET", "").strip("(")
                               .strip(")").replace("'", "").replace(",", ""), str(atom_list[count_atom]).strip("<>")
                               .strip("Atom").strip("'").strip(","))
                    cursor.execute(sqlatom, valatom)
                    count_atom += 1
                    count_atom_id += 1
                if count_atom == len(atom_list):
                    count_atom = 0
            #Binding_Protein_Ligand
            for value in pdbPocketFile:
                get_id.append(value.split("/"))
            get_str = ''.join(str(e) for e in get_id)
            get_str.strip("[]").strip("" "")
            id_list_filler = get_str.split(",")
            count_id = -1
            count_list = 0
            result_Three_Letter_final = []
            while counter < len(pdbPocketFile):
                id_list.append(id_list_filler[countli])
                countli += 6
                counter += 1
                count_while = 0
            while count_while < len(id_list):
                count_while += 1
                count_id += 1
                for key, value in dict_three_letter.items():
                    if id_list[count_id].__str__().strip("' '") in key.lower():
                        valln = (str(id_list[count_id].replace("'", "")).upper(), value[0])
                        count_list += 1
                        cursor.execute(sqll, valln)
                        print(count_list)
        #Ligand
        countsmiles = -1
        for value in smilesList:
            if countsmiles != len(smilesList):
                countsmiles += 1
                vals = (mol2substructurefiles[countsmiles], smilesList[countsmiles])  # Values for SQL statement.
            cursor.execute(sqls, vals)  # Execute SQL statement.
        cursor.execute("DELETE FROM end_concept.ligand a USING ( SELECT MIN(ctid) as ctid, ligand_id FROM end_concept.ligand GROUP BY ligand_id HAVING COUNT(*) > 1) b WHERE a.ligand_id = b.ligand_id AND a.ctid <> b.ctid")
        connection.commit()

    except (Exception, psycopg2.Error) as error:
        print("Error while connecting to PostgreSQL", error)
    finally:
        # closing database connection.
        if (connection):
            cursor.close()
            connection.close()
            print("PostgreSQL connection is closed")




def main():
    openfolder()
    Protein_LigandData()
    Mol2Splitter()
    PDBProteinParser()
    Interactions()
    # ProteinData()



main()
