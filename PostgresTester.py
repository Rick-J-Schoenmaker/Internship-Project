# Author: Rick Schoenmaker
# Date: 12-12-2019
# Function: This script extracts data from the PROLI database and calculates numbers such as
# the amount of interactions and averages.
import os

import numpy
import psycopg2

try:
    connection = psycopg2.connect(user="postgres",
                                  password="Gerbils1",
                                  host="localhost",
                                  port="5433",
                                  database="postgres")
    cursor = connection.cursor()
    list6117 = []
    pdbids = []  # Initializing  a list for pdb ids
    interactionspdb = []
    donorLigand = []
    donorProtein = []
    acceptorLigand = []
    acceptorProtein = []
    hydrophobeLigand = []
    hydrophobeProtein = []
    hydrophobeLumpedLigand = []
    hydrophobeLumpedProtein = []
    aromaticLigand = []
    aromaticProtein = []
    posionizableLigand = []
    posionizableProtein = []
    negionizableLigand = []
    negionizableProtein = []
    listPosionizableLigandSortDup = []
    listPosionizableProteinSortDup = []
    listAromaticLigandSortDup = []
    listAromaticProteinSortDup = []


    def Statements():

        # Extract pdb ids from the database of the protein groups that interact joined on atom ids.
        cursor.execute(
            "SELECT end_concept.pro_conform.pdb_id FROM end_concept.group_atoms INNER JOIN end_concept.interactions on group_id_pro = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id ")
        pdbid = cursor.fetchall()  # fetch the resulting pdb ids from the SQL statement
        for i in pdbid:
            pdbids.append(i)

        # Extract pdb ids and interaction ids from the database of the ligand groups that interact joined on atom ids.
        cursor.execute(
            "SELECT end_concept.lig_conform.pdb_id, end_concept.interactions.inter_id FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON atom_id_fk = atom_id")
        interli = cursor.fetchall()
        for i in interli:
            interactionspdb.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.lig_conform.x, y, z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id WHERE group_type_id_fk = 1")
        donorli = cursor.fetchall()
        for i in donorli:
            donorLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.pro_conform.x, y, z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id WHERE group_type_id_fk = 1 ")
        donorpro = cursor.fetchall()
        for i in donorpro:
            donorProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id,end_concept.lig_conform.x, y, z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id WHERE group_type_id_fk = 2")
        acceptorli = cursor.fetchall()
        for i in acceptorli:
            acceptorLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id,end_concept.pro_conform.x, y, z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id WHERE group_type_id_fk =  2")
        acceptorpro = cursor.fetchall()
        for i in acceptorpro:
            acceptorProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id WHERE group_type_id_fk =  3")
        aromaticli = cursor.fetchall()
        for i in aromaticli:
            aromaticLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id  WHERE group_type_id_fk =  3")
        aromaticpro = cursor.fetchall()
        for i in aromaticpro:
            aromaticProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id WHERE group_type_id_fk =  4 ")
        hydrophobeli = cursor.fetchall()
        for i in hydrophobeli:
            hydrophobeLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id,end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id WHERE group_type_id_fk =  4 ")
        hydrophobepro = cursor.fetchall()
        for i in hydrophobepro:
            hydrophobeProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id WHERE group_type_id_fk =  5 ")
        hydrophobelumpedli = cursor.fetchall()
        for i in hydrophobelumpedli:
            hydrophobeLumpedLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id,end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id  WHERE group_type_id_fk =  5 ")
        hydrophobelumpedpro = cursor.fetchall()
        for i in hydrophobelumpedpro:
            hydrophobeLumpedProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id WHERE group_type_id_fk =  6")
        posionizableli = cursor.fetchall()
        for i in posionizableli:
            posionizableLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id  WHERE group_type_id_fk =  6")
        posionizablepro = cursor.fetchall()
        for i in posionizablepro:
            posionizableProtein.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id  WHERE group_type_id_fk = 8")
        negionizableli = cursor.fetchall()
        for i in negionizableli:
            negionizableLigand.append(i)

        cursor.execute(
            "SELECT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id WHERE group_type_id_fk = 8")
        negionizablepro = cursor.fetchall()
        for i in negionizablepro:
            negionizableProtein.append(i)


    # noinspection PySimplifyBooleanCheck,PySimplifyBooleanCheck
    def Proteinwithinteraction():
        folders = []  # list with all the folders within the PDBBind 2018version folder.
        for entry in os.scandir("/home/rick/Documenten/Rick'spdb"):
            if entry.is_dir():
                folders.append(str(entry.path)[32:40].upper())  # fill list with folderpaths
        folders.sort()  # Sort the list with folders
        listpdbidsort = []  # Initializing  a list for sorted pdb ids without duplicates
        count = 0  # Initializing  a count
        countpdbid = -1
        for i in pdbids:
            countpdbid += 1
            if listpdbidsort.__contains__(pdbids[countpdbid]) == False:
                listpdbidsort.append(pdbids[countpdbid])
                count += 1
        listpdbidsort.sort()
        for i in listpdbidsort:
            list6117.append(str(i).strip('(').strip(')').strip("'").replace("'", "").replace(",", ""))
        list6117.sort()
        set1 = set(folders)
        set2 = set(list6117)
        missing = list(sorted(set1 - set2))
        print(missing)


    def Calcaverageinteraction():
        strinteractionspdb = ''.join(str(e) for e in interactionspdb)
        strinteractionspdb = strinteractionspdb.replace(")", ',')
        listinteractionspdb = strinteractionspdb.split(",")
        listinteractionspdb = [s.replace("'", '') for s in listinteractionspdb]
        listinteractionspdb = [s.replace("(", '') for s in listinteractionspdb]
        listinteractionspdb = [s.strip("''") for s in listinteractionspdb]

        countlistinter = len(listinteractionspdb) - 2
        listcountperstructure = []
        countaverageinter = 0
        templist = []
        countpdb = -2
        countinter = -1
        interlist = []
        while countinter != countlistinter:
            countpdb += 2
            countinter += 2
            if templist.__contains__(listinteractionspdb.__getitem__(countinter)) == False:
                templist.append(listinteractionspdb.__getitem__(countinter))
                interlist.append(listinteractionspdb.__getitem__(countpdb))
        print(len(templist))
        print(interlist[0:20])
        count = -1
        while count != len(list6117) - 1:
            count += 1
            print(count)
            listcountperstructure.append(list6117.__getitem__(count))
            listcountperstructure.append(interlist.count(list6117.__getitem__(count)))
            countaverageinter += interlist.count(list6117.__getitem__(count))
        print(countaverageinter)
        countaverageinter = countaverageinter / len(list6117)
        print(countaverageinter)
        print(len(listcountperstructure))

        with open("/home/rick/countinterperstuc.csv", "w") as outfile:
            for entries in listcountperstructure:
                outfile.write(str(entries))
                outfile.write("\n")

    def DonorAcceptor():
        # donor Ligand - acceptor Protein interaction calculating function that prepares the list for extracting distances.
        donorLigand.sort()
        acceptorProtein.sort()
        strdonorligand = ''.join(str(e) for e in donorLigand)
        strdonorligand = strdonorligand.replace(")", ',')
        listDonorLigandSort = strdonorligand.split(",")
        listDonorLigandSort = [s.replace("'", '') for s in listDonorLigandSort]
        listDonorLigandSort = [s.replace("(", '') for s in listDonorLigandSort]

        strAcceptorProtein = ''.join(str(e) for e in acceptorProtein)
        strAcceptorProtein = strAcceptorProtein.replace(")", ',')
        listAcceptorProteinSort = strAcceptorProtein.split(",")
        listAcceptorProteinSort = [s.replace("'", '') for s in listAcceptorProteinSort]
        listAcceptorProteinSort = [s.replace("(", '') for s in listAcceptorProteinSort]
        countx = 1
        county = 2
        countz = 3
        count = 0
        countdistancedalp = 0
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        while count != len(donorLigand):
            a = numpy.array((float(listDonorLigandSort[countx]), float(listDonorLigandSort[county]),
                             float(listDonorLigandSort[countz])))
            b = numpy.array((float(listAcceptorProteinSort[countx]), float(listAcceptorProteinSort[county]),
                             float(listAcceptorProteinSort[countz])))
            dist = numpy.linalg.norm(a - b)
            countdistancedalp += dist
            countx += 4
            county += 4
            countz += 4
            count += 1
            if 0 <= dist <= 0.5:
                count0 += 1
            if 0.5 <= dist <= 1:
                count1 += 1
            if 1 <= dist <= 1.5:
                count2 += 1
            if 1.5 <= dist <= 2:
                count3 += 1
            if 2 <= dist <= 2.5:
                count4 += 1
            if 2.5 <= dist <= 3:
                count5 += 1
            if 3 <= dist <= 3.5:
                count6 += 1
            if 3.5 <= dist <= 4:
                count7 += 1
            if 4 <= dist <= 4.5:
                count8 += 1
            if 4.5 <= dist <= 5:
                count9 += 1
        print("Average distance between LigandDonor and ProteinAcceptor:")
        print(countdistancedalp / len(donorLigand))
        print("Distances LigandDonor and ProteinAcceptor")
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        # acceptor Ligand - donor Protein interaction calculating function that prepares the list for extracting distances.
        acceptorLigand.sort()
        donorProtein.sort()
        strAcceptorLigand = ''.join(str(e) for e in acceptorLigand)
        strAcceptorLigand = strAcceptorLigand.replace(")", ',')
        listAcceptorLigandSort = strAcceptorLigand.split(",")
        listAcceptorLigandSort = [s.replace("'", '') for s in listAcceptorLigandSort]
        listAcceptorLigandSort = [s.replace("(", '') for s in listAcceptorLigandSort]

        strDonorProtein = ''.join(str(e) for e in donorProtein)
        strDonorProtein = strDonorProtein.replace(")", ',')
        listDonorProteinSort = strDonorProtein.split(",")
        listDonorProteinSort = [s.replace("'", '') for s in listDonorProteinSort]
        listDonorProteinSort = [s.replace("(", '') for s in listDonorProteinSort]
        countx = 1
        county = 2
        countz = 3
        count = 0
        countdistanceadlp = 0
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        while count != len(donorProtein):
            a = numpy.array((float(listAcceptorLigandSort[countx]), float(listAcceptorLigandSort[county]),
                             float(listAcceptorLigandSort[countz])))
            b = numpy.array((float(listDonorProteinSort[countx]), float(listDonorProteinSort[county]),
                             float(listDonorProteinSort[countz])))
            dist = numpy.linalg.norm(a - b)
            countdistanceadlp += dist
            countx += 4
            county += 4
            countz += 4
            count += 1
            if 0 <= dist <= 0.5:
                count0 += 1
            if 0.5 <= dist <= 1:
                count1 += 1
            if 1 <= dist <= 1.5:
                count2 += 1
            if 1.5 <= dist <= 2:
                count3 += 1
            if 2 <= dist <= 2.5:
                count4 += 1
            if 2.5 <= dist <= 3:
                count5 += 1
            if 3 <= dist <= 3.5:
                count6 += 1
            if 3.5 <= dist <= 4:
                count7 += 1
            if 4 <= dist <= 4.5:
                count8 += 1
            if 4.5 <= dist <= 5:
                count9 += 1
        print("Average distance between LigandAcceptor and ProteinDonor:")
        print(countdistanceadlp / len(donorProtein))
        print("Average distance between all Acceptor and Donor interactions:")
        countdis = countdistanceadlp / len(donorProtein) + countdistancedalp / len(donorLigand)
        print(countdis / 2)

        print("Distances LigandAcceptor and ProteinDonor")
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        print("Total number of interactions between Donors and Acceptors:")
        print(len(donorProtein) + len(donorLigand))


    def Hydrophobic():
        ## Interaction between Hydrophobic atoms.
        hydrophobeLigand.sort()
        hydrophobeProtein.sort()
        strHydrophobeLigand = ''.join(str(e) for e in hydrophobeLigand)
        strHydrophobeLigand = strHydrophobeLigand.replace(")", ',')
        listHydrophobeLigandSort = strHydrophobeLigand.split(",")
        listHydrophobeLigandSort = [s.replace("'", '') for s in listHydrophobeLigandSort]
        listHydrophobeLigandSort = [s.replace("(", '') for s in listHydrophobeLigandSort]
        listHydrophobeLigandSort = [s.replace(")", '') for s in listHydrophobeLigandSort]

        strHydrophobeProtein = ''.join(str(e) for e in hydrophobeProtein)
        strHydrophobeProtein = strHydrophobeProtein.replace(")", ',')
        listHydrophobeProteinSort = strHydrophobeProtein.split(",")
        listHydrophobeProteinSort = [s.replace("'", '') for s in listHydrophobeProteinSort]
        listHydrophobeProteinSort = [s.replace("(", '') for s in listHydrophobeProteinSort]
        listHydrophobeProteinSort = [s.replace(")", '') for s in listHydrophobeProteinSort]

        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistancehydrohydro = 0
        countidlig = 0
        countidpro = 0
        listHydrophobeLigandSortDup = []
        listHydrophobeProteinSortDup = []

        while countidlig + 3 <= len(listHydrophobeLigandSort):
            if listHydrophobeLigandSortDup.__contains__(listHydrophobeLigandSort[countidlig]) == False:
                listHydrophobeLigandSortDup.append(listHydrophobeLigandSort[countidlig])
                listHydrophobeLigandSortDup.append(listHydrophobeLigandSort[countidlig + 1])
                listHydrophobeLigandSortDup.append(listHydrophobeLigandSort[countidlig + 2])
                listHydrophobeLigandSortDup.append(listHydrophobeLigandSort[countidlig + 3])
            countidlig += 4
        while countidpro + 3 <= len(listHydrophobeProteinSort):
            if listHydrophobeProteinSortDup.__contains__(listHydrophobeProteinSort[countidpro]) == False:
                listHydrophobeProteinSortDup.append(listHydrophobeProteinSort[countidpro])
                listHydrophobeProteinSortDup.append(listHydrophobeProteinSort[countidpro + 1])
                listHydrophobeProteinSortDup.append(listHydrophobeProteinSort[countidpro + 2])
                listHydrophobeProteinSortDup.append(listHydrophobeProteinSort[countidpro + 3])
            countidpro += 4
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countlengthli = len(listHydrophobeLigandSortDup) - 4
        countlengthpro = len(listHydrophobeProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listHydrophobeLigandSortDup[countinterl] == listHydrophobeProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listHydrophobeLigandSortDup[countlx]), float(listHydrophobeLigandSortDup[countly]),
                         float(listHydrophobeLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listHydrophobeProteinSortDup[countpx]), float(listHydrophobeProteinSortDup[countpy]),
                         float(listHydrophobeProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancehydrohydro += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between Hydrophobe atoms on the ligand and Hydrophobe atoms on the protein")
        print(countdistancehydrohydro / countnumberinter)
        print("Number of interactions between Hydrophobe atoms on the ligand and Hydrophobe atoms on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        ## Interaction between LumpedHydrophobic atoms.
        hydrophobeLumpedLigand.sort()
        hydrophobeLumpedProtein.sort()
        strHydrophobeLumpedLigand = ''.join(str(e) for e in hydrophobeLumpedLigand)
        strHydrophobeLumpedLigand = strHydrophobeLumpedLigand.replace(")", ',')
        listHydrophobeLumpedLigandSort = strHydrophobeLumpedLigand.split(",")
        listHydrophobeLumpedLigandSort = [s.replace("'", '') for s in listHydrophobeLumpedLigandSort]
        listHydrophobeLumpedLigandSort = [s.replace("(", '') for s in listHydrophobeLumpedLigandSort]
        listHydrophobeLumpedLigandSort = [s.replace(")", '') for s in listHydrophobeLumpedLigandSort]

        strHydrophobeLumpedProtein = ''.join(str(e) for e in hydrophobeLumpedProtein)
        strHydrophobeLumpedProtein = strHydrophobeLumpedProtein.replace(")", ',')
        listHydrophobeLumpedProteinSort = strHydrophobeLumpedProtein.split(",")
        listHydrophobeLumpedProteinSort = [s.replace("'", '') for s in listHydrophobeLumpedProteinSort]
        listHydrophobeLumpedProteinSort = [s.replace("(", '') for s in listHydrophobeLumpedProteinSort]
        listHydrophobeLumpedProteinSort = [s.replace(")", '') for s in listHydrophobeLumpedProteinSort]

        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistancehydroluhydrolu = 0
        countidlig = 0
        countidpro = 0
        listHydrophobeLumpedLigandSortDup = []
        listHydrophobeLumpedProteinSortDup = []

        while countidlig + 3 <= len(listHydrophobeLumpedLigandSort):
            if listHydrophobeLumpedLigandSortDup.__contains__(listHydrophobeLumpedLigandSort[countidlig]) == False:
                listHydrophobeLumpedLigandSortDup.append(listHydrophobeLumpedLigandSort[countidlig])
                listHydrophobeLumpedLigandSortDup.append(listHydrophobeLumpedLigandSort[countidlig + 1])
                listHydrophobeLumpedLigandSortDup.append(listHydrophobeLumpedLigandSort[countidlig + 2])
                listHydrophobeLumpedLigandSortDup.append(listHydrophobeLumpedLigandSort[countidlig + 3])
            countidlig += 4
        while countidpro + 3 <= len(listHydrophobeLumpedProteinSort):
            if listHydrophobeLumpedProteinSortDup.__contains__(listHydrophobeLumpedProteinSort[countidpro]) == False:
                listHydrophobeLumpedProteinSortDup.append(listHydrophobeLumpedProteinSort[countidpro])
                listHydrophobeLumpedProteinSortDup.append(listHydrophobeLumpedProteinSort[countidpro + 1])
                listHydrophobeLumpedProteinSortDup.append(listHydrophobeLumpedProteinSort[countidpro + 2])
                listHydrophobeLumpedProteinSortDup.append(listHydrophobeLumpedProteinSort[countidpro + 3])
            countidpro += 4
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countlengthli = len(listHydrophobeLumpedLigandSortDup) - 4
        countlengthpro = len(listHydrophobeLumpedProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listHydrophobeLumpedLigandSortDup[countinterl] == listHydrophobeLumpedProteinSortDup[countinterp]:
                    a = numpy.array((float(listHydrophobeLumpedLigandSortDup[countlx]),
                                     float(listHydrophobeLumpedLigandSortDup[countly]),
                                     float(listHydrophobeLumpedLigandSortDup[countlz])))
                    b = numpy.array((float(listHydrophobeLumpedProteinSortDup[countpx]),
                                     float(listHydrophobeLumpedProteinSortDup[countpy]),
                                     float(listHydrophobeLumpedProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancehydroluhydrolu += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between LumpedHydrophobe atoms on the ligand and LumpedHydrophobe atoms on the protein")
        print(countdistancehydroluhydrolu / countnumberinter)
        print(
            "Number of interactions between LumpedHydrophobe atoms on the ligand and LumpedHydrophobe atoms on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        ## Interaction between a hydrophobic atom on ligand with LumpedHydrophobic atoms on the protein.
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countdistancehydrohydrolu = 0
        countnumberinter = 0
        countlengthli = len(listHydrophobeLigandSortDup) - 4
        countlengthpro = len(listHydrophobeLumpedProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listHydrophobeLigandSortDup[countinterl] == listHydrophobeLumpedProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listHydrophobeLigandSortDup[countlx]), float(listHydrophobeLigandSortDup[countly]),
                         float(listHydrophobeLigandSortDup[countlz])))
                    b = numpy.array((float(listHydrophobeLumpedProteinSortDup[countpx]),
                                     float(listHydrophobeLumpedProteinSortDup[countpy]),
                                     float(listHydrophobeLumpedProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancehydrohydrolu += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between a Hydrophobe atom on the ligand and LumpedHydrophobe atoms on the protein")
        print(countdistancehydrohydrolu / countnumberinter)
        print(
            "Number of interactions between a Hydrophobe atom on the ligand and LumpedHydrophobe atoms on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        # Interactions between LumpedHydrophobic atoms on the ligand and a Hydrophobe atom on the protein.
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countdistancehydroluhydro = 0
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countlengthli = len(listHydrophobeLumpedLigandSortDup) - 4
        countlengthpro = len(listHydrophobeProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listHydrophobeLumpedLigandSortDup[countinterl] == listHydrophobeProteinSortDup[countinterp]:
                    a = numpy.array((float(listHydrophobeLumpedLigandSortDup[countlx]),
                                     float(listHydrophobeLumpedLigandSortDup[countly]),
                                     float(listHydrophobeLumpedLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listHydrophobeProteinSortDup[countpx]), float(listHydrophobeProteinSortDup[countpy]),
                         float(listHydrophobeProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancehydroluhydro += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between LumpedHydrophobe atoms on the ligand and a Hydrophobe atom on the protein")
        print(countdistancehydroluhydro / countnumberinter)
        print(
            "Number of interactions between LumpedHydrophobe atoms on the ligand and a Hydrophobe atom on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)


    def Aromatic():
        # Interactions between aromatic rings parts:
        aromaticLigand.sort()
        aromaticProtein.sort()
        strAromaticLigand = ''.join(str(e) for e in aromaticLigand)
        strAromaticLigand = strAromaticLigand.replace(")", ',')
        listAromaticLigandSort = strAromaticLigand.split(",")
        listAromaticLigandSort = [s.replace("'", '') for s in listAromaticLigandSort]
        listAromaticLigandSort = [s.replace("(", '') for s in listAromaticLigandSort]
        listAromaticLigandSort = [s.replace(")", '') for s in listAromaticLigandSort]

        strAromaticProtein = ''.join(str(e) for e in aromaticProtein)
        strAromaticProtein = strAromaticProtein.replace(")", ',')
        listAromaticProteinSort = strAromaticProtein.split(",")
        listAromaticProteinSort = [s.replace("'", '') for s in listAromaticProteinSort]
        listAromaticProteinSort = [s.replace("(", '') for s in listAromaticProteinSort]
        listAromaticProteinSort = [s.replace(")", '') for s in listAromaticProteinSort]

        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistancearoaro = 0
        listcountaromatic = []
        countidlig = 0
        countidpro = 0
        while countidlig + 3 <= len(listAromaticLigandSort):
            if listAromaticLigandSortDup.__contains__(listAromaticLigandSort[countidlig]) == False:
                listAromaticLigandSortDup.append(listAromaticLigandSort[countidlig])
                listAromaticLigandSortDup.append(listAromaticLigandSort[countidlig + 1])
                listAromaticLigandSortDup.append(listAromaticLigandSort[countidlig + 2])
                listAromaticLigandSortDup.append(listAromaticLigandSort[countidlig + 3])
            countidlig += 4
        while countidpro + 3 <= len(listAromaticProteinSort):
            if listAromaticProteinSortDup.__contains__(listAromaticProteinSort[countidpro]) == False:
                listAromaticProteinSortDup.append(listAromaticProteinSort[countidpro])
                listAromaticProteinSortDup.append(listAromaticProteinSort[countidpro + 1])
                listAromaticProteinSortDup.append(listAromaticProteinSort[countidpro + 2])
                listAromaticProteinSortDup.append(listAromaticProteinSort[countidpro + 3])
            countidpro += 4
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countlengthli = len(listAromaticLigandSortDup) - 4
        countlengthpro = len(listAromaticProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listAromaticLigandSortDup[countinterl] == listAromaticProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listAromaticLigandSortDup[countlx]), float(listAromaticLigandSortDup[countly]),
                         float(listAromaticLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listAromaticProteinSortDup[countpx]), float(listAromaticProteinSortDup[countpy]),
                         float(listAromaticProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countdistancearoaro += dist
                    countnumberinter += 1
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between Aromatic rings on the ligand and Aromatic rings on the protein:")
        print(countdistancearoaro / countnumberinter)
        print("Number of interactions between Aromatic rings on the ligand and Aromatic rings on the protein:")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

    def CationAnion():
        # # Interactions between cation and anion #######################################################
        posionizableLigand.sort()
        negionizableProtein.sort()
        negionizableLigand.sort()
        posionizableProtein.sort()
        strPosionizableLigand = ''.join(str(e) for e in posionizableLigand)
        strPosionizableLigand = strPosionizableLigand.replace(")", ',')
        listPosionizableLigandSort = strPosionizableLigand.split(",")
        listPosionizableLigandSort = [s.replace("'", '') for s in listPosionizableLigandSort]
        listPosionizableLigandSort = [s.replace("(", '') for s in listPosionizableLigandSort]
        listPosionizableLigandSort = [s.replace(")", '') for s in listPosionizableLigandSort]

        strNegionizableProtein = ''.join(str(e) for e in negionizableProtein)
        strNegionizableProtein = strNegionizableProtein.replace(")", ',')
        listNegionizableProteinSort = strNegionizableProtein.split(",")
        listNegionizableProteinSort = [s.replace("'", '') for s in listNegionizableProteinSort]
        listNegionizableProteinSort = [s.replace("(", '') for s in listNegionizableProteinSort]
        listNegionizableProteinSort = [s.replace(")", '') for s in listNegionizableProteinSort]

        strNegionizableLigand = ''.join(str(e) for e in negionizableLigand)
        strNegionizableLigand = strNegionizableLigand.replace(")", ',')
        listNegionizableLigandSort = strNegionizableLigand.split(",")
        listNegionizableLigandSort = [s.replace("'", '') for s in listNegionizableLigandSort]
        listNegionizableLigandSort = [s.replace("(", '') for s in listNegionizableLigandSort]
        listNegionizableLigandSort = [s.replace(")", '') for s in listNegionizableLigandSort]

        strPosionizableProtein = ''.join(str(e) for e in posionizableProtein)
        strPosionizableProtein = strPosionizableProtein.replace(")", ',')
        listPosionizableProteinSort = strPosionizableProtein.split(",")
        listPosionizableProteinSort = [s.replace("'", '') for s in listPosionizableProteinSort]
        listPosionizableProteinSort = [s.replace("(", '') for s in listPosionizableProteinSort]
        listPosionizableProteinSort = [s.replace(")", '') for s in listPosionizableProteinSort]
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistanceposlnegp = 0
        listNegionizableProteinSortDup = []
        countidlig = 0
        countidpro = 0
        while countidlig + 3 <= len(listPosionizableLigandSort):
            if listPosionizableLigandSortDup.__contains__(listPosionizableLigandSort[countidlig]) == False:
                listPosionizableLigandSortDup.append(listPosionizableLigandSort[countidlig])
                listPosionizableLigandSortDup.append(listPosionizableLigandSort[countidlig + 1])
                listPosionizableLigandSortDup.append(listPosionizableLigandSort[countidlig + 2])
                listPosionizableLigandSortDup.append(listPosionizableLigandSort[countidlig + 3])
            countidlig += 4
        while countidpro + 3 <= len(listNegionizableProteinSort):
            if listNegionizableProteinSortDup.__contains__(listNegionizableProteinSort[countidpro]) == False:
                listNegionizableProteinSortDup.append(listNegionizableProteinSort[countidpro])
                listNegionizableProteinSortDup.append(listNegionizableProteinSort[countidpro + 1])
                listNegionizableProteinSortDup.append(listNegionizableProteinSort[countidpro + 2])
                listNegionizableProteinSortDup.append(listNegionizableProteinSort[countidpro + 3])
            countidpro += 4

        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countlengthli = len(listPosionizableLigandSortDup) - 4
        countlengthpro = len(listNegionizableProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listPosionizableLigandSortDup[countinterl] == listNegionizableProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listPosionizableLigandSortDup[countlx]), float(listPosionizableLigandSortDup[countly]),
                         float(listPosionizableLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listNegionizableProteinSortDup[countpx]), float(listNegionizableProteinSortDup[countpy]),
                         float(listNegionizableProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countdistanceposlnegp += dist
                    countnumberinter += 1
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between Posionizable atoms on the ligand and Negionizable atoms on the protein")
        print(countdistanceposlnegp / countnumberinter)
        print("Number of interactions between Posionizable atoms on the ligand and Negionizable atoms on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        ###### anion on ligand with cation on protein
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistanceneglposp = 0
        listNegionizableLigandSortDup = []
        countidlig = 0
        countidpro = 0

        while countidlig + 3 <= len(listNegionizableLigandSort):
            if listNegionizableLigandSortDup.__contains__(listNegionizableLigandSort[countidlig]) == False:
                listNegionizableLigandSortDup.append(listNegionizableLigandSort[countidlig])
                listNegionizableLigandSortDup.append(listNegionizableLigandSort[countidlig + 1])
                listNegionizableLigandSortDup.append(listNegionizableLigandSort[countidlig + 2])
                listNegionizableLigandSortDup.append(listNegionizableLigandSort[countidlig + 3])
            countidlig += 4
        while countidpro + 3 <= len(listPosionizableProteinSort):
            if listPosionizableProteinSortDup.__contains__(listPosionizableProteinSort[countidpro]) == False:
                listPosionizableProteinSortDup.append(listPosionizableProteinSort[countidpro])
                listPosionizableProteinSortDup.append(listPosionizableProteinSort[countidpro + 1])
                listPosionizableProteinSortDup.append(listPosionizableProteinSort[countidpro + 2])
                listPosionizableProteinSortDup.append(listPosionizableProteinSort[countidpro + 3])
            countidpro += 4

        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countlengthlig = len(listNegionizableLigandSortDup) - 4
        countlengthpro = len(listPosionizableProteinSortDup) - 4
        while countinterl != countlengthlig:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listNegionizableLigandSortDup[countinterl] == listPosionizableProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listNegionizableLigandSortDup[countlx]), float(listNegionizableLigandSortDup[countly]),
                         float(listNegionizableLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listPosionizableProteinSortDup[countpx]), float(listPosionizableProteinSortDup[countpy]),
                         float(listPosionizableProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countdistanceneglposp += dist
                    countnumberinter += 1
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between Negionizable atoms on the ligand and Posionizable atoms on the protein")
        print(countdistanceneglposp / countnumberinter)
        print("Number of interactions between Negionizable atoms on the ligand and Posionizable atoms on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)


    def CationAromatic():
        # Interactions between cation on ligand and aromatic rings on the protein.
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countdistancecatlarop = 0
        countinterl = 0
        countinterp = 0
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countlengthli = len(listPosionizableLigandSortDup) - 4
        countlengthpro = len(listAromaticProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listPosionizableLigandSortDup[countinterl] == listAromaticProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listPosionizableLigandSortDup[countlx]), float(listPosionizableLigandSortDup[countly]),
                         float(listPosionizableLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listAromaticProteinSortDup[countpx]), float(listAromaticProteinSortDup[countpy]),
                         float(listAromaticProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancecatlarop += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between cation atoms on the ligand and aromatic rings on the protein")
        print(countdistancecatlarop / countnumberinter)
        print("Number of interactions between cation atoms on the ligand and a aromatic rings on the protein")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

        # Interactions between aromatic rings on the ligand and cation on the protein.
        count0 = 0
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8 = 0
        count9 = 0
        countdistancearolcatp = 0
        countinterl = 0
        countinterp = 0
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countlengthli = len(listAromaticLigandSortDup) - 4
        countlengthpro = len(listPosionizableProteinSortDup) - 4
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            while countinterp != countlengthpro:
                countinterp += 4
                countpx += 4
                countpy += 4
                countpz += 4
                if listAromaticLigandSortDup[countinterl] == listPosionizableProteinSortDup[countinterp]:
                    a = numpy.array(
                        (float(listAromaticLigandSortDup[countlx]), float(listAromaticLigandSortDup[countly]),
                         float(listAromaticLigandSortDup[countlz])))
                    b = numpy.array(
                        (float(listPosionizableProteinSortDup[countpx]), float(listPosionizableProteinSortDup[countpy]),
                         float(listPosionizableProteinSortDup[countpz])))
                    dist = numpy.linalg.norm(a - b)
                    countnumberinter += 1
                    countdistancearolcatp += dist
                    if 0 <= dist <= 0.5:
                        count0 += 1
                    if 0.5 <= dist <= 1:
                        count1 += 1
                    if 1 <= dist <= 1.5:
                        count2 += 1
                    if 1.5 <= dist <= 2:
                        count3 += 1
                    if 2 <= dist <= 2.5:
                        count4 += 1
                    if 2.5 <= dist <= 3:
                        count5 += 1
                    if 3 <= dist <= 3.5:
                        count6 += 1
                    if 3.5 <= dist <= 4:
                        count7 += 1
                    if 4 <= dist <= 4.5:
                        count8 += 1
                    if 4.5 <= dist <= 5:
                        count9 += 1
            countinterp = -4
            countpx = -3
            countpy = -2
            countpz = -1
        print("Average distance between aromatic rings on the ligand and cation atoms on the protein")
        print(countdistancearolcatp / countnumberinter)
        print("Number of interactions between aromatic rings on the ligand and cation atoms on the protein ")
        print(countnumberinter)
        print(count0)
        print(count1)
        print(count2)
        print(count3)
        print(count4)
        print(count5)
        print(count6)
        print(count7)
        print(count8)
        print(count9)

    connection.commit()

    def main():
        Statements()
        Proteinwithinteraction()
        Calcaverageinteraction()
        DonorAcceptor()
        Hydrophobic()
        Aromatic()
        CationAnion()
        CationAromatic()

    main()
except (Exception, psycopg2.Error) as error:
    print("Error while connecting to PostgreSQL", error)
finally:
    # closing database connection.
    if (connection):
        cursor.close()
        connection.close()
        print("PostgreSQL connection is closed")
