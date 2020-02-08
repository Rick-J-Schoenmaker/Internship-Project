# Author: Rick Schoenmaker
# Date: 29-01-2020
# Function: This script extracts data from the PROLI database and calculates numbers such as
# the amount of interactions and averages.
import os

import numpy
import psycopg2
import numpy as np
import matplotlib.pyplot as plt

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

    listPosionizableProteinSort = []
    listPosionizableLigandSort = []

    listAromaticProteinSort = []
    listAromaticLigandSort = []


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
            "SELECT DISTINCT end_concept.interactions.inter_id, end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id WHERE group_type_id_fk =  4 ")
        hydrophobeli = cursor.fetchall()
        for i in hydrophobeli:
            hydrophobeLigand.append(i)

        cursor.execute(
            "SELECT DISTINCT end_concept.interactions.inter_id,end_concept.group_atoms.pseudo_x, pseudo_y, pseudo_z FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id WHERE group_type_id_fk =  4 ")
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
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listDonorLigandSort) - 5
        countlengthpro = len(listAcceptorProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listDonorLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listDonorLigandSort[countlx]), float(listDonorLigandSort[countly]),
                 float(listDonorLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listAcceptorProteinSort[countinterp])
            distanceProteinlist.append((float(listAcceptorProteinSort[countpx]), float(listAcceptorProteinSort[countpy]),
                 float(listAcceptorProteinSort[countpz])))
        dict_protein = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein.items():
            if dict_ligand.__contains__(kp):
                a = numpy.array(dict_protein.get(kp))
                b = numpy.array(dict_ligand.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancedalp += dist
            if 0 <= dist <= 0.25:
                count0 += 1
            if 0.25 <= dist <= 0.5:
                count1 += 1
            if 0.5 <= dist <= 0.75:
                count2 += 1
            if 0.75 <= dist <= 1:
                count3 += 1
            if 1 <= dist <= 1.25:
                count4 += 1
            if 1.25 <= dist <= 1.5:
                count5 += 1
            if 1.5 <= dist <= 1.75:
                count6 += 1
            if 1.75 <= dist <= 2:
                count7 += 1
            if 2 <= dist <= 2.25:
                count8 += 1
            if 2.25 <= dist <= 2.5:
                count9 += 1
            if 2.5 <= dist <= 2.75:
                count10 += 1
            if 2.75 <= dist <= 3:
                count11 += 1
            if 3 <= dist <= 3.25:
                count12 += 1
            if 3.25 <= dist <= 3.5:
                count13 += 1
            if 3.5 <= dist <= 3.75:
                count14 += 1
            if 3.75 <= dist <= 4:
                count15 += 1
            if 4 <= dist <= 4.25:
                count16 += 1
            if 4.25 <= dist <= 4.5:
                count17 += 1
            if 4.5 <= dist <= 4.75:
                count18 += 1
            if 4.75 <= dist <= 5:
                count19 += 1
        print("Average distance between LigandDonor and ProteinAcceptor:")
        print(countdistancedalp / countnumberinter)
        print("Distances LigandDonor and ProteinAcceptor")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances LigandDonor - ProteinAcceptor')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()

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
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listAcceptorLigandSort) - 5
        countlengthpro = len(listDonorProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listAcceptorLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listAcceptorLigandSort[countlx]), float(listAcceptorLigandSort[countly]),
                 float(listAcceptorLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listDonorProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listDonorProteinSort[countpx]), float(listDonorProteinSort[countpy]),
                 float(listDonorProteinSort[countpz])))
        dict_protein = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein.items():
            if dict_ligand.__contains__(kp):
                a = numpy.array(dict_protein.get(kp))
                b = numpy.array(dict_ligand.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistanceadlp += dist
            if 0 <= dist <= 0.25:
                count0 += 1
            if 0.25 <= dist <= 0.5:
                count1 += 1
            if 0.5 <= dist <= 0.75:
                count2 += 1
            if 0.75 <= dist <= 1:
                count3 += 1
            if 1 <= dist <= 1.25:
                count4 += 1
            if 1.25 <= dist <= 1.5:
                count5 += 1
            if 1.5 <= dist <= 1.75:
                count6 += 1
            if 1.75 <= dist <= 2:
                count7 += 1
            if 2 <= dist <= 2.25:
                count8 += 1
            if 2.25 <= dist <= 2.5:
                count9 += 1
            if 2.5 <= dist <= 2.75:
                count10 += 1
            if 2.75 <= dist <= 3:
                count11 += 1
            if 3 <= dist <= 3.25:
                count12 += 1
            if 3.25 <= dist <= 3.5:
                count13 += 1
            if 3.5 <= dist <= 3.75:
                count14 += 1
            if 3.75 <= dist <= 4:
                count15 += 1
            if 4 <= dist <= 4.25:
                count16 += 1
            if 4.25 <= dist <= 4.5:
                count17 += 1
            if 4.5 <= dist <= 4.75:
                count18 += 1
            if 4.75 <= dist <= 5:
                count19 += 1
        print("Average distance between LigandAcceptor and ProteinDonor:")
        print(countdistanceadlp / countnumberinter)
        print("Distances LigandAcceptor and ProteinDonor")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances LigandAcceptor - ProteinDonor')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()
        print("Average distance between LigandAcceptor and ProteinDonor:")
        print(countdistanceadlp / len(donorProtein))
        print("Average distance between all Acceptor and Donor interactions:")
        countdis = countdistanceadlp / len(donorProtein) + countdistancedalp / len(donorLigand)
        print(countdis / 2)

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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listHydrophobeLigandSort) -5
        countlengthpro = len(listHydrophobeProteinSort) -5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listHydrophobeLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listHydrophobeLigandSort[countlx]), float(listHydrophobeLigandSort[countly]),
                 float(listHydrophobeLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listHydrophobeProteinSort[countinterp])
            distanceProteinlist.append((float(listHydrophobeProteinSort[countpx]), float(listHydrophobeProteinSort[countpy]),
                         float(listHydrophobeProteinSort[countpz])))
        dict_protein = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein.items():
            if dict_ligand.__contains__(kp):
                a = numpy.array(dict_protein.get(kp))
                b = numpy.array(dict_ligand.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancehydrohydro += dist
            if 0 <= dist <= 0.25:
                count0 += 1
            if 0.25 <= dist <= 0.5:
                count1 += 1
            if 0.5 <= dist <= 0.75:
                count2 += 1
            if 0.75 <= dist <= 1:
                count3 += 1
            if 1 <= dist <= 1.25:
                count4 += 1
            if 1.25 <= dist <= 1.5:
                count5 += 1
            if 1.5 <= dist <= 1.75:
                count6 += 1
            if 1.75 <= dist <= 2:
                count7 += 1
            if 2 <= dist <= 2.25:
                count8 += 1
            if 2.25 <= dist <= 2.5:
                count9 += 1
            if 2.5 <= dist <= 2.75:
                count10 += 1
            if 2.75 <= dist <= 3:
                count11 += 1
            if 3 <= dist <= 3.25:
                count12 += 1
            if 3.25 <= dist <= 3.5:
                count13 += 1
            if 3.5 <= dist <= 3.75:
                count14 += 1
            if 3.75 <= dist <= 4:
                count15 += 1
            if 4 <= dist <= 4.25:
                count16 += 1
            if 4.25 <= dist <= 4.5:
                count17 += 1
            if 4.5 <= dist <= 4.75:
                count18 += 1
            if 4.75 <= dist <= 5:
                count19 += 1
        print("Average distance between Hydrophobe atoms on the ligand and Hydrophobe atoms on the protein")
        print(countdistancehydrohydro / countnumberinter)
        print("Number of interactions between Hydrophobe atoms on the ligand and Hydrophobe atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Hydrophobic Atoms - Hydrophobic Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()

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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listHydrophobeLumpedLigandSort) - 5
        countlengthpro = len(listHydrophobeLumpedProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listHydrophobeLumpedLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listHydrophobeLumpedLigandSort[countlx]), float(listHydrophobeLumpedLigandSort[countly]),
                 float(listHydrophobeLumpedLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listHydrophobeLumpedProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listHydrophobeLumpedProteinSort[countpx]), float(listHydrophobeLumpedProteinSort[countpy]),
                 float(listHydrophobeLumpedProteinSort[countpz])))

        dict_protein_hydro_lumped = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_hydro_lumped = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_hydro_lumped.items():
            if dict_ligand_hydro_lumped.__contains__(kp):
                a = numpy.array(dict_protein_hydro_lumped.get(kp))
                b = numpy.array(dict_ligand_hydro_lumped.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancehydroluhydrolu += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
        print("Average distance between LumpedHydrophobe atoms on the ligand and LumpedHydrophobe atoms on the protein")
        print(countdistancehydroluhydrolu / countnumberinter)
        print("Number of interactions between LumpedHydrophobe atoms on the ligand and LumpedHydrophobe atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances LumpedHydrophobic Atoms - LumpedHydrophobic Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()


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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countdistancehydrohydrolu = 0
        countnumberinter = 0
        countlengthli = len(listHydrophobeLigandSort) - 5
        countlengthpro = len(listHydrophobeLumpedProteinSort) - 5
        dist = 0
        for kp, vp in dict_protein_hydro_lumped.items():
            if dict_ligand.__contains__(kp):
                a = numpy.array(dict_protein_hydro_lumped.get(kp))
                b = numpy.array(dict_ligand.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancehydrohydrolu += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1

        print("Average distance between a Hydrophobe atom on the ligand and LumpedHydrophobe atoms on the protein")
        print(countdistancehydrohydrolu / countnumberinter)
        print("Number of interactions between a Hydrophobe atom on the ligand and LumpedHydrophobe atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Hydrophobic Atoms - LumpedHydrophobic Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()

        # # Interactions between LumpedHydrophobic atoms on the ligand and a Hydrophobe atom on the protein.
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countdistancehydroluhydro = 0
        countnumberinter = 0
        countlengthli = len(listHydrophobeLumpedLigandSort) - 5
        countlengthpro = len(listHydrophobeProteinSort) - 5
        dist = 0
        for kp, vp in dict_protein.items():
            if dict_ligand_hydro_lumped.__contains__(kp):
                a = numpy.array(dict_protein.get(kp))
                b = numpy.array(dict_ligand_hydro_lumped.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancehydroluhydro += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1

        print("Average distance between a LumpedHydrophobe atom on the ligand and Hydrophobe atoms on the protein")
        print(countdistancehydroluhydro / countnumberinter)
        print("Number of interactions between a LumpedHydrophobe atom on the ligand and Hydrophobe atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances LumpedHydrophobic Atoms - Hydrophobic Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()



    def Aromatic_CationAnion():
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        count20 = 0
        count21 = 0
        countlengthli = len(listAromaticLigandSort) - 5
        countlengthpro = len(listAromaticProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listAromaticLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listAromaticLigandSort[countlx]), float(listAromaticLigandSort[countly]),
                 float(listAromaticLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listAromaticProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listAromaticProteinSort[countpx]), float(listAromaticProteinSort[countpy]),
                 float(listAromaticProteinSort[countpz])))

        dict_protein_aro = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_aro = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_aro.items():
            if dict_ligand_aro.__contains__(kp):
                a = numpy.array(dict_protein_aro.get(kp))
                b = numpy.array(dict_ligand_aro.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancearoaro += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
                if 5 <= dist <= 5.25:
                    count20 += 1
                if 5.25 <= dist <= 5.5:
                    count21 += 1
        print("Average distance between Aromatic atoms on the ligand and Aromatic atoms on the protein")
        print(countdistancearoaro / countnumberinter)
        print("Number of interactions between Aromatic atoms on the ligand and Aromatic atoms on the protein")
        print(countnumberinter)
        x = np.arange(22)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19, count20, count21]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Aromatic Pseudo-Atoms - Aromatic Pseudo-Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5', '5 | 5.25',
                       '5.25 | 5.5'], rotation=90)

        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)

        plt.tight_layout()
        plt.show()

        # Interactions between cation and anion #######################################################
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listPosionizableLigandSort) - 5
        countlengthpro = len(listNegionizableProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listPosionizableLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listPosionizableLigandSort[countlx]), float(listPosionizableLigandSort[countly]),
                 float(listPosionizableLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listNegionizableProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listNegionizableProteinSort[countpx]), float(listNegionizableProteinSort[countpy]),
                 float(listNegionizableProteinSort[countpz])))

        dict_protein_neg = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_pos = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_neg.items():
            if dict_ligand_pos.__contains__(kp):
                a = numpy.array(dict_protein_neg.get(kp))
                b = numpy.array(dict_ligand_pos.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistanceposlnegp += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
        print("Average distance between Posionizable atoms on the ligand and Negionizable atoms on the protein")
        print(countdistanceposlnegp / countnumberinter)
        print("Number of interactions between Posionizable atoms on the ligand and Negionizable atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Cation Atoms - Anion Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()


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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listNegionizableLigandSort) - 5
        countlengthpro = len(listPosionizableProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listNegionizableLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listNegionizableLigandSort[countlx]), float(listNegionizableLigandSort[countly]),
                 float(listNegionizableLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listPosionizableProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listPosionizableProteinSort[countpx]), float(listPosionizableProteinSort[countpy]),
                 float(listPosionizableProteinSort[countpz])))

        dict_protein_pos = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_neg = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_pos.items():
            if dict_ligand_neg.__contains__(kp):
                a = numpy.array(dict_protein_pos.get(kp))
                b = numpy.array(dict_ligand_neg.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistanceneglposp += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
        print("Average distance between Negionizable atoms on the ligand and Posionizable atoms on the protein")
        print(countdistanceneglposp / countnumberinter)
        print("Number of interactions between Negionizable atoms on the ligand and Posionizable atoms on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Anion Atoms - Cation Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()



        # Interactions between cation on ligand and aromatic rings on the protein.
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistancecatlarop = 0
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        print(len(listPosionizableLigandSort))
        print(len(listAromaticProteinSort))
        countlengthli = len(listPosionizableLigandSort) - 5
        countlengthpro = len(listAromaticProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listPosionizableLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listPosionizableLigandSort[countlx]), float(listPosionizableLigandSort[countly]),
                 float(listPosionizableLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listAromaticProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listAromaticProteinSort[countpx]), float(listAromaticProteinSort[countpy]),
                 float(listAromaticProteinSort[countpz])))

        dict_protein_aro = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_pos = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_aro.items():
            if dict_ligand_pos.__contains__(kp):
                a = numpy.array(dict_protein_aro.get(kp))
                b = numpy.array(dict_ligand_pos.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancecatlarop += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
        print("Average distance between cation atoms on the ligand and aromatic rings on the protein")
        print(countdistancecatlarop / countnumberinter)
        print("Number of interactions between cation atoms on the ligand and a aromatic rings on the protein")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Cation Atoms - Aromatic Pseudo-Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()

        # Interactions between aromatic rings on the ligand and cation on the protein.
        countinterl = -4
        countinterp = -4
        countlx = -3
        countly = -2
        countlz = -1
        countpx = -3
        countpy = -2
        countpz = -1
        countnumberinter = 0
        countdistancearolcatp = 0
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
        count10 = 0
        count11 = 0
        count12 = 0
        count13 = 0
        count14 = 0
        count15 = 0
        count16 = 0
        count17 = 0
        count18 = 0
        count19 = 0
        countlengthli = len(listAromaticLigandSort) - 5
        countlengthpro = len(listPosionizableProteinSort) - 5
        interLigandList = []
        distanceLigandlist = []
        interProteinList = []
        distanceProteinlist = []
        while countinterl != countlengthli:
            countinterl += 4
            countlx += 4
            countly += 4
            countlz += 4
            interLigandList.append(listAromaticLigandSort[countinterl])
            distanceLigandlist.append(
                (float(listAromaticLigandSort[countlx]), float(listAromaticLigandSort[countly]),
                 float(listAromaticLigandSort[countlz])))
        while countinterp != countlengthpro:
            countinterp += 4
            countpx += 4
            countpy += 4
            countpz += 4
            interProteinList.append(listPosionizableProteinSort[countinterp])
            distanceProteinlist.append(
                (float(listPosionizableProteinSort[countpx]), float(listPosionizableProteinSort[countpy]),
                 float(listPosionizableProteinSort[countpz])))

        dict_protein_aro = {kp: vp for kp, vp in zip(interProteinList, distanceProteinlist)}
        dict_ligand_aro = {kl: vl for kl, vl in zip(interLigandList, distanceLigandlist)}
        for kp, vp in dict_protein_aro.items():
            if dict_ligand_aro.__contains__(kp):
                a = numpy.array(dict_protein_aro.get(kp))
                b = numpy.array(dict_ligand_aro.get(kp))
                dist = numpy.linalg.norm(a - b)
                countnumberinter += 1
                countdistancearolcatp += dist
                if 0 <= dist <= 0.25:
                    count0 += 1
                if 0.25 <= dist <= 0.5:
                    count1 += 1
                if 0.5 <= dist <= 0.75:
                    count2 += 1
                if 0.75 <= dist <= 1:
                    count3 += 1
                if 1 <= dist <= 1.25:
                    count4 += 1
                if 1.25 <= dist <= 1.5:
                    count5 += 1
                if 1.5 <= dist <= 1.75:
                    count6 += 1
                if 1.75 <= dist <= 2:
                    count7 += 1
                if 2 <= dist <= 2.25:
                    count8 += 1
                if 2.25 <= dist <= 2.5:
                    count9 += 1
                if 2.5 <= dist <= 2.75:
                    count10 += 1
                if 2.75 <= dist <= 3:
                    count11 += 1
                if 3 <= dist <= 3.25:
                    count12 += 1
                if 3.25 <= dist <= 3.5:
                    count13 += 1
                if 3.5 <= dist <= 3.75:
                    count14 += 1
                if 3.75 <= dist <= 4:
                    count15 += 1
                if 4 <= dist <= 4.25:
                    count16 += 1
                if 4.25 <= dist <= 4.5:
                    count17 += 1
                if 4.5 <= dist <= 4.75:
                    count18 += 1
                if 4.75 <= dist <= 5:
                    count19 += 1
        print("Average distance between aromatic rings on the ligand and cation atoms on the protein")
        print(countdistancearolcatp / countnumberinter)
        print("Number of interactions between aromatic rings on the ligand and cation atoms on the protein ")
        print(countnumberinter)
        x = np.arange(20)
        y = [count0, count1, count2, count3, count4, count5, count6, count7, count8, count9, count10, count11, count12,
             count13, count14, count15, count16, count17, count18, count19]

        plt.figure()
        # plt.plot(x, y)
        plt.bar(x, y, alpha=0.2, color=('blue'), edgecolor='cyan')
        plt.xlabel("Distance in (Å)")
        plt.ylabel("Count of interactions")
        plt.title('Distances Aromatic Pseudo-Atoms - Cation Atoms ')
        plt.xticks(x, ['0 | 0.25', '0.25 | 0.5', '0.5 | 0.75', '0.75 | 1', '1 | 1.25', '1.25 | 1.5', '1.5 | 1.75',
                       '1.75 | 2', '2 | 2.25', '2.25 | 2.5', '2.5 | 2.25', '2.25 | 2.5', '2.5 | 3', '3.25 | 3.5',
                       '3.5 | 3.75', '3.75 | 4', '4 | 4.25', '4.25 | 4.5', '4.5 | 4.75', '4.75 | 5'], rotation=90)
        # set parameters for tick labels
        plt.tick_params(axis='x', which='major', labelsize=9)
        plt.tight_layout()
        plt.show()
    connection.commit()

    def main():
        Statements()
        DonorAcceptor()
        Hydrophobic()
        Aromatic_CationAnion()


    main()
except (Exception, psycopg2.Error) as error:
    print("Error while connecting to PostgreSQL", error)
finally:
    # closing database connection.
    if (connection):
        cursor.close()
        connection.close()
        print("PostgreSQL connection is closed")
