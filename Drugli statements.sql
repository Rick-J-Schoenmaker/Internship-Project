/* Ligand interacting atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id;
/* Ligand interacting Donor atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 1;
/* Ligand interacting Acceptor atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 2;
/* Ligand interacting Aromatic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 3;
/* Ligand interacting Hydrophobic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 4;
/* Ligand interacting LumpedHydrophobic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 5;
/* Ligand interacting Posionizable atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 6;
/* Ligand interacting ZNBinder atoms (should be 0) */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 7;
/* Ligand interacting atoms NegIonizable */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 8;

/* Protein interacting atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id;
/* Protein interacting Donor atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 1;
/* Protein interacting Acceptor atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 2;
/* Protein interacting Aromatic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 3;
/* Protein interacting Hydrophobic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 4;
/* Protein interacting LumpedHydrophobic atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 5;
/* Protein interacting Posionizable atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 6;
/* Protein interacting ZNBinder atoms (should be 0) */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 7;
/* Protein interacting NegIonizable atoms */
SELECT count(distinct (atom_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id WHERE group_type_id = 8;

/* Structures with an interactions*/
SELECT count(distinct (pdb_id)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id;

/* Total number of interactions */
SELECT count(distinct (interactions)) FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id;
/* Total number of interactions per ligand structure */
SELECT count(distinct (interactions)),pdb_id,lig_id FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id GROUP BY pdb_id,lig_id;

/* Total number of interactions per protein structure */
SELECT count(distinct (interactions)),pdb_id FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id GROUP BY pdb_id;


/* Total number of interacting  protein atoms per structure */
SELECT count(distinct (atom_id)),pdb_id  FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_pro  = group_id INNER JOIN end_concept.pro_conform ON group_atoms.atom_id_fk = pro_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id GROUP BY pdb_id;

/* Total number of interacting  ligand atoms per structure */
SELECT count(distinct (atom_id)),pdb_id,lig_id  FROM end_concept.group_atoms INNER JOIN end_concept.interactions ON group_id_li  = group_id INNER JOIN end_concept.lig_conform ON group_atoms.atom_id_fk = lig_conform.atom_id INNER JOIN end_concept.group_type on group_atoms.group_type_id_fk = group_type.group_type_id GROUP BY pdb_id, lig_id;