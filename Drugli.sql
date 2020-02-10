create table protein
(
    protein_id text,
    seq        text
);

alter table protein
    owner to postgres;

create table ligand
(
    ligand_id text,
    smiles    text
);

alter table ligand
    owner to postgres;

create table atoms
(
    atom_id integer,
    atom    text
);

alter table atoms
    owner to postgres;

create table group_type
(
    group_type_id serial,
    type          text
);

alter table group_type
    owner to postgres;

create table group_atoms
(
    group_type_id_fk integer,
    atom_id_fk       integer,
    group_id         integer,
    pseudo_x         text,
    pseudo_y         text,
    pseudo_z         text
);

alter table group_atoms
    owner to postgres;

create table interactions
(
    group_id_li  integer,
    group_id_pro integer,
    inter_id     serial not null
);

alter table interactions
    owner to postgres;

create table interacting_ligand_protein
(
    protein_id text,
    ligand_id  text
);

alter table interacting_ligand_protein
    owner to postgres;

create table pro_conform
(
    pdb_id   text,
    atom_id  integer,
    chain_id text,
    x        text,
    y        text,
    z        text
);

alter table pro_conform
    owner to postgres;

create table lig_conform
(
    pdb_id  text,
    lig_id  text,
    atom_id serial not null,
    x       text,
    y       text,
    z       text
);

alter table lig_conform
    owner to postgres;


