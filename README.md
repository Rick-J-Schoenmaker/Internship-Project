# The DrugLi Project

Python code to create database of protein-ligand
interactions with structural data derived from the PDBBind database.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites
```
Retrieve the structural data from PDBBind: http://www.pdbbind.org.cn/\
```
```
Python (3.7)
```
```
Conda(4.7.12)
```
```
Numpy(1.16.4)
```
```
RDKit (Q3 2019)
```



## Deployment
To create the database it is neccesary to run the scripts in the following steps:
1. Run the SQL script Drugli.sql to create the tables for the database.
2. Run the DruglikeFiltering.py script to create a druglike dataset from the structural data obtained from the PDBBind database.
3. Run the DatabasePopulator.py script to fill the tables in the database with data from the druglike dataset.
4. Run the PostgresTester.py script to retrieve histograms and other information about the interactions from the database.
5. (Optional) Run the SQL scripts from Drugli statements.sql to query the database

## Authors
* **Rick J Schoenmaker**



