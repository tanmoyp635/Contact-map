# Contact-map
This Python code computes a contact map for any system size, excluding i, i±1 residue pairs. It uses a 4.5 Å cutoff between any heavy atoms of residue pairs. From an equilibrated PDB file, it generates a pseudo list of potential contacts, retaining only those present in at least 75% of the frames.
The input files are hard coded and are named as follows: I) equilibrated PDB: reference_contactmap.pdb, II) DCD file: concatenated_NP_large.dcd, iii) topology file: stripped.pic-med-non-phosphorylated-final.psf, output: conmap_test.dat and image: test.png. You can change them whenever you need them. 
You can run the code simply by running: python contact_map.py
you need to have numpy, matplotlib and mdtraj installed.
