# DNA-Search
'Searching' using mutation in DNA

Taking the analogy too far, this code simulates small numbers of base pairs (9,18)
and mutates them at a rate commensurate with bacterial cell replication (20-40 min).

Based on the mutation rate per bp, the output is analyzed to determine how long it would
take to generate complete coverage of the search space. 

The results are fairly grim in that it would take years to search even small numbers of bp
in this manner without large populations of cells/plasmids to start with. 

![9 bp 1e5 gens](https://user-images.githubusercontent.com/30641156/226087875-fceafaa6-6a28-468d-ab9b-ed57526dafcf.png)

![18 bp 1e3 gens 1e3 pop cap mu is 0 1](https://user-images.githubusercontent.com/30641156/226087881-0b533f7c-59c6-4343-8bb8-d509e1ea3afd.png)
