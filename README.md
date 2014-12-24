Lisp-Land
=========

Artificial life simulation in Lisp

###Project Goals:
* Create a realistic environment that life can interact with, and be influenced by.
* Possibility for many singlecellular organisms, as well as a range of viruses.
* Realistic genetics system, support for mutation, inheritance, and evolution.

Idea/Method
===========

###Standard Cell Class:
  Every cell in the simulation is represented as a property list. Each list contains the following elements in order: POS (Position, represented as an x and y coordinate consed together ) ATP (Adenosine Triphosphate, regular integer) NA (Nucleic Acids, regular integer) AA (Amino Acids, regular integer) FA (Fatty Acids, regular integer) G (Glucose, regular integer) DNA (Deoxyribonucleic acid, represented as a list of genes, each one containing a sequence of binary numbers)

###Digital DNA Structure
  In order to reach goal number three, each cell will have a unique Digital DNA (dDNA) sequence. The dDNA will be made up of many separate genes, each encoding for a specific trait. Each gene will be a number 0-31 represented in binary, 
ex. ((0 1 0 1 1) (0 0 0 1 1) (0 1 1 1 0)) Translating to (11 3 14) Not all 31 cases will be used in every gene. In the rare case one trait requires than 31 values, two genes will be used. If data needs to become a percentage, values 1-20 are used and multiplied by 5. See “Genes Translated ” for more in depth info.

What "Genes" code for what?
===========================

This is one of the most important features in Lisp-Life, the dynamic genetics system.

#####Gene 1 - Primary Membrane Composition:
Determines permeability of the membrane (Integer 1-10) Slots 30-31 unused.
* Phospholipids (purely saturated) [0-9] = Perm:5
* Phospholipids (mixed) [10-19] = Perm:4
* Free fatty acids [20-29] = Perm:3

#####Gene 2 - Mechanisms For Photosynthesis:
Determines the presence of chlorophyll and packaged form (Integer 1-3) Slots 30-31 unused.
* No Chlorophyll [0-9] = Chlorop:1
* Free Chlorophyll [10-19] = Chlorop:2
* Chloroplasts [20-29] = Chlorop:3
