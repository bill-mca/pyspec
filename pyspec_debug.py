from pyteomics.mass import *
from pyteomics.parser import *
from pyspec3 import *
import random

def new_sequence():
  return ''.join([random.choice(std_amino_acids) for aa  \
    in range(random.randint(1, 500))])
    
sequences = [new_sequence() for x in range(100)]

peps = [Peptide(sequence) for sequence in sequences]

#print ([pep.mass() for pep in peps])

#print ([[pep.mass(x) for x in range(10)] for pep in peps[:10]])

def test_digest():
  proteases = [random.choice(expasy_rules.keys()) for x in range(random.randint(1, 3))]
  print digest(Peptide(new_sequence()), proteases).describe()
  
#for x in range(10):
#  print repr(digest(Peptide(new_sequence()), 'trypsin'))
  
#for x in range(10): test_digest()

x_linker = (7, 3, 'C')
c = cross_link(x_linker, 'trypsin', Peptide(new_sequence()), Peptide(new_sequence()))
c.skipped_cuts = 4
#print c.find(4000, 20)

h = cross_link(x_linker, 'trypsin', Peptide('RCCCCCRCCCCCCCKCCCCCCCCCCCCCCCKCCCCCCCR'),\
	Peptide('CCCCCCCCCCCCCCCCCCCCCCCCCRCCCCCCCCCCKCCCCCCCCCCRCCCCCCCCCCCCCCCCKCCCCCCR'))
h.skipped_cuts = 5
#print h.find(3200, 120)

p = Peptide('KFPIGALCDSSLIWMCPLFGGHSHRCCLIPFLSIHHRFIILPSACVMW')
#set([(1, 0, FIILPSACVMW, FPIGALCDSSLIWMCPLFGGHSHR)]) # 3921 + CL
a = cross_link(x_linker, 'trypsin', p)
a.skipped_cuts = 4
print a.find(3931, 4)