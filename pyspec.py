#! /usr/bin/python2

from pyteomics.parser import *
from pyteomics.mass import *
from pyteomics.auxiliary import PyteomicsError
from re import finditer, findall

"""
Pyspec a python library for analysing crosslink protein digest mass spectrometry
Copyright (C) 2020  Bill McAlister

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.txt

"""
# dev ideas:
# Make an environment object?
# it could spawn new peptides env.new()
# it could be called instead of isodistro etc. 
# with ends -NH2 -OH

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def carry_on(mass, target, increment):
	"""Carry on is a little function to help handling using either -ve or 
        +ve values for the mass of a crosslinker (for disulfides). It saves 
        defining 6 or more lambdas. only returns true if the target could be 
        reached by adding more 'increment' to mass."""
	if increment < 0 and target > mass: 
		return True

	elif increment > 0 and target < mass: 
		return True

	else: return False

class inrange:
	def __init__(self, bottom, top):
		assert all([is_number(top), is_number(bottom)])
		self.top = top
		self.bottom = bottom

	def __repr__(self):
		return "inRange: {0}-{1}".format(self.bottom, self.top)

	def __call__(value):
		return value < self.top and value > self.bottom

class environment:
        
	def __init__(isodistro, aacomp):
		self.aacomp = kwargs.get(aacomp, std_aa_comp)
		self.isodistro = kwargs.get(isodistro, nist_mass)

	def new(sequence, name=''):
		return peptide(sequence, name, 
			aacomp=self.aacomp,
			isodisto=self.isodistro)

class peptide:	
	def __init__(self, sequence, name='', **kwargs):
		"""Most of the functionality of this module. The an amino acid 
                sequence and turns it into a peptide object. The advantage of 
                using the peptide object is that it is a more pythonic 
                representation of a protein than just a naked sequence. One 
                important feature is that the peptide can be fed an "isodistro"
                or "aacomp" as kwargs at initialisation. these are dictionaries 
                that allow the user to specify the propensity of different 
                isotopes and the chemical composition of the amino acids. one 
                can even make up whole new amino acids by adding a new entry to 
                the aacomp dictionary. The modx system allows the description of
                modified / unnatural amino acids by prepending lowercase 
                characters to an amino acids one letter name eg M is methionine,
                oM could be oxidised methionine and ooM double oxidised.
		"""

		if isinstance(sequence, str):
			self.sequence = parse(sequence,
                                              show_unmodified_termini=False, 
				              split=False,
				              allow_unknown_modifications=False,
				              **kwargs)

		elif isinstance(sequence, list):
			self.sequence = sequence
			
		else: raise TypeError("This is not a valid sequence.")

		self.name = name

		self.isodistro = kwargs.get('isodistro', nist_mass)#<<<<<--
		# Will these two lines cost heaps of RAM or is a 
		# reference to the existing data created?
		# It saves many lines of "try", "except" code.
		self.aacomp = kwargs.get('aacomp', std_aa_comp)

	def __float__(self):
		return self.mass()

	def __nonzero__(self):
		return bool(self.sequence)

	def __key(self):
		return repr(self)

	def __hash__(self):
		return hash(self.__key())

	def __eq__(self, other):

		if isinstance(other, peptide):
			return self.mass() == other.mass()
		else: return self.mass() == other

	def __gt__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, peptide):
			return self.mass > other.mass
		else: return self.mass() > other

	def __lt__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, peptide):
			return self.mass < other.mass
		else: return self.mass() < other

	def __ge__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, peptide):
			return self.mass >= other.mass
		else: return self .mass()>= other

	def __le__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, peptide):
			return self.mass <= other.mass
		else: return self.mass() <= other
	
	def __cmp__(self, other):
		"""Compares peptides to each other. Those with greater mass score higher.
		This is the function called when a list is sorted. So a list of peptide
		objects will be sorted from lowest mass to highest.
		""" 
		try: assert isinstance(other, peptide)
		except AssertionError: raise TypeError("Can't" 
			"compare a peptide to a non peptide.")
		return cmp(self.mass(), other.mass())

	def __repr__(self):
		"""print(peptide) --> the peptide's sequence as a string"""
		return ''.join(self.sequence)

	def __str__(self, *args):
		return ''.join(self.sequence)
		
	def describe(self):
		"""Prints verbose details about the peptide instance onto screen. Useful if you
		want to quickly summarise what a list of peptides contains."""
		sequence = ''.join(self.sequence) 
		if len(sequence) > 25: sequence = sequence[:25] + '...'
		return "Peptide: {0} residues | {1} Da | sequence:  {2}".format(
			len(self), int(self.mass()*100)/100.00, sequence)

	def __getitem__(self, value):
		"""peptide['PEPTIDE'] --> The frequency of occurance of 'PEPTIDE' in the 
		amino acid sequence.

		peptide['trypsin'] --> a digest of peptide by trypsin.

		peptide[3:7] --> ''.join(peptide.sequence)[3:7]

		peptide[7] --> ''.join(peptide.sequence)[7]
		"""
		if isinstance(value, str): 
		 	if value in amino_acids:
				return [g.end() for g in
                                        finditer(expasy_rules[protease],
				                 self.sequence)]
			elif expasy_rules.has_key(value):
				digests.get(value, self.digest(value))
			else: raise TypeError

		elif isinstance(value, slice) and self.sequence:
			return self.sequence[value.start:value.stop:value.step]

		elif self.sequence and isinstance(value, int): 
			return self.sequence[value]
		
		else: raise TypeError(value)
	
	def __len__(self):
		"""len(peptide) --> an integer indicating the sequence length"""
		return length(self.sequence, labels=self.aacomp.keys())

	def __contains__(self, x):
		"""true if x is a substring of the peptide's sequence false otherwise"""
		x = x.upper()
		try:
			assert all([a in amino_acids for a in x])
		except AssertionError:
			print("The provided subsequence was not a protein sequence.")
		if self.find(x) != -1:
			return True
		else: return False

	def mass(self, charge=0, *args, **kwargs):
		"""
		"""
		# Need a good system for calculating charge m/z
		# preferably with Na+ support.
		isodistro = kwargs.get('isodistro', self.isodistro)

		aacomp = kwargs.get('aacomp', self.aacomp)
		amino_acids = aacomp.keys()

		for key in set(amino_acids).intersection(kwargs.keys()):
			if isinstance(kwargs[key], Composition):
				aacomp[key] = kwargs[key]
		
			elif isinstance(kwargs[key], dict):
				isodistro[key] = kwargs[key]

		for key in set(aacomp).intersection(kwargs.keys()):
			aacomp[key] = kwargs[key]

		if {'av', 'average', 'Average', 'AV'}.intersection(args):
			return calculate_mass(self.sequence, 
						average = True,
						aa_comp=aacomp,
						mass_data=isodistro)

		elif {'mp','mpm','MPM','MP','probable','most probable',
			'most probable mass'}.intersection(args):
			mp = most_probable_isotopic_composition(self.sequence,
							aa_comp=aacomp, 
							mass_data=isodistro)
			return calculate_mass(mp)

		else:  # mono-isotopic
			return calculate_mass(self.sequence, 
				aa_comp = aacomp, 
				mass_data = isodistro,)

	def count(self, x):
		"""An easy way of counting the frequency of an amino acid."""
		return self.sequence.count(x)

	def startswith(self, prefix, start=None, end=None):
		if not self.sequence: print("Unknown Sequence.")
		return self.sequence.startswith(prefix, start, end)

	def endswith(self, suffix, start=None, end=None):
		if not self.sequence: print("Unknown Sequence.")
		return self.sequence.endswith(suffix, start, end)


	def find(self, subseq, tolerance=1):
		if isinstance(subseq, string):
			return self.sequence.find(subseq)
		else:
			x = 0
			while x < len(self):
				i = x + 1

				while i < len(self):
					mass = self.mass()

					if  mass > target_mass + tolerance: break
					elif mass < target_mass - tolerance: pass
					else: return True
					i += 1
				x += 1

			return False


	def find_fragment(self, target_mass, tolerance=1):
		"""Takes the mass of an unknown fragment originating from the peptide's sequence
		and returns all feasible substrings that correspond to a fragment of the given
		mass. the data is a list of (mass,sequence) tuples sorted in terms of mass;
		lowest to highest."""

		# Used a set here to stop non-unique results clogging up the RAM
		fragments = set([])

		# x is the index of the start of a potential fragment. 
		# i is the index of the end of a potential fragment. 
		for x in range(len(self.sequence)):
			i = x + 1

			while i < len(self.sequence):
				mass = self.mass()

				if  mass > target_mass + tolerance: break
				elif mass < target_mass - tolerance: pass
				else: fragments.add(self.sequence[x:i])

				i += 1
		
		return sorted([peptide(fragment) for fragment in fragments])	
		
	def iupac(self):
		return modX2IUPAC(self)

def modX2IUPAC(pep):
	"""changes Xmod into the old-fashioned IUPAC amino acid alphabet.
	This destroys any information about chemical modifications or
	unnatural amino acids that the sequence holds."""
	if isinstance(pep, peptide):
		pep = ''.join(pep.sequence)
	elif isinstance(pep, list):
		pep = ''.join(pep)
	if '-' in pep:
		pep = pep[pep.find('-'):pep.rfind('-')]
	return ''.join([x for x in pep if x.isupper()])			


class digest:
	# What about those nifty -OH ends? the digested peptides won't have them
	# probably need a system to include them.

	def __gen_cuts(self):
		# I'm using reg-ex via the 're' module. because it is far more
		# powerful than str.find() 
		cuts = []
		for protease in self.proteases:
			[g.end() for g in finditer(expasy_rules[protease],
				self.sequence)]
		return cuts


	def gen_frags(self, skip=1):
	# skip = 0 returns an empty list.
	# skip = 1 means every potential site is cut.
	# skip = 2 means every second site is cut ect.
		solutions =[]
		for x in range(skip):
			frags = []
			last = 0
			for i in range(x, len(self.cuts), skip):
				frags.append(self.sequence[last:self.cuts[i]])
				last = self.cuts[i]
			frags.append(self.sequence[last:])
			solutions.extend(frags)
		return list(solutions)
	
	def vfrags(self, skip=2):
		"""vfrags is to be used to return all the fragments from a digest (up to a 
		given skiped cut count) as peptide objects with the parent aacomp and 
		isodistro."""
		frags = []
		for x in range(1, skip+1):
			frags.extend(gen_frags(x))
		return [peptide(frag, isodistro=self.protein.isodistro,
			aacomp=self.protein.aacomp) for frag in frags] 
			


	def __init__(self, peptide, *args):
		data = []
		for x in args: 
			if isinstance(x, list) or isinstance(x, tuple):
				data.extend(x)
			elif isinstance(x, dict):
				data.extend(x.values())
		
		proteases = [s for s in args if s in expasy_rules.keys()]
		del data

		try: assert len(peptides) == 1
		except AssertionError: raise PyteomicsError("You need to specify exactly 1 "
			"peptide to perform a digest.")

		try: assert proteases
		except AssertionError: raise PyteomicsError("You must specify at least 1 "
			"protease to perform a digest.")

		# Some weird stuff trying to make the __init__ take diverse input.
		# doesn't need to take multiple protein sequences. if you want to do that
		# then make two digests with different sequences.
		self.protein = peptides[0]; del peptides
		if self.protein.name: self.proteiname = protein.name
		else: self.proteiname = 'a generic peptide'

		self.sequence = ''.join(self.protein.sequence)

		self.cuts = self.__gen_cuts()
		
		print('cuts:',  self.cuts)
		self.fragments = ([peptide(frag, 
			aacomp = self.protein.aacomp, 
			isodistro = self.protein.isodistro)
			for frag in self.gen_frags()])

		self.av_len = str(int(sum([len(frag) for frag in 
			self.fragments]))/len(self.fragments))
	
	def __getitem__(self, value):
		if isinstance(value, int):
			return peptide(self.fragments[value])
		# Next for if the value was a sequence.
		else: return self.fragments.index(value)

	def __nonzero__(self):
		"""Checks if the instance is True or False."""
		return bool(self.fragments)

	def __repr__(self):
		"""Returns the ModX sequence."""
		return '|'.join(self.fragments)

	def __str__(self):
		if len(self.__repr__()) > 63:
			seq = self.__repr__()[:60] +'...'
		else: seq = self.__repr__()

		return seq

	def __contains__(self, x):
		if isinstance(x, int):
			return x in [int(mi_peptide_mass(frag)) for frag in self.fragments]
		if isinstance(x, str):
			return x in self.fragments
		
	def decribe(self):
		if len(self.proteases) == 1:
			prot = proteases[0]
		elif len(self.proteases) > 1:
			prot = '{0} and {1}'.format(', '.join(proteases[:-1]),
				proteases[-1])
		
		output = ([
			"Digest of {0} by {1}:".format(self.protein, prot),
			"Number of Fragments: {0}".format(str(len(self.fragments))),
			"Average Length of Fragments: {0} Residues".format(self.av_len),
			"sequence: {0}\n".format(self.__str__()) 
			])
		if self.protein == prot:
			output[0] = ("{0} self digest:".format(self.protein))
					
		return '\n'.join(output)

	def show_frags(self):
		for frag in sorted([peptide(frag) for frag in self.fragments]):
			print(frag)


class cross_link:

	def __init__(xlinker, *args, **kwargs):

		# make two digests then when you want to look for a mass:
		# first check if it is a digested peptide
		# then check if it is a peptide with crosslinker attached
		# ^ to do this you will need to be able to search for 
		# cross-linker attachment sites.
		# a crosslinker should probably have a regex rule for it's
		# attachment site. The rule needs to be modified at execution
		# so that the attachment site isn't next to a cleaved residue 
		# then check if it is a pair of crosslinked peptides 
		peptides = [arg for arg in args if (isinstance(arg, peptide))]
		proteases = [arg for arg in args if (isinstance(arg, protease))]

		# This is the way I will allow the switching on and off of the 
		# deeper search functions. should allow turning off the shallow\
		# options too.
		self.skipped_cuts = kwargs.get('skipped_cuts', 1)
		self.tandem_xlinks = kwargs.get('tandem_xlinks', False)
		self.pure_peptides = kwargs.get('pure_peptides', True)
		self.crowding = kwargs.get('crowding', True)
		
		self.precision = 0.5 # in Da
		
		self.mass = xlinker[0]
		self.caps = xlinker[1]
		self.rule = xlinker[2]	

		self.digests = [digest(peptide, proteases)
                                for peptide in peptides]
		
	def find(self, target):
		#frags = [[frag for frag in digest] for digest in self.digests]
		#frags = [item for sublist in frags for item in sublist]
		# ^ Incomprehensible list comprehension. flattens a nested list.

                ### 2020.01.05 Debug this threw an error. I tried correcting
                ### it with little memory of the purpose of the code (see below)
                
		# frags = [frag for digest.vfrags(self.skipped_cuts)
                #          in self.digests]
		# 	#for frag in digest.vfrags(self.skipped_cuts)] #hybrid
		# print(frags)

                frags = [frag for frag in digest.vfrags(self.skipped_cuts)]
		print(frags)

		p = self.precision / 2.0000
		pure = inrange(target - p, target + p)
		solutions = {(0, 0, pep) for pep in frags if pure(pep)}

		frags = [frag for frag in frags if re.search(self.xlinker.rule, 
			frag.__repr__())] 
		# ^ This is a list of lists with peptides that have a crosslink site in them.

		frags = list({(pep.mass(), pep) for pep in frags if pep < target})
		
		# I Could rewrite this bit below to use carry_on but since i need to
		# reverse/forward sort i might as well leave ok() in tact as it is
		# easier to follow what is going on with this code.
 		if self.mass < 0:
			ok = lambda x: x < target
			frags.sort(None, None, True) # reverse sorting
		else: 
			ok = lambda x: x > target
			frags.sort(None, None, False) # forward sorting

		xlinker = (self.mass, self.caps, self.rule)

		while frags:
			pep1 = frags.pop(0)
                        # using pop so it never gets checked again.
			solutions.update(x_link(target, xlinker,
                                                pep1[1], pep1[1]))
			for pep2 in frags:
				if ok(pep1[0] + pep2[0]):
					solutions.update(x_link(target, xlinker,
						                self.precision,
                                                                pep1, pep2))
			# Remember that at this point pep1/2 are (mass, peptide)
                        # tuples to minimise clock time.
				else: break
		return solutions


def x_link(target, xlinker, p, *args):
	mass = xlinker[0]
	caps = xlinker[1]
	rule = xlinker[2]
	solutions = []

	peps, pm = [x[1] for x in args], sum([x[0] for x in args])

	is_solution = inrange(target - p, target + p)

	# Total X-Link Mass
	txlm = lambda x, y: sum([mass*x, caps*y, pm])

	# Now work out how many places there are in the peptides to put x_linker
	sites = sum([len([re.findall(rule, ''.join(x.sequence)) for x in peps]
        )])

	for i in range(1, sites+1):
		m = txlm(i, 0)

		if is_solution(m): 
			solutions.append(tuple([i, 0].extend(peps)))

		for x in range(1, i+1):
			z = txlm(i, x)
				
			if is_solution(z):
				solutions.append(tuple([i, 0].extend(peps)))

			elif not carry_on(z, target, caps): break

		        elif not carry_on(m, target, mass): break

	return solutions

print calculate_mass('GGGGGGGGGGGGGGGGGGG')
p = peptide('GGGGFSDALKPMCVVSTRGGGGGGGGGG')
print p
print p.mass(0, 'av')
d = digest(p, 'trypsin')
print d
