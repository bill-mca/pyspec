# pyspec version 3: pyteomics mod.
from pyteomics.parser import *
from pyteomics.mass import *
from pyteomics.auxiliary import PyteomicsError
from re import finditer, findall
from data import *

# If you are here to understand this program this is a link to the
# pyteomics tutorial: http://pythonhosted.org/pyteomics/


   ## TO DO ##

# A cool structure would be nist_mass defind as 
# a class with __getitem__ so that it could automatically
# return masses of charges
# e.g. if '+' in 'Mg++' -> masses['Mg++'.removeall('+')] - 'e*' * 'Mg'.count('+') 
# it could then count cahrges to add into calculate_mass
# and would automatically support Na, Mg ect..

# ISOFORMS

# ways of saving aacomp and isodistro
# ways of saving peptides
# open fastas
# re-do what-weighs
# make it so that xlink.find() checks a list of masses all at once.
# make X-Link take just a list of peptides.
# write out a bunch of data: modified aa, masses of elements, 
# refine x-link so that the toggle settings work. esp. crowding

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def carry_on(mass, target, increment):
	"""Carry on is a little function to help handling using either -ve or +ve values 
	for the mass of a crosslinker (for disulfides). It saves defining 4 or more 
	lambdas. It only returns true if the target could be reached by adding more
	'increment' to mass."""
	if increment > 0 and mass < target: 
	      return True

	elif increment < 0 and mass > target: 
		return True

	else: return False

class inrange:
	def __init__(self, bottom, top):
		assert all([is_number(top), is_number(bottom)])
		self.top = top
		self.bottom = bottom

	def __repr__(self):
		return "inrange({0}, {1})".format(self.bottom, self.top)

	def __call__(self, value):
		return value < self.top and value >= self.bottom

class Environment:
	def __init__(self, isodistro=None, aacomp=None):
		if aacomp:
		    self.aacomp = aacomp
		else: aacopm = std_aa_comp
		
		if isodistro:
		    self.isodistro = isodistro
		else: isodisto = nist_mass 

	def new(self, sequence, name=''):
		return Peptide(sequence, name, 
			aacomp=self.aacomp,
			isodisto=self.isodistro)

class Peptide:
	
	def __init__(self, sequence, name='', **kwargs):
		"""Most of the functionality of this module. The an amino acid sequence and 
		turns it into a	peptide object. The advantage of using the peptide object 
		is that it is a more pythonic representation of a protein than just a naked
		sequence. One important feature is that the peptide can be fed an 
		"isodistro" or "aacomp" as kwargs at initialisation. these are dictionaries
		that allow the user to specify the propensity of different isotopes and the
		chemical composition of the amino acids. one can even make up whole new 
		amino acids by adding a new entry to the aacomp dictionary. The modx system
		allows the description of modified / unnatural amino acids by prepending
		lowercase characters to an amino acids one letter name eg M is methionine,
		oM could be oxidised methionine and ooM double oxidised.
		"""

		if isinstance(sequence, str):
			self.sequence = parse(sequence, 					  					show_unmodified_termini=False, 
				split=False,
				allow_unknown_modifications=True,
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

	def __nonzero__(self):
		return bool(self.sequence)

	def __key(self):
		return (repr(self), self.name)

	def __hash__(self):
		return hash(self.__key())

	def __eq__(self, other):
		"""Because this program is mass spec oriented I have decided to make __eq__
		compare the mass of the two peptides rather than checking if their sequences are
		identical this means that different sequences can be indistinguishable (as in an
		MS machine). 
		e.g. PEPTIDE == EPPEDTI -> True or LLLL == IIII -> True
		"""
		if isinstance(other, Peptide):
			return self.mass() == other.mass()
		else: return self.mass() == other

	def __gt__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, Peptide):
			return self.mass > other.mass
		else: return self.mass() > other

	def __lt__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, Peptide):
			return self.mass < other.mass
		else: return self.mass() < other

	def __ge__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, Peptide):
			return self.mass >= other.mass
		else: return self .mass()>= other

	def __le__(self, other):
		"""compares peptide's mass to another peptide or a number."""
		if isinstance(other, Peptide):
			return self.mass <= other.mass
		else: return self.mass() <= other
	
	def __cmp__(self, other):
		"""Compares peptides to each other. Those with greater mass score higher.
		This is the function called when a list is sorted. So a list of peptide
		objects will be sorted from lowest mass to highest.
		""" 
		try: assert isinstance(other, Peptide)
		except AssertionError: raise TypeError("Can't" 
			"compare a peptide to a non peptide.")
		return cmp(self.mass(), other.mass())

	def __repr__(self):
		return ''.join(self.sequence)

	def __str__(self):
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

		peptide[3:7] --> ''.join(peptide.sequence)[3:7]

		peptide[7] --> ''.join(peptide.sequence)[7]
		"""
		if isinstance(value, str): 
		 	if value in self.aacomp.keys():
				return self.sequence.count(value)
			else: raise TypeError

		else: return self.sequence[value]
	
	def __len__(self):
		"""len(peptide) --> an integer indicating the sequence length"""
		if self.sequence:
			return length(self.sequence, labels=self.aacomp.keys())
		else: return 0

	def __contains__(self, x):
		"""true if x is a substring of the peptide's sequence false otherwise"""
		if self.find(x) != -1:
			return True
		else: return False

	def mass(self, charge=0, *args, **kwargs):
		"""A funtion for finding the mass of the sequence. Charge defults to
		zero however specifying an integer here simulates addition of H+
		on a mass spectrum. 
			i.e. p.mass(3) == (p.mass()+3)/3
		
		This function draws from the isotopic distribution specified for the
		peptide instance. this way if many peptides were created with
		different conditions (e.g. N15 tag on one, another has heavy leucine)
		they could easily be kept track of and passed to cross_link together.

		Passing "av" or "average" to args will return the average mass

		Passing "mp" or "most probable" to args will return th most probable 
		mass for the peptide (the highest peak in the mass spectrum).
		
		Note that for the above the charge argument must be specified and
		still works.
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

		for key in set(aacomp.keys()).intersection(kwargs.keys()):
			aacomp[key] = kwargs[key]

		if {'av', 'average', 'Average', 'AV'}.intersection(args):
			return calculate_mass(self.sequence, charge=charge,
						average = True,
						aa_comp=self.aacomp,
						mass_data=isodistro)

		elif {'mp','mpm','MPM','MP','probable','most probable',
			'most probable mass'}.intersection(args):
			mp = most_probable_isotopic_composition(self.sequence,
							aa_comp=aacomp, 
							mass_data=isodistro)
			return calculate_mass(mp, charge=charge, 
				aa_comp=self.aacomp, mass_data=isodistro)

		else:  # mono-isotopic
			return calculate_mass(self.sequence, charge=charge,
				aa_comp = aacomp, 
				mass_data = isodistro,)

	def count(self, x):
		"""An easy way of counting the frequency of an amino acid."""
		return self.sequence.count(x)

	def modify(self, old, new):
		"""Returns a modified copy of the peptide with all occurances of 
		The substring "old" swapped for the string "new". This is the 
		easiest way of applying a fixed modification to a peptide. 
		To perform a variable modification: peptide[index] = mod_aa"""
		assert isinstance(old, str) and isinstance(new, string)
		return Peptide((''.join(self.sequence)).replace(old, new))

	def startswith(self, prefix, start=None, end=None):
		return ''.join(self.sequence).startswith(prefix, start, end)
	# These two need to be modernised for the modX stuff.
	def endswith(self, suffix, start=None, end=None):
		return ''.join(self.sequence).endswith(suffix, start, end)


	def find(self, subseq):
		assert isinstance(subseq, string)
		return ''.join(self.sequence).find(subseq)



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
		
		return sorted([Peptide(fragment) for fragment in fragments])	
		
	def iupac(self):
		return modX2IUPAC(self)

def modX2IUPAC(pep):
	"""changes Xmod into the old-fashioned IUPAC amino acid alphabet.
	This destroys any information about chemical modifications or
	unnatural amino acids that the sequence holds."""
	if isinstance(pep, Peptide):
		pep = ''.join(pep.sequence)
	elif isinstance(pep, list):
		pep = ''.join(pep)
	if '-' in pep:
		pep = pep[pep.find('-'):pep.rfind('-')]
	return ''.join([x for x in pep if x.isupper()])			


class digest:

	def __gen_cuts(self):
		# I'm using reg-ex via the 're' module. because it is far more
		# powerful than str.find() 
		cuts = set([])
		
		for protease in self.proteases:
			cuts.update({g.end() for g in finditer(expasy_rules[protease],
				self.sequence)})
		return sorted(list(cuts))


	def gen_frags(self, step=1):
	# skip = 0 returns an empty list.
	# skip = 1 means every potential site is cut.
	# skip = 2 means every second site is cut ect.
		solutions = []
		for x in range(step):
		  # If we are stepping by 3 then there are 3 reading frames
		  # this loop cycles through the reading frames.
			frags = []
			last = 0
			for i in range(x, len(self.cuts), step):

				frags.append(self.sequence[last:self.cuts[i]])
				last = self.cuts[i]
			frags.append(self.sequence[last:])
			solutions.extend(frags)
		return solutions
	
	def vfrags(self, skip=2):
		"""vfrags is to be used to return all the fragments from a digest (up to a 
		given skipped cut count) as peptide objects with the parent aacomp and 
		isodistro."""
		frags = []
		for x in range(1, skip+1):
			frags.extend(self.gen_frags(x))
		return [self.environment.new(frag, 'A fragment of' + self.name)\
		  for frag in frags] 
			

	
	def __init__(self, peptide, *args):
		"""Digest returns lists of peptide instance fragments that 
		inherit their isotopic information and modifications from
		their parent. 
		"""
		data = []
		for x in args: 
		    if isinstance(x, set):
			data.extend(list(x))

		# it would be nice if you could digest by giving a regex rule as an arg
		for x in args: 
			if isinstance(x, list) or isinstance(x, tuple):
				data.extend(x)
			elif isinstance(x, dict):
				data.extend(x.values())
		data.extend(args)
		self.proteases = [x for x in data if x in expasy_rules.keys()]
		del data

		
		try: assert self.proteases
		except AssertionError: raise PyteomicsError("You must specify at least 1 "
			"protease to perform a digest.")
		
		self.environment = Environment(peptide.isodistro, peptide.aacomp)
		
		if peptide.name: self.name = peptide.name
		else: self.name = 'Generic Peptide'

		self.sequence = ''.join(peptide.sequence)

		self.cuts = self.__gen_cuts()
		
		self.fragments = ([self.environment.new(frag, 'Fragment of' + self.name) 
				  for frag in self.gen_frags()])

		self.av_len = str(int(sum([len(frag) for frag in 
			self.fragments]))/len(self.fragments))
	
	def __iter__(self):
		return self.fragments
	      
	def __getitem__(self, value):
		if isinstance(value, int):
			return Peptide(self.fragments[value])
		# Next for if the value was a sequence.
		else: return self.fragments.index(value)

	def __nonzero__(self):
		"""Checks if the instance is True or False."""
		return bool(self.fragments)

	def __repr__(self):
		"""Returns the ModX sequence."""
		return '|'.join([str(frag) for frag in self.fragments])

	def __str__(self):
		if len(self.__repr__()) > 63:
			seq = self.__repr__()[:60] +'...'
		else: seq = self.__repr__()

		return seq

	def __contains__(self, x):
		if isinstance(x, int):
			return x in [frag.mass() for frag in self.fragments]
		if isinstance(x, str):
			return x in [str(frag) for frag in self.fragments]
		
	def describe(self):
		if len(self.proteases) == 1:
			prot = self.proteases[0]
		elif len(self.proteases) > 1:
			prot = '{0} and {1}'.format(', '.join(self.proteases[:-1]),
				self.proteases[-1])
		
		output = ([
			"Digest of {0} by {1}:".format(self.name, prot),
			"Number of Fragments: {0}".format(str(len(self.fragments))),
			"Average Length of Fragments: {0} Residues".format(self.av_len),
			"sequence: {0}\n".format(self.__repr__()) 
			])
		if self.name == prot:
			output[0] = ("{0} self digest:".format(prot))
					
		return '\n'.join(output)

	def show_frags(self):
		for frag in sorted(self.fragments):
			print(frag.describe())


class cross_link:

	def __init__(self, xlinker, *args, **kwargs):
		"""The cross_link class is for analysing the mass spectrum of
		a peptide cross_linking experiment. Since there are too many 
		feasible combinations of bridge-fragment the class is just a 
		grouping of conditions that then generate the alternatives to 
		test if any match a mass of interest.

		The "xlinker" is a tuple of (mass, caps, rule) where:
			- Mass is the mass of the crosslinker when it spans
			  two adjacent peptides.
			- Caps is any mass that the crosslinker has added to 
			  it if one end is not attached (e.g. BS3 is hydrolysed)
			- Rule is a regex rule that describes the sites to which 
			  the crosslink bridge is attached.

		WARNING! be careful when using this code because it produces false
		positives. Critically assess the results. A known issue is that 
		solutions are returned in which a protease cuts a site that a 
		bridge is attached to. There are also false negatives because the
		function does not allow for tandem crosslinking such as long daisy 
		chains of peptide-bridge-peptide-bridge-peptide-bridge...
		"""
		assert isinstance(tuple, xlinker)
		peptides = [arg for arg in args if (isinstance(arg, Peptide))]
		proteases = set(expasy_rules.keys()).intersection(args)
		assert peptides
		assert proteases
		# This is the way I will allow the switching on and off of the 
		# deeper search functions. should allow turning off the shallow\
		# options too.
		self.skipped_cuts = kwargs.get('skipped_cuts', 1)
#		self.tandem_xlinks = kwargs.get('tandem_xlinks', False)
		self.pure_peptides = kwargs.get('pure_peptides', True)
#		self.crowding = kwargs.get('crowding', True)
		
		self.mass = xlinker[0]
		self.caps = xlinker[1]
		self.rule = xlinker[2]	

		self.digests = [digest(peptide, proteases) for peptide in peptides]

	def __call__(self, result):
		"""This was quickly cobbled together to fascilitate checking results
		that an instance turns out from find.

		XL(XL.find(2500, 0)[0]) --> 2500"""

		pm = sum([pep.mass() for pep in result [2:]])
		return sum([self.mass*result[0], self.caps*result[1], pm])
	      
	def find(self, target, precision=1):
		"""A solution returned by find is a tuple with 2 integers at [0] and [1].
		the integers denote the stociometries of cross-linker and unattached ends of
		crosslinker respectivly. Thus [1] should never be higher than [0] and for the
		peptides to be joined [0] should be >= 1."""
		#frags = [[frag for frag in digest] for digest in self.digests]
		#frags = [item for sublist in frags for item in sublist]


		frags = [frag for digest in self.digests
			for frag in digest.vfrags(self.skipped_cuts)] #hybrid
		# ^ Incomprehensible list comprehension. flattens a nested list.

		p = precision/2.000
		pure = inrange(target - p, target + p)
		solutions = {(0, 0, pep) for pep in frags if pure(pep)}

		frags = [frag for frag in frags if re.search(self.rule, 
			frag.__repr__())] 
		# ^ This is a list of lists with peptides that have a crosslink site in them.

					# The +3 ensures you catch all disulfides
		frags = list({(pep.mass(), pep) for pep in frags if pep < target+3})

		# I Could rewrite this bit below to use carry_on but since i need to
		# reverse/forward sort i might as well leave ok() in tact as it is
		# easier to follow what is going on with this code.
 		if self.mass < 0:
			ok = lambda x: x > target
			frags.sort(None, None, True) # reverse sorting
		else: 
			ok = lambda x: x < target
			frags.sort(None, None, False) # forward sorting

		xlinker = (self.mass, self.caps, self.rule)

		while frags:
			pep1 = frags.pop(0) # using pop so it never gets checked again.
			solutions.update(x_link(target, xlinker, p, pep1, pep1))
			for pep2 in frags:
				if ok(pep1[0] + pep2[0]):
					solutions.update(x_link(target, xlinker, 
						p, pep1, pep2))
			# Remember that at this point pep1/2 are (mass, peptide) tuples.
			# To minimise clock time.
				else: break
				
		return list(solutions)


def x_link(target, xlinker, p, *args):
	"""This function is just a backend and isn't really useful to humans. 
	args must be a (mass, peptide) tuple"""
	mass = xlinker[0]
	caps = xlinker[1]
	rule = xlinker[2]
	solutions = []

	# pm is the mass of all the peptides, peps is all the peptides in a list.
	# Henceforth they need not be associated so: 
	peps, pm = [x[1] for x in args], sum([x[0] for x in args])

	is_solution = inrange(target - p, target + p)

	# Total X-Link Mass
	txlm = lambda x, y: sum([mass*x, caps*y, pm])

	# Now work out how many places there are in the peptides to put x_linker 
	sites = sum([len(re.findall(rule, ''.join(x.sequence))) for x in peps])

	for i in range(1, sites+1):
		m = txlm(i, 0)
		#print (i, 0), sites, args
		#print m
		if is_solution(m): 
			print m
			solutions.append(tuple([i, 0] + peps))
			
		if caps:
		    if len(args) > 1: g = i
		    else: g = i+1

		    for x in range(1, g):
			z = txlm(i, x)
			#print (i, x), args
			#print z
			if is_solution(z):
				#print txlm(i,x)
				solutions.append(tuple([i, x] + peps))

			elif not carry_on(z, target, caps): break

		if not carry_on(m, target, mass): break

	return solutions
