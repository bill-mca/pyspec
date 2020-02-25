from pyteomics.mass import Composition

hha = "MIHLFSAGIKMAKTKQEWLYQLRRCSSVNTLEKII"\
	"HKNRDSLLNSERESFNSAADHRLAELITGKLYD"\
	"RIPKEIWKYVR"

mtypekinase = "MPFGNTHNKYKLNYKSEEEYPDLSKHN"\
		"NHMAKVLTPDLYKKLRDKETPSGFT"\
		"LDDVIQTGVDNPGHPFIMTVGCVAG"\
		"DEESYTVFKDLFDPIIQDRHGGFKP"\
		"TDKHKTDLNHENLKGGDDLDPHYVL"\
		"SSRVRTGRSIKGYTLPPHCSRGERR"\
		"AVEKLSVEALNSLTGEFKGKYYPLK"\
		"SMTEQEQQQLIDDHFLFDKPVSPLL"\
		"LASGMARDWPDARGIWHNDNKSFLV"\
		"WVNEEDHLRVISMEKGGNMKEVFRR"\
		"FCVGLQKIEEIFKKAGHPFMWNEHL"\
		"GYVLTCPSNLGTGLRGGVHVKLAHL"\
		"SKHPKFEEILTRLRLQKRGTGGVDT"\
		"AAVGSVFDISNADRLGSSEVEQVQL"\
		"VVDGVKLMVEMEKKLEKGQSIDDMI"\
		"PAQK"

tomb = "MHHHHHHMDEYSPKRHDIAQLKFLCETLYHDCLANLEESNHG"\
	"WVNDPTSAINLQLNELIEHIATFALNYKIKYNEDNKL"\
	"IEQIDEYLDDTFMLFSSYGINMQDLQKWRKSGNRLFRCF"\
	"VNATKENPASLSC"

# CROSS LINKERS

# cross linkers are tuples of: 
# [0]mass, when it is a bridge.
# [1]cap, any additional mass present when an end is bot bound to a peptide.
# [2]rule, a regex rule that describes the residues to which this CL binds.
DiSulf = (-2, 0, 'C')

bs3 = (572.43, 1, '[RK]')

nist_mass = {
	'H': {1: (1.0078250320710, 0.99988570),
		2: (2.01410177784, 0.00011570),
		3: (3.016049277725, 0.0),
		0: (1.0078250320710, 1.0)},
    
	'H+': {1: (1.00727646677, 1.0),
		0: (1.00727646677, 1.0)},

	'Neu':{1: (1.00727646677, 1.0),
		0: (1.00727646677, 1.0)},
	# Use Neu if you want to mak just one 
	# atom in one aa heavier. eg N15 leu tag.

	'e*': {0: (0.00054857990943, 1.0)}, 

	'C': {12: (12.0000000, 0.98938),
		13: (13.0033548378, 0.01078),
		14: (14.0032419894, 0.0),
		0: (12.0000000, 1.0)},

	'N': {14: (14.00307400486, 0.9963620),
		15: (15.00010889827, 0.0036420),
		0: (14.00307400486, 1.0)},

	'O': {16: (15.9949146195616, 0.9975716),
		17: (16.9991317012, 0.000381),
		18: (17.99916107, 0.0020514),
		0: (15.9949146195616, 1.0)},
    
	'P': {31: (30.9737616320, 1.0000),
		0: (30.9737616320, 1.0000)},

	'S': {32: (31.9720710015, 0.949926),
	          33: (32.9714587615, 0.00752),
	          34: (33.9678669012, 0.042524),
	          36: (35.9670807620, 0.00011),
	           0: (31.9720710015, 1.0)},

	'Se': { 74:   (73.9224764,  0.0089),
		76:   (75.9192136,  0.0937), 
		77:   (76.9199140,  0.0763), 
		78:   (77.9173091,  0.2377), 
		80:   (79.9165213,  0.4961),
		82:   (81.9166994,  0.0873),
		80:   (79.9165213,  0.4961)}
    }


std_aa_comp = {
    'A':   Composition({'H': 5, 'C': 3, 'O': 1, 'N': 1}),
    'C':   Composition({'H': 5, 'C': 3, 'S': 1, 'O': 1, 'N': 1}),
    'D':   Composition({'H': 5, 'C': 4, 'O': 3, 'N': 1}),
    'E':   Composition({'H': 7, 'C': 5, 'O': 3, 'N': 1}),
    'F':   Composition({'H': 9, 'C': 9, 'O': 1, 'N': 1}),
    'G':   Composition({'H': 3, 'C': 2, 'O': 1, 'N': 1}),
    'H':   Composition({'H': 7, 'C': 6, 'N': 3, 'O': 1}),
    'I':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'K':   Composition({'H': 12, 'C': 6, 'N': 2, 'O': 1}),
    'L':   Composition({'H': 11, 'C': 6, 'O': 1, 'N': 1}),
    'M':   Composition({'H': 9, 'C': 5, 'S': 1, 'O': 1, 'N': 1}),
    'seM': Composition({'H': 9, 'C': 5, 'Se': 1, 'O': 1, 'N': 1}),
    'N':   Composition({'H': 6, 'C': 4, 'O': 2, 'N': 2}),
    'P':   Composition({'H': 7, 'C': 5, 'O': 1, 'N': 1}),
    'Q':   Composition({'H': 8, 'C': 5, 'O': 2, 'N': 2}),
    'R':   Composition({'H': 12, 'C': 6, 'N': 4, 'O': 1}),
    'S':   Composition({'H': 5, 'C': 3, 'O': 2, 'N': 1}),
    'T':   Composition({'H': 7, 'C': 4, 'O': 2, 'N': 1}),
    'V':   Composition({'H': 9, 'C': 5, 'O': 1, 'N': 1}),
    'W':   Composition({'C': 11, 'H': 10, 'N': 2, 'O': 1}),
    'Y':   Composition({'H': 9, 'C': 9, 'O': 2, 'N': 1}),
    'H-':  Composition({'H': 1}),
    '-OH': Composition({'O': 1, 'H': 1}),
    }
