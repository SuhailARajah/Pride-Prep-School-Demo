######################################################
# Exploring the Central Dogma of Biology with Python
#
# Look over this code, what do you think it does?
# Make a hypothesis about what some of the commands do.
# Talk about it with the person(s) near you.
# It's ok if you make a wrong guess or hypothesis!
# The idea is to learn by tinkering!
#
# Requires: a file called dna.txt that has on the first
# line an organism and some dna as follows:
# human ACATGCTAGAATAGCCGCATGTACTAGTTAA
######################################################

def reverse(s):
    return s[::-1]

def reverse_complement(s):
  # Make a dictionary so we can lookup complements easily
  lookup_complement = { 'A':'T', 'C':'G', 'G':'C', 'T':'A'}
  rc = ""
  for base in reverse(s):
    rc = rc + lookup_complement[base]
  return rc

def transcribe(s) :
  # Make a dictionary for transcribing the strands
  transcription = { 'T':'U', 'A':'A', 'C':'C', 'G':'G' }
  rna = ""
  for base in s:
    rna = rna + transcription[base]
  return rna

def translate(s) :
  translation = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"STOP", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  rna = s
  start = rna.find('AUG')  # Why is this looking for AUG? Any Ideas?
  proteins = []
  while start != -1:
      # found ORF at location start (What the heck is an ORF?)
      protein = ""
      loc = start
      while start+2 < len(rna):
          codon = rna[start:start+3]
          start += 3
          # check if this is a stop codon (Stop what?)
          if translation[codon] == "STOP": break;
          protein = protein + translation[codon]
      # add location and protein to list (a list of what?)
      proteins.append( (loc, protein) )
      # search for next ORF in this rna because of frame shifts (I don't see a picture frame here?)
      # start searching just one location past our last ORF
      start = rna.find('AUG',loc+1)
  return proteins

# open a file called dna.txt must be in the same
# folder as this python script
with open ("./dna.txt", "r") as myfile:
    data = myfile.readlines()
    seqs = [ seq.split() for seq in data ]

# get organism and dna from first line of the file
organism = seqs[0][0].strip()
dna = seqs[0][1].strip()

# Create the reverse complement of the dna string
rc_dna = reverse_complement(dna)

print( organism, "  forward dna: ", dna )
print( organism, " rev comp dna: ", rc_dna )

# Transcribe both dna strings into rna strings
forward_rna = transcribe(dna)
reverse_rna = transcribe(rc_dna)

print( organism, "  forward rna: ", forward_rna )
print( organism, "  reverse rna: ", reverse_rna )

# Translate both rna strings into lists of potential proteins
forward_proteins = translate(forward_rna)
reverse_proteins = translate(reverse_rna)

print( organism, " possible forward proteins: ", forward_proteins )
print( organism, " possible reverse proteins: ", reverse_proteins )
