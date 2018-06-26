import argparse
import bioc
import re
import itertools

#pattern = re.compile("[A-Z0-9]* [ia]s a tumor antigen")
pattern1 = re.compile("(?P<term>\S+) [ia]s a (\S*\s)?(\S*\s)?(\S*\s)?(\S*\s)?tumor antigen")
#pattern2 = re.compile("the tumor antigen (?P<term>\S+)")
#pattern3 = re.compile("the tumor antigen gene (?P<term>\S+)")
#pattern4 = re.compile("the tumor antigen protein (?P<term>\S+)")
#patterns = [pattern1,pattern2,pattern3,pattern4]
patterns = [pattern1]

stoplist = {'gene','protein','rna','is','may','cell','cells'}

def searchForTumorAntigens(document,acceptableTerms):
	assert isinstance(document,bioc.BioCDocument)
	pmid = document.infons['pmid']
	title = document.infons['title']
	journal = document.infons['journal']
	year = document.infons['year']
	for passage in document.passages:
		text = passage.text
		if "tumor antigen" in text:
			sentences = text.split('.')
			sentences = [ s.strip() for s in sentences if "tumor antigen" in s ]
			for s,pattern in itertools.product(sentences,patterns):			
				search = pattern.search(s)
				if search:
					term = search.groupdict()['term']
					termLower = term.lower()
					if termLower in acceptableTerms and not termLower in stoplist:
						out = [term,pmid,title,journal,year,s]
						outTxt = "\t".join(out)
						print(outTxt)
						#print("%%s\t%s" % (pmid,term,s))
	#print(document.infons)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Search for mentions of tumor antigens')
	parser.add_argument('--biocFile',required=True,type=str,help='Filename of BioC file to search')
	parser.add_argument('--acceptableTermsFile',required=True,type=str,help='List of terms to help filter')
	args = parser.parse_args()

	with open(args.acceptableTermsFile) as f:
		acceptableTerms = set([ line.strip().lower() for line in f ])

	with bioc.iterparse(args.biocFile) as parser:
		for document in parser:
			searchForTumorAntigens(document,acceptableTerms)

