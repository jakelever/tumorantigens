import sys
import itertools
import kindred
import pickle
import argparse
import codecs
import time
import re
import string
from collections import defaultdict,Counter
import json

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

def filterCorpus(corpus,filterTerms):
	filtered = kindred.Corpus()
	for doc in corpus.documents:
		termsFound = any( ft in doc.text.lower() for ft in filterTerms )
		if termsFound:
			filtered.addDocument(doc)
	return filtered

def parseAndFindEntities(biocFile,wordlistPickle,outSentencesFilename):
	print("%s : start" % now())

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)

	#with open(filterTermsFile,'r') as f:
	#	filterTerms = [ line.strip().lower() for line in f ]
	strictFilterTerms = ['tumor antigen','tumour antigen','tumor-antigen','tumour-antigen']
	weakFilterTerms = ['antigen']

	timers = Counter()

	outSentences = []

	currentID = None
	duplicateCheck = set()

	print("%s : processing..." % now())
	parser = kindred.Parser()
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=True,detectMicroRNA=False,acronymDetectionForAmbiguity=True,mergeTerms=True)
	for corpusno,corpus in enumerate(kindred.iterLoadDataFromBioc(biocFile)):
		startTime = time.time()
		corpus = filterCorpus(corpus,weakFilterTerms)
		timers['filter'] += time.time() - startTime

		startTime = time.time()
		parser.parse(corpus)
		timers['parser'] += time.time() - startTime
		print("%s : parsed" % now())

		startTime = time.time()
		ner.annotate(corpus)
		timers['ner'] += time.time() - startTime
		print("%s : ner" % now())

		startTime = time.time()

		for doc in corpus.documents:

			# Reset the duplicate check set for each new PMID
			if doc.metadata['id'] != currentID:
				currentID = doc.metadata['id']
				duplicateCheck = set()

			for sentence in doc.sentences:
				sentenceTextLower = sentence.text.lower()



				#print(sentence.text)
				#print(sentence.entitiesWithLocations)
				entityTypesInSentence = set([ entity.entityType for entity,tokenIndices in sentence.entityAnnotations ])
				gotKeyword = 'keyword' in entityTypesInSentence
				gotGene = 'gene' in entityTypesInSentence
				gotProtein = 'protein' in entityTypesInSentence
				gotCancer = 'cancer' in entityTypesInSentence

				containsStrictTerm = any( ft in sentenceTextLower for ft in strictFilterTerms )
				containsWeakTerm = any( ft in sentenceTextLower for ft in weakFilterTerms )

				topicMatch = containsStrictTerm or (containsWeakTerm and gotCancer)

				if topicMatch and gotKeyword and (gotGene or gotProtein):
					sentenceText = sentence.text.strip(string.whitespace + ',')

					if not sentenceText in duplicateCheck:
						tmpData = dict(doc.metadata)
						tmpData['sentence'] = sentenceText
						outSentences.append(tmpData)
						duplicateCheck.add(sentenceText)

		timers['entitiesAdded'] += time.time() - startTime

		print("%s : entities added" % now())
		sys.stdout.flush()

	with open(outSentencesFilename,'w') as f:
		json.dump(outSentences,f,indent=2)

	print("%s : done" % now())
	
	for section,sectiontime in timers.items():
		print("%s\t%f" % (section,sectiontime))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	#parser.add_argument('--filterTerms',required=True)
	parser.add_argument('--wordlistPickle',required=True)
	parser.add_argument('--outSentencesFilename',required=True)

	args = parser.parse_args()

	parseAndFindEntities(args.biocFile,args.wordlistPickle,args.outSentencesFilename)

