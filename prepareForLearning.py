import kindred
import argparse
import os
import json
import time
import pickle
from collections import defaultdict

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Prepare data for active learning')
	parser.add_argument('--sentenceDir',required=True,type=str,help='Directory containing JSON files')
	parser.add_argument('--wordlistPickle',required=True,type=str,help='Pickle on entity names')
	parser.add_argument('--outRelationsVectorsAndMetadata',required=True,type=str,help='Output stuff for quick loading')
	args = parser.parse_args()

	print("%s : started" % now())
	with open(args.wordlistPickle,'rb') as f:
		oldTermLookup = pickle.load(f)
	print("%s : wordlist loaded" % now())

	termLookup = {}
	for term,details in oldTermLookup.items():
		detailsDict = { termtype:termid for termtype,termid in sorted(list(details)) }
		if 'gene' in detailsDict:
			newDetails = set([('geneOrProtein', detailsDict['gene'])])
			termLookup[term] = newDetails
		elif 'protein' in detailsDict:
			newDetails = set([('geneOrProtein', detailsDict['protein'])])
			termLookup[term] = newDetails
		else:
			termLookup[term] = details
	print("%s : simplifying wordlist" % now())

	corpus = kindred.Corpus()
	alreadySeenSentences = set()
	for i,f in enumerate(sorted(os.listdir(args.sentenceDir))):
		if not f.endswith('.json'):
			continue

		sentenceFile = os.path.join(args.sentenceDir,f)
		with open(sentenceFile) as f:
			sentenceData = json.load(f)

		for sentence in sentenceData:
			if sentence["sentence"] in alreadySeenSentences:
				continue
			alreadySeenSentences.add(sentence["sentence"])

			# Check for probable lists
			commaCount = sentence["sentence"].count(',')
			semiColonCount = sentence["sentence"].count(';')
			if commaCount >= 5 or semiColonCount >= 3:
				continue

			metadata = dict(sentence)
			del metadata["sentence"]
			doc = kindred.Document(sentence["sentence"],metadata=metadata)
			corpus.addDocument(doc)

		#if i > 100:
		#	break
	print("%s : corpus loaded" % now())

	parser = kindred.Parser()
	parser.parse(corpus)
	print("%s : parsed" % now())

	ner = kindred.EntityRecognizer(lookup=termLookup,detectVariants=False,detectFusionGenes=False,detectMicroRNA=False,acronymDetectionForAmbiguity=True,mergeTerms=True,removePathways=True)
	ner.annotate(corpus)
	print("%s : ner" % now())

	candidateBuilder = kindred.CandidateBuilder(entityCount=2,acceptedEntityTypes=[('keyword','geneOrProtein')])
	unfilteredCandidateRelations = candidateBuilder.build(corpus)
	print("%s : candidateBuilder (%d)" % (now(),len(unfilteredCandidateRelations)))

	entityIDToDoc = {}
	for doc in corpus.documents:
		for entity in doc.entities:
			entityIDToDoc[entity.entityID] = doc

	candidateRelations,metadata = [],[]
	for cr in unfilteredCandidateRelations:
		entityIDToEntity = { entity.entityID:entity for entity,tokenIndices in cr.sentence.entityAnnotations }
		entityIDToTokenLocs = { entity.entityID:tokenIndices for entity,tokenIndices in cr.sentence.entityAnnotations }
		keyword,geneOrProtein = [ entityIDToEntity[entityID] for entityID in cr.entityIDs ]
		keywordLoc,geneOrProteinLoc = [ entityIDToTokenLocs[entityID] for entityID in cr.entityIDs ]

		docs = list(set([ entityIDToDoc[entityID] for entityID in cr.entityIDs ]))
		assert len(docs) == 1
		doc = docs[0]

		if len(geneOrProtein.text) < 3:
			continue

		#print('keywordLoc',keywordLoc)
		#print('geneOrProteinLoc',geneOrProteinLoc)
		distBetweenEntities = min( abs(min(keywordLoc)-max(geneOrProteinLoc)), abs(min(geneOrProteinLoc)-max(keywordLoc)) )
		#print('distBetweenEntities',distBetweenEntities)
		#assert False
		if distBetweenEntities > 10:
			continue

		candidateRelations.append(cr)
		metadata.append(doc.metadata)
	print("%s : candidate filtering (%d)" % (now(),len(candidateRelations)))

	vectorizer = kindred.Vectorizer(entityCount=2)
	vectors = vectorizer.fit_transform(candidateRelations)
	print("%s : vectorizer" % now())

	with open(args.outRelationsVectorsAndMetadata,'wb') as outF:
		pickle.dump((candidateRelations,vectors,metadata),outF)
	print("%s : saved" % now())
