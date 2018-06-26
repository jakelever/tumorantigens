import sys
import kindred
import argparse
import os
import json
import time
import pickle
from sklearn.linear_model import LogisticRegression
import random
from collections import defaultdict
import numpy as np
import math

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

class RESPONSE:
	POSITIVE = 1
	NEGATIVE = 0
	ENTITYERROR = -1
		
	TABLE = {'y':POSITIVE,'n':NEGATIVE,'x':ENTITYERROR}

def promptRelation(cr):
	assert isinstance(cr, kindred.Relation)
	entityIDToEntity = { entity.entityID:entity for entity,tokenIndices in cr.sentence.entityAnnotations }

	keyword,geneOrProtein = [ entityIDToEntity[entityID] for entityID in cr.entityIDs ]

	keywordStart,keywordEnd = keyword.position[0]
	geneOrProteinStart,geneOrProteinEnd = geneOrProtein.position[0]

	charByChar = list(cr.sentence.text)
	charByChar[keywordStart] = bcolors.FAIL + charByChar[keywordStart]
	charByChar[keywordEnd-1] += bcolors.ENDC
	charByChar[geneOrProteinStart] = bcolors.HEADER + charByChar[geneOrProteinStart]
	charByChar[geneOrProteinEnd-1] += bcolors.ENDC
	
	sentence = "".join(charByChar)

	print()
	print('#'*30)
	print(sentence)

	response = None
	while not response in RESPONSE.TABLE:
		response = input('Antigen? Positive=y,Negative=n,EntityError=x: ')

	return RESPONSE.TABLE[response]
	
def classifyForProbs(responses, vectors):
	trainIndices = sorted(list(responses.keys()))
	#testIndices = [ i for i in range(vectors.shape[0]) if not i in trainIndices ]

	trainIndices = [ i for i in trainIndices if not responses[i] == RESPONSE.ENTITYERROR ]

	trainX = vectors[trainIndices,:]
	trainY = [ 1 if responses[i] == RESPONSE.POSITIVE else 0 for i in trainIndices ]
	#testX = vectors[testIndices,:]
	#print(trainX.shape,len(trainY))
	#print(testX.shape)

	clf = LogisticRegression(class_weight='balanced',random_state=1)
	clf.fit(trainX,trainY)
	assert list(clf.classes_) == [0,1]
	#probs = clf.predict_proba(testX)
	probs = clf.predict_proba(vectors)
	#print(probs.shape, clf.classes_)
	#probsWithIndex = { testIndex:probs[i,1] for i,testIndex in enumerate(testIndices) }
	return probs[:,1].tolist()

def bootstrapProbs(allResponses, vectors):
	bootstrapCount = 10
	bootstrapFraction = 0.8
	countRequirement = 5

	probs = defaultdict(list)
	for _ in range(bootstrapCount):
		#print("bootstrapping...")
		subsetSize = int(round(bootstrapFraction*len(allResponses)))
		subsetKeys = random.sample(list(allResponses.keys()), subsetSize)
		subset = { k:allResponses[k] for k in subsetKeys }
		tmpProbs = classifyForProbs(subset, vectors)
		for index,prob in enumerate(tmpProbs):
			probs[index].append(prob)

	#standardDevs = [ (np.std(vals),index) for index,vals in probs.items() if len(vals) >= countRequirement ]
	#standardDevs = sorted(standardDevs,reverse=True)
	#return standardDevs

	#scores = [ (abs(0.5 - (sum( 1 for v in vals if v > 0.5 ) / len(vals))),index) for index,vals in probs.items() if len(vals) >= countRequirement ]
	scores = [ (np.mean( [ abs(v-0.5) for v in vals ] ),index) for index,vals in probs.items() if len(vals) >= countRequirement ]

	ambig = [ 1 if (max(vals) > 0.5 and min(vals) < 0.5) else 0 for index,vals in probs.items() ]
	ambigCount = sum(ambig)
	ambigPerc = round(100*ambigCount/len(ambig),1)
	print()
	print("Ambig: %f%% (%d/%d)" % (ambigPerc,ambigCount,len(ambig)))

	#for (index,vals),score in zip(probs.items(), scores):
	#	print(index,vals, score)
	#	moo = input('x:')

	scores = sorted(scores)
	return scores

oldPredictionProbs = None
def buildKnowledgebase(allResponses, vectors, candidateRelations, metadata, ID2Term, outFile='antigens.tsv', threshold=0.5):
	global oldPredictionProbs
	print()

	predictionProbs = classifyForProbs(allResponses, vectors)
	if not oldPredictionProbs is None:
		norm = [ (x-y)**2 for x,y in zip(predictionProbs,oldPredictionProbs) ]
		change = math.sqrt(sum(norm))
		print ("Prob change: %f" % change)
	oldPredictionProbs = predictionProbs

	known = [ i for i,response in allResponses.items() if response == RESPONSE.POSITIVE ]
	predicted = [ i for i,prob in enumerate(predictionProbs) if prob >= threshold ]
	combined = sorted(list(set(known + predicted)))

	combined = [ i for i in combined if not (i in allResponses and allResponses[i] == RESPONSE.ENTITYERROR) ]

	with open(outFile,'w') as outF:
		headers = ['HUGOorUniprotID','normalizedTerm','termInText','pmid','journal','year','title','probability','sentence']
		outF.write("\t".join(headers) + "\n")
		for i in combined:
			cr = candidateRelations[i]
			m = metadata[i]
			prob = predictionProbs[i]

			entityIDToEntity = { entity.entityID:entity for entity,tokenIndices in cr.sentence.entityAnnotations }
			keyword,geneOrProtein = [ entityIDToEntity[entityID] for entityID in cr.entityIDs ]

			geneOrProteinID = geneOrProtein.externalID
			if geneOrProteinID in ID2Term:
				isGeneOrProtein,normalized = ID2Term[geneOrProteinID]
			else:
				isGeneOrProtein,normalized = "ambigious","ambigious"

			outData = [geneOrProtein.externalID,normalized,geneOrProtein.text,m['pmid'],m['journal'],m['year'],m['title'],prob,cr.sentence.text]
			outF.write("\t".join(map(str,outData)) + "\n")


	print("Written %d lines to %s" % (len(combined),outFile))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Actively learn a model')
	parser.add_argument('--relationsVectorsAndMetadata',required=True,type=str,help='Output stuff for quick loading')
	parser.add_argument('--genes',required=True,type=str,help='Gene wordlist')
	parser.add_argument('--proteins',required=True,type=str,help='Protein wordlist')
	args = parser.parse_args()

	ID2Term = {}
	with open(args.genes) as f:
		for line in f:
			termid,singleterm,terms,entrezid = line.strip('\n').split('\t')
			assert not termid in ID2Term
			ID2Term[termid] = ('gene',singleterm)
	with open(args.proteins) as f:
		for line in f:
			termid,singleterm,terms = line.strip('\n').split('\t')
			assert not termid in ID2Term
			ID2Term[termid] = ('protein',singleterm)

	with open(args.relationsVectorsAndMetadata,'rb') as f:
		candidateRelations,vectors,metadata = pickle.load(f)

	vectors = vectors.tocsr()
	print(vectors.shape)

	allResponses = {}
	if os.path.isfile('allResponses.pickle'):
		with open('allResponses.pickle','rb') as f:
			tmpVectors,allResponses = pickle.load(f)
			#assert (tmpVectors.nonzero() == vectors.nonzero()).all()
			assert tmpVectors.shape == vectors.shape
			assert np.sum(tmpVectors) == np.sum(vectors)

	#buildKnowledgebase(allResponses, vectors, candidateRelations, metadata, ID2Term, outFile='antigen.tsv')
	#sys.exit(0)

	while True:
		i = random.randint(0,len(candidateRelations)-1)
		cr = candidateRelations[i]
		allResponses[i] = promptRelation(cr)


		if len(allResponses) > 10:
			uncertainty = bootstrapProbs(allResponses,vectors)
			#print(uncertainty[:10])
			for std_dev,i in uncertainty[:1]:
				cr = candidateRelations[i]
				allResponses[i] = promptRelation(cr)
				#allResponses[i] = random.choice([RESPONSE.POSITIVE,RESPONSE.NEGATIVE])


			buildKnowledgebase(allResponses, vectors, candidateRelations, metadata, ID2Term, outFile='antigen.tsv')
			with open('allResponses.pickle','wb') as f:
				pickle.dump((vectors,allResponses),f)


