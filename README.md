# Tumor Antigens

This project used [Kindred](https://github.com/jakelever/kindred) to extract mentions of tumor antigens in PubMed abstracts and PubMed Central papers.

The relevant files are described below:

- [wordlistLoader.py](https://github.com/jakelever/tumorantigens/blob/master/wordlistLoader.py) - Load up the list of gene names and prepare it for quick searching against sentences
- [findSentences.py](https://github.com/jakelever/tumorantigens/blob/master/findSentences.py) - Find sentences that mention "tumor antigen" (plus more spellings) and a gene name
- [prepareForLearning.py](https://github.com/jakelever/tumorantigens/blob/master/prepareForLearning.py) - Parse the sentences and vectorize them using Kindred
- [activelyLearn.py](https://github.com/jakelever/tumorantigens/blob/master/activelyLearn.py) - Use active learning to find ambigious sentences and request annotation. Continously builds the final knowledgebase as more sentences are annotated.

