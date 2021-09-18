#!/usr/bin/env python3

# unpackEF.py

# Based on
# parsefeaturejsons.py

# but with a basically different goal. This version
# of the script does not normalize features, and it
# writes the files to disk as a pseudo-text with
# words separated by spaces.

# it expects to have the following files living in the same
# folder:

# SonicScrewdriver.py
# CorrectionRules.txt
# VariantSpellings.txt

# Note that if you want to run this on your own machine,
# you will need to change the "rootpath," below --
# circa line 470 -- which is the
# folder where we expect to find Extracted Features living.

# Also, you'll need to create a metadata file that pairs
# docids with paths to the features.

# Example of USAGE:

# run parsefeaturejsons ids2pathlist.tsv

import csv, os, sys, bz2, random, json, math
from collections import Counter

import numpy as np
import pandas as pd

# import utils
import SonicScrewdriver as utils

# By default, this script does spelling normalization according to
# rules that update archaic spelling and normalize to British where
# practice differs across the Atlantic.

# It also corrects ocr errors.

translator = dict()

with open('CorrectionRules.txt', encoding = 'utf-8') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        if len(row) < 2:
            continue
        translator[row[0]] = row[1]

with open('VariantSpellings.txt', encoding = 'utf-8') as f:
    reader = csv.reader(f, delimiter = '\t')
    for row in reader:
        if len(row) < 2:
            continue
        translator[row[0]] = row[1]

stopwords = set()

allowedvocab = set()

with open('modelvocab.tsv', encoding = 'utf-8') as f:
    lines = f.readlines()
    for l in lines:
        l = l.strip()
        fields = l.split('\t')
        if int(fields[1]) < 30 or (len(fields[0]) < 2 and fields[0] != "i"):
            stopwords.add(fields[0])
        elif len(fields[0]) > 0 and fields[0][0].isdigit():
            stopwords.add(fields[0])
        elif '«' in fields[0] or '»' in fields[0]:
            continue
        else:
            allowedvocab.add(fields[0])

inlist = "‘’"
outlist = "''"
removallist = '“”"(){}[]<>_/■^*\t«»'
punctzapper = str.maketrans(inlist, outlist, removallist)

def normalize_token(token):
    ''' Normalizes a token by lowercasing it and and splitting
    on the hyphen.
    '''
    global stopwords, allowedvocab, translator, punctzapper

    token = token.lower()

    if token in stopwords:
        return []
    elif token in allowedvocab:
        return [token]

    token = token.translate(punctzapper)

    token = token.strip("'.*,-——")

    if len(token) < 2:
        return []
    elif token[0].isdigit() and token[-1].isdigit():
        return []

    if token in translator and translator[token] in allowedvocab:
        return [translator[token]]

    if 'ﬁ' in token or 'ﬂ' in token:
        token = token.replace('ﬁ', 'fi')
        token = token.replace('ﬂ', 'fl')
    elif '&quot;' in token:
        token = token.replace('&quot', '')
    elif token.endswith("'s"):
        token = token.replace("'s", "")


    if token in allowedvocab:
        return  [token]
    elif '-' in token:
        return [x for x in token.split('-') if x in allowedvocab]
        # I never want to treat hyphenated words as distinct
        # features; for modeling it's preferable to count the
        # pieces
    elif '—' in token:
        return [x for x in token.split('—') if  x in allowedvocab]
    elif '—' in token:
        return [x for x in token.split('—') if x in allowedvocab]
    elif '_' in token:
        return [x for x in token.split('_') if x in allowedvocab]
    elif ',' in token:
        return [x for x in token.split(',') if x in allowedvocab]
    else:
        return []

def add_feature(feature, count, features):
    ''' Adds feature-count to a dictionary
    if feature is already in the dictionary.
    '''

    if feature in features:
        features[feature] = count


class VolumeFromJson:

    # A data object that contains page-level wordcounts read from
    # json.

    # It normalizes wordcounts by lower-casing, and by folding certain
    # categories of word together; see normalize_token above.

    # It also includes functions that allow a volume to divide itself
    # according to instructions. E.g.: "volume, cut yourself into
    # three parts, after leaving out certain pages!"

    def __init__(self, volumepath, volumeid):
        '''Initializes a LoadedVolume by reading wordcounts from
        a json file. By default it reads all the pages. But if
        skip-front and skip-back are set to positive values,
        it will skip n pages.'''


        if volumepath.endswith('bz2'):
            with bz2.open(volumepath, mode = 'rt', encoding = 'utf-8') as f:
                thestring = f.read()
        else:
            with open(volumepath, encoding = 'utf-8') as f:
                thestring = f.read()

        thejson = json.loads(thestring)

        self.volumeid = thejson['id']

        pagedata = thejson['features']['pages']

        self.numpages = len(pagedata)
        self.pagecounts = []
        self.totalcounts = Counter()
        self.totaltokens = 0
        self.bodytokens = 0

        chunktokens = 0
        typesinthischunk = set()
        # a set of types in the current 10k-word chunk; progress
        # toward which is tracked by chunktokens

        self.integerless_pages = 0
        self.skipped_pages = 0
        compromise_pg = 0

        capitalizedbodytokens = 0

        for i in range(self.numpages):
            thispagecounts = Counter()
            thisbodytokens = 0
            thisheadertokens = 0
            thispage = pagedata[i]

            # There are really two ways of numbering pages. They come in an order,
            # which gives them an inherent ordinality (this is the *first* page). But
            # they also have cardinal *labels* attached, in the "seq" field. These labels
            # are usually, but not necessarily, convertible to integers. (Usually "00000001",
            # but could be "notes.") *Usually* they are == to the ordinal number,
            # but again, not necessarily! The world is full of fun edge cases!

            # In this loop, i is the ordinal page number, and cardinal_page is the cardinal
            # label; its value will be -1 if it can't be converted to an integer.

            # compromise_pg skips pages that have no integer seq, but otherwise
            # proceeds ordinally

            try:
                cardinal_page = int(thispage['seq'])
            except:
                cardinal_page = -1

            if cardinal_page > 0:
                compromise_pg += 1
            elif cardinal_page < 0:
                self.integerless_pages += 1

            if cardinal_page >= 0:

                bodywords = thispage['body']['tokenPosCount']
                for token, partsofspeech in bodywords.items():

                    normaltokenlist = normalize_token(token)

                    # Notice that we treat each word as a list, to permit
                    # counting both parts of a hyphenated word.
                    # But usually this will be a list of one.

                    for normaltoken in normaltokenlist:

                        for part, count in partsofspeech.items():
                            thisbodytokens += count
                            chunktokens += count
                            thispagecounts[normaltoken] += count

                headerwords = thispage['header']['tokenPosCount']
                for token, partsofspeech in headerwords.items():
                    normaltokenlist = normalize_token(token)

                    for normaltoken in normaltokenlist:
                        normaltoken = "#header" + normaltoken

                        for part, count in partsofspeech.items():
                            thisheadertokens += count
                            thispagecounts[normaltoken] += count

                # You will notice that I treat footers (mostly) as part of the body
                # Footers are rare, and rarely interesting.

                footerwords = thispage['footer']['tokenPosCount']
                for token, partsofspeech in footerwords.items():

                    normaltokenlist = normalize_token(token)

                    for normaltoken in normaltokenlist:

                        for part, count in partsofspeech.items():
                            thisbodytokens += count
                            chunktokens += count
                            thispagecounts[normaltoken] += count

                self.pagecounts.append(thispagecounts)

                for key, value in thispagecounts.items():
                    self.totalcounts[key] += value

                self.totaltokens += thisbodytokens
                self.totaltokens += thisheadertokens
                self.bodytokens += thisbodytokens

            else:
                # print(i, cardinal_page, compromise_pg)
                self.skipped_pages += 1

        # We are done with the __init__ method for this volume.

    def write_volume(self, outpath, override = False, translator = dict(), use_headers = False, skip_front = 0, skip_back = 0, docid = 'docid', tokensperpage = 186):

        global stopwords

        ''' This writes volume text after using a translation table to normalize,
        e.g., British or archaic spelling.

        It can be instructed to skip a certain fraction of pages at the front or back of the volume.
        '''

        startposition = int(skip_front * len(self.pagecounts))
        endposition = len(self.pagecounts) - int(skip_back * len(self.pagecounts))

        if startposition > endposition:
            print('Error in page trimming')
            sys.exit(0)

        pagedata = self.pagecounts[startposition: endposition]

        numpages = len(pagedata)
        numtokens = numpages * tokensperpage

        chunks = math.ceil(numtokens / 10000)

        if chunks > 40:
            chunks = 40

        pagesperchunk = int(numpages / chunks) + 1

        chunklist = []
        thischunk = []
        pagesinthischunk = 0

        metarows = []
        metarow = dict()

        for idx, page in enumerate(pagedata):
            for token, count in page.items():
                if token.startswith('#header') and not use_headers:
                    continue
                elif len(token) < 1:
                    continue

                thischunk.extend([token] * count)

            pagesinthischunk += 1

            if pagesinthischunk >= pagesperchunk:
                metarow['id'] = docid + '_' + str(len(chunklist))
                metarow['tokens'] = len(thischunk)
                metarow['skipped_pages'] = self.skipped_pages
                metarow['trimmed_pages'] = int(skip_front * len(self.pagecounts)) + int(skip_back * len(self.pagecounts))
                metarow['pagesinchunk'] = pagesinthischunk
                metarows.append(metarow)
                metarow = dict()
                chunklist.append(' '.join(thischunk))
                pagesinthischunk = 0
                thischunk = []

        metarow['id'] = docid + '_' + str(len(chunklist))
        metarow['tokens'] = len(thischunk)
        metarow['skipped_pages'] = self.skipped_pages
        metarow['trimmed_pages'] = int(skip_front * len(self.pagecounts)) + int(skip_back * len(self.pagecounts))
        metarow['pagesinchunk'] = pagesinthischunk
        metarows.append(metarow)
        chunklist.append(' '.join(thischunk))

        with open(outpath, mode = 'a', encoding = 'utf-8') as f:
            for idx, chunk in enumerate(chunklist):
                f.write(docid + '_' + str(idx) + ' lbl ' + chunk + '\n')

        return metarows

    def get_wordfreqs(self, translator = dict(), use_headers = False, skip_front = 0, skip_back = 0):

        global stopwords

        ''' This returns word frequencies for a volume.
        '''

        startposition = int(skip_front * len(self.pagecounts))
        endposition = len(self.pagecounts) - int(skip_back * len(self.pagecounts))

        if startposition > endposition:
            print('Error in page trimming')
            sys.exit(0)

        pagedata = self.pagecounts[startposition: endposition]

        tokenctr = Counter()

        for page in pagedata:
            for token, count in page.items():

                if token.startswith('#header') and not use_headers:
                    continue
                elif len(token) < 1:
                    continue
                elif len(token) < 2 and not token[0].isalpha():
                    continue
                elif token in translator:
                    token = translator[token]

                # if token in stopwords:
                    #continue

                tokenctr[token] += count

        return tokenctr

def write_doc_freqs(wordfreqs, docfreqs, authorsets):
    allwords = set()
    worddict = dict()
    for word, count in wordfreqs.items():
        lowword = word.lower()
        allwords.add(lowword)
        if lowword not in worddict:
            worddict[lowword] = set()
        worddict[lowword].add(word)

    with open('docfreqs.tsv', mode = 'w', encoding = 'utf-8') as f:
        f.write('word\tallcount\talldocs\tcapcount\tcapdocs\tburstiness\tcapratio\tcapburstiness\tauthorcount\n')
        for word in allwords:

            allcount = 0
            alldocs = 0
            capcount = 0
            capdocs = 0

            for spelling in worddict[word]:
                allcount += wordfreqs[spelling]
                alldocs += docfreqs[spelling]
                if len(spelling) > 1 and spelling[0].isupper() and spelling[1].islower():
                    capcount += wordfreqs[spelling]
                    capdocs += docfreqs[spelling]

            burstiness = allcount / (alldocs + .1)
            capratio = capcount / (allcount + .1)
            capburstiness = capcount / (capdocs + .1)
            authorcount = len(authorsets[word])

            row = [word, allcount, alldocs, capcount, capdocs, burstiness, capratio, capburstiness, authorcount]
            row = '\t'.join([str(x) for x in row]) + '\n'
            f.write(row)


if __name__ == "__main__":

    args = sys.argv


    # The root path where volumes are stored is hard-coded here:

    rootpath = '/Volumes/TARDIS/work/ef/fic/'
    outpath = '/Users/tunder/work/tomallet/cohort4.txt'

    # You will need to change that if you're using it on your own machine.


    if len(args) < 2:
        print('This script requires a command-line argument:')
        print('a metadata file to use.')
        sys.exit(0)

    else:
        missing = 0
        path_to_meta = args[1]
        pathpairs = []
        themissing = []

        meta = pd.read_csv(path_to_meta, dtype = 'object', sep = '\t')

        ctr = 0
        all_outrows = []


        for idx, row in meta.iterrows():

            ctr += 1

            if ctr % 100 == 1:
                print(ctr)

            docid = row['docid']
            path = row['path']
            tokensperpage = float(row['tokensperpage'])
            if pd.isnull(tokensperpage) or tokensperpage < 1:
                tokensperpage = 100
                print('tokensperpage error ', docid)

            inpath = rootpath + path

            if not os.path.isfile(inpath):
                missing += 1
                print('missing')

            else:
                vol = VolumeFromJson(inpath, docid)
                metarows = vol.write_volume(outpath, override = True, translator = translator, use_headers = False, skip_front = .15, skip_back = 0.05, docid = docid,
                    tokensperpage = tokensperpage)
                all_outrows.extend(metarows)

        print(missing)


        columns = ['id', 'tokens', 'skipped_pages', 'trimmed_pages', 'pagesinchunk']

        with open('parsing_metadata3.tsv', mode = 'w', encoding = 'utf-8') as f:
            scribe = csv.DictWriter(f, fieldnames = columns, delimiter = '\t')
            scribe.writeheader()

            for row in all_outrows:
                if 'pages' in row:
                    row['pagesinchunk'] = row.pop('pages')
                if 'pagesinthischunk' in row:
                    row['pagesinchunk'] = row.pop('pagesinthischunk')
                scribe.writerow(row)











