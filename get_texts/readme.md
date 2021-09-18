get_texts
==========

Once metadata is constructed, we need to download HathiTrust Extracted Features, unpack them, discard the first 15% and last 5% of pages, and divide them into roughly 10,000-word chunks.

This is done in ```unpackEF3.py```. The "3" indicates that it's the latest in a long line of similar scripts.

Part of the process is also defining a restricted vocabulary for the model. This can be found in ```modelvocab.tsv```. Each word is paired with the number of authors found using it. The first few words, with zeroes, are actually very common words that we excluded from the model.

Very rare words are excluded. The process behind that is documented (although I admit, not fully or clearly) in ```MakeVocab.ipynb```.