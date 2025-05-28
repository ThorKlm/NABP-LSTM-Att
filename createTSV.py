import pickle
import random
import os
from tqdm import tqdm
random.seed(123)

CDR_kmer = 3
Ag_kmer = 3

def sequence2tokens(seq, kmer):
    tokens = []
    for i in range(len(seq) - (kmer-1)):
        tokens.append(seq[i:i + kmer])
    return tokens

def CDR_Ag_create_tsv(CDR_Ag_pairs, mode, CDR_kmer, Ag_kmer, CDR_TSV, Ag_TSV):
    import os
    os.makedirs(os.path.dirname(CDR_TSV), exist_ok=True)
    os.makedirs(os.path.dirname(Ag_TSV), exist_ok=True)
    f1 = open(CDR_TSV, "w")
    f2 = open(Ag_TSV, "w")

    for pair in tqdm(CDR_Ag_pairs):
        ID = pair[0]
        CDR_seq = pair[1]
        Ag_seq = pair[2]
        label = pair[3]

        label_str = str(label)
        CDR = str(CDR_seq).lower()
        Ag = str(Ag_seq).lower()

        CDR_tokens = " ".join(sequence2tokens(CDR, CDR_kmer))
        Ag_tokens = " ".join(sequence2tokens(Ag, Ag_kmer))

        f1.write(f"{mode}\t{ID}\t{label_str}\t{CDR_tokens}\n")
        f2.write(f"{mode}\t{ID}\t{label_str}\t{Ag_tokens}\n")

    f1.close()
    f2.close()

def CDR_Ag_createTrainValTestTSV():
    with open('data/asPICKLE/train_CDR123_antigen.pickle', 'rb') as f:
        train = pickle.load(f)
    with open('data/asPICKLE/val_CDR123_antigen.pickle', 'rb') as f:
        val = pickle.load(f)
    with open('data/asPICKLE/test_CDR123_antigen.pickle', 'rb') as f:
        test = pickle.load(f)

    base = f"data_prev/asTSV/cdr_kmer{CDR_kmer}_ag_kmer{Ag_kmer}/"
    CDR_Ag_create_tsv(train, "train", CDR_kmer, Ag_kmer, base + "CDR_tr.tsv", base + "Ag_tr.tsv")
    CDR_Ag_create_tsv(val, "val", CDR_kmer, Ag_kmer, base + "CDR_val.tsv", base + "Ag_val.tsv")
    CDR_Ag_create_tsv(test, "test", CDR_kmer, Ag_kmer, base + "CDR_te.tsv", base + "Ag_te.tsv")

if not os.path.exists('./data/asTSV'):
    os.makedirs('./data/asTSV')

print("started")
CDR_Ag_createTrainValTestTSV()
print("finshed")