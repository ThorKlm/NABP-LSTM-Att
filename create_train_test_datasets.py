import pickle
import os
import glob
import pandas as pd
from sklearn.model_selection import train_test_split
import random
from tqdm import tqdm


random.seed(123)

def merge_cdrs(cdr1, cdr2, cdr3, linker="XXX"):
    return (cdr1 or "") + linker + (cdr2 or "") + linker + (cdr3 or "")

# Load all CSVs from data/nano into one DataFrame
nano_files = glob.glob("data/nano/*.csv")
df_list = [pd.read_csv(f) for f in nano_files]
all_data = pd.concat(df_list, ignore_index=True)
all_data.dropna(subset=["pdb"], inplace=True)
# Drop duplicates by keeping the first occurrence (or you can aggregate if needed)
all_data_unique = all_data.drop_duplicates(subset="pdb", keep="first")
all_data_unique.set_index("pdb", inplace=True)
record_dict = all_data_unique.to_dict(orient="index")

# Load cluster and binding info
with open('data/asPICKLE/clusters.pickle', 'rb') as binary_reader:
    clusters_raw = pickle.load(binary_reader)

with open('data/asPICKLE/intra_group_binding.pickle', 'rb') as binary_reader:
    intra_pairs = pickle.load(binary_reader)

with open('data/asPICKLE/inter_group_binding.pickle', 'rb') as binary_reader:
    inter_pairs = pickle.load(binary_reader)

# Resolve clusters to full 6-element data
clusters = []
for cluster in clusters_raw:
    resolved_cluster = []
    for pdb1 in cluster:
        try:
            r1 = record_dict[pdb1]
            item = (
                r1["Hchain_sequence"], r1["CDRH1"], r1["CDRH2"], r1["CDRH3"],
                r1["Antigen sequence_1"], pdb1
            )
            resolved_cluster.append(item)
        except KeyError:
            continue
    clusters.append(resolved_cluster)

# Resolve intra-group pairs
resolved_intra = []
for pdb1, pdb2 in intra_pairs:
    try:
        r1 = record_dict[pdb1]
        item = (
            r1["Hchain_sequence"], r1["CDRH1"], r1["CDRH2"], r1["CDRH3"],
            r1["Antigen sequence_1"], pdb1
        )
        resolved_intra.append(item)
    except KeyError:
        continue

# Resolve inter-group pairs
resolved_inter = []
for pdb1, pdb2 in inter_pairs:
    try:
        r1 = record_dict[pdb1]
        item = (
            r1["Hchain_sequence"], r1["CDRH1"], r1["CDRH2"], r1["CDRH3"],
            r1["Antigen sequence_1"], pdb1
        )
        resolved_inter.append(item)
    except KeyError:
        continue

# Split into train/test
train_pos, test_pos = [], []
for cluster in tqdm(clusters, desc="Splitting clusters"):
    if len(cluster) == 0:
        continue  # avoid empty clusters
    tr, te = train_test_split(cluster, test_size=0.2, random_state=123)
    train_pos.extend(tr)
    test_pos.extend(te)

tr, te = train_test_split(resolved_intra, test_size=0.2, random_state=123)
train_pos.extend(tr)
test_pos.extend(te)

# Sample negative pairs equal to total positives
neg_data = random.sample(resolved_inter, len(train_pos) + len(test_pos))

# Label the data
train_data_pos = [(i[0], i[1], i[2], i[3], i[4], 1) for i in train_pos]
test_data_pos  = [(i[0], i[1], i[2], i[3], i[4], 1) for i in test_pos]
train_data_neg, test_data_neg = train_test_split(neg_data, test_size=0.2, random_state=123)
train_data_neg = [(i[0], i[1], i[2], i[3], i[4], 0) for i in train_data_neg]
test_data_neg  = [(i[0], i[1], i[2], i[3], i[4], 0) for i in test_data_neg]

# Split validation from training
train_data_pos, val_data_pos = train_test_split(train_data_pos, test_size=0.05, random_state=123)
train_data_neg, val_data_neg = train_test_split(train_data_neg, test_size=0.05, random_state=123)

# Save labeled sets
os.makedirs("data/asPICKLE", exist_ok=True)
pickle.dump(train_data_pos, open("data/asPICKLE/train_data_pos.pickle", "wb"))
pickle.dump(val_data_pos, open("data/asPICKLE/val_data_pos.pickle", "wb"))
pickle.dump(test_data_pos, open("data/asPICKLE/test_data_pos.pickle", "wb"))

pickle.dump(train_data_neg, open("data/asPICKLE/train_data_neg.pickle", "wb"))
pickle.dump(val_data_neg, open("data/asPICKLE/val_data_neg.pickle", "wb"))
pickle.dump(test_data_neg, open("data/asPICKLE/test_data_neg.pickle", "wb"))

# Combine + shuffle
train_data_all = train_data_pos + train_data_neg
val_data_all = val_data_pos + val_data_neg
test_data_all = test_data_pos + test_data_neg

random.shuffle(train_data_all)
random.shuffle(val_data_all)
random.shuffle(test_data_all)

pickle.dump(train_data_all, open("data/asPICKLE/train_data_all.pickle", "wb"))
pickle.dump(val_data_all, open("data/asPICKLE/val_data_all.pickle", "wb"))
pickle.dump(test_data_all, open("data/asPICKLE/test_data_all.pickle", "wb"))

# Format for model input: (ID, merged_CDRs, Ag_seq, label)
def format_model_data(data):
    return [
        (item[0], merge_cdrs(item[1], item[2], item[3]), item[4], item[5])
        for item in data
    ]

train_CDR_antigen = format_model_data(train_data_all)
val_CDR_antigen   = format_model_data(val_data_all)
test_CDR_antigen  = format_model_data(test_data_all)

pickle.dump(train_CDR_antigen, open("data/asPICKLE/train_CDR123_antigen.pickle", "wb"))
pickle.dump(val_CDR_antigen, open("data/asPICKLE/val_CDR123_antigen.pickle", "wb"))
pickle.dump(test_CDR_antigen, open("data/asPICKLE/test_CDR123_antigen.pickle", "wb"))

print("Dataset creation complete")
print("Train:", len(train_CDR_antigen))
print("Val:", len(val_CDR_antigen))
print("Test:", len(test_CDR_antigen))
