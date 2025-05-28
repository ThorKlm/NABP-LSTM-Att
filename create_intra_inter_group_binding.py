from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer
import pandas as pd
import os
import pickle
from itertools import combinations, product
from tqdm import tqdm

# === PARAMETERS ===
input_dir = "data/nano"
output_csv = "clusters.csv"
similarity_threshold = 0.9  # Equivalent to CD-HIT 90% identity

# === Step 1: Read and Merge All CSV Files ===
print("Reading CSV files...")
all_data = []
for file in os.listdir(input_dir):
    if file.endswith(".csv"):
        df = pd.read_csv(os.path.join(input_dir, file))
        all_data.append(df)

full_df = pd.concat(all_data, ignore_index=True)
print(f"Loaded {len(full_df)} entries.")

# === Step 2: Merge Antigen Sequences ===
def merge_sequences(row):
    return ''.join([
        str(row.get('Antigen sequence_1', '') or ''),
        str(row.get('Antigen sequence_2', '') or ''),
        str(row.get('Antigen sequence_3', '') or '')
    ])

full_df["merged_antigen_seq"] = full_df.apply(merge_sequences, axis=1)

# Drop rows with empty antigen sequences
full_df = full_df[full_df["merged_antigen_seq"].str.strip() != ""]

# === Step 3: Sequence Vectorization (k-mer based) ===
print("Vectorizing antigen sequences using 3-mers...")
vectorizer = CountVectorizer(analyzer='char', ngram_range=(3, 3))
X = vectorizer.fit_transform(full_df["merged_antigen_seq"])

# === Step 4: Similarity Matrix ===
print("Calculating similarity matrix...")
similarity_matrix = cosine_similarity(X)

# === Step 5: Clustering ===
print("Clustering sequences...")
# Convert threshold into distance
distance_matrix = 1 - similarity_matrix
clustering = AgglomerativeClustering(
    affinity="precomputed",
    linkage="average",
    distance_threshold=1 - similarity_threshold,
    n_clusters=None
)
clusters = clustering.fit_predict(distance_matrix)

# === Step 6: Save Clusters ===
print("Saving clusters...")
full_df["cluster"] = clusters
full_df[["pdb", "cluster"]].to_csv(output_csv, index=False)
print(f"Clusters saved to {output_csv}")

clusters_df = pd.read_csv("clusters.csv")
input_dir = "data/nano"

# Merge all .csv data
print("Loading nanobody-antigen data...")
all_data = []
for file in os.listdir(input_dir):
    if file.endswith(".csv"):
        df = pd.read_csv(os.path.join(input_dir, file))
        all_data.append(df)
df_full = pd.concat(all_data, ignore_index=True)

# Add cluster info
df_full = df_full.merge(clusters_df, on="pdb", how="inner")
print(f"Entries with clusters: {len(df_full)}")

# === Generate intra-group (positive) and inter-group (negative) pairs ===
intra_pairs = []
inter_pairs = []

print("Generating pairs...")
grouped = df_full.groupby("cluster")

for cluster_id, group in tqdm(grouped):
    pdbs = group["pdb"].tolist()
    intra_pairs.extend(list(combinations(pdbs, 2)))  # All unique intra-cluster pairs

# Inter-group: use product of entries from distinct clusters
cluster_ids = sorted(df_full["cluster"].unique())
for i in tqdm(range(len(cluster_ids))):
    for j in range(i + 1, len(cluster_ids)):
        cluster_i = df_full[df_full["cluster"] == cluster_ids[i]]["pdb"]
        cluster_j = df_full[df_full["cluster"] == cluster_ids[j]]["pdb"]
        inter_pairs.extend(list(product(cluster_i, cluster_j)))

print(f"Intra-group pairs: {len(intra_pairs)}")
print(f"Inter-group pairs: {len(inter_pairs)}")

# Save as pickles
os.makedirs("data/asPICKLE", exist_ok=True)
with open("data/asPICKLE/intra_group_binding.pickle", "wb") as f:
    pickle.dump(intra_pairs, f)

with open("data/asPICKLE/inter_group_binding.pickle", "wb") as f:
    pickle.dump(inter_pairs, f)

print("Pickles saved to data/asPICKLE/")