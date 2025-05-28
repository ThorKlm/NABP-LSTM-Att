import pickle
from features import CDR_Ag_Processor
import tokenization
from features import convert_examples_to_features
import os
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)


def get_Feature(tsv_file_path, kmer, max_seq_length, chunk_size=100000, prefix="chunk", output_dir="./data/features", skip_existing=True):
    vocab_file = f"vocab/vocab_{kmer}kmer.txt"
    tokenizer = tokenization.FullTokenizer(vocab_file=vocab_file, do_lower_case=True)
    processor = CDR_Ag_Processor()
    label_list = processor.get_labels()
    cdr_number_list = ["1", "2", "3"]

    lines = processor._read_tsv(tsv_file_path)
    total = len(lines)

    for i in range(0, total, chunk_size):
        chunk_file = os.path.join(output_dir, f"{prefix}_features_{i}_{min(i+chunk_size, total)}.pickle")
        if skip_existing and os.path.exists(chunk_file):
            print(f"Skipping existing chunk: {chunk_file}")
            continue

        examples = processor.get_examples(tsv_file_path, start=i, end=min(i + chunk_size, total))
        print(f"Processing chunk {i}-{min(i+chunk_size, total)}")
        chunk_features = convert_examples_to_features(examples, label_list, cdr_number_list, max_seq_length, tokenizer)

        with open(chunk_file, 'wb') as writer:
            pickle.dump(chunk_features, writer)

        del chunk_features  # Free memory
        tf.keras.backend.clear_session()  # Optional but helpful



# the max len of the CDRs is 24, so when k = 1 the max tokens is 24
# the max len of the AGs is 2373 , so when k = 3 the max tokens is (2373 - (3-1)) = 2371
cdr_kmer = 3
ag_kmer = 3
if not os.path.exists('./data/features'):
    os.makedirs('./data/features')

tsv_path = "data_prev/asTSV/cdr_kmer" + str(cdr_kmer) + "_ag_kmer" + str(ag_kmer) + "/"
ag_max_seq_length = 2371
cdr_max_seq_length = 24
cdr_tsv_file_name = "CDR_te.tsv"
ag_tsv_file_name = "Ag_te.tsv"
cdr_tsv_file_path = tsv_path + cdr_tsv_file_name
ag_tsv_file_path = tsv_path + ag_tsv_file_name
'''
cdr_features = get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length)
ag_features = get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/cdr_features_te.pickle', 'wb') as binary_writer:
    pickle.dump(cdr_features, binary_writer)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/ag_features_te.pickle', 'wb') as binary_writer:
    pickle.dump(ag_features, binary_writer)
'''
tsv_path = "data_prev/asTSV/cdr_kmer" + str(cdr_kmer) + "_ag_kmer" + str(ag_kmer) + "/"
output_dir = f"data_prev/features/cdr_kmer{cdr_kmer}_ag_kmer{ag_kmer}"
os.makedirs(output_dir, exist_ok=True)
cdr_tsv_file_path = tsv_path + "CDR_te.tsv"
ag_tsv_file_path = tsv_path + "Ag_te.tsv"
get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length, chunk_size=100000,
            prefix="cdr_te", output_dir=output_dir)
get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length, chunk_size=100000,
            prefix="ag_te", output_dir=output_dir)
##################################################

'''
cdr_tsv_file_name = "CDR_tr.tsv"
ag_tsv_file_name = "Ag_tr.tsv"
cdr_tsv_file_path = tsv_path + cdr_tsv_file_name
ag_tsv_file_path = tsv_path + ag_tsv_file_name
cdr_features = get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length)
ag_features = get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/cdr_features_tr.pickle', 'wb') as binary_writer:
    pickle.dump(cdr_features, binary_writer)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/ag_features_tr.pickle', 'wb') as binary_writer:
    pickle.dump(ag_features, binary_writer)
##################################################
'''
cdr_tsv_file_path = tsv_path + "CDR_tr.tsv"
ag_tsv_file_path = tsv_path + "Ag_tr.tsv"
get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length, chunk_size=100000,
            prefix="cdr_tr", output_dir=output_dir)
get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length, chunk_size=100000,
            prefix="ag_tr", output_dir=output_dir)

'''
cdr_tsv_file_name = "CDR_val.tsv"
ag_tsv_file_name = "Ag_val.tsv"
cdr_tsv_file_path = tsv_path + cdr_tsv_file_name
ag_tsv_file_path = tsv_path + ag_tsv_file_name
cdr_features = get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length)
ag_features = get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/cdr_features_val.pickle', 'wb') as binary_writer:
    pickle.dump(cdr_features, binary_writer)
with open('data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/ag_features_val.pickle', 'wb') as binary_writer:
    pickle.dump(ag_features, binary_writer)
'''
cdr_tsv_file_path = tsv_path + "CDR_val.tsv"
ag_tsv_file_path = tsv_path + "Ag_val.tsv"

get_Feature(cdr_tsv_file_path, cdr_kmer, cdr_max_seq_length, chunk_size=100000,
            prefix="cdr_val", output_dir=output_dir)
get_Feature(ag_tsv_file_path, ag_kmer, ag_max_seq_length, chunk_size=100000,
            prefix="ag_val", output_dir=output_dir)
print("done with the main operation")

import os
import pickle

def load_all_features(directory, prefix):
    """Load and concatenate all pickle files starting with a given prefix."""
    all_features = []
    for file in sorted(os.listdir(directory)):
        if file.startswith(prefix) and file.endswith(".pickle"):
            path = os.path.join(directory, file)
            print(f"Loading {path}")
            with open(path, 'rb') as f:
                features = pickle.load(f)
                all_features.extend(features)
    return all_features

def merge_and_save_features(base_output_dir, prefixes):
    """
    Merges and saves all chunked features for the given prefixes.

    Args:
        base_output_dir: directory where chunked features are stored
        prefixes: list of prefixes like ['cdr_te', 'ag_te', ...]
    """
    for prefix in prefixes:
        print(f"\nProcessing split: {prefix}")
        all_features = load_all_features(base_output_dir, prefix)
        merged_path = os.path.join(base_output_dir, f"{prefix}_combined.pickle")
        print(f"Saving merged features to {merged_path}")
        with open(merged_path, 'wb') as f:
            pickle.dump(all_features, f)

# === CONFIGURATION ===
cdr_kmer = 3
ag_kmer = 3
base_output_dir = f'data_prev/features/cdr_kmer{cdr_kmer}_ag_kmer{ag_kmer}'

# Make sure this matches the prefixes used during chunked saving
prefixes = [
    'cdr_te', 'ag_te',
    'cdr_tr', 'ag_tr',
    'cdr_val', 'ag_val'
]

# Run the merge and save
merge_and_save_features(base_output_dir, prefixes)