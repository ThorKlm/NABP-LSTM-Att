from keras.preprocessing.sequence import pad_sequences
from model import get_model
from sklearn.metrics import roc_auc_score,average_precision_score
import pickle

cdr_kmer = 3
ag_kmer = 1

with open('data_prev/data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/cdr_features_te.pickle', 'rb') as binary_reader:
    cdr_features_te = pickle.load(binary_reader)
with open('data_prev/data_prev/features/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/ag_features_te.pickle', 'rb') as binary_reader:
    ag_features_te = pickle.load(binary_reader)


model=None
model=get_model()
directory_path = './model/cdr_kmer' + str(cdr_kmer) + '_ag_kmer' + str(ag_kmer) + '/'
print("directory_path", directory_path)
model.load_weights(directory_path + "Model99.h5")

# === Load vocab and tokenizer ===
cdr_kmer = 3
ag_kmer = 1

def load_vocab(vocab_file):
    with open(vocab_file, 'r') as f:
        return {line.strip(): idx for idx, line in enumerate(f)}

cdr_vocab_file = f'vocab/vocab_{cdr_kmer}kmer.txt'
ag_vocab_file = f'vocab/vocab_{ag_kmer}kmer.txt'

cdr_vocab = load_vocab(cdr_vocab_file)
ag_vocab = load_vocab(ag_vocab_file)

def sequence_to_kmers(seq, k):
    return [seq[i:i + k] for i in range(len(seq) - k + 1)]

def kmers_to_ids(kmers, vocab):
    return [vocab.get(kmer, 0) for kmer in kmers]  # 0 = unknown token

# === User-provided amino acid sequences ===
nanobody_samples = [["nbGFP_6xzf", "GFPVNRY", "SSAGDR", "NVGFEY", "QVQLVESGGALVQPGGSLRLSCAASGFPVNRYSMRWYRQADTNNDGWIEGDELKEREWVAGMSSAGDRSSYEDSVKGRFTISRDDARNTVYLQMNSLKPEDTAVYYCNVNVGFEYWGQGTQVTVSS"],
                     ["mCherry_8im0", "GRPFSEY", "RSSGT", "SRVDTDSPAFYDY", "QVQLVESGGGLVQAGGSLRLSCAVSGRPFSEYNLGWFRQAPGKEREFVARIRSSGTTVYTDSVKGRFSASRDNAKNMGYLQLNSLEPEDTAVYYCAMSRVDTDSPAFYDYWGQGTQVTVST"],
                     ["nbALB_8y9t", "GFTFSRY", "NSGGTY", "NSGDGKRYCSGGYCYRS", "EVQLQESGGGLVQPGGSLRLSCAASGFTFSRYWMFWVRQAPGKGLEWISDINSGGTYTRYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCATNSGDGKRYCSGGYCYRSRGQGTLVTVSS"],
                     ["nblys_1mel", "GYTIGPY", "NMGGGI", "DSTIYASYYECGHGLSTGGYGYDS", "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS"],
                     ["nbSARS_7f5h", "GSDFSSS", "SEGS", "VDRWYDY", "QVQLQESGGGLVQAGGSLRLSCAASGSDFSSSTMGWYRQAPGKQREFVAISSEGSTSYAGSVKGRFTISRDNAKNTVYLQMNSLEPEDTAVYYCNVVDRWYDYWGQGTQVTVSA"],
                     ["nbNAT_8zoy", "GFPVTNF", "YSTGI", "KDNGAWRQNYDY", "EVQLVESGGGLVQAGGSLRLSCAASGFPVTNFEMYWYRQAPGKEREWVAAIYSTGITEYADSVKGRFTISRDNSKNTVYLQMNSLKPEDTAVYYCNVKDNGAWRQNYDYWGQGTQVTVSS"]]


antigen_samples = [["aeGFP", "VSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"],
                   ["mCherry", "GMVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"],
                   ["albumin", "RGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"],
                   ["lysozyme", "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"],
                   ["sars-cov-2 rbd", "AGSPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTGTLEVLFQ"],
                   ["NAT", "APRDGDAQPRETWGKKIDFLLSVVGFAVDLANVWRFPYLCYKNGGGAFLIPYTLFLIIAGMPLFYMELALGQYNREGAATVWKICPFFKGVGYAVILIALYVGFYYNVIIAWSLYYLFSSFTLNLPWTDCGHTWNSPNCTDPKLLNGSVLGNHTKYSKYKFTPAAEFYERGVLHLHESSGIHDIGLPQWQLLLCLMVVVIVLYFSLWKGVKTSGKVVWITATLPYFVLFVLLVHGVTLPGASNGINAYLHIDFYRLKEATVWIDAATQIFFSLGAGFGVLIAFASYNKFDNNCYRDALLTSSINCITSFVSGFAIFSILGYMAHEHKVNIEDVATEGAGLVFILYPEAISTLSGSTFWAVVFFVMLLALGLDSSMGGMEAVITGLADDFQVLKRHRKLFTFGVTFSTFLLALFCITKGGIYVLTLLDTFAAGTSILFAVLMEAIGVSWFYGVDRFSNDIQQMMGFRPGLYWRLCWKFVSPAFLLFVVVVSIINFKPLTYDDYIFPPWANWVGWGIALSSMVLVPIYVIYKFLSTQGSLWERLAYGITPENEHHLVAQRDIRQFQLQHWLAI"]]

map = [[[] for _ in range(len(nanobody_samples))] for _ in range(len(antigen_samples))]
for CDR_index in [1, 2, 3]:
    for nb_index, nanobody_sample in enumerate(nanobody_samples):
        cdr_input_ids_list = []
        cdr_number_ids_list = []
        ag_input_ids_list = []

        for ag_index, antigen_sample in enumerate(antigen_samples):
            nanobody_seq_raw = nanobody_sample[CDR_index]
            antigen_seq_raw = antigen_sample[1]
            # === Tokenize and encode ===
            cdr_kmers = sequence_to_kmers(nanobody_seq_raw.lower(), cdr_kmer)
            ag_kmers = sequence_to_kmers(antigen_seq_raw.lower(), ag_kmer)

            cdr_input_ids = kmers_to_ids(cdr_kmers, cdr_vocab)
            ag_input_ids = kmers_to_ids(ag_kmers, ag_vocab)
            # print("cdr_input_ids", cdr_input_ids)
            # print("ag_input_ids", ag_input_ids)

            cdr_number_ids = [CDR_index-1] * len(cdr_input_ids)  # Use 0 for unknown CDR number, or customize as needed

            cdr_input_ids_list.append(cdr_input_ids)
            cdr_number_ids_list.append(cdr_number_ids)
            ag_input_ids_list.append(ag_input_ids)

        # === Pad sequences ===
        cdr_input_ids_padded = pad_sequences(cdr_input_ids_list, maxlen=24, padding='post', truncating='post')
        cdr_number_ids_padded = pad_sequences(cdr_number_ids_list, maxlen=24, padding='post', truncating='post')
        ag_input_ids_padded = pad_sequences(ag_input_ids_list, maxlen=2371, padding='post', truncating='post')

        # === Run prediction ===
        prediction = model.predict([cdr_input_ids_padded, cdr_number_ids_padded, ag_input_ids_padded])
        # print("prediction", prediction)
        for ag_index, antigen_sample in enumerate(antigen_samples):
            binding_prob = prediction[ag_index][0]
            # print("binding_prob", binding_prob, nanobody_sample[0], antigen_sample[0])
            map[nb_index][ag_index].append(binding_prob)

    import numpy as np
    import matplotlib.pyplot as plt

    # Print the results
    print(map)
    # Compute the average where data_prev is present, else np.nan
    avg_scores = np.array([[cell[CDR_index-1] if cell else np.nan for cell in row] for row in map])

    # Plot
    plt.figure(figsize=(8, 6))
    im = plt.imshow(avg_scores, cmap='viridis', interpolation='nearest')

    plt.colorbar(im, label="Average PPI Score (E)")
    plt.title("Nanobody–Antigen Interaction Heatmap, CDR:" + str(CDR_index))

    # Label axes
    plt.xlabel("Nanobody Index")
    plt.ylabel("Antigen Index")
    plt.xticks(ticks=np.arange(len(map[0])), labels=[nanobody_samples[i][0] for i in range(len(map[0]))])
    plt.xticks(rotation=45)
    plt.yticks(ticks=np.arange(len(map)), labels=[antigen_samples[i][0] for i in range(len(map))])

    # Annotate with values if desired
    for i in range(avg_scores.shape[0]):
        for j in range(avg_scores.shape[1]):
            val = avg_scores[i, j]
            if not np.isnan(val):
                plt.text(j, i, f"{val:.2f}", ha="center", va="center", color="white", fontsize=8)

    plt.tight_layout()
    plt.show()
    # map = np.array(map)
    # print("map", map)

import numpy as np
import matplotlib.pyplot as plt

# Print the results
print(map)
# Compute the average where data_prev is present, else np.nan
avg_scores = np.array([[np.mean(cell) if cell else np.nan for cell in row] for row in map])

# Plot
plt.figure(figsize=(8, 6))
im = plt.imshow(avg_scores, cmap='viridis', interpolation='nearest')

plt.colorbar(im, label="Average PPI Score (E)")
plt.title("Nanobody–Antigen Interaction Heatmap")

# Label axes
plt.xlabel("Nanobody Index")
plt.ylabel("Antigen Index")
plt.xticks(ticks=np.arange(len(map[0])), labels=[nanobody_samples[i][0] for i in range(len(map[0]))])
plt.xticks(rotation=45)
plt.yticks(ticks=np.arange(len(map)), labels=[antigen_samples[i][0] for i in range(len(map))])

# Annotate with values if desired
for i in range(avg_scores.shape[0]):
    for j in range(avg_scores.shape[1]):
        val = avg_scores[i, j]
        if not np.isnan(val):
            plt.text(j, i, f"{val:.2f}", ha="center", va="center", color="white", fontsize=8)

plt.tight_layout()
plt.show()