from keras.layers import *
from keras.models import *
from keras import backend as K
from keras.engine.topology import Layer
from keras import initializers


MAX_LEN_cdr = 24
NB_WORDS_cdr = 8001  # kmer1 -> 21 , kmer2 -> 401, kmer3 -> 8001
NB_cdr_number_ids = 3

MAX_LEN_ag = 2371
NB_WORDS_ag = 21  # kmer1 -> 21 , kmer2 -> 401, kmer3 -> 8001

EMBEDDING_DIM = 100

filters = 256

cdr_kernel_size = 6
cdr_pool_size = cdr_strides = 4

ag_kernel_size = 60
ag_pool_size = ag_strides = 20

lstm_size = 50

att_size = 50

dt_ratio = 0.5

class AttLayer(Layer):
    def __init__(self, attention_dim):
        self.init = initializers.RandomNormal(seed=10)
        self.supports_masking = True
        self.attention_dim = attention_dim
        super(AttLayer, self).__init__()

    def build(self, input_shape):
        assert len(input_shape) == 3
        self.W = K.variable(self.init((input_shape[-1], self.attention_dim)))
        self.b = K.variable(self.init((self.attention_dim, )))
        self.u = K.variable(self.init((self.attention_dim, 1)))
        self.trainable_weights = [self.W, self.b, self.u]
        super(AttLayer, self).build(input_shape)

    def compute_mask(self, inputs, mask=None):
        return mask

    def call(self, x, mask=None):
        # size of x :[batch_size, sel_len, attention_dim]
        # size of u :[batch_size, attention_dim]
        # uit = tanh(xW+b)
        uit = K.tanh(K.bias_add(K.dot(x, self.W), self.b))
        ait = K.dot(uit, self.u)
        ait = K.squeeze(ait, -1)

        ait = K.exp(ait)

        if mask is not None:
            # Cast the mask to floatX to avoid float64 upcasting in theano
            ait *= K.cast(mask, K.floatx())
        ait /= K.cast(K.sum(ait, axis=1, keepdims=True) + K.epsilon(), K.floatx())
        ait = K.expand_dims(ait)
        weighted_input = x * ait
        output = K.sum(weighted_input, axis=1)

        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[-1])


def get_model():
    cdrs_ids = Input(shape=(MAX_LEN_cdr,))
    cdrs_number_ids = Input(shape=(MAX_LEN_cdr,))

    ags_ids = Input(shape=(MAX_LEN_ag,))

    emb_cdr_ids = Embedding(NB_WORDS_cdr, EMBEDDING_DIM, trainable=True)(cdrs_ids)
    emb_cdr_number_ids = Embedding(NB_cdr_number_ids, EMBEDDING_DIM, trainable=True)(cdrs_number_ids)

    emb_cdr = Add()([emb_cdr_ids, emb_cdr_number_ids ])
    emb_cdr_bn = BatchNormalization()(emb_cdr)
    emb_cdr_dt = Dropout(dt_ratio)(emb_cdr_bn)

    emb_ag_ids = Embedding(NB_WORDS_ag, EMBEDDING_DIM, trainable=True)(ags_ids)
    emb_ag_bn = BatchNormalization()(emb_ag_ids)
    emb_ag_dt = Dropout(dt_ratio)(emb_ag_bn)

    cdr_conv_layer = Conv1D(filters = filters, kernel_size = cdr_kernel_size,padding = "valid",activation='relu')(emb_cdr_dt)
    cdr_max_pool_layer = MaxPooling1D(pool_size = cdr_pool_size, strides = cdr_strides)(cdr_conv_layer)

    ag_conv_layer = Conv1D(filters = filters, kernel_size = ag_kernel_size,padding = "valid",activation='relu')(emb_ag_dt)
    ag_max_pool_layer = MaxPooling1D(pool_size = ag_pool_size, strides = ag_strides)(ag_conv_layer)

    merge_layer=Concatenate(axis=1)([cdr_max_pool_layer, ag_max_pool_layer])
    bn=BatchNormalization()(merge_layer)
    dt=Dropout(dt_ratio)(bn)

    l_lstm = Bidirectional(LSTM(lstm_size, return_sequences=True))(dt)
    l_att = AttLayer(att_size)(l_lstm)

    preds = Dense(1, activation='sigmoid')(l_att)

    model = Model(inputs=[cdrs_ids, cdrs_number_ids, ags_ids],outputs= [preds])

    model.compile(loss='binary_crossentropy',optimizer='adam')

    return model