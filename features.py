# coding=utf-8
# Copyright 2018 The Google AI Language Team Authors.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import csv
import tokenization
import tensorflow as tf
from tqdm import tqdm


class InputExample(object):
  def __init__(self, guid, text_a, cdr_number, label=None):
    """Constructs a InputExample.

    Args:
      guid: Unique id for the example.
      text_a: string. The untokenized text of the first sequence. For single
        sequence tasks, only this sequence must be specified.
      CDR_number: string. The untokenized text of the CDR number
      label: (Optional) string. The label of the example. This should be
        specified for train and dev examples, but not for test examples.
    """
    self.guid = guid
    self.text_a = text_a
    self.cdr_number = cdr_number
    self.label = label

class InputFeatures(object):
  """A single set of features of data_prev."""

  def __init__(self,
               input_ids,
               label_id,
               cdr_number_ids,
               is_real_example=True):
    self.input_ids = input_ids
    self.label_id = label_id
    self.cdr_number_ids = cdr_number_ids
    self.is_real_example = is_real_example

class DataProcessor(object):
  """Base class for data_prev converters for sequence classification data_prev sets."""

  def get_train_examples(self, data_dir):
    """Gets a collection of `InputExample`s for the train set."""
    raise NotImplementedError()

  def get_dev_examples(self, data_dir):
    """Gets a collection of `InputExample`s for the dev set."""
    raise NotImplementedError()

  def get_test_examples(self, data_dir):
    """Gets a collection of `InputExample`s for prediction."""
    raise NotImplementedError()

  def get_labels(self):
    """Gets the list of labels for this data_prev set."""
    raise NotImplementedError()

  @classmethod
  def _read_tsv(cls, input_file, quotechar=None):
    """Reads a tab separated value file."""
    with tf.gfile.Open(input_file, "r") as f:
      reader = csv.reader(f, delimiter="\t", quotechar=quotechar)
      lines = []
      for line in reader:
        lines.append(line)
      return lines


class CDR_Ag_Processor(DataProcessor):
  def get_examples(self, file_path, start=0, end=None):
    lines = self._read_tsv(file_path)
    lines = lines[start:end] if end is not None else lines[start:]
    examples = []
    for line in lines:
      guid = line[0] + "-" + line[1]
      label = line[2]
      text_a = line[3]
      examples.append(InputExample(guid=guid, text_a=text_a, cdr_number=None, label=label))
    return examples

  def get_labels(self):
    return ["0", "1"]

  def _create_examples(self, lines):
    """Creates examples for the training and dev sets."""
    examples = []
    for (i, line) in enumerate(lines):
      set_type = tokenization.convert_to_unicode(line[0])
      ID = tokenization.convert_to_unicode(line[1])
      guid = "%s-%s" % (set_type, ID)
      label = tokenization.convert_to_unicode(line[2])
      text_a = tokenization.convert_to_unicode(line[3])
      cdr_number = tokenization.convert_to_unicode(line[4])

      examples.append(
          InputExample(guid=guid, text_a=text_a, cdr_number=cdr_number, label=label))
    return examples

def convert_single_example(ex_index, example, label_list, max_seq_length, tokenizer):
  tokens = tokenizer.tokenize(example.text_a)
  tokens = tokens[:max_seq_length]
  input_ids = tokenizer.convert_tokens_to_ids(tokens)
  label_id = int(example.label)

  while len(input_ids) < max_seq_length:
    input_ids.append(0)

  return InputFeatures(input_ids=input_ids, label_id=label_id, cdr_number_ids=None)

def convert_examples_to_features(examples, label_list, cdr_number_list, max_seq_length,
                                 tokenizer):
  """Convert a set of `InputExample`s to a list of `InputFeatures`."""

  features = []
  for (ex_index, example) in tqdm(enumerate(examples), total = len(examples)):
    if ex_index % 10000 == 0:
      tf.logging.info("Writing example %d of %d" % (ex_index, len(examples)))

    feature = convert_single_example(ex_index, example, label_list,
                                     max_seq_length, tokenizer)
    features.append(feature)
  return features

def main(_):
  tf.logging.set_verbosity(tf.logging.INFO)

if __name__ == "__main__":
  tf.app.run()
