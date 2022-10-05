import argparse
import os
import re
import pandas as pd
import spacy
from spacy.lang.char_classes import ALPHA, ALPHA_LOWER, ALPHA_UPPER
from spacy.lang.char_classes import CONCAT_QUOTES, LIST_ELLIPSES, LIST_ICONS
from spacy.util import compile_infix_regex

# Modify tokenizer infix patterns
infixes = (
    LIST_ELLIPSES
    + LIST_ICONS
    + [
        r"(?<=[{a}\/\|])(?:[\(])(?=[{a}/\|])".format(a=ALPHA),
        r"(?<=[{a}])(?:[\/])(?=[{a}\(])".format(a=ALPHA)
    ]
)

infix_re = compile_infix_regex(infixes)
nlp = spacy.load("en_core_web_sm")
nlp.tokenizer.infix_finditer = infix_re.finditer


def read_a1(a1_file):
    a1_df=pd.read_csv(a1_file, sep='\t', header=None)
    a1_df.columns = ['TermID', 'EntityType', 'Entity']
    entities = a1_df[a1_df['TermID'].str.startswith('T')].copy()
    entities[['EntityType', 'Start End']] = entities.EntityType.str.split(' ', n=1, expand=True)
    continuous_entities = entities[~entities['Start End'].str.contains(';')].copy()
    continuous_entities[['Start', 'End']] = continuous_entities['Start End'].str.split(' ', expand=True)
    return continuous_entities


def read_txt(txt_file):
    with open(txt_file, 'r') as rd:
        content = rd.read()
    doc = nlp(content)
    sents = doc.sents
    start=0
    end=0
    i=0
    sent_dict = {}
    for s in sents:
        #print("Senetence" + str(i))
        idx = content.index(s.text)
        #print("idx:" + str(idx))
        if i==0:
            end=len(s.text)
        else:
            start = idx
            end = start + len(s.text)
        sent_dict[i] = {'Start': start, 'End': end, 'Text': s.text}
        i=i+1
    return sent_dict, doc


def check_entity_boundary(iob, token, entity_flag, entity_idx, continuous_entities, entity_map, pmid, recursion_flag):
    #print("Token:" + token.text + ", Entity idx: " + str(entity_idx) + ", recursion flag: " + str(recursion_flag))
    #print("token idx: " + str(token.idx) + ", strt: " + str(continuous_entities.iloc[entity_idx]['Start']) + ", End: " + str(continuous_entities.iloc[entity_idx]['End']))
    if token.idx>=continuous_entities.iloc[entity_idx]['Start']:
        if token.idx<continuous_entities.iloc[entity_idx]['End']:
            if entity_flag == 0:
                iob.write(token.text + "\t" + "B-" + entity_map[continuous_entities.iloc[entity_idx]['EntityType']] + "\n")
            else:
                iob.write(token.text + "\t" + "I-" + entity_map[continuous_entities.iloc[entity_idx]['EntityType']] + "\n")
            #print(t.text)
            entity_flag = 1
            if recursion_flag == 1:
                print("in recursion, Consecutive Entity: " + token.text)
            return entity_idx, entity_flag
        elif recursion_flag == 1:
            # This function was called in recurssion to check if there is any consecutive entities. So we looked ahead and could not find consecutive entities.
            # Hence, returning previous entity index
            return entity_idx-1, entity_flag
        else:
            print("Recursion called, current entity_idx: " + str(entity_idx))
            entity_flag = 0
            entity_idx_recurssion = entity_idx
            if entity_idx+1<continuous_entities.shape[0]:
                entity_idx_recurssion, entity_flag = check_entity_boundary(iob, token, entity_flag, entity_idx+1, continuous_entities, entity_map, pmid, 1)
            if entity_idx_recurssion == entity_idx:
                # this means there is no consecutive entities
                if token.text.strip()!="" and token.text.strip()!=pmid:
                    iob.write(token.text + "\tO\n")
                elif token.text.strip()=="":
                    iob.write("\n")
                if entity_flag==1:
                    entity_idx=entity_idx+1
                entity_flag = 0
            else:
                # there is a consecutive entity, the token has been written to the file, we now return incremented entity_idx and entity_flag
                return entity_idx_recurssion, entity_flag
    else:
        # We arrive here for the tokens that are before any entities.
        if token.text.strip()!="" and token.text.strip()!=pmid:
            iob.write(token.text + "\tO\n")
        elif token.text.strip()=="":
            iob.write("\n")
    return entity_idx, entity_flag


def convert_to_iob(iob_file, continuous_entities, sent_dict, doc):
    sent=0
    entity_idx = 0
    entity_flag = 0
    entity_map = {"Pharmacological_substance": "PS", "Disorder": "DS", "Subject": "S"}
    pmid = os.path.basename(iob_file).replace('.tsv','')
    #print("continuous_entities.shape:")
    #print(continuous_entities.shape[0])
    with open(iob_file, 'w') as iob:
        for t in doc:
            if entity_idx>=continuous_entities.shape[0]:
                if t.text.strip()!="" and t.text.strip()!=pmid:
                    iob.write(t.text + "\tO\n")
                elif t.text.strip()=="":
                    iob.write("\n")
            else:
                entity_idx, entity_flag = check_entity_boundary(iob,
                                                                t,
                                                                entity_flag,
                                                                entity_idx,
                                                                continuous_entities,
                                                                entity_map,
                                                                pmid,
                                                                0)
            if t.idx + len(t.text) == sent_dict[sent]['End']:
                iob.write("\n")
                sent = sent+1


def main_process(data_dir):
    a1_files = [f for f in os.listdir(data_dir) if f.endswith('.a1')]
    entity_types = ["Pharmacological_substance", "Disorder", "Subject"]
    for a1 in a1_files:
        only_entities = pd.DataFrame()
        pmid = os.path.basename(a1).replace('.txt', '')
        print("PMID:" + str(pmid))
        # print(a1)
        a1_file = os.path.join(data_dir, a1)
        a2_file = os.path.join(data_dir, a1.replace('a1', 'a2'))
        txt_file = os.path.join(data_dir, a1.replace('a1', 'txt'))
        iob_file = os.path.join(data_dir, a1.replace('a1', 'tsv'))
        if os.path.getsize(a1_file) > 0:
            continuous_entities = read_a1(a1_file)
            print(continuous_entities['EntityType'].unique())
            only_entities = continuous_entities[continuous_entities['EntityType'].isin(entity_types)].copy()
            only_entities[['Start', 'End']] = only_entities[['Start', 'End']].astype('int')
            only_entities = only_entities.sort_values(by=['Start'])
            # print(only_entities)
            # print(only_entities.shape)
        if os.path.getsize(txt_file) > 0:
            sent_dict, doc = read_txt(txt_file)
        convert_to_iob(iob_file, only_entities, sent_dict, doc)
        # print(sent_dict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='This script convert BioNLP 2013 annotation files and create IOB for \
    the NERs annotation')
    parser.add_argument("-f", "--folder", nargs=1, required=True, help="BioNLP annotation folder", metavar="PATH")
    args = parser.parse_args()
    # data_dir="/Users/saha/Workspace/PHAEDRA_corpus/dev/"
    main_process(args.folder[0])
    # Due to inconsistent line breaks in the original text files, the IOB file will may have multiple line breaks \
    # between sentences. From command line either run following command on individual files or merge training folder \
    # together and create train.tsv
    ## Folder wise data merge and replacing multiple line breaks with single like break
    # cat -s train/*.tsv > train.tsv
    ## Replacing multiple line breaks with single like break
    ## cat -s file_name.tsv > file_name.tsv