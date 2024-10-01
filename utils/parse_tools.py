__all__ = ["clean_variant_id", "is_synonymous", "aggregate_scores"]
from Bio.PDB.Polypeptide import protein_letters_3to1
import numpy as np
import pandas as pd

protein_letters_1to3 = {v:k.title() for k,v in protein_letters_3to1.items()}
protein_letters_1to3['Z'] = 'Glx'
protein_letters_1to3['X'] = 'Xaa'
protein_letters_1to3['B'] = 'Asx'
protein_letters_1to3['J'] = 'Xle'
protein_letters_1to3['*'] = "Ter"

def clean_variant_id(s : str) -> str:
    """
    s : one letter amino acid substitution string

    Returns a cleaned variant id in the format of p.Leu123Ser
    """
    if not len(s):
        return ""
    hgvs_pro = protein_letters_1to3[s[0]] + s[1:-1]
    if s[-1] == "*":
        hgvs_pro += "Ter"
    elif s[-1] == "~":
        hgvs_pro += "del"
    else:
        hgvs_pro += protein_letters_1to3[s[-1]]
    return "p."+ hgvs_pro

def is_synonymous(hgvs_pro : str) -> bool:
    """
    hgvs_pro : str

    Returns True if the variant is synonymous
    """
    try:
        if not len(hgvs_pro):
            return False
    except TypeError:
        print(hgvs_pro)
        raise TypeError()
    if hgvs_pro[-1] == "=":
        return True
    if hgvs_pro[2:5] == hgvs_pro[-3:]:
        return True
    return False

def aggregate_scores(row, score_columns,score_agg_method=None):
    """
    row : pandas.Series
    score_columns : list of str

    Returns the average of the scores in the score_columns
    """
    scores = np.array(pd.to_numeric(row[score_columns],errors='coerce')).reshape(-1)
    scores = scores[~np.isnan(scores)]
    if score_agg_method is None:
        return scores
    if score_agg_method == "median":
        return [np.median(scores),]

def format_GRCh38_coordinates(scoreset,scoreset_config):
    """
    scoreset : pandas.DataFrame
    scoreset_config : dict

    Returns a DataFrame with the GRCh38 coordinates
    """
    scoreset = scoreset.assign(CHROM=str(scoreset_config['CHROM']),
                               POS=scoreset.loc[:,scoreset_config['pos_col']].astype(str),
                               REF=scoreset.loc[:,scoreset_config['ref_col']],
                               ALT=scoreset.loc[:,scoreset_config['alt_col']])
    return scoreset

