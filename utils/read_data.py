__all__ = ['load_config', 'load_scoreset', 'detects_splicing']

import yaml
import pandas as pd

def load_config(config_path):
    with open(config_path, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            return None

def load_scoreset_config(scoreset_name, config_filepath='/home/dzeiberg/mave_preprocessing/datasets.yaml'):
    config = load_config(config_filepath)
    return config[scoreset_name]

def load_scoreset(scoreset_name,
                    scoreset_config=None,
                    config_filepath='/home/dzeiberg/mave_preprocessing/datasets.yaml',
                    curation_filepath='/home/dzeiberg/mave_preprocessing/MAVE_Curation.xlsx'):
    if scoreset_config is None:
        config = load_config(config_filepath)
        scoreset_config = config[scoreset_name]
    header_row = scoreset_config.get('header_row', 0)
    author_transcript = scoreset_config.get('author_transcript', None)
    if author_transcript is None:
        curation_df = pd.read_excel(curation_filepath, sheet_name='MAVE Curation (new)',header=1)
        curation_df = curation_df.assign(curation_name=curation_df.loc[:,'Dataset Name'].str.strip()).set_index("curation_name")
        author_transcript = curation_df.loc[scoreset_config.get("curation_name"),'Transcript ID']
    if scoreset_config['extension'] == "xlsx":
        sheet_name = scoreset_config.get('sheet_name', 0)
        scoreset = pd.read_excel(scoreset_config['filepath'], sheet_name=sheet_name, header=header_row)
    elif scoreset_config['extension'] == "csv":
        scoreset = pd.read_csv(scoreset_config['filepath'], header=header_row)
    else:
        raise ValueError("Invalid extension")
    if "." in author_transcript:
        author_transcript = author_transcript.split(".")[0]
    scoreset = scoreset.assign(author_transcript=author_transcript)

    score_columns = scoreset_config.get('score_columns', 'score')
    if isinstance(score_columns, str):
        score_columns = [score_columns,]
    for score_column in score_columns:
        if score_column not in scoreset.columns:
            raise ValueError(f"Score column {score_column} not found in scoreset")

    amino_acid_substitution_column = scoreset_config['amino_acid_substitution_column']
    scoreset.loc[:,]
    scoreset = scoreset.assign(**{amino_acid_substitution_column: scoreset.loc[:,amino_acid_substitution_column].astype(str).str.replace("nan","")})
    return scoreset
    
def detects_splicing(scoreset_name,scoreset_config=None,config_filepath='/home/dzeiberg/mave_preprocessing/datasets.yaml',
                    curation_filepath='/home/dzeiberg/mave_preprocessing/MAVE_Curation.xlsx'):
    if scoreset_config is None:
        config = load_config(config_filepath)
        scoreset_config = config[scoreset_name].get('curation_name')
    curation_df = pd.read_excel(curation_filepath, sheet_name='MAVE Curation (new)',header=1)
    curation_df = curation_df.assign(curation_name=curation_df.loc[:,'Dataset Name'].str.strip()).set_index("curation_name")
    detects = curation_df.loc[scoreset_config,'Detects splicing variants?']
    if isinstance(detects, bool):
        return detects
    if detects == "TRUE":
        return True
    return False

def load_score_range(scoreset_name,scoreset_config=None,config_filepath='/home/dzeiberg/mave_preprocessing/datasets.yaml',
                    curation_filepath='/home/dzeiberg/mave_preprocessing/MAVE_Curation.xlsx'):
    """
    Load the curated score thresholds for a given scoreset

    Returns:
    -------
    flipped: bool
        Whether the score thresholds are flipped (i.e. higher scores are more likely to be pathogenic)
    lower_threshold: float
        The lower threshold for a score
        if flipped, this is the maximum value of the normal range, otherwise the maximum value of the abnormal range
    upper_threshold: float
        The upper threshold for a score
        if flipped, this is the minimum value of the abnormal range, otherwise the minimum value of the normal range
    """
    if scoreset_config is None:
        config = load_config(config_filepath)
        scoreset_config = config[scoreset_name]
        curated_name = scoreset_config.get('curation_name')
    curation_df = pd.read_excel(curation_filepath, sheet_name='MAVE Curation (new)',header=1)
    curation_df = curation_df.assign(curation_name=curation_df.loc[:,'Dataset Name'].str.strip()).set_index("curation_name")
    if curation_df.loc[curated_name,'Score ranges reported in text?'] == "Reported":
        normal_min = float(curation_df.loc[curated_name,"Normal min (published)"])
        normal_max = float(curation_df.loc[curated_name,"Normal max (published)"])
        abnormal_min = float(curation_df.loc[curated_name,"Abnormal min (published)"])
        abnormal_max = float(curation_df.loc[curated_name,"Abnormal max (published)"])
    elif scoreset_config.get('score_threshold_source',None) is not None:
        score_threshold_source = scoreset_config['score_threshold_source']
        normal_min = float(curation_df.loc[score_threshold_source,"Normal min (published)"])
        normal_max = float(curation_df.loc[score_threshold_source,"Normal max (published)"])
        abnormal_min = float(curation_df.loc[score_threshold_source,"Abnormal min (published)"])
        abnormal_max = float(curation_df.loc[score_threshold_source,"Abnormal max (published)"])
    else:
        return None, None, None
    if normal_max <= abnormal_min or abnormal_max > normal_min:
        flipped = True
        return flipped, normal_max, abnormal_min
    return False, abnormal_max, normal_min
        

def format_hgvs_nuc(scoreset,scoreset_config):
    nucleotide_substitution_column = scoreset_config.get('nucleotide_substitution_column',None)
    if nucleotide_substitution_column is None:
        return scoreset
    nucleotide_substitution_format = scoreset_config['nucleotide_substitution_format']
    if nucleotide_substitution_format == 'hgvs_nuc_with_transcript':
        scoreset = scoreset.assign(hgvs_nuc=scoreset.loc[:,nucleotide_substitution_column].str.split(":").str[1])
    elif nucleotide_substitution_format == 'Erwood_custom':
        scoreset = scoreset.assign(hgvs_nuc=scoreset.apply(lambda row: f"c.{row['Wild type Base']}{row.CDS}>{row['Edited Base']}", axis=1))
    else:
        assert nucleotide_substitution_format == 'hgvs_nuc'
        scoreset = scoreset.assign(hgvs_nuc=scoreset.loc[:,nucleotide_substitution_column])
    assert 'hgvs_nuc' in scoreset.columns
    return scoreset