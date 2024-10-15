from utils.parse_tools import clean_variant_id, is_synonymous, aggregate_scores, format_GRCh38_coordinates
from utils.read_data import load_scoreset,detects_splicing,load_scoreset_config,load_score_range, format_hgvs_nuc
from utils.mapping_utils import merge_with_clinvar, merge_with_gnomAD, assign_sample_labels, assign_author_labels, merge_with_clinvar_by_hgvs_nuc,merge_with_gnomAD_by_hgvs_nuc, postprocess_score_range, merge_with_clinvar_by_GRCh38_coordinates, merge_with_gnomAD_by_GRCh38_coordinates,validate_hgvs_pro
from utils.external_data_utils import load_external_data
from pathlib import Path
import joblib

def assign_hgvs_pro(scoreset, scoreset_config):
    amino_acid_substitution_format = scoreset_config['amino_acid_substitution_format']
    amino_acid_substitution_column = scoreset_config['amino_acid_substitution_column']
    assert amino_acid_substitution_column in scoreset.columns
    [isinstance(x,str) for x in scoreset.loc[:,amino_acid_substitution_column]], scoreset.loc[:,amino_acid_substitution_column]
    if amino_acid_substitution_format == 'single_letter':
        scoreset = scoreset.assign(hgvs_pro=scoreset.loc[:,amino_acid_substitution_column].apply(clean_variant_id)) # hgvs_pro
    elif amino_acid_substitution_format == "p.single_letter":
        scoreset = scoreset.assign(hgvs_pro=scoreset.loc[:,amino_acid_substitution_column].apply(lambda s: clean_variant_id(s[2:]))) # hgvs_pro
    elif amino_acid_substitution_format == 'hgvs_pro_with_transcript':
        scoreset = scoreset.assign(hgvs_pro=scoreset.loc[:,amino_acid_substitution_column].str.split(":").str[1])
    else:
        assert amino_acid_substitution_format == 'hgvs_pro'
        scoreset = scoreset.assign(hgvs_pro=scoreset.loc[:,amino_acid_substitution_column])
        # remove any rows that aren't a valid protein substitution
        scoreset.loc[scoreset.hgvs_pro.str.slice(0,2) != "p.", 'hgvs_pro'] = ""
    scoreset = scoreset.assign(hgvs_pro=scoreset.hgvs_pro.fillna(""))
    return scoreset

def pipeline_HEAD(scoreset_name):
    scoreset_config = load_scoreset_config(scoreset_name)
    scoreset = load_scoreset(scoreset_name,scoreset_config=scoreset_config)
    scoreset = assign_hgvs_pro(scoreset, scoreset_config)
    assert 'hgvs_pro' in scoreset.columns
    assert [isinstance(x,str) for x in scoreset.hgvs_pro], scoreset.hgvs_pro
    scoreset = scoreset.assign(is_synonymous=scoreset.hgvs_pro.apply(is_synonymous)) # is_synonymous
    scoreset = scoreset.assign(scores = scoreset.apply(lambda row: aggregate_scores(row, scoreset_config['score_columns'],scoreset_config.get("score_agg_method",None)),axis=1)) # scores
    scoreset = scoreset[scoreset.scores.apply(len) > 0]
    assay_detects_splicing = detects_splicing(scoreset_name)
    cache_dir = Path(scoreset_config.get('cache_dir',"/data/dzeiberg/mave_calibration/cache/external_data"))
    cache_dir.mkdir(exist_ok=True,parents=True)
    uniprot_acc = scoreset_config['uniprot_acc']
    pkl_file = cache_dir / f"{uniprot_acc}.pkl"
    if pkl_file.exists():
        gene_info, gnomAD_df, clinvar_df = joblib.load(pkl_file)
    else:
        gene_info, gnomAD_df, clinvar_df = load_external_data(uniprot_acc)
        joblib.dump((gene_info, gnomAD_df, clinvar_df), pkl_file)
    return scoreset_config, scoreset,assay_detects_splicing, gene_info, gnomAD_df, clinvar_df


def pipeline_TAIL(scoreset, scoreset_name, assay_detects_splicing, scoreset_config):
    scoreset = validate_hgvs_pro(scoreset)
    scoreset = assign_sample_labels(scoreset,assay_detects_splicing)
    scores_flipped, lower_threshold, upper_threshold = load_score_range(scoreset_name)
    scoreset = assign_author_labels(scoreset, scores_flipped, lower_threshold, upper_threshold,scoreset_config)
    scoreset = postprocess_score_range(scoreset, scoreset_config)
    return scoreset

def run_pipeline_A(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df):
    scoreset = merge_with_clinvar(scoreset, clinvar_df)
    scoreset = merge_with_gnomAD(scoreset, gnomAD_df)
    scoreset = pipeline_TAIL(scoreset, scoreset_name, assay_detects_splicing, scoreset_config)
    return scoreset

def pipeline_A(scoreset_name):
    """
    Pipeline A: scoreset only has the variant in protein substitution format
    """
    scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df = pipeline_HEAD(scoreset_name)
    return run_pipeline_A(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df)

def run_pipeline_B(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df):
    scoreset = format_hgvs_nuc(scoreset, scoreset_config)
    scoreset = merge_with_clinvar_by_hgvs_nuc(scoreset, clinvar_df)
    scoreset = merge_with_gnomAD_by_hgvs_nuc(scoreset, gnomAD_df)
    scoreset = pipeline_TAIL(scoreset, scoreset_name, assay_detects_splicing, scoreset_config)
    return scoreset

def pipeline_B(scoreset_name):
    """
    Pipeline B: scoreset has the variant in nucleotide substitution format
    """
    scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df = pipeline_HEAD(scoreset_name)
    return run_pipeline_B(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df)

def run_pipeline_C(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df):
    scoreset = format_GRCh38_coordinates(scoreset, scoreset_config)
    scoreset = merge_with_clinvar_by_GRCh38_coordinates(scoreset, clinvar_df)
    scoreset = merge_with_gnomAD_by_GRCh38_coordinates(scoreset, gnomAD_df)
    scoreset = pipeline_TAIL(scoreset, scoreset_name, assay_detects_splicing, scoreset_config)
    return scoreset

def pipeline_C(scoreset_name):
    """
    Pipeline C: scoreset has the variant with GRCh38 coordinates
    """

    scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df = pipeline_HEAD(scoreset_name)
    return run_pipeline_C(scoreset_name, scoreset_config,scoreset,assay_detects_splicing,gene_info, gnomAD_df, clinvar_df)
