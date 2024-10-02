__all__ = ["load_mapping_data","summarize_clinvar_group","get_clinvar_summaries","translate_refseq_to_ensembl","translate_ensembl_to_refseq","gather_gnomAD_info","merge_with_clinvar","merge_with_gnomAD","assign_sample_labels", "validate_hgvs_pro"]

import joblib
import pandas as pd
import tqdm
from typing import List
import numpy as np
import logging
from Bio.PDB.Polypeptide import protein_letters_3to1
protein_letters_1to3 = {v:k.title() for k,v in protein_letters_3to1.items()}
tqdm.tqdm.pandas()

def load_mapping_data(**kwargs):
    """
    Load mapping data from a pickle file

    Required Kwargs
    ----------
    mapping_filepath : str
        Path to the mapping data pickle file
        default : "/data/dzeiberg/mave_calibration/cache/mapping_data.pkl"

    Returns
    -------
    metadata : pd.DataFrame
        Dataset to uniprot_acc dataframe
    gene_info : pd.DataFrame
        uniprot_acc to gene info (transcripts, HGNC, symbols, etc.)
    gnomAD_df : pd.DataFrame
        gnomAD data frame
    spliceAI_df : pd.DataFrame
        spliceAI data frame
    clinvar_df : pd.DataFrame
        clinvar data frame
    """
    mapping_data = joblib.load(kwargs.get('mapping_filepath',"/data/dzeiberg/mave_calibration/cache/mapping_data.pkl"))
    metadata = mapping_data['metadata']
    gene_info = mapping_data['gene_info']
    gnomAD_df = mapping_data['gnomad_df']
    spliceAI_df = mapping_data['spliceAI_df']
    clinvar_df = mapping_data['clinvar_df']
    return metadata, gene_info, gnomAD_df, spliceAI_df, clinvar_df

def summarize_clinvar_group(grp):
    p_lp = grp.is_pathogenic.sum()
    b_lb = grp.is_benign.sum()
    conflicting = grp.is_conflicting.sum()#(grp.CLNSIG == "Conflicting_classifications_of_pathogenicity").sum()
    VUS = grp.is_VUS.sum()#(grp.CLNSIG == "Uncertain_significance").sum()
    alleleIDs = "|".join(grp.ID_x.astype(str).values)
    spliceAI_max = grp.spliceAI_score.max()
    return pd.Series(dict(num_p_lp=p_lp,num_b_lb=b_lb,num_conflicting=conflicting,num_VUS=VUS,clinvar_alleleIDs=alleleIDs,clinvar_records=len(grp),clinvar_spliceAI_max=spliceAI_max))

def get_clinvar_summaries(clinvar_df,refseq_transcript):
    """
    ClinVar summaries for a given refseq transcript

    Parameters
    ----------
    clinvar_df : pd.DataFrame
        ClinVar data frame
    
    refseq_transcript : str
        RefSeq transcript ID (e.g. NM_XXXXXX)
    
    Returns
    -------
        hgvs_pro_summaries : pd.DataFrame
            Summary of ClinVar data for each protein HGVS variant
            num_p_lp, num_b_lb, num_conflicting, num_VUS, clinvar_alleleIDs, clinvar_records, clinvar_spliceAI_max
    """
    assert refseq_transcript[:2] == "NM" and ("." not in refseq_transcript), "Refseq transcript must be in the form NM_XXXXXX"
    group_summaries = clinvar_df[(clinvar_df.transcript_base == refseq_transcript)].groupby("hgvs_pro").progress_apply(summarize_clinvar_group)
    hgvs_pro_summaries = pd.DataFrame(group_summaries)
    hgvs_pro_summaries = hgvs_pro_summaries[hgvs_pro_summaries.index.str.len() > 0]
    return hgvs_pro_summaries

def translate_refseq_to_ensembl(refseq_transcript,**kwarg):
    """
    Translate a RefSeq transcript to an Ensembl transcript using Ensembl v.112

    Parameters
    ----------
    refseq_transcript : str
        RefSeq transcript ID (e.g. NM_XXXXXX or NM_XXXXXX.X)
    
    Required kwargs
    ---------------
    transcript_table_filepath : str
        Path to the Ensembl transcript table
        default "/data/dzeiberg/mave_calibration/cache/transcript_mapping_table.tsv"

    Returns
    -------
    ensembl_transcripts : list
        List of Ensembl transcript Stable IDs
    """
    with_version = "." in refseq_transcript
    table = pd.read_csv(kwarg.get("transcript_table_filepath",
        "/data/dzeiberg/mave_calibration/transcript_mapping_table.tsv"),sep="\t")
    if not with_version:
        Ensembl_transcript_stable_ids = table.loc[table.display_label.str.contains(refseq_transcript+".",regex=False),'stable_id'].values
    else:
        Ensembl_transcript_stable_ids = table.loc[table.display_label == refseq_transcript,'stable_id'].values
    return Ensembl_transcript_stable_ids

def translate_ensembl_to_refseq(ensembl_transcript,**kwarg):
    """
    Translate an Ensembl transcript to a RefSeq transcript using Ensembl v.112

    Parameters
    ----------
    ensembl_transcript : str
        Ensembl transcript Stable ID (e.g. ENSTXXXXXX)

    Required kwargs
    ---------------
    transcript_table_filepath : str
        Path to the Ensembl transcript table
        default "/data/dzeiberg/mave_calibration/cache/transcript_mapping_table.tsv"
    
    Returns
    -------
    refseq_transcripts : list
        List of RefSeq transcript IDs
    """
    with_version = "." in ensembl_transcript
    table = pd.read_csv(kwarg.get("transcript_table_filepath",
        "/data/dzeiberg/mave_calibration/cache/transcript_mapping_table.tsv"),sep="\t")
    if not with_version:
        refseq = table.loc[table.stable_id == ensembl_transcript,'display_label']
    else:
        refseq = table.loc[table.stable_id == ensembl_transcript.split(".")[0],'display_label']
    return refseq.values


def merge_with_clinvar(scoreset,clinvar_df):
    """
    Merge scoreset with clinvar data

    Parameters
    ----------
    scoreset : pd.DataFrame
        Scoreset data frame
    
    clinvar_df : pd.DataFrame
        ClinVar data frame
    
    Returns
    -------
    scoreset : pd.DataFrame
        Scoreset data frame with clinvar data merged
    """
    clinvar_hgvs_pro_summaries = clinvar_df.groupby("hgvs_pro").apply(summarize_clinvar_group)
    scoreset_w_clinvar = pd.merge(scoreset.set_index("hgvs_pro"),clinvar_hgvs_pro_summaries,
                            left_index=True,
                            right_index=True,
                            how="left",
                            validate="one_to_one",
                            suffixes=("","_clinvar"))
    return scoreset_w_clinvar # clinvar_spliceAI_max, clinvar_records, num_p_lp, num_b_lb, num_conflicting, num_VUS, clinvar_alleleIDs

def merge_with_clinvar_by_hgvs_nuc(scoreset, clinvar_df):
    author_transcript = scoreset.author_transcript.iloc[0]
    if "NM_" in author_transcript:
        Ensembl_transcript_stable_ids = translate_refseq_to_ensembl(author_transcript)
    else:
        Ensembl_transcript_stable_ids = [author_transcript,]
    clinvar_hgvs_nuc_summaries = clinvar_df[(clinvar_df.Feature.isin(Ensembl_transcript_stable_ids)) & \
                                        (clinvar_df.hgvs_nuc.isin(scoreset.hgvs_nuc))].groupby("hgvs_nuc").apply(summarize_clinvar_group)
    scoreset_w_clinvar = pd.merge(scoreset.dropna(subset='hgvs_nuc').set_index("hgvs_nuc"),clinvar_hgvs_nuc_summaries,
                            left_index=True,
                            right_index=True,
                            how="left",
                            validate="one_to_one",
                            suffixes=("","_clinvar"))
    return scoreset_w_clinvar

def merge_with_clinvar_by_GRCh38_coordinates(scoreset, clinvar_df):
    scoreset = scoreset.reset_index().set_index(['CHROM','POS','REF','ALT'])
    clinvar_df.set_index(['CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    clinvar_summaries = clinvar_df[clinvar_df.index.isin(scoreset.index)].groupby(level=list(range(4)),).apply(summarize_clinvar_group)
    num_dups = scoreset.index.duplicated().sum()
    if num_dups > 0:
        scoreset = scoreset[~scoreset.index.duplicated(keep='first')]
        logging.warning(f"Removed {num_dups} duplicate GRCh38 coordinates, not expecting any; check input file")
    scoreset_w_clinvar = pd.merge(scoreset,clinvar_summaries,
                            left_index=True,
                            right_index=True,
                            how="left",
                            validate="one_to_one",
                            suffixes=("","_clinvar"))
    return scoreset_w_clinvar

def get_transcript_ids(scoreset):
    author_transcript = scoreset.author_transcript.iloc[0]
    if "NM_" in author_transcript:
        refseq_transcript = [author_transcript,]
        Ensembl_transcript_stable_ids = translate_refseq_to_ensembl(refseq_transcript[0])
    else:
        Ensembl_transcript_stable_ids = [author_transcript,]
        refseq_transcript = translate_ensembl_to_refseq(Ensembl_transcript_stable_ids)
    return refseq_transcript, Ensembl_transcript_stable_ids

def gather_gnomAD_info(gnomAD_df : pd.DataFrame,
                        Ensembl_transcript_stable_ids : List[str] = [],
                        RefSeq_transcript_ids : List[str] = [],):
    """
    Gather gnomAD information for a given Ensembl transcript

    Parameters
    ----------
    gnomAD_df : pd.DataFrame
        gnomAD data frame
    
    Ensembl_transcript_stable_id : List[str]
        Ensembl transcript Stable ID

    Returns
    -------
    gnomAD_info : pd.DataFrame
        gnomAD v4.1 genomes + exomes summarized at the amino acid substitution level
        gnomAD_variants_maxAC_AF, gnomAD_variants_max_spliceAI_score, gnomAD_variants_VCF_INFO
         - gnomAD_variants_maxAC_AF : allele frequency of the record with the highest allele count
         - gnomAD_variants_max_spliceAI_score : maximum spliceAI score
         - gnomAD_variants_VCF_INFO : GRCh38 genomic coordinates of the variants (e.g. "1:12345:A:T|2:34567:G:C")
    """
    Ensembl_transcript_stable_ids = [x.split(".")[0] for x in Ensembl_transcript_stable_ids]
    RefSeq_transcript_ids = [x.split(".")[0] for x in RefSeq_transcript_ids]
    gnomAD_variants = gnomAD_df[(gnomAD_df.Feature.isin(set(Ensembl_transcript_stable_ids))) | \
                                (gnomAD_df.Feature.isin(RefSeq_transcript_ids))].reset_index().groupby('hgvs_pro')
    gnomAD_variants_maxAC_AF = gnomAD_variants.apply(lambda grp: grp.loc[grp.AC.idxmax(),'AF'],include_groups=False)
    gnomAD_variants_max_spliceAI_score = gnomAD_variants.spliceAI_score.max()
    gnomAD_variants_VCF_INFO = gnomAD_variants.apply(lambda grp: "|".join(grp.loc[:,['CHROM',"POS",'REF','ALT']].astype(str).apply(":".join, axis=1).values))
    return pd.DataFrame(dict(AF=gnomAD_variants_maxAC_AF,
                            gnomAD_spliceAI_score=gnomAD_variants_max_spliceAI_score,
                            gnomAD_variants_VCF_INFO=gnomAD_variants_VCF_INFO))

def merge_with_gnomAD(scoreset, gnomAD_df):
    """
    Merge scoreset with gnomAD data

    Parameters
    ----------
    scoreset : pd.DataFrame
        Scoreset data frame
    
    gnomAD_df : pd.DataFrame
        gnomAD data frame
    
    refseq_transcript : str
        RefSeq transcript ID (e.g. NM_XXXXXX)
    
    Returns
    -------
    scoreset : pd.DataFrame
        Scoreset data frame with gnomAD data merged
    """
    
    refseq_transcript, Ensembl_transcript_stable_ids = get_transcript_ids(scoreset)
    gnomAD_info = gather_gnomAD_info(gnomAD_df,Ensembl_transcript_stable_ids,refseq_transcript)
    scoreset_processed = pd.merge(scoreset,gnomAD_info,
                                        left_index=True,right_index=True,how="left",validate="one_to_one").reset_index()
    return scoreset_processed

def merge_with_gnomAD_by_hgvs_nuc(scoreset, gnomAD_df):
    if scoreset.index.name != "hgvs_nuc":
        scoreset.set_index("hgvs_nuc",inplace=True)
    if gnomAD_df.index.name != "hgvs_nuc":
        gnomAD_df.set_index("hgvs_nuc",inplace=True)
    refseq_transcript, Ensembl_transcript_stable_ids = get_transcript_ids(scoreset)
    gnomAD_matches = gnomAD_df[(gnomAD_df.index.isin(scoreset.index)) & \
                            (gnomAD_df.Feature.isin(set(Ensembl_transcript_stable_ids).union(set(refseq_transcript))))].sort_values(by='AC',ascending=False)
    gnomAD_matches_unique = gnomAD_matches[~gnomAD_matches.index.duplicated(keep='first')]
    # do not use hgvs_pro column from gnomAD, prefer that from scoreset
    if "hgvs_pro" in scoreset.columns:
        gnomAD_matches_unique = gnomAD_matches_unique.drop("hgvs_pro",axis=1)
    scoreset_processed = pd.merge(scoreset,gnomAD_matches_unique,
                                        left_index=True,right_index=True,how="left",validate="one_to_one")
    scoreset_processed.rename(columns=dict(spliceAI_score='gnomAD_spliceAI_score'),inplace=True)
    return scoreset_processed
        
def merge_with_gnomAD_by_GRCh38_coordinates(scoreset, gnomAD_df):
    gnomAD_unique = gnomAD_df.reset_index().drop_duplicates(subset=['CHROM',"POS",'REF','ALT'],keep='first').set_index(['CHROM',"POS",'REF','ALT'])
    scoreset_w_gnomad = pd.merge(scoreset,gnomAD_unique,
                            left_index=True,
                            right_index=True,
                            how="left",
                            validate="one_to_one",
                            suffixes=("","_gnomAD"))
    return scoreset_w_gnomad

def validate_hgvs_pro(scoreset):
    if "hgvs_pro" not in scoreset.columns:
        scoreset = scoreset.assign(hgvs_pro=scoreset.apply(check_row_hgvs,axis=1))
    return scoreset

def check_row_hgvs(row):
    hgvs_pro = ""
    if "hgvs_pro" in row.index:
        hgvs_pro = row.hgvs_pro
    else:
        if "hgvs_pro_x" in row.index and not pd.isna(row.hgvs_pro_x) and "?" not in row.hgvs_pro_x:
            hgvs_pro = row.hgvs_pro_x
        if len(hgvs_pro) and "hgvs_pro_y" in row.index and not pd.isna(row.hgvs_pro_y) and "?" not in row.hgvs_pro_y:
            assert row.hgvs_pro_y == hgvs_pro
    return hgvs_pro

def assign_sample_labels(scoreset, detects_splicing):
    """
    Assign labels to the scoreset, used to determine which of the samples this instance will be used to calibrate (i.e. P/LP, B/LB, gnomAD, synonymous)

    Parameters
    ----------
    scoreset : pd.DataFrame
        Scoreset data frame
    
    detects_splicing : bool
        Whether the assay detects splicing
    
    Returns
    -------
    scoreset : pd.DataFrame
        Scoreset data frame with labels assigned
    """
    scoreset = scoreset.assign(labels = scoreset.progress_apply(lambda row: _assign_labels(row,detects_splicing),axis=1))
    return scoreset

def _assign_labels(row : pd.Series, detects_splicing : bool) -> List[str]:
    """
    Assign labels to the row 

    Parameters
    ----------
    row : pd.Series
        record
    
    Returns
    -------
    labels : list
        list of labels for the row subset of {'P/LP','B/LB','gnomAD' ,'synonymous','splicing'}
    """
    labels = []
    if (not detects_splicing) and (('clinvar_spliceAI_max' in row.index and row.clinvar_spliceAI_max > 0.5) or \
        ('gnomAD_spliceAI_score' in row.index and row.gnomAD_spliceAI_score > 0.5) or \
            ('spliceAI_score' in row.index and row.spliceAI_score > 0.5)):
        labels.append("splicing")
        return labels
    # hgvs_pro = ""
    # if "hgvs_pro" in row.index:
    #     hgvs_pro = row.hgvs_pro
    # else:
    #     if "hgvs_pro_x" in row.index and not pd.isna(row.hgvs_pro_x) and "?" not in row.hgvs_pro_x:
    #         hgvs_pro = row.hgvs_pro_x
    #     if len(hgvs_pro) and "hgvs_pro_y" in row.index and not pd.isna(row.hgvs_pro_y) and "?" not in row.hgvs_pro_y:
    #         assert row.hgvs_pro_y == hgvs_pro
    hgvs_pro = row.hgvs_pro
    
    if hgvs_pro[-3:] == "Ter":
        labels.append("nonsense")
        return labels
    if hgvs_pro[-3:] == "del":
        labels.append("deletion")
        return labels
    if row.AF > 0:
        labels.append("gnomAD")
    if row.num_p_lp > 0:
        labels.append("P/LP")
    if row.num_b_lb > 0:
        labels.append("B/LB")
    if row.num_VUS > 0:
        labels.append("VUS")
    if row.num_conflicting > 0:
        labels.append("conflicting")
    if row.is_synonymous:
        labels.append("synonymous")
    else:
        labels.append("non-synonymous")
    
    return labels

def assign_author_labels(scoreset, scores_flipped, lower_threshold, upper_threshold,scoreset_config):
    if scoreset_config.get("author_label_column",None) is not None:
        return scoreset.assign(author_labels=scoreset.loc[:,scoreset_config['author_label_column']])
    elif scores_flipped is None:
        return scoreset.assign(author_labels='N/A')
    labels = []
    for score in list(map(np.mean, scoreset.scores)):
        if scores_flipped:
            if score > upper_threshold:
                labels.append("Functionally_Abnormal")
            elif score < lower_threshold:
                labels.append("Functionally_Normal")
            else:
                labels.append("Itermediate")
        else:
            if score < lower_threshold:
                labels.append("Functionally_Abnormal")
            elif score > upper_threshold:
                labels.append("Functionally_Normal")
            else:
                labels.append("Intermediate")
    scoreset = scoreset.assign(author_labels=labels)
    return scoreset

def postprocess_score_range(scoreset, scoreset_config):
    postprocess_steps = scoreset_config.get('postprocess_scores',[])
    for step in postprocess_steps: 
        if step == "logit":
            scoreset = scoreset.assign(scores_pre_logit=scoreset.scores)
            scoreset = scoreset.assign(scores=scoreset.scores.apply(lambda x: np.log(x) - np.log(1-x)))
        elif step == "normalize_benign":
            benign_scores = np.concatenate(scoreset.loc[scoreset.labels.apply(lambda x: "B/LB" in x),"scores"].values).reshape(-1,)
            benign_mean = np.mean(benign_scores)
            benign_std = np.std(benign_scores)
            scoreset = scoreset.assign(scores_pre_normalize=scoreset.scores)
            scoreset= scoreset.assign(scores=scoreset.scores.apply(lambda x: (x - benign_mean)/benign_std))
        else:
            raise NotImplementedError("Postprocessing method {} not implemented".format(scoreset_config['postprocess_scores']))
    return scoreset

