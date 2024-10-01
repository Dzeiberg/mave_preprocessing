import requests
import pandas as pd
from pathlib import Path
import urllib.request
import subprocess
from datetime import datetime
from pysam import VariantFile
import json
from io import StringIO

def get_gene_info(uniprot_acc, **kwargs):
    url = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_acc}?"
    payload = {}
    headers = {}

    response = requests.request("GET", url, headers=headers, data=payload)
    json = response.json()
    seq = json['sequence']['sequence']
    mane_record = [rec for rec in json['dbReferences'] if rec['type'] == "MANE-Select"][0]
    hgnc_id = [rec for rec in json['dbReferences'] if rec['type'] == "HGNC"][0]['id']
    gene_id = [rec for rec in json['dbReferences'] if rec['type'] == "GeneID"][0]['id']
    ensembl_nuc = mane_record['id']
    ensembl_prot = mane_record['properties']['protein sequence ID']
    refseq_nuc = mane_record['properties']['RefSeq nucleotide sequence ID']
    refseq_prot = mane_record['properties']['RefSeq protein sequence ID']
    gene_name = json['gene'][0]['name']['value']
    cache_dir = Path(kwargs.get("cache_dir","/tmp"))
    cache_dir.mkdir(exist_ok=True)
    local_pth = cache_dir / "MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz"
    if not local_pth.exists():
        urllib.request.urlretrieve("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz",local_pth)
    mane_gtf = pd.read_csv(local_pth, sep="\t",compression='gzip',header=None)
    mane_record = mane_gtf[(mane_gtf[2] == "gene") & \
        (mane_gtf[8].str.contains(f'gene_name "{gene_name}"',regex=False))].iloc[0]
    CHROM = mane_record[0].replace("chr","")
    START,STOP = mane_record[3],mane_record[4]
    STRAND = mane_record[6]

    return dict(seq=seq,
                MANE_RefSeq_nuc=refseq_nuc,
                gene_name=gene_name,
                MANE_RefSeq_prot=refseq_prot,
                CHROM=CHROM,
                START=START,
                STOP=STOP,
                HGNC_ID=hgnc_id,
                STRAND=STRAND,
                GeneID=gene_id,
                ENSEMBL_nuc=ensembl_nuc,
                ENSEMBL_prot=ensembl_prot)

def queryClinVarVCF(CHROM,START,STOP,gene_name=None, gene_id=None,**kwargs):
    filepath = Path(kwargs.get("filepath"))
    assert filepath.exists(), "filepath {} does not exist".format(filepath)
    picard_filepath = Path(kwargs.get("picard_filepath"))
    assert picard_filepath.exists(), "picard_filepath does not exist"
    write_dir = Path(kwargs.get("write_dir","/tmp"))
    write_dir.mkdir(exist_ok=True)

    output_file = write_dir / f"ClinVar_selectvariants_chr{CHROM}:{START}-{STOP}_{str(datetime.now()).replace(' ','_')}.vcf"
    cmd = f"gatk SelectVariants -V {filepath} -L {CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output {output_file}"
    subprocess.run(cmd.split(" "))


    tsvout = str(output_file).replace('.vcf','.tsv')
    variants2table = f"gatk VariantsToTable -V {output_file} -O {tsvout}"
    subprocess.run(variants2table.split(" "))

    clinVar_df = pd.read_csv(tsvout,delimiter='\t')
    if gene_name is not None and gene_id is not None:
        clinVar_df = clinVar_df[clinVar_df.GENEINFO == f"{gene_name}:{gene_id}"]

    clinVar_df = annotate_variants(clinVar_df)
    clinVar_df.rename(columns={"#CHROM":'CHROM'},inplace=True)
    return clinVar_df

def annotate_variants(df):
    low_qual_rev_stat = {'no_assertion_criteria_provided','no_classification_provided'}
    df = df.assign(CHROM=df.CHROM.astype(str),
                   POS=df.POS.astype(str),
                   REF=df.REF.astype(str),
                   ALT=df.ALT.astype(str),
                   is_pathogenic=(df.CLNSIG.isin({"Pathogenic","Likely_pathogenic","Pathogenic/Likely_pathogenic"})) & (~df.CLNREVSTAT.isin(low_qual_rev_stat)),
                   is_benign=(df.CLNSIG.isin({"Benign","Likely_benign","Benign/Likely_benign"})) & (~df.CLNREVSTAT.isin(low_qual_rev_stat)),
                   is_conflicting=df.CLNSIG == "Conflicting_classifications_of_pathogenicity",
                   is_VUS=df.CLNSIG == "Uncertain_significance")
    return df
    

def queryGnomAD(CHROM,START,STOP,HGNC_ID,**kwargs):
    """
    Query gnomAD for missense variants in a gene

    Steps:
    1) Get the chromosomal coordinates of the gene from the MANE GTF file
    2) Use GATK SelectVariants to extract variants in the gene from the gnomAD exomes and genomes VCF files
        2A) Filter for SNPs
        2B) Exclude filtered variants
    3) Use GATK MergeVcfs to combine the exomes and genomes VCF files
    4) Use GATK VariantsToTable to convert the combined VCF file to a TSV file
    5) Manually parse the VEP annotations in the TSV file
    6) Filter for missense variants

    Required args:
    - CHROM: str : The chromosome for which to query gnomAD
    - START: int : The minimum position in the chromosome for which to query gnomAD
    - STOP: int : The maximum position in the chromosome for which to query gnomAD

    Required kwargs:
    - gnomad_vcf_root: str : Path to the root directory of the gnomAD vcf directory
    - picard_filepath: str : Path to the picard.jar file

    Optional kwargs:
    - write_dir: str : Path to the directory where the output files will be written : default "/tmp"

    Returns:
    - missense_df: pd.DataFrame : A DataFrame containing parsed VEP annotations for matched missense variants in gnomAD exomes and genomes
    """
    gnomad_vcf_root = Path(kwargs.get("gnomad_vcf_root"))
    assert gnomad_vcf_root.exists(), "gnomad_vcf_root does not exist"
    write_dir = Path(kwargs.get("write_dir","/tmp"))
    write_dir.mkdir(exist_ok=True)
    java = Path(kwargs.get("java"))
    picard_filepath = Path(kwargs.get("picard_filepath"))
    assert picard_filepath.exists(), "picard_filepath does not exist"

    gnomAD_exomes_filepath = gnomad_vcf_root / f"exomes/gnomad.exomes.v4.1.sites.chr{CHROM}.vcf.bgz"
    gnomAD_genomes_filepath = gnomad_vcf_root / f"genomes/gnomad.genomes.v4.1.sites.chr{CHROM}.vcf.bgz"
    exomes_output_File = write_dir / f"selectvariants_{str(datetime.now()).replace(' ','_')}.exomes.vcf"
    genomes_output_File = write_dir / f"selectvariants_{str(datetime.now()).replace(' ','_')}.genomes.vcf"
    cmd = f"gatk SelectVariants -V {gnomAD_exomes_filepath} -L chr{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output {exomes_output_File}"
    subprocess.run(cmd.split(" "))
    cmd = f"gatk SelectVariants -V {gnomAD_genomes_filepath} -L chr{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output {genomes_output_File}"
    subprocess.run(cmd.split(" "))
    output_File = write_dir / f"combinevariants_{str(datetime.now()).replace(' ','_')}.vcf"
    cmd = f'{java} -jar {picard_filepath} MergeVcfs I={exomes_output_File} I={genomes_output_File} O={output_File}'
    subprocess.run(cmd.split(" "))
    tsvout = str(output_File).replace('.vcf','.tsv')
    variants2table = f"gatk VariantsToTable -V {output_File} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -ASF AC -ASF AF -ASF vep -O {tsvout}"
    subprocess.run(variants2table.split(" "))
    gnomAD_df = pd.read_csv(tsvout,delimiter='\t')
    vep_columns = get_vep_columns_from_vcf_header(output_File)
    vep_df = parse_vep(gnomAD_df,columns=vep_columns)
    gnomAD_df = pd.merge(gnomAD_df,vep_df,left_index=True,right_on='index',validate='one_to_many')
    gene_df = gnomAD_df[gnomAD_df.HGNC_ID == HGNC_ID]
    gene_df = gene_df.assign(CHROM=gene_df.CHROM.astype(str).str.replace("chr",""),
                                POS=gene_df.POS.astype(str),
                                REF=gene_df.REF.astype(str),
                                ALT=gene_df.ALT.astype(str))
    return gene_df

def get_vep_columns_from_vcf_header(vcf_file):
    vcf = VariantFile(vcf_file)
    return vcf.header.info['vep'].description.split("Format: ")[1].split("|")
    
def parse_vep(df,columns):
    vep_series = df.vep.apply(lambda r: list(map(lambda s: dict(zip(columns,s.split('|'))),r.split(","))))
    vep_df = pd.DataFrame(vep_series,index=df.index).explode('vep')
    vep_df = pd.DataFrame.from_records(vep_df.vep.values,index=vep_df.index).reset_index()
    return vep_df

def run_vep(df,**kwargs):
    write_dir = Path(kwargs.get("write_dir","/data/dzeiberg/mave_calibration/vep_data"))
    write_dir.mkdir(exist_ok=True)
    vep_input_file = write_dir / f"vep_input_{str(datetime.now()).replace(' ','_')}.vcf"
    vep_output_file = write_dir / f"vep_output_{str(datetime.now()).replace(' ','_')}.vcf"
    df.loc[:,['CHROM','POS','ID','REF','ALT']].to_csv(vep_input_file,sep='\t',index=False)
    cmd = f"sudo docker run -v /data/dbs/vep_data/:/data -v {write_dir}:/inputs ensemblorg/ensembl-vep vep --input_file /inputs/{vep_input_file.name} --output_file STDOUT --format vcf --vcf --offline --assembly GRCh38 --everything --fork 4 --cache --dir_cache /data --fasta /data/fasta/Homo_sapiens.GRCh38.cdna.all.fa.bgz --force_overwrite --no_stats"
    result = subprocess.run(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with open(vep_output_file,'w') as f:
        f.write(result.stdout)
    vcf = VariantFile(vep_output_file)
    vep_cols = vcf.header.info['CSQ'].description.split("Format: ")[1].split("|")
    
    vep_results = [line for line in result.stdout.split("\n") if line[:2] != "##"]
    vep_df = pd.read_csv(StringIO('\n'.join(vep_results)),delimiter='\t').dropna(subset='INFO')

    vep_series = vep_df.INFO.apply(lambda r: list(map(lambda s: dict(zip(vep_cols,s.split('|'))),r.split(","))))
    info_df = pd.DataFrame(vep_series,index=vep_df.index).explode('INFO')

    info_df = pd.DataFrame.from_records(info_df.INFO.values,index=info_df.index)
    vep_df = pd.merge(vep_df,info_df,left_index=True,right_index=True,validate='one_to_many')
    vep_csv_output_file = write_dir / f"vep_output_{str(datetime.now()).replace(' ','_')}.vcf"
    vep_df.to_csv(vep_csv_output_file)
    return vep_df

def load_external_data(uniprot_acc,**kwargs):
    gene_info = get_gene_info(uniprot_acc)
    with open(kwargs.get("external_config_path","/home/dzeiberg/mave_preprocessing/external_data.json")) as f:
        external_config = json.load(f)
    gnomAD_df = queryGnomAD(gene_info['CHROM'],
                            gene_info['START'],
                            gene_info['STOP'],
                            gene_info['HGNC_ID'],
                            **external_config['gnomAD'],
                            **external_config['external_tools'])
    clinvar_df = queryClinVarVCF(gene_info['CHROM'],
                                 gene_info['START'],
                                 gene_info['STOP'],
                                 gene_name=gene_info['gene_name'],
                                 gene_id=gene_info['GeneID'],
                                 **external_config['clinvar']['GRCh38'],
                                 **external_config['external_tools'])
    vep_df = run_vep(clinvar_df)
    clinvar_df = join_clinvar(clinvar_df,vep_df)
    clinvar_df = clinvar_df.assign(hgvs_pro=clinvar_df.HGVSp.str.split(":").str[1],
                                   hgvs_nuc=clinvar_df.apply(make_hgvs_nuc,axis=1))
    gnomAD_df = gnomAD_df.assign(hgvs_pro=gnomAD_df.HGVSp.str.split(":").str[1],
                                    hgvs_nuc=gnomAD_df.HGVSc.str.split(":").str[1])

    spliceAI_scores = querySpliceAI(gene_info['CHROM'],
                                    gene_info['START'],
                                    gene_info['STOP'],
                                    **external_config['splice_ai'],
                                    **external_config['external_tools'],
                                    assembly='hg38')
    clinvar_df = pd.merge(clinvar_df,spliceAI_scores,on=['CHROM','POS','REF','ALT'],validate='many_to_one',how='left')
    gnomAD_df = pd.merge(gnomAD_df,spliceAI_scores,on=['CHROM','POS','REF','ALT'],validate='many_to_one',how='left')

    return gene_info, gnomAD_df, clinvar_df

def extract_two_uppercase(s):
    # Extract only uppercase characters from the string
    uppercase_chars = [char for char in s if char.isupper()]
    
    # Return the first two uppercase characters or an empty list if less than two
    return (uppercase_chars[0], uppercase_chars[1]) if len(uppercase_chars) >= 2 else ("", "")

def make_hgvs_nuc(row):
    CDS_pos = row.CDS_position
    ref,alt = extract_two_uppercase(row.Codons)
    return f"c.{CDS_pos}{ref}>{alt}"

def querySpliceAI(chrom, position_min, position_max,**kwargs):
        """
        Query SpliceAI for spliceAI scores in a region of the genome

        Required args:
        - chrom: str : The chromosome for which to query SpliceAI
        - position_min: int : The minimum position in the chromosome for which to query SpliceAI
        - position_max: int : The maximum position in the chromosome for which to query SpliceAI

        Required kwargs:
        - spliceAIFilePath: str : Path to the SpliceAI VCF file

        Optional kwargs:
        - write_dir: str : Path to the directory where the output files will be written : default "/tmp"
        
        """
        write_dir = Path(kwargs.get('write_dir',"/tmp"))
        write_dir.mkdir(exist_ok=True)
        spliceAIFilePath = Path(kwargs.get('spliceAIFilePath'))
        assembly = kwargs.get("assembly")
        assert assembly in ['hg38','hg19'], "assembly must be one of ['hg38','hg19']"
        assert spliceAIFilePath.exists(), "spliceAIFilePath does not exist"
        output_filepath = write_dir / f'splice_ai_query_result.{str(datetime.now()).replace(" ","_")}.vcf'
        cmd = f"gatk SelectVariants -V {spliceAIFilePath} -L {chrom}:{max(position_min,1)}-{position_max} --output {output_filepath}"
        subprocess.run(cmd.split(" "))
        result_df = pd.read_csv(output_filepath,comment='#',delimiter='\t',header=None,
                        names='CHROM POS ID REF ALT QUAL FILTER INFO'.split(" "),
                            dtype={k : str for k in 'CHROM POS REF ALT'.split(" ")})
        result_df = result_df.assign(spliceAI_score=result_df.INFO.apply(lambda s: max(list(map(float,
                                                                                            s.split("|")[2:6])))))
        result_df.to_csv(write_dir / f"splice_ai_match_{assembly}_CHR{chrom}_{position_min}_{position_max}.tsv",sep='\t',index=False)
        return result_df


def join_clinvar(clinvar_df,vep_df):
    vep_df = vep_df.assign(CHROM=vep_df.loc[:,'#CHROM'].astype(str),
                            POS=vep_df.loc[:,'POS'].astype(str),
                            REF=vep_df.loc[:,'REF'].astype(str),
                            ALT=vep_df.loc[:,'ALT'].astype(str))
    return pd.merge(clinvar_df,vep_df,on=['CHROM','POS','REF','ALT'],validate='one_to_many')

if __name__ == "__main__":
    gene_info, gnomAD_df, clinvar_df = load_external_data("Q92560")

