# Methods to run MMseqs2 and calculate number of proteome matches

import os
import subprocess
import pandas as pd


def isoleucine_to_leucine(sequence):
    return sequence.replace("I", "L") if isinstance(sequence, str) else sequence


def setup_mmseqs_dirs(search_tmp_dir="./mmseqs2_tmp"):
    """Create directories for MMseqs2 proteome matches search."""
    # dir for MMseqs2 files
    os.makedirs(search_tmp_dir, exist_ok=True) 
    # dir for MMseqs2 runs tmp files
    tmp_files_dir = os.path.join(search_tmp_dir, "tmp")
    os.makedirs(tmp_files_dir, exist_ok=True)
    # dir target DB (multiple files)
    target_db_dir = os.path.join(search_tmp_dir, "targetDB")
    os.makedirs(target_db_dir, exist_ok=True)
    # dir query DB (multiple files)
    query_db_dir = os.path.join(search_tmp_dir, "queryDB")
    os.makedirs(query_db_dir, exist_ok=True)
    # dir raw search results (multiple files)
    search_result_dir = os.path.join(search_tmp_dir, "results")
    os.makedirs(search_result_dir, exist_ok=True)
    # path to search results in .m8 format
    search_result_path = os.path.join(search_tmp_dir, "search_result.m8") 
    return search_tmp_dir, tmp_files_dir, target_db_dir, query_db_dir, search_result_dir, search_result_path


def create_query_fasta(output_data, query_fasta_path):
    """
    Create query database of de novo predicted peptides 
    (with I replaced by L)
    """
    
    with open(query_fasta_path, "w") as f:
        f.write(
            "\n".join([
                f">{idx}\n{isoleucine_to_leucine(peptide)}" # I -> L
                for idx, peptide 
                in output_data.sequence_no_ptm.reset_index().values
            ])
        )
    print("queryDB (fasta):", query_fasta_path)


def create_target_fasta(reference_proteome_path, contam_path, target_fasta_path):
    """
    Create target database for MMseqs2 search:
    - join reference proteome with common contaminants
    - replace I with L (indistinguishable by mass)
    """
    
    cmd = [
        "cat",
        reference_proteome_path,
        contam_path,
        "| sed -e 's/I/L/g'",
        ">",
        target_fasta_path,
    ]
    subprocess.run(" ".join(cmd), shell=True, check=True)
    print("targetDB (fasta):", target_fasta_path)


def run_mmseqs(
    target_fasta_path,
    query_fasta_path,
    target_db_dir,
    query_db_dir,
    search_result_dir,
    search_result_path,
    tmp_files_dir,
    args=[],
):
    SEARCH_COLS = "query,target,qaln,taln,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"
    QUERY_KEY = "target"
    
    print("\n CLEAN EXISTING RESULTS")
    cmd = ["rm -rf", search_result_dir + "/*"]
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n CREATE TARGET DB")
    target_db_path = os.path.join(target_db_dir, "targetDB")
    cmd = [
        "mmseqs",
        "createdb",
        target_fasta_path,
        target_db_path,
        "-v 1",
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n CREATE QUERY DB")
    query_db_path = os.path.join(query_db_dir, "queryDB")
    cmd = [
        "mmseqs",
        "createdb",
        query_fasta_path,
        query_db_path,
        "-v 1",
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    print("\n SEARCH/MAP SEQUENCES")
    # mmseqs map <i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir> [options]
    cmd = [
        "mmseqs",
        "map",
        target_db_path,
        query_db_path,
        os.path.join(search_result_dir, "search_result"),
        tmp_files_dir,
        "-a", # backtrace alignments 
        "-e inf",
        "--cov-mode 1",
#         "-v 1",
    ] + args
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)

    print("\n CONVERT RESULTS TO .m8")
    # mmseqs convertalis queryDB targetDB resultDB resultDB.m8
    cmd = [
        "mmseqs",
        "convertalis",
        target_db_path,
        query_db_path,
        os.path.join(search_result_dir, "search_result"),
        search_result_path,
        f"--format-output {SEARCH_COLS}",
        "-v 1"
    ]
    print(" ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True, check=True)
    
    search_df = pd.read_csv(search_result_path, sep="\t", names=SEARCH_COLS.split(","))
    search_df = search_df.drop_duplicates(QUERY_KEY)
    search_df = search_df.set_index(search_df[QUERY_KEY].apply(int))
    search_df = search_df[search_df.mismatch <= 2] # filter s.t. n_mismatches <= 2
    return search_df
