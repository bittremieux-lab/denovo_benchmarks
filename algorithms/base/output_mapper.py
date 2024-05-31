"""Base class with methods for OutputMapper."""

from pyteomics import proforma

class OutputMapperBase:
    def _format_scores(self, scores):
        """
        Write a list of float per-token scores
        into a string of float scores separated by ','.
        """
        return ",".join(map(str, scores))

    def format_scan(self, scan):
        """
        Represent spectrum scan id as F{file_i}:{scan_id} string, where
        - `file_i` - number of the .mgf file in a dataset 
            (lexicographically sorted)
        - `scan_id` - scan id of the spectrum 
            (SCANS= in an .mgf file).
        """
        return scan

    def format_scan_index(self, scan_index):
        """
        Represent spectrum scan id as F{file_i}:{scan_index} string, where
        - `file_i` - number of the .mgf file in a dataset 
            (lexicographically sorted)
        - `scan_index` - number of the spectrum in an .mgf file.
        """
        return scan_index

    def format_sequence(self, sequence):
        """
        Convert peptide sequence to the common output data format 
        (ProForma with modifications represented 
        in the delta mass notation).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.  
        """
        return sequence

    def format_sequence_and_scores(self, sequence, aa_scores):
        """
        Convert peptide sequence to the common output data format
        (ProForma with modifications represented 
        in the delta mass notation) and modify per-token scores 
        if needed.
        (Only needed if per-token scores have to be modified 
        to correspond the transformed sequence in ProForma format.
        Otherwise use `format_sequence` method).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original algorithm output format.
        aa_scores: str
            String of per-token scores for each token in the sequence.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the common output data format.
        transformed_aa_scores: str
            String of per-token scores corresponding to each token
            in the transformed sequence.
        """
        sequence = self.format_sequence(sequence)
        return sequence, aa_scores

    def simulate_token_scores(self, pep_score, sequence):
        """
        Define proxy per-token scores from the peptide score
        if per-token scores are not provided by the model.
        Expects the sequence to be already in 
        the ProForma delta mass notation!
        """
        try:
            seq = proforma.parse(sequence)
        except:
            print(sequence)
        n_tokens = len(seq[0])
        if seq[1]["n_term"]:
            n_tokens += len(seq[1]["n_term"])
        if seq[1]["c_term"]:
            n_tokens += len(seq[1]["c_term"])
            
        scores = [str(pep_score),] * n_tokens
        return self._format_scores(scores)

    def format_output(self, output_data):
        """
        Transform ['scans', 'sequence', 'score', 'aa_scores'] columns 
        of `output_data` dataframe to the common outout format.
        Assumes that predicted sequences are provided 
        for all dataframe entries (no NaNs).

        Parameters
        ----------
        output_data : pd.DataFrame
            Dataframe with algorithm outputs. Must contain columns:
            - 'sequence' - predicted peptide sequence;
            - 'score' - confidence score for the predicted sequence;
            - 'aa_scores' - per-amino acid scores, if available. 
                Otherwise, the whole peptide `score` will be used 
                as a score for each amino acid.
            - 'scans' or 'scan_indices' - information to match 
                each prediction with its ground truth sequence. 
                One of two columns must be provided.

        Returns
        -------
        transformed_output_data : pd.DataFrame
            Dataframe with algorithm predictions
            in the common output data format. 
        """
        
        if "aa_scores" in output_data:
            output_data[["sequence", "aa_scores"]] = output_data.apply(
                lambda row: self.format_sequence_and_scores(row["sequence"], row["aa_scores"]),
                axis=1,
                result_type="expand",
            )
        
        else:
            output_data["sequence"] = output_data["sequence"].apply(
                self.format_sequence,
            )
            output_data["aa_scores"] = output_data.apply(
                lambda row: self.simulate_token_scores(row["score"], row["sequence"]), 
                axis=1,
            )
            
        if "scans" in output_data:
            output_data["scans"] = output_data["scans"].apply(self.format_scan)

        if "scan_indices" in output_data:
            output_data["scan_indices"] = output_data["scan_indices"].apply(
                self.format_scan_index
            )

        return output_data
