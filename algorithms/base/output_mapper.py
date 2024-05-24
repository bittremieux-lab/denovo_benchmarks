from pyteomics import proforma

class OutputMapperBase:
    def _format_scores(self, scores):
        """
        Write a list of float per-token scores
        into a string of float scores separated by ','.
        """
        return ",".join(map(str, scores))

    def format_scan(self, scan):
        return scan

    def format_scan_index(self, scan_index):
        return scan_index

    def format_sequence(self, sequence):
        return sequence

    def format_sequence_and_scores(self, sequence, aa_scores):
        sequence = self.format_sequence(sequence)
        return sequence, aa_scores

    def simulate_token_scores(self, pep_score, sequence):
        """
        TODO.
        (if per-token scores are not provided by the model,
        define proxy per-token scores from the peptide score)
        """
        seq = proforma.parse(sequence)
        n_tokens = len(seq[0])
        if seq[1]["n_term"]:
            n_tokens += len(seq[1]["n_term"])
        if seq[1]["c_term"]:
            n_tokens += len(seq[1]["c_term"])
            
        scores = [str(pep_score),] * n_tokens
        return self._format_scores(scores)

    def format_output(self, output_data):
        """TODO."""
        # TODO: assumes that all the predicted sequences are provided (no NaNs)
        
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
