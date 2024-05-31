"""Base class with methods for InputMapper."""

class InputMapperBase:
    def __init__(self,):
        pass

    def format_sequence(self, sequence):
        """
        Convert peptide sequence to the algorithm input format.
        Only required if the algorithm uses sequence information 
        at runtime (e.g., for internal performance evaluation).

        Parameters
        ----------
        sequence : str
            Peptide sequence in the original format.

        Returns
        -------
        transformed_sequence : str
            Peptide sequence in the algorithm input format.
        """
        return sequence

    def format_input(self, spectrum, file_i=0):
        """
        Convert the spectrum (annotation sequence and params) to the
        input format expected by the algorithm.

        Parameters
        ----------
        spectrum : dict
            Peptide sequence in the original format.
        file_i: int
            Number of .mgf file being processed. Used to ensure a unique
            scan_id for each spectrum.

        Returns
        -------
        transformed_spectrum : dict
            Peptide sequence in the algorithm input format.
        """
        spectrum["params"]["seq"] = self.format_sequence(
            spectrum["params"]["seq"]
        )

        # add file id to the spectrum scan data
        spectrum["params"]["scans"] = "F{}:{}".format(
            file_i, spectrum["params"]["scans"]
        )
        return spectrum
