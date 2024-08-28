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

    def format_input(self, spectrum, filename):
        """
        Convert the spectrum (annotation sequence and params) to the
        input format expected by the algorithm.

        Parameters
        ----------
        spectrum : dict
            Peptide sequence in the original format.
        filename: int
            Name of .mgf file being processed. Used to ensure a unique
            scan_id for each spectrum.

        Returns
        -------
        transformed_spectrum : dict
            Peptide sequence in the algorithm input format.
        """
        
        spectrum["params"]["seq"] = "PEPTIDE"

        scan_id = spectrum["params"]["scans"]
        spectrum["params"]["scans"] = filename + ":" + scan_id
        return spectrum
