"""Base class with methods for InputMapper."""

class InputMapperBase:
    def __init__(self,):
        pass

    def format_input(self, spectrum):
        """
        Convert the spectrum (annotation sequence and params) to the
        input format expected by the algorithm.

        Parameters
        ----------
        spectrum : dict
            Peptide sequence in the original format.

        Returns
        -------
        transformed_spectrum : dict
            Peptide sequence in the algorithm input format.
        """
        # Any input format changes
        
        # Dummy annotation if expected by the algorithm
        spectrum["params"]["seq"] = "PEPTIDE"

        return spectrum
