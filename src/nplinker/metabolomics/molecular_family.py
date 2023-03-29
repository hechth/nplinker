from typing_extensions import Self

from nplinker.metabolomics.spectrum import Spectrum
from nplinker.strain_collection import StrainCollection

class MolecularFamily():

    def __init__(self, family_id: int):
        """Class to model molecular families.

        Args:
            family_id(int): Id for the molecular family.
        """
        self.id: int = -1
        self.family_id: int = family_id
        self.spectra: list[Spectrum] = []
        self.family = None
        self.spectra_ids: set[int] = set()

    # def has_strain(self, strain):
    #     for spectrum in self.spectra:
    #         if spectrum.has_strain(strain):
    #             return True

    #     return False

    @property
    def strains(self) -> StrainCollection:
        """Get strains of spectra in the molecular family.

        Returns:
            set[StrainCollection]: StrainCollection of strains from which the spectra in the molecular family are coming.
        """
        strains: StrainCollection = StrainCollection()
        for spectrum in self.spectra:
            for strain in spectrum.strains:
                strains.add(strain)
        return strains

    def add_spectrum(self, spectrum: Spectrum):
        """Add a spectrum to the spectra list.

        Args:
            spectrum(Spectrum): Spectrum to add to the molecular family.
        """
        self.spectra.append(spectrum)

    def __str__(self) -> str:
        return 'MolFam(family_id={}, spectra={})'.format(
            self.family_id, len(self.spectra))

    def __eq__(self, other: Self) -> bool:
        return bool(self.id == other.id)

    def __hash__(self) -> int:
        return hash(self.id)


def map_spectra_to_families(spec_dict: dict[int, Spectrum], molecular_families: list[MolecularFamily]):
    """Map the spectra to the molecular families.

    Args:
        spec_dict(dict[int, Spectrum]): Dictionary mapping integer to Spectra.
        molecular_families(list[MolecularFamily]): _description_

    Examples:
        >>> 
        """
    for i, x in enumerate(molecular_families):
        x.id = i
        for spec_id in x.spectra_ids:
            x.add_spectrum(spec_dict[spec_id])
            spec_dict[spec_id].family = x
            spec_dict[spec_id].family_id = x.id
            
