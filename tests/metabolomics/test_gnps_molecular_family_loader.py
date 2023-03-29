import os
import numpy
import pytest
from nplinker.metabolomics.gnps.gnps_molecular_family_loader import \
    GNPSMolecularFamilyLoader
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.molecular_family import MolecularFamily, map_spectra_to_families
from nplinker.metabolomics.spectrum import Spectrum
from .. import DATA_DIR


@pytest.fixture
def molecular_families(spec_dict) -> list[MolecularFamily]:
    return make_families(spec_dict.values())

@pytest.mark.parametrize("filename", [
    os.path.join(DATA_DIR, "edges.pairsinfo"),
    DATA_DIR / "edges.pairsinfo"
])
def test_has_molecular_families(filename):
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual) == 25769
    assert len(actual[0].spectra_ids) == 19


def test_families_are_identical(spec_dict, molecular_families):
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    actual = GNPSMolecularFamilyLoader(filename).families()

    actual.sort(key= lambda x: min(x.spectra_ids))

    map_spectra_to_families(spec_dict, actual)


    for x in molecular_families:
        for spec in x.spectra:
            x.spectra_ids.add(spec.spectrum_id)

    numpy.testing.assert_array_equal(actual, molecular_families)

