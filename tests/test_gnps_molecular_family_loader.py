import os
import numpy

import pytest
from nplinker.metabolomics.gnps.gnps_molecular_family_loader import GNPSMolecularFamilyLoader
from nplinker.metabolomics.molecular_family import MolecularFamily

from .test_metabolomics import molecular_families, spec_dict
from . import DATA_DIR

@pytest.fixture
def molecular_families_gnps():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    return sut.families()


@pytest.mark.parametrize("filename", [
    os.path.join(DATA_DIR, "edges.pairsinfo"),
    DATA_DIR / "edges.pairsinfo"
])
def test_has_molecular_families(filename):
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual) == 25769
    assert len(actual[0].spectra_ids) == 19


def test_families_are_identical(spec_dict, molecular_families: list[MolecularFamily]):
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    actual = GNPSMolecularFamilyLoader(filename).families()

    actual.sort(key= lambda x: min(x.spectra_ids))

    for i, x in enumerate(actual):
        x.id = i
        for spec_id in x.spectra_ids:
            x.add_spectrum(spec_dict[spec_id])
            

    for x in molecular_families:
        for spec in x.spectra:
            x.spectra_ids.add(spec.spectrum_id)

    numpy.testing.assert_array_equal(actual, molecular_families)