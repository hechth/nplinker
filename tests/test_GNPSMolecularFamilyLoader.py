import os
import pytest
from nplinker.metabolomics.GNPSMolecularFamilyLoader import \
    GNPSMolecularFamilyLoader
from nplinker.metabolomics.singleton_family import SingletonFamily
from . import DATA_DIR
from .test_metabolomics import spec_dict, molecular_families


@pytest.fixture
def molecular_families_gnps():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    return sut.families()


def test_has_molecular_families():
    filename = os.path.join(DATA_DIR, "edges.pairsinfo")
    sut = GNPSMolecularFamilyLoader(filename)
    actual = sut.families()
    assert len(actual) == 25769
    assert len(actual[0].spectra_ids) == 19


def test_families_equal(spec_dict, molecular_families_gnps, molecular_families):

    for fam in molecular_families_gnps:
        for spec_id  in fam.spectra_ids:
            fam.add_spectrum(spec_dict[spec_id])
            spec_dict[spec_id].family = fam
            spec_dict[spec_id].family_id = fam.family_id
    
    def sort_method(x):
        return len(x.spectra)
    
    molecular_families.sort(key = sort_method, reverse = True)
    molecular_families_gnps.sort(key = sort_method, reverse = True)

    assert molecular_families == molecular_families_gnps
    

def test_count_non_singleton_families(molecular_families_gnps):
    actual = list(filter(lambda x: not isinstance(x, SingletonFamily), molecular_families_gnps))
    assert len(actual) == 29