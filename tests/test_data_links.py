from nplinker.scoring.linking.data_linking import DataLinks
from nplinker.strain_collection import StrainCollection
from nplinker.strains import Strain
from .test_metabolomics import spec_dict, spec_with_families
from .test_GNPSMolecularFamilyLoader import molecular_families_gnps


class GCFDummy():
    def __init__(self):
        self.id = 100

class MolecularFamilyDummy():
    def __init__(self):
        self.id = 200


def test_collect_mappings_spec(spec_with_families):
    sut = DataLinks()
    sut.collect_mappings_spec(spec_with_families.values())
    actual = sut.mapping_spec.shape

    assert actual == (25935,3)


def test_collect_mappings_spec_v2(molecular_families_gnps):
    sut = DataLinks()
    sut.collect_mappings_spec_v2(molecular_families_gnps)
    actual = sut.mapping_spec.shape

    assert actual == (25935,3)

def test_common_strains(spec_with_families):
    sut = DataLinks()
    strains = StrainCollection()
    strains.add(Strain(13))
    sut.matrix_strain_spec(spec_with_families.values(), strains)
    sut.collect_mappings_spec(spec_with_families.values())
    sut.data_family_mapping()
    actual = sut.common_strains([GCFDummy()], [MolecularFamilyDummy()], True)

    assert actual is not None


def test_data_family_mapping(spec_with_families):
    sut = DataLinks()
    strains = StrainCollection()
    strains.add(Strain(13))
    sut.matrix_strain_spec(spec_with_families.values(), strains)
    sut.collect_mappings_spec(spec_with_families.values())
    sut.data_family_mapping()

def test_data_family_mapping_v2(spec_dict, molecular_families_gnps):
    sut = DataLinks()
    strains = StrainCollection()
    strains.add(Strain(13))
    sut.matrix_strain_spec(spec_dict.values(), strains)
    sut.collect_mappings_spec_v2(molecular_families_gnps)
    sut.data_family_mapping()