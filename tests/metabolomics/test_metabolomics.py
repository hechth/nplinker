from nplinker.metabolomics.gnps.gnps_molecular_family_loader import GNPSMolecularFamilyLoader
from nplinker.metabolomics.gnps.gnps_spectrum_loader import GNPSSpectrumLoader
from nplinker.metabolomics.metabolomics import make_families
from nplinker.metabolomics.metabolomics import load_dataset
from nplinker.strain_collection import StrainCollection
from .. import DATA_DIR


def test_load_spectra(spec_dict):
    assert len(spec_dict.keys()) > 0


def test_load_dataset():
    strains = StrainCollection()
    strains.add_from_file(DATA_DIR / "strain_mappings.csv")

    mgf_file = DATA_DIR / "spectra.mgf"
    edges_file = DATA_DIR /  "edges.pairsinfo"
    nodes_file = DATA_DIR / "nodes.tsv"

    spec_dict, spectra, molecular_families, unknown_strains = load_dataset(
        strains,
        mgf_file,
        edges_file,
        nodes_file
    )

    assert isinstance(spec_dict, dict)
    assert len(spectra) > 1

    expected_spectra = GNPSSpectrumLoader(mgf_file).spectra()
    assert spectra == expected_spectra

    expected_families = GNPSMolecularFamilyLoader(edges_file).families()
    assert molecular_families == expected_families

    


def test_make_families(spec_dict):
    families = make_families(spec_dict.values())
    assert len(families) == 25769
