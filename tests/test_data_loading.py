import pathlib
import h5py
import pytest
import data_loading

repo_root = pathlib.Path(__file__).resolve().parents[1]

def test__detect_tristan_data_version_verify_v1():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v1'
    file_names = ['flds.tot.041', 'param.041', 'prtl.tot.041', 'spect.041']

    for name in file_names:
        with h5py.File(data_dir / name, 'r') as file:
            version = data_loading.__detect_tristan_data_version(file)
            assert version == 1

def test__detect_tristan_data_version_verify_v2():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070', 'spec.tot.00070']

    for name in file_names:
        with h5py.File(data_dir / name, 'r') as file:
            version = data_loading.__detect_tristan_data_version(file)
            assert version == 2


def test__detect_tristan_data_version_incorrect_data():
    data_dir = repo_root / 'tests' / 'data' / 'unsupported_data'
    file_names = ['flds.tot.041', 'param.041', 'prtl.tot.041', 'spect.041']

    for name in file_names:
        with h5py.File(data_dir / name, 'r') as file:
            with pytest.raises(ValueError):
                version = data_loading.__detect_tristan_data_version(file)

def test__insert_directory_default():
    pass

def test__insert_directory_specified_location():
    pass

def test__verify_file_path():
    pass

def test___handle_tristan_v2():
    pass

def test_load_dataset():
    pass