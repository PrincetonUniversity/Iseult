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
    test_path = pathlib.Path('/dir1/dir2/dir3/file.txt')
    test_path = data_loading.__insert_directory(test_path, 'new_dir')

    fiducial_path = pathlib.Path('/dir1/dir2/dir3/new_dir/file.txt')
    assert test_path == fiducial_path

def test__insert_directory_specified_location():
    test_path = pathlib.Path('/dir1/dir2/dir3/file.txt')
    test_path = data_loading.__insert_directory(test_path, 'new_dir', 3)

    fiducial_path = pathlib.Path('/dir1/dir2/new_dir/dir3/file.txt')
    assert test_path == fiducial_path

def test__verify_file_path_v1_data():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v1'
    file_names = ['flds.tot.041', 'param.041', 'prtl.tot.041', 'spect.041']

    for name in file_names:
        original_path = data_dir / name
        found_path = data_loading.__verify_file_path(original_path, 'dataset_name')
        assert found_path == original_path

def test__verify_file_path_v2_data_single_dir():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070', 'spec.tot.00070']

    for name in file_names:
        original_path = data_dir / name
        found_path = data_loading.__verify_file_path(original_path, 'dataset_name')
        assert found_path == original_path

def test__verify_file_path_v2_data_multi_dir():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'standard_structure'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070', 'spec.tot.00070']

    for name in file_names:
        found_path = data_loading.__verify_file_path(data_dir / name, 'dataset_name')

        if 'flds' in name:
            fiducial_path = data_dir / 'flds' / name
        elif 'prtl' in name:
            fiducial_path = data_dir / 'prtl' / name
        elif 'spec' in name:
            fiducial_path = data_dir / 'spec' / name
        else:
            fiducial_path = data_dir / name

        assert found_path == fiducial_path

def test___handle_tristan_v2():
    pass

def test_load_dataset():
    pass