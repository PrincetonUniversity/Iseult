import pathlib
import h5py
import numpy as np
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
                data_loading.__detect_tristan_data_version(file)

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
        found_path = data_loading.__verify_file_path(original_path)
        assert found_path == original_path

def test__verify_file_path_v2_data_single_dir():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070', 'spec.tot.00070']

    for name in file_names:
        original_path = data_dir / name
        found_path = data_loading.__verify_file_path(original_path)
        assert found_path == original_path

def test__verify_file_path_v2_data_multi_dir():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'standard_structure'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070', 'spec.tot.00070']

    for name in file_names:
        found_path = data_loading.__verify_file_path(data_dir / name)

        if 'flds' in name:
            fiducial_path = data_dir / 'flds' / name
        elif 'prtl' in name:
            fiducial_path = data_dir / 'prtl' / name
        elif 'spec' in name:
            fiducial_path = data_dir / 'spec' / name
        else:
            fiducial_path = data_dir / name

        assert found_path == fiducial_path

def test___handle_tristan_v2_datasets():
    data_dir = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory'
    file_names = ['flds.tot.00070', 'params.00070', 'prtl.tot.00070']#, 'spec.tot.00070']
    dataset_names = {'params.00070':['c','c_omp','gamma0','interval','istep',
                                     'me','mi','mx0','my0','ppc0','sigma',
                                     'sizex','sizey','qi','stride','time'],
                     'prtl.tot.00070':['ue','ve','we','ui','vi','wi','xe','ye',
                                       'ze','xi','yi','zi','inde','indi','proce',
                                       'proci'],
                     'flds.tot.00070':['bx','by','bz','ex','ey','ez','jx','jy',
                                       'jz','dens','densi']}
    fiducial_data = {'params.00070':[np.array((0.45)),np.array((25.0)),1,
                                     np.array((5.)),np.array((1)),np.array((1.)),
                                     np.array((1.)),np.array((6000)),
                                     np.array((12000)),np.array((16.)),
                                     np.array((0.3)),np.array((2.)),
                                     np.array((240.)),np.array((1.)),
                                     np.array((50.)),np.array((325025.))],
                     'prtl.tot.00070':[np.array((18.)),np.array((22.)),
                                       np.array((26.)),np.array((19.)),
                                       np.array((23.)),np.array((27.)),
                                       np.array((34.)),np.array((38.)),
                                       np.array((42.)),np.array((35.)),
                                       np.array((39.)),np.array((43.)),
                                       np.full((5), 10),np.full((5), 11),
                                       np.full((5), 14),np.full((5), 15),
                                       np.full((5), 10)],
                     'flds.tot.00070':[np.full((1,5,10), 0),np.full((1,5,10), 1),
                                       np.full((1,5,10), 2),np.full((1,5,10), 11),
                                       np.full((1,5,10), 12),np.full((1,5,10), 13),
                                       np.full((1,5,10), 14),np.full((1,5,10), 15),
                                       np.full((1,5,10), 16),np.full((1,5,10), 7),
                                       np.full((1,5,10), 4)]}

    for name in file_names:
        with h5py.File(data_dir / name, 'r') as file:
            for i, dataset_name in enumerate(dataset_names[name]):
                # catch the cases with warnings
                if dataset_name in ['gamma0']:
                    with pytest.warns(UserWarning):
                        result = data_loading.__handle_tristan_v2(data_dir, file, dataset_name, slice(None))
                else:
                    result = data_loading.__handle_tristan_v2(data_dir, file, dataset_name, slice(None))
                    assert np.array_equiv(result, fiducial_data[name][i]), \
                        f'Datasets {dataset_name} in file {name} do not match expectations.' \
                        f'Expected {fiducial_data[name][i] = } got {result = }'

def test___handle_tristan_v2_KeyError():
    data_path = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory' / 'spec.tot.00070'
    with h5py.File(data_path, 'r') as file:
        with pytest.raises(KeyError):
            data_loading.__handle_tristan_v2(data_path, file, 'bad_name', slice(None))

def test___handle_tristan_v2_slicing():
    data_path = repo_root / 'tests' / 'data' / 'tristan_v1' / 'flds.tot.041'
    test_slice = (slice(0,1), slice(1,5,1), slice(2,8,2))
    fiducial_arr = np.zeros((1,4,3))
    with h5py.File(data_path, 'r') as file:
        result = data_loading.__handle_tristan_v2(data_path, file, 'bx', test_slice)
        assert np.array_equiv(fiducial_arr, result)

def test_load_dataset_v1_data():
    file_path = repo_root / 'tests' / 'data' / 'tristan_v1' / 'flds.tot.041'
    dataset = 'by'
    fiducial_dataset = np.full((1,5,10), 1)

    test_dataset = data_loading.load_dataset(file_path, dataset)
    assert np.array_equiv(fiducial_dataset, test_dataset)

def test_load_dataset_v2_data():
    file_path = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory' / 'flds.tot.00070'
    dataset = 'by'
    fiducial_dataset = np.full((1,5,10), 1)

    test_dataset = data_loading.load_dataset(file_path, dataset)
    assert np.array_equiv(fiducial_dataset, test_dataset)

def test_load_dataset_unsupported_data():
    file_path = repo_root / 'tests' / 'data' / 'unsupported_data' / 'flds.tot.041'
    dataset = 'by'

    with pytest.raises(ValueError):
        data_loading.load_dataset(file_path, dataset)

def test_load_dataset_slicing():
    file_path = repo_root / 'tests' / 'data' / 'tristan_v2' / 'single_directory' / 'flds.tot.00070'
    dataset = 'by'
    test_slice = (slice(0,1), slice(1,5,1), slice(2,8,2))
    fiducial_dataset = np.full((1,4,3), 1)

    test_dataset = data_loading.load_dataset(file_path, dataset, test_slice)
    assert np.array_equiv(fiducial_dataset, test_dataset)

def test_load_dataset_scalar_return():
    file_path = repo_root / 'tests' / 'data' / 'tristan_v1' / 'param.041'
    dataset = 'c_omp'
    fiducial_dataset = 8

    test_dataset = data_loading.load_dataset(file_path, dataset)
    assert test_dataset == fiducial_dataset
