#!/usr/bin/env python3
import numpy as np
import h5py
import pathlib
import warnings

# define mapping between tristan v1 and v2 data, v1 data is considered the default
__v2_map = {
            # Parameters file
            'c':'algorithm:c',
            'c_omp':'plasma:c_omp',
            'gamma0':'set-to-one',
            'interval':'output:interval',
            'istep':'output:istep',
            'me':'particles:m1',
            'mi':'particles:m2',
            'mx0':'grid:mx0',
            'my0':'grid:my0',
            'ppc0':'plasma:ppc0',
            'sigma':'plasma:sigma',
            'sizex':'node_configuration:sizex',
            'sizey':'node_configuration:sizey',
            'qi':'particles:ch2',
            'stride':'output:stride',
            'time':'timestep', # or maybe timestep/plasma:c_omp or possibly time:last
            # Particles file
            'ue':'u_1',
            've':'v_1',
            'we':'w_1',
            'ui':'u_2',
            'vi':'v_2',
            'wi':'w_2',
            'xe':'x_1',
            'ye':'y_1',
            'ze':'z_1',
            'xi':'x_2',
            'yi':'y_2',
            'zi':'z_2',
            'inde':'ind1',
            'indi':'ind2',
            'proce':'proc1',
            'proci':'proc1',
            # Fields file
            'bx':'bx',
            'by':'by',
            'bz':'bz',
            'ex':'ex',
            'ey':'ey',
            'ez':'ez',
            'jx':'jx',
            'jy':'jy',
            'jz':'jz',
            'dens':'compute_dens', # dens1 + dens2
            'densi':'dens2',
            }

# =============================================================================
def __detect_tristan_data_version(file: h5py.File) -> int:

    # The fields to query were chosen only because they don't exist in the other version of Tristan
    #  Files:  particles fields   spectra  paramater
    v1_keys = ('che',    'densi', 'gamma', 'acool')
    v2_keys = ('ind_1',   'dens1', 'ebins', 'algorithm:c')

    if any(key in file for key in v1_keys):
        return 1
    elif any(key in file for key in v2_keys):
        return 2
    else:
        raise ValueError('Data file format is not supported. If this is a data file from Tristan v1 or v2 please submit a bug report.')
# =============================================================================

# =============================================================================
def __verify_file_path(file_path: pathlib.Path, dataset_name: str):
    # Check that the file exists. If not, try to handle the special cases before raising exception
    if not file_path.exists():
        if 'param.' in file_path.name:
            file_path = file_path.with_name(file_path.name.replace('param', 'params'))
        elif 'spect.' in file_path.name:
            file_path = file_path.with_name(file_path.name.replace('spect', 'spec.tot'))
            warnings.warn('Spectra not yet supported with Tristan v2 data. Spectra plots will show dummy data with a value of 1.')
            if dataset_name in ['gmin','spece','specerest','specp','specprest','umean']:
                return np.ones((10,10))
            elif dataset_name in ['xsl', 'gamma']:
                return np.ones(10)
            else:
                return 1
        else:
            raise FileNotFoundError(f'File not found at path: {file_path}')

    return file_path
# =============================================================================

# =============================================================================
def __handle_tristan_v2(file_path: pathlib.Path, file: h5py.File, dataset_name: str, dataset_slice: tuple | slice):
    # Make sure the dataset is properly mapped
    try:
        dataset_name = __v2_map[dataset_name]
    except KeyError:
        raise KeyError(f'Requested dataset "{dataset_name}" is not mapped between Tristan v1 and v2 data ' \
                        f'when attempting to read from the file at {file_path}. Please add appropriate mapping')

    # Check if this dataset requires additional handling, if not then return it and exit early
    if dataset_name not in [__v2_map['dens'], __v2_map['gamma0']]:
        # Check that the dataset exists, return zero data and print warning if it doesn't.
        if dataset_name not in file:
            warnings.warn(f'{file_path} does not contain the dataset "{dataset_name}". Returning zero valued data.')
            return np.zeros(10)
        else:
            return file[dataset_name][dataset_slice]

    # Datasets that need special handling
    if dataset_name == __v2_map['dens']:
        return file['dens1'][dataset_slice] + file['dens2'][dataset_slice]
    elif dataset_name == __v2_map['gamma0']:
        warnings.warn('"gamma0" is not present in Tristan v2 datasets. Setting gamma0=1')
        return 1
    # elif dataset_name == __v2_map['?????']:
    #     return
    else:
        raise ValueError(f'Dataset "{dataset_name}" was indicated to require special handling but no clause was supplied to do that handling.')
# =============================================================================

# =============================================================================
def load_dataset(file_path: str | pathlib.Path, dataset_name: str, dataset_slice: tuple | slice) -> np.array:

    # First check argument types
    assert isinstance(file_path, pathlib.Path) or isinstance(file_path, str)
    assert isinstance(dataset_name, str)
    assert isinstance(dataset_slice, tuple) or isinstance(dataset_slice, slice)
    if isinstance(dataset_slice, tuple):
        assert len(dataset_slice) <= 3
        for element in dataset_slice:
            assert isinstance(element, slice)

    # Convert file_path to pathlib for easier usage later
    file_path = pathlib.Path(file_path)

    # Verify the file_path and perform any required fixes for Tristan v2 data
    file_path = __verify_file_path(file_path, dataset_name)
    if not isinstance(file_path, pathlib.Path):
        # in this case the file path handling has returned dummy data instead of
        # a file_path since reading that data is not supported yet. Given that we
        # will just return the dummy data directly
        return file_path

    # open file
    with h5py.File(file_path, 'r') as file:
        # detect v1 or v2 data
        data_version = __detect_tristan_data_version(file)

        # Raise exception and exit if the version is wrong
        if data_version != 1 and data_version != 2:
            raise ValueError(f'Improper Tristan data version detected. Version detected is "{data_version}", should be "1" or "2".')

        # If the data is from Tristan v1 then return it and exit early
        if data_version == 1:
            # might need to add a try/except here to catch KeyErrors
            return file[dataset_name][dataset_slice]

        # If the data is from Tristan v2 then perform whatever handling is needed. Note that at this point `data_version` must be 2
        return __handle_tristan_v2(file_path, file, dataset_name, dataset_slice)
# =============================================================================
