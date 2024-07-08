#!/usr/bin/env python3
import numpy as np
import h5py
import pathlib

# define mapping between tristan v1 and v2 data, v1 data is considered the default
__v2_map = {
            # Parameters file
            'c_omp':'plasma:c_omp',
            'istep':'output:istep',
            'me':'particles:m1',
            'mi':'particles:m2',
            'ppc0':'plasma:ppc0',
            'qi':'particles:ch2',
            'stride':'output:stride',
            'time':'timestep', # or maybe timestep/plasma:c_omp
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
    v2_keys = ('bx_1',   'dens1', 'ebins', 'algorithm:c')

    if any(key in file for key in v1_keys):
        return 1
    elif any(key in file for key in v2_keys):
        return 2
    else:
        raise ValueError('Data file format is not supported. If this is a data file from Tristan v1 or v2 please submit a bug report.')
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

    # open file
    with h5py.File(file_path, 'r') as file:
        # detect v1 or v2 data
        data_version = __detect_tristan_data_version(file)

        # Raise exception and exit if the version is wrong
        if data_version != 1 or data_version != 2:
            raise ValueError(f'Improper Tristan data version detected. Version detected is "{data_version}", should be "1" or "2".')

        # If the data is from Tristan v1 then return it and exit early
        if data_version == 1:
            return file[dataset_name][dataset_slice]

        # =====
        # If the data is from Tristan v2 then perform whatever handling is needed. Note that at this point `data_version` must be 2
        # =====

        # Make sure the dataset is properly mapped
        try:
            dataset_name = __v2_map[dataset_name]
        except KeyError:
            raise KeyError(f'Requested dataset "{dataset_name}" is not mapped between Tristan v1 and v2 data. Please add appropriate mapping')

        # check if this dataset requires additional handling, if not then return it and exit early
        if dataset_name not in [__v2_map['dens']]:
            return file[dataset_name][dataset_slice]

        # Datasets that need special handling
        if dataset_name == __v2_map['dens']:
            return file['dens1'][dataset_slice] + file['dens2'][dataset_slice]
        # elif dataset_name == __v2_map['?????']:
        #     pass
        else:
            raise ValueError(f'Dataset "{dataset_name}" was indicated to require special handling but no clause was supplied to do that handling.')
# =============================================================================
