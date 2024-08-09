#!/usr/bin/env python3
import numpy as np
import h5py
import pathlib
import warnings

# =============================================================================
def __detect_tristan_data_version(file: h5py.File) -> int:
    """Determine which version of Tristan the data file belongs to.

    Parameters
    ----------
    file : h5py.File
        The HDF5 file to identify the version of.

    Returns
    -------
    int
        An integer indicating which version of Tristan the file was generated with. Should be 1 or 2

    Raises
    ------
    ValueError
        If the version of Tristan cannot be identified.
    """
    # The fields to query were chosen only because they don't exist in the other version of Tristan
    #  Files:  particles fields   spectra  paramater
    v1_keys = ('che',   'densi', 'gamma', 'acool')
    v2_keys = ('ind_1', 'dens1', 'ebins', 'algorithm:c')

    if any(key in file for key in v1_keys):
        return 1
    elif any(key in file for key in v2_keys):
        return 2
    else:
        raise ValueError('Data file format is not supported. If this is a data file from Tristan v1 or v2 please submit a bug report.')
# =============================================================================

# =============================================================================
def __insert_directory(path: pathlib.Path, inserted_directory: str, location: int = -1) -> pathlib.Path:
    """Insert a directory at a specific spot in the path

    Parameters
    ----------
    path : pathlib.Path
        The path to insert a directory into.
    inserted_directory : str
        The name of the directory to insert
    location : int, optional
        The location to insert the directory, by default -1 to be right before the file name

    Returns
    -------
    pathlib.Path
        The pathlib.Path with the inserted directory in the path
    """
    path_list = list(pathlib.PurePath(path).parts)
    path_list.insert(location, inserted_directory)
    return pathlib.Path(pathlib.PurePath('').joinpath(*path_list))
# =============================================================================

# =============================================================================
def __verify_file_path(file_path: pathlib.Path) -> pathlib.Path:
    """Verify that a HDF5 file exists at `file_path` and if it doesn't then try to find it based on known directory structures of Tristan v1 and v2 outputs

    Parameters
    ----------
    file_path : pathlib.Path
        The path to the file

    Returns
    -------
    pathlib.Path
        The confirmed or corrected path to the file

    Raises
    ------
    FileNotFoundError
        If the file cannot be found then this error is raised.
    """
    unmodified_path = file_path
    # Check that the file exists. If not, try to handle the special cases before raising exception
    if file_path.exists():
        return file_path

    # Now that we've verified the file doesn't exist let's try some fixes

    # Correct the name for Tristan v2 data
    if 'param.' in file_path.name:
        file_path = file_path.with_name(file_path.name.replace('param', 'params'))
    elif 'spect.' in file_path.name:
        file_path = file_path.with_name(file_path.name.replace('spect', 'spec.tot'))

    # Check if the corrected path exists. If it does then it's Tristan v2 data
    # that is all in one directory. If not then we need to handle the case where
    # the data is in a subdirectory
    if file_path.exists():
        return file_path

    # Data is in a subdirectory so let's add that to the path
    if 'flds' in file_path.name:
        return __insert_directory(file_path, 'flds')
    elif 'prtl' in file_path.name:
        return __insert_directory(file_path, 'prtl')
    elif 'spec' in file_path.name:
        return __insert_directory(file_path, 'spec')

    raise FileNotFoundError(f'File not found at path: {unmodified_path} and could not be found in any standard places.')
# =============================================================================

# =============================================================================
def __handle_tristan_v2_spectra(file_path: pathlib.Path, file: h5py.File, dataset_name: str) -> np.ndarray:
    spectral_data = file[dataset_name]
    spectral_data = np.sum(spectral_data, axis=(1,2))
    return spectral_data
# =============================================================================

# =============================================================================
def __handle_tristan_v2(file_path: pathlib.Path, file: h5py.File, dataset_name: str, dataset_slice: tuple | slice) -> np.ndarray | int | np.float64:
    """Load Tristan v2 data and perform any necessary transformation to convert it to the same format as Tristan v1 data

    Parameters
    ----------
    file_path : pathlib.Path
        The path to the file.
    file : h5py.File
        The HDF5 file to open
    dataset_name : str
        The name of the dataset to load, should be the Tristan v1 name
    dataset_slice : tuple | slice
        How to slice the data

    Returns
    -------
    np.ndarray | int
        The loaded data. Can either be an array or a scalar

    Raises
    ------
    KeyError
        Raised if there isn't a mapping between the Tristan v1 and v2 data for the provided dataset
    ValueError
        Raised if the dataset requires special handling but no clause for the handling exists
    """
    # define mapping between tristan v1 and v2 data, v1 data is considered the default
    v2_map = {# Parameters file
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
              'inde':'ind_1',
              'indi':'ind_2',
              'proce':'proc_1',
              'proci':'proc_2',
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
              # Spectra
              # HACK: These mappings are all first pass guesses, they still need to be verified
              'gamma':'ebins',
            #   'gmax':'',
            #   'gmin':'',
              'spece':'n1',
              'specp':'n3',
              'specerest':'restframe_unsupported',
              'specprest':'restframe_unsupported',
              'xsl':'xbins',
              }

    # Make sure the dataset is properly mapped
    try:
        dataset_name = v2_map[dataset_name]
    except KeyError:
        raise KeyError(f'Requested dataset "{dataset_name}" is not mapped between Tristan v1 and v2 data ' \
                        f'when attempting to read from the file at {file_path}. Please add appropriate mapping')

    # Check if this dataset requires additional handling, if not then return it and exit early
    special_handling_list = [v2_map['dens'], v2_map['gamma0'],
                             v2_map['spece'], v2_map['specerest'],
                             v2_map['specp'], v2_map['specprest']]
    if dataset_name not in special_handling_list:
        # Check that the dataset exists, return zero data and print warning if it doesn't.
        if dataset_name not in file:
            warnings.warn(f'{file_path} does not contain the dataset "{dataset_name}". Returning zero valued data.')
            return np.zeros(10)
        else:
            return file[dataset_name][dataset_slice]

    # Datasets that need special handling
    if dataset_name == v2_map['dens']:
        return file['dens1'][dataset_slice] + file['dens2'][dataset_slice]
    elif dataset_name == v2_map['gamma0']:
        warnings.warn('"gamma0" is not present in Tristan v2 datasets. Setting gamma0=1')
        return np.array([1])
    elif dataset_name in [v2_map['spece'], v2_map['specp']]:
        return __handle_tristan_v2_spectra(file_path, file, dataset_name)[dataset_slice]
    elif dataset_name in [v2_map['specerest'], v2_map['specprest']]:
        warnings.warn('Rest frame spectra are not supported with Tristan v2 Data. Using dummy data.')
        shape = list(file['n1'].shape)
        return np.ones((shape[0],shape[-1])) # only use the energy and x dimensions since we would sum over the others normally
    else:
        raise ValueError(f'Dataset "{dataset_name}" was indicated to require special handling but no clause was supplied to do that handling.')
# =============================================================================

# =============================================================================
def load_dataset(file_path: str | pathlib.Path, dataset_name: str, dataset_slice: tuple | slice = slice(None)) -> np.ndarray | int | np.float64:
    """Load a dataset from a datafile. Automatically detects the format and loads+transforms the data to match Tristan v1 formats

    Parameters
    ----------
    file_path : str | pathlib.Path
        The path to the file
    dataset_name : str
        The Tristan v1 name of the dataset to load
    dataset_slice : tuple | slice, optional
        How to slice the dataset. Can either be a Slice object or a tuple of
        Slice object where the nth element is the slice for the nth dimension,
        by default slice(None) which loads all the data

    Returns
    -------
    np.ndarray | int | np.float64
        The loaded dataset

    Raises
    ------
    ValueError
        The file is not from Tristan v1 or v2
    """
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
    file_path = __verify_file_path(file_path)

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
            loaded_data = file[dataset_name][dataset_slice]
        elif data_version == 2:
            # If the data is from Tristan v2 then perform whatever handling is needed.
            loaded_data = __handle_tristan_v2(file_path, file, dataset_name, dataset_slice)

        # Reduce to a scalar if the array is size 1
        if loaded_data.size == 1:
            return loaded_data.item()

        return loaded_data
# =============================================================================
