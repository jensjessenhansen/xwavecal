from configparser import ConfigParser
import tempfile
import os
import mock
from glob import glob
import pytest

from xwavecal.main import reduce_data

@pytest.mark.e2e
@mock.patch('xwavecal.fibers.IdentifyFibers.get_calibration_filename',
            return_value='xwavecal/tests/data/nres_test_data/cpt_nres03_20190405_0014_fibers_011.fits')
def test_reduce_data(mock_arc_template):
    # TODO this test should be replaced with a proper pytest fixture which makes the tempfile
    #  once per session. Then we can make independent tests for lampflat and arc creation.
    with tempfile.TemporaryDirectory() as temp_directory:
        args = type('test_args', (), {'input_dir': 'xwavecal/tests/data/nres_test_data/',
                                      'output_dir': temp_directory, 'fpack': False})
        config = ConfigParser()
        config.read('xwavecal/tests/data/test_config.ini')
        config.set('reduction', 'database_path', '"' + os.path.join(temp_directory, 'test.db') + '"')
        data_paths = glob('xwavecal/tests/data/nres_test_data/*w00*.fits*')
        reduce_data(data_paths, args=args, config=config)
        # check that the lampflat, traces and blaze are in the database via query for match
        # check that the correct number of traces exists.
        data_paths = glob('xwavecal/tests/data/nres_test_data/*a00*.fits*')
        reduce_data(data_paths, args=args, config=config)
        # check that the wavecals are in the database and query for match.
        # TODO implement benchmarks to check:
        # check that the wavelength solution for fiber 1 has 32 overlaps with peaks>6 and 16 marked as good.
        # check that the wavelength solution for fiber 2 has 30 overlaps with peaks>6 and 17 marked as good.
        # check that the mad for the wcs for fiber 1 (fiber 2) is near 0.0036 Ang (0.0054 Ang).
        assert True
