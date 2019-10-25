import numpy as np
from astropy.io import fits
from astropy.table import table


def spectrum_table_generator(order_fluxes, order_stderrors, extension_name):
    """Creates astropy table with the column format expected by xwavecal and converts it to a fits hdu"""

    table_column_dict = {'id': [], 'flux': [], 'stderr': [], 'pixel': [], 'fiber': [], 'ref_id': [], 'wavelength': []}
    for order_id, (order_flux, order_stderr) in enumerate(zip(order_fluxes, order_stderrors)):
        table_column_dict['flux'].append(order_flux)
        table_column_dict['stderr'].append(order_stderr)
        table_column_dict['pixel'].append(np.arange(len(order_flux)))
        table_column_dict['id'].append(order_id)
        table_column_dict['wavelength'].append(np.ones_like(order_flux, dtype=float) * np.nan)
        table_column_dict['ref_id'].append(np.nan)
        table_column_dict['fiber'].append(0)

    tab = table.Table(table_column_dict)
    fits_hdu = fits.table_to_hdu(tab)
    fits_hdu.name = extension_name
    return fits_hdu


def fitscreator(order_fluxes, order_stderrors, deblazed_fluxes, deblazed_stderrors, raw_image, header_keys):
    """Creates a fits file that has the format expected by xwavecal"""

    spectrum_hdu = spectrum_table_generator(order_fluxes, order_stderrors, 'SPECBOX')
    blzcorr_hdu = spectrum_table_generator(deblazed_fluxes, deblazed_stderrors, 'BLZCORR')

    primary_hdu = fits.PrimaryHDU(raw_image)
    out_frame = fits.HDUList()
    out_frame.append(primary_hdu)
    out_frame.append(spectrum_hdu)
    out_frame.append(blzcorr_hdu)

    for key, value in header_keys.items():
        out_frame[0].header[key] = value

    return out_frame


if __name__ == '__main__':
    thar_file = fits.open('thar_s1_2018-10-27T17-48-01.fits')
    blaze_file = fits.open('mastercont_20180322_slit6.fits')

    output_fits_name = 'test.fits'

    #  Minimal set of header keywords needed for xwavecal to run
    header_keys = {'OBSTYPE': 'wavecal',
             'GAIN': 1.,
             'RDNOISE': 9.5,
             'OBJECTS': 'thar&none&none',
             'DATE-OBS': '2018-10-27T17:48:01.461717',
             'TELESCOP': 'SONG1',
             'INSTRUME': 'Spectrograph',
             'SITEID': 'Teide',
             'FRAMENUM': 'some_unique_id'}

    raw_image = np.zeros((2048, 2048))

    #  Create list containing flux and stderror arrays (orders need to be increasingly bluer)
    order_fluxes, order_stderrors, deblazed_fluxes, deblazed_stderrors = [], [], [], []
    for i in range(20, 0, -1):  # reversed order as the spectral orders in SONG extracted spectrum files run from blue to red
        flux = thar_file[0].data[1, i, :]
        deblazed_flux = flux / blaze_file[0].data[i, :]

        order_fluxes.append(flux)
        order_stderrors.append(np.sqrt(flux))
        deblazed_fluxes.append(deblazed_flux)
        deblazed_stderrors.append(np.sqrt(flux) / blaze_file[0].data[i, :])

    thar_file.close()
    blaze_file.close()

    fitsfile = fitscreator(order_fluxes, order_stderrors, deblazed_fluxes, deblazed_stderrors, raw_image, header_keys)
    fitsfile.writeto(output_fits_name)
    fitsfile.close()
