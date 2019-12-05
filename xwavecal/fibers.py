import numpy as np
from copy import deepcopy

from xwavecal.utils.misc_utils import normalize_by_brightest, n_largest_elements
from xwavecal.utils.correlate import correlate2d
from xwavecal.utils.fiber_utils import lit_wavecal_fibers, lit_fibers
from xwavecal.stages import ApplyCalibration, Stage
from xwavecal.images import DataProduct

import logging
logger = logging.getLogger(__name__)


class MakeFiberTemplate(Stage):
    def __init__(self, runtime_context=None):
        super(MakeFiberTemplate, self).__init__(runtime_context=runtime_context)

    @property
    def calibration_type(self):
        return 'FIBERS'

    def do_stage(self, image):
        order = self.runtime_context.template_trace_id
        logger.info('Generating a template with center order {0}'.format(order))
        spec = image.data_tables[self.runtime_context.main_spectrum_name]
        num_fibers = len(lit_fibers(image))
        orders = [order - num_fibers, order, order + num_fibers]
        template_region = np.logical_or.reduce(tuple(spec['id'] == i for i in orders))
        template = DataProduct(data=spec[template_region], header=deepcopy(image.header),
                               data_name=self.calibration_type.lower(), translator=image.translator)
        template.set_header_val('type', self.calibration_type.lower())
        template.fiber0_lit, template.fiber1_lit, template.fiber2_lit = image.fiber0_lit, image.fiber1_lit, image.fiber2_lit
        return image, template


class IdentifyFibers(ApplyCalibration):
    def __init__(self, runtime_context=None):
        super(IdentifyFibers, self).__init__(runtime_context=runtime_context)

    @property
    def calibration_type(self):
        return 'FIBERS'

    def apply_master_calibration(self, image, template_path):
        """
        :param image: DataProduct
        :param template_path: path to the arc lamp template which will be used to identify the fibers
        :return: image. The spectrum data_table has two columns appended to it: fiber and ref_id.

        NOTE: We follow the convention that the lowest lying ThAr fiber (i.e. has smallest x and y coordinate
        on the detector) corresponds to the smallest lit fiber index. E.g. if 0, 1 are lit with ThAr, then the lowest
        lying ThAr fiber is denoted 0. If 1,2 are lit, then the lowest lying ThAr fiber is denoted 1.
        ref_ids must be such that the mean wavelength of an order always increases with increasing ref_id

        NOTE: This stage relies on the fact that only 2 fibers are lit at one time.
        """
        if len(image.data_tables[self.runtime_context.main_spectrum_name]['flux']) == 0:
            logger.error('Image has length 0 spectrum. Skipping fiber identification')
        elif image.num_wavecal_fibers() == 0:
            logger.info('Image does not have any fibers lit with ThAr, skipping fiber identification.', )
        else:
            logger.info('Identifying fibers via cross correlation.')
            spectrum = image.data_tables[self.runtime_context.main_spectrum_name]

            read_noise = image.get_header_val('read_noise')
            signal_to_noise = spectrum['flux'] / (np.sqrt(np.abs(spectrum['flux']) + read_noise ** 2))
            spectrum_to_search = normalize_by_brightest(signal_to_noise)

            template = self.construct_single_fiber_template(template_path,
                                                            self.calibration_type.lower(),
                                                            num_lit_fibers=image.num_lit_fibers())
            matched_ids = self.identify_matching_orders(spectrum_to_search, template,
                                                        num_arc_lamp_fibers=image.num_wavecal_fibers())

            fiber_ids = self.build_fiber_column(matched_ids, image, spectrum)
            ref_ids = self.build_ref_id_column(matched_ids, fiber_ids, self.runtime_context.ref_id)
            # TODO refactor storing information about which spectra exist.
            image.set_header_val('IDTEMPL', (template_path, 'ID of the fiber template.'))
            for key in [self.runtime_context.main_spectrum_name, self.runtime_context.blaze_corrected_spectrum_name]:
                if image.data_tables.get(key) is not None:
                    image.data_tables[key]['fiber'] = fiber_ids
                    image.data_tables[key]['ref_id'] = ref_ids
        return image

    @staticmethod
    def build_fiber_column(matched_ids, image, spectrum, low_fiber_first=True):
        """
        :param matched_ids: array: the indices of spectrum who are spectrally matched to ref_id
        :param image: xwavecal.images.DataProduct
        :param spectrum: astropy.tables.Table: the extracted spectrum.
        :param low_fiber_first: whether the lowest lying fiber (lowest in terms of row index in spectrum)
        should be associated with the first lit fiber. For example:
           if lit_wavecal_fibers(image) = [1, 2] (fibers 1 and 2 are lit with a wavelength calibration lamp)
           and matched_ids = [10, 11] (the rows 10 and 11 in the spectrum have been matched to the template)
           then if label_low_fiber:
              spectrum[10]['fiber_id'] will equal 1, and spectrum[11]['fiber_id'] will equal 2
           If not label_low_fiber:
              spectrum[10]['fiber_id'] will equal 2, and spectrum[11]['fiber_id'] will equal 1

        :return: array:
                 array of fiber identification numbers such that spectrum[fiber_ids == 1] would yield a complete
                 extracted spectrum for the 1 fiber on the echelle spectrograph. spectrum[fiber_ids == 2] would
                 be the spectrum for fiber 2.
        """
        fiber_sequence = lit_fibers(image) if low_fiber_first else lit_fibers(image)[::-1]
        match_fiber = lit_wavecal_fibers(image)[0] if low_fiber_first else lit_wavecal_fibers(image)[-1]
        num_lit = len(lit_fibers(image))
        sequence_ids = np.arange(len(spectrum['id'])) % num_lit
        sequence_ids += np.where(fiber_sequence == match_fiber)[0] - sequence_ids[matched_ids[0]]
        return fiber_sequence[sequence_ids % num_lit]

    @staticmethod
    def build_ref_id_column(matched_ids, fiber_ids, ref_id):
        """
        :param matched_ids: array: the indices of spectrum who are spectrally matched to ref_id
        :param ref_id: the reference id which is to be assigned to those orders which have been spectrally
        matched with the template and whose row id is given in matched_ids.
        :param fiber_ids: array: integers which designate to each row of spectrum, a fiber id. See
        IdentifyFibers.build_fiber_column().
        :return: array:
                 array of reference ids (order indices), each entry corresponding
                 to each row of the spectrum, to be used in 1/(m0+i) of the grating equation.
                 We assume that spectrum is such that spectrum[i] is bluer than spectrum[j] if i > j.
                 Thus, ref_ids[i] > ref_ids[j] (where ref_ids[i] is the reference id of the single order spectrum[i]).
                 This array is such that ref_ids[matched_ids[k]] = ref_id for all valid indices k
                 of matched_ids.
        """
        ref_ids = np.ones_like(fiber_ids)
        return ref_ids

    @staticmethod
    def construct_single_fiber_template(template_path, ext_name, num_lit_fibers=2):
        template = DataProduct.load(template_path, extension_name=ext_name)
        single_fiber = template.data['flux'].data
        template = np.zeros((num_lit_fibers * (single_fiber.shape[0] - 1) + 1, single_fiber.shape[1]))
        template[::num_lit_fibers] = single_fiber
        return template

    @staticmethod
    def identify_matching_orders(two_d_spectrum, template, num_arc_lamp_fibers=2):
        order_likelyhood = correlate2d(two_d_spectrum, template, max_lag=100)
        one_d_likelyhood = np.max(order_likelyhood, axis=1)
        matched_ids = n_largest_elements(one_d_likelyhood, n=num_arc_lamp_fibers)
        return np.array(matched_ids, dtype=np.int)
