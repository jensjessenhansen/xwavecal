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
        elif image.num_wavecal_fibers() >= 1:
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
            ref_ids = self.build_ref_id_column(matched_ids, image, spectrum, self.runtime_context.ref_id)
            # TODO refactor storing information about which spectra exist.
            image.set_header_val('IDTEMPL', (template_path, 'ID of the fiber template.'))
            for key in [self.runtime_context.main_spectrum_name, self.runtime_context.blaze_corrected_spectrum_name]:
                if image.data_tables.get(key) is not None:
                    image.data_tables[key]['fiber'] = fiber_ids
                    image.data_tables[key]['ref_id'] = ref_ids
        else:
            logger.info('Image does not have any fibers lit with ThAr, skipping fiber identification.', )
        return image

    @staticmethod
    def build_fiber_column(matched_ids, image, spectrum):
        # note this assumes that 2 fibers are lit. other_fiber = lit... will fail if only one is lit.
        fiber_ids = np.ones_like((spectrum['id']), dtype=np.int) * (-1)
        fiber_ids[min(matched_ids) % 2:: 2] = min(lit_wavecal_fibers(image))
        other_fiber = lit_fibers(image)[lit_fibers(image) != min(lit_wavecal_fibers(image))][0]
        fiber_ids[(min(matched_ids) + 1) % 2:: 2] = other_fiber
        return fiber_ids

    @staticmethod
    def build_ref_id_column(matched_ids, image, spectrum, ref_id):
        # note this assumes that 2 fibers are lit. other_fiber = lit... will fail if only one is lit.
        lowest_lying_arc = min(lit_wavecal_fibers(image))
        ref_ids = np.ones_like((spectrum['id']))
        ref_ids[1::2], ref_ids[::2] = np.arange(len(ref_ids[1::2])), np.arange(len(ref_ids[::2]))
        ref_ids[min(matched_ids) % 2:: 2] += ref_id - ref_ids[min(matched_ids)]
        other_fiber = lit_fibers(image)[lit_fibers(image) != lowest_lying_arc][0]
        other_fiber_match = min(matched_ids) + other_fiber - lowest_lying_arc
        ref_ids[other_fiber_match % 2:: 2] += ref_id - ref_ids[other_fiber_match]
        return ref_ids

    @staticmethod
    def construct_single_fiber_template(template_path, ext_name, num_lit_fibers=2):
        template = DataProduct.load(template_path, extension_name=ext_name)
        single_fiber = template.data['flux'].data
        template = np.zeros((num_lit_fibers * single_fiber.shape[0] - 1, single_fiber.shape[1]))
        template[::num_lit_fibers] = single_fiber
        return template

    @staticmethod
    def identify_matching_orders(two_d_spectrum, template, num_arc_lamp_fibers=2):
        order_likelyhood = correlate2d(two_d_spectrum, template, max_lag=100)
        one_d_likelyhood = np.max(order_likelyhood, axis=1)
        matched_ids = n_largest_elements(one_d_likelyhood, n=num_arc_lamp_fibers)
        return np.array(matched_ids, dtype=np.int)