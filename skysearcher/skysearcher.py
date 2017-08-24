'''
===========================================================================
skysearcher v0.1.0
August 2017
Sol W. Courtney
swc2124@columbia.edu
Columbia University NYC, NY
===========================================================================
python skysearcher/skysearcher.py

This is the main program and the only script intended to be run. 
All settings are in the rc.cfg file.  If there is not an rc.cfg 
file, then run this program and one will be written.

All values that are set by the configuration file are marked 
like this: <<ss_lib.variable_value>>.  All values are defined 
in the ss_lib so this script stays clean and readable.

GitHub:
https://github.com/swc2124/skysearcher.git

rc.cfg values --> <<ss_lib.variable_value>> 
'''

from __future__ import division, absolute_import, print_function

import skysearch_lib as ss_lib

# Make sure all the directories are in place.
ss_lib.sortout_directories()

# Load the halo's grid & table (grid, table).
grid, table = ss_lib.load_halo()

# Make an empty dictionary to keep region information in.
feature_dict = {}

# An integer value that seres as an ID number for each feature
# in the "region_dict."
feature_id = 0

# --------------------------------------------------------------------------
# Begin moving outward from the halo projection's center.
#
# "annulus_id"       - An integer increasing from 0 for each radial step
#                     outward from center of halo projection.

# "annulus"          - An integer value representing the distance from
#                     the halo projection's center in unites set by
#                     the rc.cfg file.

# <<ss_lib.r_start>> - An integer value representing the starting
#                     radius (same units).

# <<ss_lib.r_stop>>  - An integer value representing the ending
#                     radius(same units).

# "annuli"           - An array of annular segments radially separated.
annuli = xrange(ss_lib.r_start, ss_lib.r_stop, 1)

for annulus_id, annulus in enumerate(annuli):

    # -----------------------------------------------------------------------
    # Determine the evaluation region's RADIAL EXTENT based on the values of:

    # <<ss_lib.annulus_scale>>        - A float representing the percentage
    #                                   of the current radius used to
    #                                   determine the evaluation area
    #                                   at this radius.

    # The percentage value multiplied with the current radius (region).
    region = annulus * ss_lib.annulus_scale

    # This is the evaluation region's starting radius (r_in).
    r_in = annulus - region

    # This is the evaluation region's ending radius (r_out).
    r_out = annulus + region

    # -----------------------------------------------------------------------
    # Determine the evaluation region's ANNULAR EXTENT based on the values of:

    # <<ss_lib.units_annulus>>        - The working unit for annular movement.
    #                                   Later on more units can be added here but
    #                                   for now there's just "degree."  This means
    #                                   that the input values best range
    #                                   from 0.0 - 360.0!

    if ss_lib.units_annulus == 'degree':

        # The starting degree of this radius' annulus (phi_start).
        phi_start = -180.0

        # The ending degree (phi_stop).
        phi_stop = 180.0

        # <<ss_lib.annulus_phi_step>> - For units of degrees, this is the number
        #                               we dived 360.0 by to get the evaluation
        #                               region's segment size.  The larger this
        #                               number, the smaller the segment.  How do
        #                               smaller or larger segment sizes effect
        #                               this algorithm? (float)

        # The size of the annular segment to evaluate (phi_step). If
        # <<ss_lib.annulus_phi_step>> = 1.0, then each evaluation regions
        # annular segment extent will be 1.0 degree.  If however,
        # <<ss_lib.annulus_phi_step>> = 2.0, the the extent will be 1/2 the
        # previous extent.  So the larger the <<ss_lib.annulus_phi_step>> the
        # smaller the evaluation region. (float)
        phi_step = 360.0 // ss_lib.annulus_phi_step

    # "annular_segments"              - Array of segment values to be iterated over
    #                                   as to complete one rotation around the given
    #                                   radius' annulus.

    # "annular_step"                  - The magnitude of a single annular step.

    annular_segments, annular_step = np.linspace(
        start=phi_start, stop=phi_stop, num=phi_step,
        endpoint=True, retstep=True, dtype=np.float16)  # np.float16 to save space.

    # "annular_segment"               - This is a single annular step segment value
    #                                   from the array "annular_segments".
    #                                   (np.float16)

    for annular_segment in annular_segments:

        # "evaluation_region_limits"  - An array of indices's.  These are the indices's
        #                               of the stars contained within the evaluation
        #                               region's limits.
        evaluation_region_limits = np.nonzero(
            np.logical_and(
                np.logical_and(
                    table['Degree'] >= annular_segment - annular_step,
                    table['Degree'] < annular_segment + annular_step),
                np.logical_and(
                    table['Rads'] >= r_in,
                    table['Rads'] < r_out)))

        # "n_boxes"                   - This is an integer value for the number of grid
        #                               spaces contained within the region.
        #                               i.e len(evaluation_region_limits)
        n_boxes = len(evaluation_region_limits[0])

        # "points"                    - An array of tuples. Each tuple consists of a
        #                               single grid space's x, y coordinates.  This
        #                               is used to be sure that this grid space is not
        #                               already part of another region.
        points = zip(
            (table['x_int'][evaluation_region_limits]).tolist(),
            (table['y_int'][evaluation_region_limits]).tolist())

        # "satids"                    - A list of satellite id numbers for this halo.
        satids = table['satids'][evaluation_region_limits]

        # "xbox_values"               - An array of xbox values for the regions
        #                               grid spaces.
        xbox_values = table['Xbox'][evaluation_region_limits]

        # if there are grid spaces within this region:
        # n_boxes > 0
        if n_boxes:

            # STEP 1.0
            # Make sure that the feature id is in the "feature_dict."
            # ------------------------------------------------------------------
            if not feature_id in feature_dict.keys():

                # make an empty dictionary for the halo's
                # satellite information.
                sats_book = {}

                # Walk through each of the halo's satids
                # and use them as keys for this dictionary.
                for satid in table.meta['satids']:

                    # set each value to 0.
                    sats_book[satid] = 0

                # Bump the feature ID index up by 1.
                feature_id += 1

                # Add a new feature object to "feature_dict."
                feature_dict[feature_id] = new_feature(
                    r_in,
                    r_out,
                    round(annular_segment, 5),
                    round(annular_segment + annular_step, 5),
                    ss_lib.halo,
                    feature_id,
                    points,
                    xbox_values,
                    count_satids(satids, sats_book))

            # STEP 2.0
            # Figure out if this feature is known or current.
            # ------------------------------------------------------------------
            # "current_feature"       - A boolean value.
            #                           True if evaluation_region_limits
            #                           belong to the current feature.
            current_feature = False

            # "known_feature"         - A boolean value.
            #                           True if evaluation_region_limits
            #                           appear in another feature.
            known_feature = False

            # Check the set of points against the points of all
            # known features one at a time.
            for known_feat_id in feature_dict.keys():

                # Designate this feature.
                _feature = feature_dict[known_feat_id]

                # check each set of (x, y) points in "points" against
                # this features points
                for point in set(points):

                    # If features share a common point, then they
                    # are the same feature.
                    if point in feature['points']:

                        # Set "known_feature" to True & break loop.
                        known_feature = True
                        break

                # Break outer loop.
                if current_feature:
                    break

            # STEP 3.0
            # Make sure that the feature id is in the "feature_dict."
            # ------------------------------------------------------------------

        # if there are no grid spaces:
        # i.e n_boxes == 0.
        else:

            pass

    # [Do something here]

    # --------------------------------------------------
    # End of "for annular_segment in annular_segments:"
    # No more annular segments.
    # End of annular evaluation.
    # --------------------------------------------------

# [Do something here]

# ------------------------------------------------------
# End of "for annulus_id, annulus in enumerate(annuli):"
# No more radial segments.
# End of evaluation.
# ------------------------------------------------------
