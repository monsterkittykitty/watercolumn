#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Inherits from class WaterColumn
# Contains Kongsberg-specific calculations

import getopt
import KMALL
import logging
import math
import sys
import watercolumn

# __name__ is module's name
logger = logging.getLogger(__name__)


class KongsbergWaterColumn(watercolumn.WaterColumn):

    def __init__(self, input_file_or_directory=None, file_handler_type="snapshot", udp_ip=None, swath_width_percent=None,
                 swath_width_m=None, num_nadir_beams=None, along_track_swaths=None, dual_swath=True,
                 bin_height_percent=None, bin_height_m=None, max_nav_gap_sec=0.04, write_directory=None):

        super().__init__(input_file_or_directory, file_handler_type, udp_ip, swath_width_percent, swath_width_m,
                       num_nadir_beams, along_track_swaths, bin_height_percent, bin_height_m, max_nav_gap_sec,
                         write_directory)

        # If dual swath mode is active, average pings from dual swath (defaults to True).
        # NOTE: along_track_swaths will take priority over dual_swath
        self.dual_swath = dual_swath

        self.print_settings()

    def extract_dg_timestamps(self, k, dg_offsets):
        """
        Extracts timestamps from datagram headers in .kmall file.
        :param: k: instance of class kmall with indexed file
        :param: dg_offsets: file offsets for desired datagram
        :return: array containing timestamps of all SKM datagrams in a file
        """

        time_dg_array = []
        datetime_dg_array = []

        for offset in dg_offsets:
            k.FID.seek(offset, 0)

            dg_header = k.read_EMdgmHeader()
            time_dg = dg_header['dgtime']
            time_dg_array.append(time_dg)

            datetime_dg = dg_header['dgdatetime']
            datetime_dg_array.append(datetime_dg)

        return time_dg_array, datetime_dg_array

    def extract_SKM_KMB_timestamps(self, k, SKM_offsets):
        """
        # TODO: Can this also work with .kmwcd? Can these have nav data?
        Extracts timestamps of SKM and KMB datagrams from .kmall file.
        :param: k: instance of class kmall with indexed file
        :param: dg_offsets: file offsets for desired datagram
        :return: array containing timestamps of all SKM datagrams in a file
        """
        time_SKM_array = []

        for offset in SKM_offsets:
            k.FID.seek(offset, 0)
            dg_SKM = k.read_EMdgmSKM()
            length_sample_array = dg_SKM['infoPart']['numSamplesArray']

            time_KMB_array = []

            for i in range(length_sample_array):
                time_KMB = (dg_SKM['sample']['KMdefault']['time_sec'][i] +
                            (dg_SKM['sample']['KMdefault']['time_nanosec'][i] / 1.0E9))
                time_KMB_array.append(time_KMB)

            # Error checking: Do timestamps of KMB datagrams always appear in ascending order or can they appear out of order?!
            sorted_time_KMB_array = sorted(time_KMB_array, key=float)
            if sorted_time_KMB_array != time_KMB_array:
                # TODO: If this prints we need to do something different with indexing/accessing KMB datagrams.
                #  This is used when comparing MWC timestamps to KMB position/attitude datagrams for interpolation...
                print("sorted_KMB_array: ", sorted_time_KMB_array)
                print("unsorted_KMB_array: ", time_KMB_array)
                #  Probably do something more robust here anyway.
                print("ERROR: KMB datagrams appear out of order.")
                exit()

            time_SKM_array.append(time_KMB_array)

        return time_SKM_array

    def basic_swath_stats(self):
        # Determine swath width and min, max, and avg depth:
        #       |         / \
        #       |       /     \
        #     y |     /         \
        #       |   /   (swath)   \
        #       |_/_________________\ _
        #                  x
        x_across_track_distance_array = []
        y_depth_array = []

        for beam in range(self.numBeams):  # 0 to self.numBeams - 1
            # Across-track angle:
            beam_point_angle_re_vertical = self.dg_MWC['beamData']['beamPointAngReVertical_deg'][beam]
            # Along-track angle:
            # Kongsberg: "Along ship steering angle of the TX beam (main lobe of transmitted pulse),
            # angle referred to transducer face. Angle as used by beamformer (includes stabilisation). Unit degree."
            # TODO: Determine tilt angle re vertical
            tilt_angle_re_tx_deg = self.tilt_angle_re_tx_deg[self.dg_MWC['beamData']['beamTxSectorNum'][beam]]



            # Index in sampleAmplitude05dB array where bottom detected:
            # Kongsberg: "Two way range in samples. Approximation to calculated distance from tx to bottom
            # detection [meters] = soundVelocity_mPerSec * detectedRangeInSamples / (sampleFreq_Hz * 2).
            # The detected range is set to zero when the beam has no bottom detection.
            # Replaced by detectedRangeInSamplesHighResolution for higher precision."
            detected_range = self.dg_MWC['beamData']['detectedRangeInSamples'][beam]

            range_to_wc_data_point = (self.soundVelocity * detected_range) / (self.sampleFreq * 2)

            # Across-track distance
            x = (range_to_wc_data_point * math.sin(math.radians(beam_point_angle_re_vertical))
                 * math.cos(math.radians(tilt_angle_re_tx_deg)))
            # Depth
            y = (range_to_wc_data_point * math.cos(math.radians(beam_point_angle_re_vertical))
                 * math.cos(math.radians(tilt_angle_re_tx_deg)))

            x_across_track_distance_array.append(x)
            y_depth_array.append(y)

        #print('minacrosstrack: ', min(x_across_track_distance_array))
        #print('maxacrosstrack: ', max(x_across_track_distance_array))
        #print('index max: ', x_across_track_distance_array.index(max(x_across_track_distance_array)))
        # TODO: If outer beam(s) have too-deep or too-shallow soundings, this could affect swath width calculation. How to make this more robust?
        swath_width = abs(min(x_across_track_distance_array)) + abs(max(x_across_track_distance_array))
        min_depth = min(y_depth_array)
        max_depth = max(y_depth_array)
        avg_depth = sum(y_depth_array) / len(y_depth_array)

        # TODO: Maybe make a PingStats or SwathStats class?
        return {'swath_width': swath_width, 'min_depth': min_depth, 'max_depth': max_depth, 'avg_depth': avg_depth}

#if __name__ == '__main__':
    # kwc = KongsbergWaterColumn(input_file_or_directory="/home/monster-kitty/Desktop/SharedFolder_ASV_III/0019_20190511_204630_ASVBEN.kmall")
    # print(kwc.average_percent_swath_width)

    def print_settings(self):
        print("Input file or directory: ", self.input_file_or_directory)
        print("Across-track averaging: ")
        print("\tPercent: ", self.swath_width_percent)
        print("\tMeters: ", self.swath_width_m)
        print("\tNadir beams: ", self.num_nadir_beams)
        print("Along-track averaging: ")
        print("\tSwaths: ", self.along_track_swaths)
        print("\tDual swath: ", self.dual_swath)
        print("Bin height: ")
        print("\tPercent depth: ", self.bin_height_percent)
        print("\tMeters: ", self.bin_height_m)
        print("Maximum acceptable navigation gap: ", self.max_nav_gap_sec, " seconds")

if __name__ == '__main__':

    # Defaults:
    _input_file_or_directory = None
    _udp_ip = None

    _swath_width_percent = None  # 10
    _swath_width_m = None
    _num_nadir_beams = None

    _along_track_swaths = None
    _dual_swath = False

    _bin_height_percent = None  # 1
    _bin_height_m = None

    _max_nav_gap_sec = None

    _write_directory = None

    # Read command line args for files/directory, num_nadir_beams:
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:p:m:n:s:d:b:a:g:",
                                   ["file=", "width_percent=", "width_m=", "num_beams=", "swaths=", "dual_swath=",
                                    "bin_percent=", "bin_m=", "nav_gap="])
    except getopt.GetoptError:
        print("kongsberg_watercolumnplot.py")
        # File or directory containing .kmall or .kmwcd file extensions
        print("-f <file_or_directory>")
        # Three options for selecting across-track width to include in watercolumn curtain:
        # 1. A percent of overall swath width, centered at nadir beam.
        # (Example: 10% of overall swath width.)
        print("-p <across_track_width_percent>")
        # 2. A specific number of meters, centered at nadir beam.
        # (Example: 10 meters gives 5 meters left and 5 meters right of nadir beam.)
        print("-m <across_track_width_m>")
        # 3. A specific number of beams, centered at nadir beam.
        # (Example: 20 beams.)
        print("-n <number_nadir_beams>")
        # TODO: This could also be done by along-track meters. Do we want to include this?
        # Two options for selecting along-track width to include in watercolumn curtain:
        # 1. A fixed number of pings.
        print("-s <number_along_track_swaths>")
        # 2. If in dual swath mode, average dual swath pings.
        print("-d <dual_swath>")
        # Two options for selecting vertical bin height:
        # TODO: Possibly this could be scaled/stretched for plotting purposes as needed and *always* use percentage of
        #  average nadir depth for each and every ping?
        # 1. Vertical bin height as a percent of average nadir depth of first ping. 1% by default.
        print("-b <bin_height_percent_depth>")
        # 2. Fixed, absolute vertical bin height. (EX: 0.25 m)
        print("-a <bin_height_m>")
        # Maximum acceptable gap in navigation to interpolate (linear) through.
        print("-g <max_nav_gap_seconds>")
        sys.exit()

    for opt, arg in opts:
        if opt == "-h":
            print("watercolumnplot.py\n")
            print("-f <file_or_directory>")
            print("\nDetermine across-track curtain width based on one of the following: ")
            print("-p <across_track_width_percent>")
            print("-m <across_track_width_m>")
            print("-n <number_nadir_beams>")
            print("\nDetermine along-track curtain width based on one of the following: ")
            print("-s <number_along_track_swaths>")
            # TODO: Is dual swath Kongsberg-specific? If not, maybe all this can go in watercolumn.py?
            print("-d <dual_swath>")
            print("\nDetermine vertical bin size based on one of the following: ")
            print("-b <bin_height_percent_depth>")
            print("-a <bin_height_m>")
            print("-g <max_nav_gap_seconds>")
            sys.exit()
        elif opt in ("-f", "--file"):
            _input_file_or_directory = arg
        elif opt in ("-p", "--width_percent"):
            _swath_width_percent = int(arg)
        elif opt in ("-m, --width_m"):
            _swath_width_m = int(arg)
        elif opt in ("-n", "--num_beams"):
            _num_nadir_beams = int(arg)
        elif opt in ("-s", "--swaths"):
            _along_track_swaths = int(arg)
        elif opt in ("-d", "--dual_swath"):
            _dual_swath = True
        elif opt in ("-b, --bin_percent"):
            _bin_height_percent = int(arg)
        elif opt in ("-a, --bin_m"):
            _bin_height_m = int(arg)
        elif opt in ("g", "--nav_gap"):
            _max_nav_gap_sec = int(arg)

    if _input_file_or_directory is None:
        print("Must enter file or directory: watercolumnplot.py -f <file_or_directory>")
        sys.exit()

    kongsberg_wc = KongsbergWaterColumn(_input_file_or_directory, "snapshot", _udp_ip, _swath_width_percent, _swath_width_m,
                                        _num_nadir_beams, _along_track_swaths, _dual_swath, _bin_height_percent,
                                        _bin_height_m, _max_nav_gap_sec, _write_directory)
