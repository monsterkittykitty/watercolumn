#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parent class of KongsbergWaterColumn

import file_handler

class WaterColumn:
    #TODO: See https://stackoverflow.com/questions/15836713/allowing-specific-values-for-an-argparse-argument
    # for possible ways to limit choices in command line.
    def __init__(self, input_file_or_directory, swath_width_percent=None, swath_width_m=None, num_nadir_beams=None,
                 along_track_swaths=None, bin_height_percent=None, bin_height_m=None,
                 max_nav_gap_sec=0.04, file_handler_type="snapshot_directory"):

        self.input_file_or_directory = input_file_or_directory

        # 3 ways to determine across-track width to average.
        #   [NOTE: "Center" is defined as vertical. Thus, the "center" beam will vary with roll.]
        # If more than one is selected, they will be prioritized in the order given below:
        #   1. Percentage of total swath width (EX: Enter 0.25 to average center 25% of beams.)
        self.swath_width_percent = swath_width_percent
        #   2. Across-track width (EX: Enter 10 to average center 10 meters of beams (+/- 5 meters from center).)
        self.swath_width_m = swath_width_m
        #   3. Number of beams (EX: Enter 50 to average center 50 beams.)
        self.num_nadir_beams = num_nadir_beams
        # If none are selected, default of 10% swath width will be used:
        if swath_width_percent is None and swath_width_m is None and num_nadir_beams is None:
            self.swath_width_percent = 10


        # Number of pings along-track to average (EX: Enter 4 to average 4 pings.)
        self.along_track_swaths = along_track_swaths

        # 2 ways to determine bin height.
        # If more than one is selected, they will be prioritized in the order given below:
        #   1. Percentage of maximum depth in swath (EX: Enter 0.1 for a bin height of 10% of swath's maximum depth.)
        self.bin_height_percent = bin_height_percent
        #   2. Height of bin in meters, starting at transducer (EX: Enter 0.25 for a bin height of 25 cm.)
        self.bin_height_m = bin_height_m
        # If none are selected, default of 1% water depth will be used:
        if bin_height_percent is None and bin_height_m is None:
            self.bin_height_percent = 1

        # Maximum acceptable gap in navigation data to interpolate through. 0.04 seconds (25 Hz) by default.
        self.max_nav_gap_sec = max_nav_gap_sec

        # Establish file handling method.
        if file_handler_type is "monitor_directory":
            self.file_handler = file_handler.MonitorDirectory(self.input_file_or_directory)
        elif file_handler_type is "write_directory":
            self.file_handler = file_handler.WriteDirectory(self.input_file_or_directory)
        # SnapshotDirectory by default:
        else:
            self.file_handler = file_handler.SnapshotDirectory(self.input_file_or_directory)

    def interpolate_linear(self, x1, y1, x2, y2, x3):
        """
        Uses linear interpolation between known points to determine coordinate of unknown point.
        :param x1: x1 of (x1, y1) pair, corresponding to first known point
        :param y1: y1 of (x1, y1) pair, corresponding to first known point
        :param x2: x2 of (x2, y2) pair, corresponding to second known point
        :param y2: y2 of (x2, y2) pair, corresponding to second known point
        :param x3: x3 of (x3, y3) pair, corresponding to independent variable of unknown point
        :return: y3 of (x3, y3) pair, corresponding to (previously) unknown point
        """
        y3 = y1 + ((x3 - x1) / (x2 - x1)) * (y2 - y1)
        return y3

    def basic_swath_stats(self):
        pass

    def data_structure(self):
        pass

    def plot(self):
        pass

    def run(self):
        pass