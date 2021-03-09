#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Parent class of KongsbergWaterColumn

import file_handler
import os

class WaterColumn:
    #TODO: See https://stackoverflow.com/questions/15836713/allowing-specific-values-for-an-argparse-argument
    # for possible ways to limit choices in command line.
    def __init__(self, input_file_or_directory, file_handler_type="snapshot_directory", udp_ip=None,
                 swath_width_percent=None, swath_width_m=None, num_nadir_beams=None, along_track_swaths=None,
                 bin_height_percent=None, bin_height_m=None, max_nav_gap_sec=0.04, write_directory=None):

        # Data input can come in one of 3 forms:
        # 1. A path to a single MBES file
        # 2. A path to a directory containing multiple MBES files
                # - A directory may be read in 2 ways:
                # a. 'snapshot' (default): A snapshot of directory contents at time this call is initiated.
                # b. 'monitor': Monitor directory for newly created files.
        # 3. A network connection to capture UDP datagrams
        if input_file_or_directory is not None:
            if os.path.isfile(input_file_or_directory):
                self.input_file_or_directory = input_file_or_directory
                self.file_handler_type = None
                self.udp_ip = None
            elif os.path.isdir(input_file_or_directory):
                self.input_file_or_directory = input_file_or_directory
                self.file_handler_type = file_handler_type
                self.udp_ip = None
        elif udp_ip is not None:
            self.input_file_or_directory = None
            self.file_handler_type = None
            self.udp_ip = udp_ip
        else:
            print("ERROR: Must enter file, directory, or UDP IP.")

        # # TODO: Determine which method of prioritization makes more sense:
        # # TODO: Method 1:
        # # 3 ways to determine across-track width to average.
        # #   [NOTE: "Center" is defined as vertical. Thus, the "center" beam will vary with roll.]
        # # If more than one is selected, they will be prioritized in the order given below:
        # #   1. Percentage of total swath width (EX: Enter 0.25 to average center 25% of beams.)
        # self.swath_width_percent = swath_width_percent
        # #   2. Across-track width (EX: Enter 10 to average center 10 meters of beams (+/- 5 meters from center).)
        # self.swath_width_m = swath_width_m
        # #   3. Number of beams (EX: Enter 50 to average center 50 beams.)
        # self.num_nadir_beams = num_nadir_beams
        # # If none are selected, default of 10% swath width will be used:
        # if swath_width_percent is None and swath_width_m is None and num_nadir_beams is None:
        #     self.swath_width_percent = 10

        # TODO: Method 2:
        if swath_width_percent is not None:
            self.swath_width_percent = swath_width_percent
            self.swath_width_m = None
            self.num_nadir_beams = None
        elif swath_width_m is not None:
            self.swath_width_percent = None
            self.swath_width_m = swath_width_m
            self.num_nadir_beams = None
        elif num_nadir_beams is not None:
            self.swath_width_percent = None
            self.swath_width_m = None
            self.num_nadir_beams = num_nadir_beams
        else:
            # TODO: Determine if 10% swath width is appropriate default.
            self.swath_width_percent = 10
            self.swath_width_m = None
            self.num_nadir_beams = None

        # Number of pings along-track to average (EX: Enter 4 to average 4 pings.)
        self.along_track_swaths = along_track_swaths

        # # TODO: Determine which method of prioritization makes more sense:
        # # TODO: Method 1:
        # # 2 ways to determine bin height.
        # # If more than one is selected, they will be prioritized in the order given below:
        # #   1. Percentage of maximum depth in swath (EX: Enter 0.1 for a bin height of 10% of swath's maximum depth.)
        # self.bin_height_percent = bin_height_percent
        # #   2. Height of bin in meters, starting at transducer (EX: Enter 0.25 for a bin height of 25 cm.)
        # self.bin_height_m = bin_height_m
        # # If none are selected, default of 1% water depth will be used:
        # if bin_height_percent is None and bin_height_m is None:
        #     self.bin_height_percent = 1

        # TODO: Method 2:
        if bin_height_percent is not None:
            self.bin_height_percent = bin_height_percent
            self.bin_height_m = None
        elif bin_height_m is not None:
            self.bin_height_percent = None
            self.bin_height_m = bin_height_m
        else:
            # TODO: Determine if 1% swath width is appropriate default.
            self.bin_height_percent = 1
            self.bin_height_m = None

        # Maximum acceptable gap in navigation data to interpolate through. 0.04 seconds (25 Hz) by default.
        self.max_nav_gap_sec = max_nav_gap_sec

        # Establish file handling method.
        if file_handler_type is "monitor":
            self.file_handler = file_handler.MonitorDirectory(self.input_file_or_directory)
        # elif file_handler_type is "write":
        #     self.file_handler = file_handler.WriteDirectory(self.input_file_or_directory)
        # SnapshotDirectory by default:
        else:
            self.file_handler = file_handler.SnapshotDirectory(self.input_file_or_directory)

        self.write_directory = write_directory

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