#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Extract and plot nadir watercolumn beams.
# When run as main, defaults to select 10% swath width for nadir curtain averaging.
# Can specify percent swath width, number of nadir beams, or across-track width to include in curtain averaging.

import copy
import getopt
import KMALL
import logging
import math
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

#sys.path.append("../PycharmProjects")

# __name__ is module's name
logger = logging.getLogger(__name__)


class WaterColumnPlot:

    def __init__(self, input_file=None, percent_swath_width=None, num_nadir_beams=None,
                 across_track_width=None, bin_height_percent_depth=None):
        self.files = self.create_file_list(input_file)
        self.percent_swath_width = percent_swath_width
        self.num_nadir_beams = num_nadir_beams
        self.across_track_width = across_track_width
        self.bin_height_percent_depth = bin_height_percent_depth

        # Maximum acceptable gap in navigation data to interpolate through (seconds), default 25 Hz (0.04 sec):
        # TODO: Make this command-line configurable:
        self.MAX_NAV_GAP_SEC = 0.04

    # Create file list if input path is directory
    def create_file_list(self, input_file):
        temp_list = []

        if os.path.isfile(input_file):
            if input_file.lower().endswith(('.kmall', '.kmwcd')):
                temp_list.append(input_file)
        elif os.path.isdir(input_file):
            for root, dirs, files in os.walk(input_file):
                for filename in files:
                    if filename.lower().endswith(('.kmall', '.kmwcd')):
                        temp_list.append(os.path.join(root, filename))

        if len(temp_list) > 0:
            return temp_list
        else:
            logger.warning("Invalid file path: %s" % input_file)
            sys.exit(1)


    def extract_dg_timestamps(self, k, dg_offsets):
        """
        Extracts timestamps of SKM datagrams in .kmall file.
        :param: k: instance of class kmall with indexed file
        :param: dg_offsets: file offsets for desired datagram
        :return: array containing timestamps of all SKM datagrams in a file
        """

        time_dg_array = []
        datetime_dg_array = []

        for offset in dg_offsets:
            k.FID.seek(offset, 0)

            dg_header = k.read_EMdgmHeader()
            #time_sec = dg_header['time_sec']
            #time_nanosec = dg_header['time_nanosec']
            #time_dg = time_sec + (time_nanosec / 1.0E9)
            time_dg = dg_header['dgtime']
            time_dg_array.append(time_dg)

            datetime_dg = dg_header['dgdatetime']
            datetime_dg_array.append(datetime_dg)

            # self.dg_SKM = k.read_EMdgmSKM()
            # time_sec_SKM = self.dg_SKM['header']['time_sec']
            # time_nanosec_SKM = self.dg_SKM['header']['time_nanosec']
            # time_SKM = time_sec_SKM + (time_nanosec_SKM * (10 ** -9))
            # time_SKM_array.append(time_SKM)

        return time_dg_array, datetime_dg_array

    def extract_SKM_KMB_timestamps(self, k, SKM_offsets):
        time_SKM_array = []

        for offset in SKM_offsets:
            k.FID.seek(offset, 0)
            dg_SKM = k.read_EMdgmSKM()
            length_sample_array = dg_SKM['infoPart']['numSamplesArray']

            time_KMB_array = []

            for i in range(length_sample_array):
                time_KMB = (dg_SKM['sample']['KMdefault']['time_sec'][i] +
                            (dg_SKM['sample']['KMdefault']['time_nanosec'][i] / 1.0E9))
                # Pair #KMB timestamp with its index incase they appear in the file out of order.
                time_KMB_array.append([time_KMB, i])



            # Error checking: Do timestamps of KMB datagrams always appear in ascending order or can they appear out of order?!
            # sorted_time_KMB_array = sorted(time_KMB_array, key=float)
            # if sorted_time_KMB_array != time_KMB_array:
            #     # TODO: If this prints we need to do something different with indexing/accessing KMB datagrams.
            #     #  This is used when comparing MWC timestamps to KMB position/attitude datagrams for interpolation...
            #     print("sorted_KMB_array: ", sorted_time_KMB_array)
            #     print("unsorted_KMB_array: ", time_KMB_array)
            #     #  Probably do something more robust here anyway.
            #     print("ERROR: KMB datagrams appear out of order.")
            #     exit()

            # Sort #KMB array based on timestamps before adding to #SKM array.
            time_SKM_array.append(time_KMB_array.sort(key=lambda x: x[0]))

        return time_SKM_array

    def interpolate(self, x1, y1, x2, y2, x3):
        """
        Linear interpolation.
        """
        y3 = y1 + ((x3 - x1) / (x2 - x1)) * (y2 - y1)
        return y3

    def sort_dataframe_by_index(self, df, datagram_type=None):
        """
        Sort pandas dataframe by index value; optionally filter by datagram type.
        :param df: Pandas dataframe.
        :param datagram_type: Kongsberg datagram type. EX: #MRZ, #MWC...
        :return: Pandas dataframe sorted by index value; optionally filtered by datagram type.
        """
        if datagram_type is None:
            return df.sort_index()
        else:
            filtered_df = df[df['MessageType'].map(lambda x: x.strip("'b")) == datagram_type]
            return filtered_df.sort_index()

    def extract_wc_from_file(self):
        # TODO: Test this with .kmwcd files

        for fp in self.files:
            k = KMALL.kmall(fp)
            k.index_file()

            # Extract only #MWC and #SKM datagrams from k.Index and sort by timestamp:
            sorted_df_MWC = self.sort_dataframe_by_index(k.Index, "#MWC")
            sorted_df_SKM = self.sort_dataframe_by_index(k.Index, "#SKM")

            # Extract #MWC and #SKM file offsets:
            MWCOffsets = sorted_df_MWC['ByteOffset'].tolist()
            SKMOffsets = sorted_df_SKM['ByteOffset'].tolist()

            timestamps_SKM_KMB = self.extract_SKM_KMB_timestamps(k, SKMOffsets)

            #count = 0
            fullFile_avergedWaterColumnArray = []
            full_file_samples_per_bin_array = []
            full_file_beam_number_array = []

            #for offset in MWCOffsets:  # For each MWC datagram:
            # for offset in MWCOffsets[0:4]:  # Only two pings to evaluate dual swath
            # TODO: FOR TESTING: (Use only a subset of MWC records:)
            for offset in MWCOffsets[0:(int((len(MWCOffsets) / 6)))]:  # Each MWCOffset corresponds to a ping
                k.FID.seek(offset, 0)  # Find start of MWC datagram
                self.dg_MWC = k.read_EMdgmMWC()  # Read MWC datagram

                # Extract datagram timestamp:
                time_MWC_array = []
                datetime_MWC_array = []

                time_MWC = self.dg_MWC['header']['dgtime']
                time_MWC_array.append(time_MWC)

                # TODO: Not sure if I need this:
                datetime_MWC = self.dg_MWC['header']['dgdatetime']
                datetime_MWC_array.append(datetime_MWC)

                # Find corresponding SKM datagram:
                SKM_index1 = None
                SKM_index2 = None
                # TODO: probably a better way to search? Binary search?
                for i in range(len(timestamps_SKM_KMB) - 1):
                    #if time_MWC >= time_SKM_array[i] and time_MWC < time_SKM_array[i + 1]:
                    if timestamps_SKM_KMB[i][0] <= time_MWC < timestamps_SKM_KMB[i + 1][0]:
                        SKM_index1 = SKM_index2 = i
                        break

                if SKM_index1 is None:
                    print("ERROR finding matching SKM datagram.")
                    exit()

                KMB_index1 = None
                KMB_index2 = None
                # TODO: probably a better way to search? Binary search?
                for i in range(len(timestamps_SKM_KMB[SKM_index1]) - 1):
                    if timestamps_SKM_KMB[SKM_index1][i] <= time_MWC < timestamps_SKM_KMB[SKM_index1][i + 1]:
                        # Ensure gap to interpolate through is not too big:
                        if ((timestamps_SKM_KMB[SKM_index1][i + 1] - timestamps_SKM_KMB[SKM_index1][i]) < self.MAX_NAV_GAP_SEC):
                            KMB_index1 = i
                            KMB_index2 = i + 1
                        else:
                            print("ERROR: Nav gap too large to interpolate!")
                        break

                if KMB_index1 is None:
                    # This is the case where the MWC timestamp falls between SKM datagrams:
                    if (timestamps_SKM_KMB[SKM_index1][len(timestamps_SKM_KMB[SKM_index1]) - 1]
                            <= time_MWC < timestamps_SKM_KMB[SKM_index1 + 1][0]):
                        SKM_index2 = SKM_index1 + 1
                        KMB_index1 = len(timestamps_SKM_KMB[SKM_index1]) - 1
                        KMB_index2 = 0

                    else:
                        # TODO: Handle this in some elegant way. Drop that ping and move on to the next one or something...
                        print("time_MWC: ", time_MWC)
                        print("timestamps_SKM_KMB[SKM_index1][len(timestamps_SKM_KMB[SKM_index1]) - 2]: ",
                              timestamps_SKM_KMB[SKM_index1][len(timestamps_SKM_KMB[SKM_index1]) - 2])
                        print("timestamps_SKM_KMB[SKM_index1][len(timestamps_SKM_KMB[SKM_index1]) - 1]: ",
                              timestamps_SKM_KMB[SKM_index1][len(timestamps_SKM_KMB[SKM_index1]) - 1])
                        print("timestamps_SKM_KMB[SKM_index1 + 1][0]: ",
                              timestamps_SKM_KMB[SKM_index1 + 1][0])
                        print("ERROR finding matching KMB datagram.")
                        exit()

                print("SKM_index1: ", SKM_index1)
                print("KMB_index1: ", KMB_index1)

                # TODO: Not sure if it's better to re-access the file or to store all pitch values in an array?
                # Interpolate pitch. (x,y) points are (time, pitch):
                k.FID.seek(SKMOffsets[SKM_index1], 0)
                dg_SKM = k.read_EMdgmSKM()
                x1 = timestamps_SKM_KMB[SKM_index1][KMB_index1]
                y1 = dg_SKM['sample']['KMdefault']['pitch_deg'][KMB_index1]
                x2 = timestamps_SKM_KMB[SKM_index2][KMB_index2]
                y2 = dg_SKM['sample']['KMdefault']['pitch_deg'][KMB_index2]
                interpolated_pitch = self.interpolate(x1, y1, x2, y2, time_MWC)

                # print("y1, pitch: ", y1)
                # print("x1, time: ", x1)
                # print("y2, pitch: ", y2)
                # print("x2, time: ", x2)
                # print("y3, pitch: ", interpolated_pitch)
                # print("x3, time: ", time_MWC)



                # Kongsberg: "Number of swaths per ping. A swath is a complete set of across track data.
                # A swath may contain several transmit sectors and RX fans"
                self.swaths_per_ping = self.dg_MWC['cmnPart']['swathsPerPing']
                # Kongsberg: "Number of transmitting sectors (Ntx). Denotes the number of times
                # the struct EMdgmMWCtxSectorData is repeated in the datagram."
                self.num_tx_sectors = self.dg_MWC['txInfo']['numTxSectors']
                # Kongsberg: Along ship steering angle of the TX beam (main lobe of transmitted pulse),
                # angle referred to transducer face. Angle as used by beamformer (includes stabilisation). Unit degree.
                # TODO: Will need determine tilt angle relative to vertical
                # TODO: This is length 3 in testing with dual swath
                self.tilt_angle_re_tx_deg = self.dg_MWC['sectorData']['tiltAngleReTx_deg']

                # TODO: Create function to interpolate attitude data
                # TODO: Create a function to use interpolated attitude data to determine attitude at time of ping.









                # Kongsberg: "Heave at vessel reference point, at time of ping, i.e. at midpoint of first tx pulse in rxfan"
                # TODO: Will need offsets: vessel reference point to transducer
                self.heave = self.dg_MWC['txInfo']['heave_m']
                # Kongsberg: "Number of beams in this datagram (Nrx)."
                self.numBeams = self.dg_MWC['rxInfo']['numBeams']

                # Extract values needed to calculate ranges to each sampleAmplitude05dB value:
                # Kongsberg: "The sample rate is normally decimated to be approximately
                # the same as the bandwidth of the transmitted pulse. Unit hertz."
                self.sampleFreq = self.dg_MWC['rxInfo']['sampleFreq_Hz']
                # Kongsberg: "Sound speed at transducer, unit m/s."
                self.soundVelocity = self.dg_MWC['rxInfo']['soundVelocity_mPerSec']

                swath_dimensions = self.swath_dimensions_calculation()

                # TODO: Figure out how to break this out into a function:
                # If choosing curtain width based on percentage of swath width, we first need to determine swath width (above).
                # Now, place samples in appropriate bins:
                if self.percent_swath_width is not None:
                    # Calculate curtain width based on percent of swath width
                    curtain_width = (self.percent_swath_width / 100) * swath_dimensions['swath_width']
                    # Calculate bin height based on percent of maximum depth
                    vertical_bin_height = (self.bin_height_percent_depth / 100) * swath_dimensions['max_depth']

                    binned_sample_array = [[] for value in range(int(100 / self.bin_height_percent_depth))]
                    binned_count_array = [0 for value in range(int(100 / self.bin_height_percent_depth))]
                    beam_number_array_single_ping = []

                    # For every beam in ping:
                    for beam in range(self.numBeams):  # 0 to self.numBeams - 1
                        # Across-track angle:
                        beam_point_angle_re_vertical = self.dg_MWC['beamData']['beamPointAngReVertical_deg'][beam]
                        # Along-track angle:
                        tilt_angle_re_tx_deg = self.tilt_angle_re_tx_deg[self.dg_MWC['beamData']['beamTxSectorNum'][beam]]
                        tilt_angle_re_vertical_deg = tilt_angle_re_tx_deg + interpolated_pitch

                        print("tilt_angle_re_tx_deg: ", tilt_angle_re_tx_deg)
                        # print("interpolated_pitch: ", interpolated_pitch)
                        # print("tilt_angle_re_vertical_deg: ", tilt_angle_re_vertical_deg)





                        detected_range = self.dg_MWC['beamData']['detectedRangeInSamples'][beam] # Index in sampleAmplitude05dB array where bottom detected

                        # For each watercolumn data point in a single beam:
                        for i in range(detected_range + 1):  # 0 to detected_range
                            range_to_wc_data_point = (self.soundVelocity * i) / (self.sampleFreq * 2)

                            test_x = range_to_wc_data_point * math.sin(math.radians(beam_point_angle_re_vertical))  # Across-track distance

                            # Across-track distance
                            x = (range_to_wc_data_point * math.sin(math.radians(beam_point_angle_re_vertical))
                                * math.cos((math.radians(tilt_angle_re_vertical_deg))))

                            test_y = (range_to_wc_data_point * math.cos(math.radians(beam_point_angle_re_vertical))) - self.heave

                            y = (range_to_wc_data_point * math.cos(math.radians(beam_point_angle_re_vertical))
                                 * math.cos(math.radians(tilt_angle_re_vertical_deg))) - self.heave

                            # if i == detected_range:
                            #     print('test_x: ', test_x)
                            #     print('test_y: ', test_y)
                            #     print('x: ', x)
                            #     print('y: ', y)
                            #     sys.exit()

                            if (-curtain_width / 2) < x < (curtain_width / 2):  # If across-track distance falls within specified curtain:
                                # Calculate depth of watercolumn data point:
                                # TODO: Something is effed up with heave?
                                y = (range_to_wc_data_point * math.cos(math.radians(beam_point_angle_re_vertical))
                                    * math.cos(math.radians(tilt_angle_re_tx_deg))) + self.heave  # Depth
                                # Determine corresponding bin based on depth:
                                # This should create bins 0 to vertical_bin_height, vertical_bin_height+ to 2*vertical_bin_height.
                                # This means we don't need a 101st bin for max_depth; max_depth will fit in 100th bin
                                bin_index = max(0, (math.floor(y / vertical_bin_height) - 1))
                                # Place watercolumn data point in appropriate bin:
                                binned_sample_array[bin_index].append(self.dg_MWC['beamData']['sampleAmplitude05dB_p'][beam][i])
                                binned_count_array[bin_index] += 1

                                if i > 50:
                                    beam_number_array_single_ping.append(beam)

                    if count is 0:
                        print(len(binned_sample_array))

                    # After going through every beam in a ping:
                    for i in range(len(binned_sample_array)):
                        # Average bins:
                        if len(binned_sample_array[i]) > 0:
                            binned_sample_array[i] = sum(binned_sample_array[i]) / len(binned_sample_array[i])
                        else:
                            binned_sample_array[i] = 0  # If there are no values in a bin, enter zero.

                    fullFile_avergedWaterColumnArray.append(binned_sample_array)
                    full_file_samples_per_bin_array.append(binned_count_array)

                    # TODO: TESTING
                    full_file_beam_number_array.append(beam_number_array_single_ping)



                    dummy_array = copy.deepcopy(fullFile_avergedWaterColumnArray)
                    for i in range(len(dummy_array)):
                        for j in range(len(dummy_array[i])):
                            if dummy_array[i][j] is 0:
                                pass
                            else:
                                dummy_array[i][j] = 100


                    # TODO: TEST
                    #print("max count per bin: ", max(binned_count_array))  # 1400-1500ish
                    #print("min count per bin: ", min(binned_count_array))  # 0

                count = count + 1

            # After all MWC offsets in file:

            # TODO: TESTING
            # For test file 0019, this is always 68 and 187--I'm a little surprised this never changes throughout the whole file?
            # (Test file 0019 has some pretty serious roll... Min -1.0 to max +9.0--not sure why so one-sided? Beginning file -1.0 to +6.0)
            # Are the 256 (of 400) beams used for watercolumn determined based on attitude--so it will always be the center 256 beams?
            min_beam = self.numBeams
            max_beam = 0
            for i in range(len(full_file_beam_number_array)):
                if max(full_file_beam_number_array[i]) > max_beam:
                    max_beam = max(full_file_beam_number_array[i])
                if min(full_file_beam_number_array[i]) < min_beam:
                    min_beam = min(full_file_beam_number_array[i])
            print("numBeams: ", self.numBeams)
            print("min_beam: ", min_beam)
            print("max_beam: ", max_beam)

            for i in range(2):
                print("Averaged WC values, index: ", i, ": ", fullFile_avergedWaterColumnArray[i])

            # Convert python list (2D matrix) containing watercolumn curtain data to numpy matrix and transpose for correct orientation:
            np_fullFile_avergedWaterColumnArray = np.transpose(np.array(fullFile_avergedWaterColumnArray))
            np_full_file_samples_per_bin_array = np.transpose(np.array(full_file_samples_per_bin_array))

            np_dummy_array = np.transpose(np.array(dummy_array))

            # Plot!
            # Plot watercolumn:
            fig1 = plt.figure()
            #plt.pcolor: can give values for x and y values, shading=False (or something like that) to remove black boxes if necessary
            #pcolor may flip image based on what it thinks should be 0,0 in top left corner
            #pcolor: give 1D array of x values and 1D array for y values
            plt.imshow(np_fullFile_avergedWaterColumnArray, interpolation='nearest')
            fig1.suptitle('Watercolumn', fontsize=20)
            #plt.yticks(np.arange(0, swath_dimensions['max_depth'], (20 * vertical_bin_height)))
            plt.colorbar()
            plt.show()

            # Plot count per bin:
            fig2 = plt.figure()
            plt.imshow(np_dummy_array, interpolation='nearest')
            fig2.suptitle('Bins With Zero Values', fontsize=20)
            plt.colorbar()
            plt.show()

    def swath_dimensions_calculation(self):
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

    def percent_swath_width_calculation(self, swath_width, min_depth, max_depth, avg_depth):
        # Calculate curtain width based on percent of swath width
        curtain_width = (self.percent_swath_width / 100) * swath_width
        # Calculate bin height based on percent of average depth
        # TODO: This has to be max_depth or will exceed size of binned array when depth is greater than average...
        vertical_bin_height = (self.bin_height_percent_depth / 100) * max_depth

        binned_sample_array = [[] for value in range(int(100 / self.bin_height_percent_depth) + 1)]

        for beam in range(self.numBeams):  # 0 to self.numBeams - 1
            # Across-track angle:
            beam_point_angle_re_vertical = self.dg_MWC['beamData']['beamPointAngReVertical_deg'][beam]
            # Along-track angle:
            tilt_angle_re_tx_deg = self.tilt_angle_re_tx_deg[self.dg_MWC['beamData']['beamTxSectorNum']]

            # Index in sampleAmplitude05dB array where bottom detected:
            detected_range = self.dg_MWC['beamData']['detectedRangeInSamples'][beam]

            # For each watercolumn data point in a single beam:
            for i in range(detected_range + 1):  # 0 to detected range
                range_to_wc_data_point = (self.soundVelocity * i) / (self.sampleFreq * 2)
                x = range_to_wc_data_point * math.sin(math.radians(beam_point_angle_re_vertical))  # Across-track distance

                # If across-track distance falls within specified curtain:
                if (-curtain_width / 2) < x < (curtain_width / 2):
                    # Calculate depth of watercolumn data point:
                    y = range_to_wc_data_point * math.cos(math.radians(beam_point_angle_re_vertical))  # Depth
                    # Determine corresponding bin based on depth:
                    bin_index = math.floor(y / vertical_bin_height)
                    # Place watercolumn data point in appropriate bin:
                    binned_sample_array[bin_index].append(self.dg_MWC['beamData']['sampleAmplitude05dB_p'][beam][i])

        # After going through every beam in a ping:
        for i in range(len(binned_sample_array)):
            # Average bins:
            if len(binned_sample_array[i]) > 0:
                # TODO: TESTING:
                for j in binned_sample_array[i]:
                    if j is 0:
                        print("ZERO")
                        print(binned_sample_array[i])
                binned_sample_array[i] = sum(binned_sample_array[i]) / len(binned_sample_array[i])
            else:
                binned_sample_array[i] = 0  # If there are no values in a bin, enter zero.

        fullFile_avergedWaterColumnArray.append(binned_sample_array)
        pass

    def run(self):
        print("Running!")
        self.extract_wc_from_file()


if __name__ == '__main__':

    # Defaults:
    input_file = None
    percent_swath_width = 10
    num_nadir_beams = None
    across_track_width = None
    bin_height_percent_depth = 1

    # Read command line args for files/directory, num_nadir_beams:
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:p:n:w:b:", ["file=", "percent=", "num_beams=", "width=", "bin_height="])
    except getopt.GetoptError:
        print("watercolumnplot.py")
        # TODO: Not sure about .kmwcd file extensions?
        # File or directory containing .kmall or .kmwcd file extensions
        print("-f <file_or_directory>")
        # Three options for selecting across-track width to include in watercolumn curtain:
        # 1. A percent of overall swath width, centered at nadir beam.
        # (Example: 10% of overall swath width.)
        print("-p <percent_swath_width>")
        # 2. A specific number of beams, centered at nadir beam.
        # (Example: 20 beams.)
        print("-n <number_nadir_beams>")
        # 3. A specific number of meters, centered at nadir beam.
        # (Example: 10 meters gives 5 meters left and 5 meters right of nadir beam.)
        print("-w <across_track_width>")
        # Vertical bin height as a percent of maximum depth. 1% by default.
        # (Example: 100 one-meter bins in 100 meters of water.)
        print("-b <bin_height_percent_depth>")
        sys.exit()

    for opt, arg in opts:
        if opt == "-h":
            print("watercolumnplot.py")
            print("-f <file_or_directory>")
            print("-p <percent_swath_width>")
            print("-n <number_nadir_beams>")
            print("-w <across_track_width>")
            print("-b <bin_height_percent_depth>")
            sys.exit()
        elif opt in ("-f", "--file"):
            input_file = arg
        elif opt in ("-p", "--percent"):
            percent_swath_width = int(arg)
        elif opt in ("-n", "--num_beams"):
            num_nadir_beams = int(arg)
            percent_swath_width = None  # Undo default value
        elif opt in ("-w, --width"):
            across_track_width = int(arg)
            percent_swath_width = None  # Undo default value
        elif opt in ("-b, --bin_height"):
            bin_height_percent_depth = int(arg)

    if input_file is None:
        print("Must enter file or directory: watercolumnplot.py -f <file_or_directory>")
        sys.exit()

    wc_plotter = WaterColumnPlot(input_file, percent_swath_width, num_nadir_beams,
                                 across_track_width, bin_height_percent_depth)
    wc_plotter.run()