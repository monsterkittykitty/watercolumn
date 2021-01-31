#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import sys

# __name__ is module's name
logger = logging.getLogger(__name__)

# TODO: I'm not sure this type of class/subclass inheritance is necessary here?

# General FileHandler class should:
# Check valid paths (create folder if doesn't already exist?)
# Check valid file extensions
class FileHandler:

    def __init__(self, input_file_or_directory_or_directory=None):
        self.input_file_or_directory_or_directory = input_file_or_directory_or_directory

        # # 3 ways to handle files:
        # # 1. Treat input_file_or_directory as a snapshot, using contents only at time of object creation. Defaults True.
        # self.snapshot_directory = snapshot_directory
        # # 2. Continually monitor input_file_or_directory for new content. Defaults False.
        # self.monitor_directory = monitor_directory
        # # 3. Treat input_file_or_directory as a directory to write captured UDP datagrams to. Defaults False.
        # self.write_directory = write_directory

    def action(self, input_file_or_directory):
        pass


# SnapshotDirectory should:
# Create and return a list of files
class SnapshotDirectory(FileHandler):
    def __init__(self, input_file_or_directory):
        super().__init__(input_file_or_directory)

    def action(self, input_file_or_directory):
        """
        Creates and returns a list of file paths to .kmall and/or .kmwcd file(s). NOTE: This does not continue to
        monitor a directory for new files, it is only a snapshot in time.
        :param input_file_or_directory: Path to file or folder containing .kmall and/or .kmwcd files from which to
        extract watercolumn data.
        :return: Returns a list containing file path(s) to .kmall or .kmwcd file(s) at the specified location given by
        argument (input_file_or_directory). If argument is a path to an appropriate file, expect the list to be of
        length 1; if argument is a path to a directory containing multiple appropriate files, expect list to be of
        length >1.
        """
        temp_list = []

        if os.path.isfile(input_file_or_directory):
            if input_file_or_directory.lower().endswith(('.kmall', '.kmwcd')):
                temp_list.append(input_file_or_directory)
        elif os.path.isdir(input_file_or_directory):
            for root, dirs, files in os.walk(input_file_or_directory):
                for filename in files:
                    if filename.lower().endswith(('.kmall', '.kmwcd')):
                        temp_list.append(os.path.join(root, filename))

        if len(temp_list) > 0:
            return temp_list
        else:
            logger.warning("Invalid file path: %s. No .kmall or .kmwcd found." % input_file_or_directory)
            sys.exit(1)


# MonitorDirectory should:
# Create and return a list of files
# Continue to monitor directory and append (or create new?) list as needed
class MonitorDirectory(FileHandler):
    def __init__(self, input_file_or_directory):
        super().__init__(input_file_or_directory)

    # TODO:
    def action(self, input_file_or_directory):
        pass


# WriteDirectory should:
# Write file to directory when requested.
class WriteDirectory(FileHandler):
    def __init__(self, input_file_or_directory):
        super().__init__(input_file_or_directory)

    # TODO:
    def action(self, input_file_or_directory):
        pass