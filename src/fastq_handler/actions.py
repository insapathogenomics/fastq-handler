
import os
from abc import ABC, ABCMeta, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Type

import pandas as pd

from fastq_handler.records import Processed
from fastq_handler.utilities import Utils


class ProcessAction(ABCMeta):
    """
    abstract class for processing action
    """

    @staticmethod
    @abstractmethod
    def process(fastq_file: str, sample_id: str, processed: Processed):
        """
        process
        """
        pass

    def estimate_output_size(fastq_file: str, sample_id: str, processed: Processed) -> int:
        """
        output_size
        """
        pass


class ProcessActionMergeWithLast(ProcessAction):
    """get_dir_merged_last
    class to merge with last"""

    @staticmethod
    def process(fastq_file: str, sample_id: str, processed: Processed):
        """
        process
        """
        last_run_file = processed.get_id_merged_last(
            sample_id)

        if last_run_file == "":
            return

        utils = Utils()

        utils.append_file_to_gz(
            last_run_file, fastq_file
        )

    @staticmethod
    def estimate_output_size(fastq_file_size: int, sample_id: str, processed: Processed) -> int:
        """
        output_size
        """
        last_run_file = processed.get_id_merged_last(
            sample_id)

        if last_run_file == "":
            return fastq_file_size

        return fastq_file_size + os.path.getsize(last_run_file)
