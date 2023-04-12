
import os
from abc import ABC, ABCMeta, abstractmethod
from dataclasses import dataclass
from typing import Iterable, List, Optional, Type

import dnaio
import pandas as pd
import xopen
from fastq_filter import fastq_records_to_file, file_to_fastq_records

from fastq_handler.records import Processed
from fastq_handler.utilities import Utils


class ProcessAction(ABC):
    """
    abstract class for processing action
    """

    __name__ = "process_action"

    @abstractmethod
    def __init__(self):
        pass

    @staticmethod
    @abstractmethod
    def process(fastq_file: str, sample_id: str, processed: Processed):
        """
        process
        """
        pass

    @abstractmethod
    def estimate_output_size(fastq_file: str, sample_id: str, processed: Processed) -> int:
        """
        output_size
        """
        pass

    @abstractmethod
    def __str__(self):
        pass


class ProcessActionMergeWithLast(ProcessAction):
    """get_dir_merged_last
    class to merge with last"""

    __name__ = "merge_with_last"

    def __init__(self):
        super().__init__()

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

    def __str__(self):
        return __name__


DEFAULT_COMPRESSION_LEVEL = 2


def fastq_records_to_file_w_max(records: Iterable[dnaio.Sequence], filepath: str,
                                compression_level: int = DEFAULT_COMPRESSION_LEVEL, max_size: int = 1000000000, margin=500000):
    """
    write fastq records to file, stop when max size is reached
    """

    with xopen.xopen(filepath, mode='wb', threads=0,
                     compresslevel=compression_level) as output_h:

        for record in records:

            record = record.fastq_bytes()

            if os.path.getsize(filepath) > (max_size - margin):
                break

            output_h.write(record)


class ProcessActionDownsize(ProcessAction):
    """
    class to downsize fastq files
    """

    __name__ = "downsize"

    def __init__(self, max_size: int):
        self.max_size = max_size

    def process(self, fastq_file: str, sample_id: str, processed: Processed):
        """
        process
        """

        utils = Utils()

        records = []

        if os.path.getsize(fastq_file) < self.max_size:
            return

        records = file_to_fastq_records(fastq_file)

        tmp_file = utils.temp_fastq_file(fastq_file=fastq_file)

        fastq_records_to_file_w_max(
            records, tmp_file, max_size=self.max_size)

        utils.move_file(tmp_file, fastq_file)

    def estimate_output_size(self, fastq_file_size: int, sample_id: str, processed: Processed) -> int:
        """
        output_size
        """
        if fastq_file_size < self.max_size:
            return fastq_file_size

        return self.max_size

    def __str__(self):
        return __name__
