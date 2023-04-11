
import os
from dataclasses import dataclass
from typing import List, Optional, Type

from fastq_handler.actions import ProcessAction


@dataclass
class InputState:
    """
    class to hold input - defaults set for inheritance reasons. check py3.10 for better solution
    """

    fastq_dir: str = ""
    name_tag: str = ""


@dataclass
class RunParams:
    """
    class to hold params with default values
    """
    actions: Optional[List[Type[ProcessAction]]] = None
    keep_name: bool = False
    sleep_time: int = 10
    max_size: int = 100000000


@dataclass
class OutputDirs:
    """
    class to hold output dirs
    """
    output_dir: str = "output"
    logs_dirname: str = "logs"


@dataclass
class RunConfig(InputState, RunParams, OutputDirs):
    """
    class to hold run config"""

    def __post_init__(self):
        """
        post init
        """
        if self.actions is None:
            self.actions = []

        self.logs_dir = os.path.join(self.output_dir, self.logs_dirname)
        self.output_dir = os.path.abspath(self.output_dir)
