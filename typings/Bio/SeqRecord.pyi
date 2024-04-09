from typing import Iterator


class SeqRecord:

    def __iter__(self) -> Iterator[str]: ...