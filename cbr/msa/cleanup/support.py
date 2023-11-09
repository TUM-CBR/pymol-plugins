from typing import Iterable, List

def normalized(values : List[float]) -> List[float]:
    min_v = min(values)
    max_v = max(values)
    spread = max_v - min_v

    if abs(spread) < 1:
        spread = 1

    return [(value - min_v)/spread for value in values]

def score_contigous(
    values : Iterable[bool]
    ) -> List[int]:

    def __run__():
        segment_size = 0

        for value in enumerate(values):
            
            if value:
                segment_size += 1
            else:
                segment_size = 0

            yield segment_size


    return list(__run__())