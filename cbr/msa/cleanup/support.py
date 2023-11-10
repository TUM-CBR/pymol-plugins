from typing import Iterable, List

def normalized(values : List[float]) -> List[float]:
    min_v = min(values)
    max_v = max(values)
    spread = max_v - min_v

    if abs(spread) < 1:
        spread = 1

    return [(value - min_v)/spread for value in values]

def score_contigous(
    values : Iterable[bool],
    min_segment_size = 1
    ) -> List[int]:
    
    segments = set()
    segment_start = -1
    values_len = 0

    for ix,value in enumerate(values):
        
        segment_size = ix - segment_start - 1

        if not value and segment_size >= min_segment_size:

            # segment start is set when a 'false' value is
            # observed. Therefore, it's index should not be
            # included in the segment
            segments.add((segment_start + 1, ix))
            segment_start = ix
        elif not value:
            segment_start = ix

        values_len += 1

    result = list(int(0) for _ in range(0, values_len))

    for (start, end) in segments:
        for ix in range(start, end):
            result[ix] = end - start

    return result