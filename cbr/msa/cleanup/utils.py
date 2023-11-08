from typing import List

def normalized(values : List[float]) -> List[float]:
    min_v = min(values)
    max_v = max(values)
    spread = max_v - min_v

    return [(value - min_v)/spread for value in values]