import numpy as np
from typing import Iterable, List, Tuple

def lse(xs: List[float], ys: List[float]) -> Tuple[float, float]:
    a = np.vstack([xs, np.ones(len(xs))]).T
    m, b = np.linalg.lstsq(a, ys, rcond=None)[0]

    return m, b

def avg(xs: Iterable[float]) -> float:

    count = 0
    total = 0
    for x in xs:
        total += x
        count += 1

    return total/count