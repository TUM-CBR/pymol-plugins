from .data import *
from .math import *

def is_baseline_run(run : KineticsRun) -> bool:
    return int(run.run_metadata.concentration * 1000) == 0

def run_to_vel_vs_conc(
    global_attributes : GlobalAttributes,
    run : KineticsRun,
    baseline : float = 0
) -> Iterable[Point2d]:
    conc_units = global_attributes.concentration_units
    intervals = global_attributes.measurement_interval
    conc = run.run_metadata.concentration
    molar_extinction = global_attributes.molar_extinction
    distance = global_attributes.distance
    prev_product_conc = 0

    for value in run.data:

        # Remove the baseline noise
        value -= baseline

        # abs = conc * molar_extinction * distance
        # => conc = abs / (molar_extinction * distance)
        product_conc = (value / (molar_extinction * distance)) / conc_units

        current_vel = abs(prev_product_conc - product_conc) / intervals
        prev_product_conc = product_conc

        # We assume that for every molecule of the product
        # we consume a molecule of the substrate
        yield Point2d(
            x = conc - product_conc,
            y = current_vel# / conc
        )

def as_conc_vs_time_series(runs: KineticsRuns) -> List['Series[RunMetadata]']:
    raise ValueError("Not implemented")

def as_vel_vs_conc_series(runs: KineticsRuns) -> List['Series[RunMetadata]']:

    baseline_run = next(
        (run.data for run in runs.runs if is_baseline_run(run)),
        [0.0]
    )
    baseline = avg(baseline_run)

    return [
        Series(
            metadata = run.run_metadata,
            values = list(run_to_vel_vs_conc(runs.global_attributes, run, baseline))
        )
        for run in runs.runs if not is_baseline_run(run)
    ]

def combined_velocity_vs_conc(velocity_vs_conc: List['Series[RunMetadata]']) -> List[Point2d]:
    result = [
        point
        for series in velocity_vs_conc
        for point in series.values
    ]
    result.sort(key=lambda point: point.x, reverse=True)
    return result