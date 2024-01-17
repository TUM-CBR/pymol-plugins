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
        product_conc = (value / (molar_extinction * distance)) * conc_units

        current_vel = abs(prev_product_conc - product_conc) / intervals
        prev_product_conc = product_conc

        # We assume that for every molecule of the product
        # we consume a molecule of the substrate
        yield Point2d(
            x = conc - product_conc,
            y = current_vel# / conc
        )

def run_to_conc_vs_time(
        global_attributes : GlobalAttributes,
        run: KineticsRun,
        baseline : float = 0
    ) -> Iterable[Point2d]:
    conc_units = global_attributes.concentration_units
    intervals = global_attributes.measurement_interval
    molar_extinction = global_attributes.molar_extinction
    distance = global_attributes.distance
    simulation_spec = SimulationSpec(
        periods = len(run.data),
        interval = intervals
    )

    for time, value in zip(simulation_spec.times(), run.data):

        # Remove the baseline noise
        value -= baseline

        # abs = conc * molar_extinction * distance
        # => conc = abs / (molar_extinction * distance)
        product_conc = (value / (molar_extinction * distance)) * conc_units

        # We assume that for every molecule of the product
        # we consume a molecule of the substrate
        yield Point2d(
            x = time,
            y = product_conc
        )


def as_conc_vs_time_series(runs: KineticsRuns) -> List['Series[RunMetadata]']:

    baseline_ix, baseline_run = next(
        ((ix, run.data) for ix,run in enumerate(runs.runs) if is_baseline_run(run)),
        (-1, [0.0])
    )
    baseline = avg(baseline_run)

    return [
        Series(
            metadata = run.run_metadata,
            values = list(run_to_conc_vs_time(runs.global_attributes, run, baseline))
        )
        for ix, run in enumerate(runs.runs) if baseline_ix != ix
    ]

def as_vel_vs_conc_series(
        runs: KineticsRuns,
        samples_to_consider: int
    ) -> 'Series[RunVelocityMetadata]':

    conc_vs_time = as_conc_vs_time_series(runs)

    def slope(input: 'Series[RunMetadata]'):
        points = input.values[0:samples_to_consider]
        m, _ = lse(
            [point.x for point in points],
            [point.y for point in points]
        )
        return m

    return Series(
        metadata = RunVelocityMetadata(),
        values =  [
            Point2d(
                series.metadata.concentration * runs.global_attributes.concentration_units,
                slope(series)
            )
            for series in conc_vs_time
        ]
    )