"""Immutable catalog snapshots supplied to state-conditioned proposals."""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ProposalContext:
    """Active unit-prior element coordinates at a transdimensional proposal.

    ``unit_elements_by_population`` is indexed by PCAT population. Each entry
    has shape ``(number_of_active_elements, parameters_per_element)``. The
    arrays are copied and made read-only so callbacks cannot alter sampler state.
    """

    population_index: int
    unit_elements_by_population: tuple[np.ndarray, ...]

    def __post_init__(self):
        population_index = int(self.population_index)
        populations = tuple(np.array(values, dtype=float, copy=True) for values in self.unit_elements_by_population)
        if population_index < 0 or population_index >= len(populations):
            raise ValueError("population_index must identify an entry in unit_elements_by_population.")
        for values in populations:
            if values.ndim != 2 or not np.isfinite(values).all():
                raise ValueError("population element coordinates must be finite two-dimensional arrays.")
            if np.any((values < 0.) | (values > 1.)):
                raise ValueError("population element coordinates must lie in the unit cube [0, 1].")
            values.setflags(write=False)
        object.__setattr__(self, "population_index", population_index)
        object.__setattr__(self, "unit_elements_by_population", populations)

    @property
    def unit_elements(self):
        """Return active unit coordinates for the population being proposed."""
        return self.unit_elements_by_population[self.population_index]

    @property
    def element_counts(self):
        """Return the active element count in each population."""
        return tuple(values.shape[0] for values in self.unit_elements_by_population)


__all__ = ["ProposalContext"]