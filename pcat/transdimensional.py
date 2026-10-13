"""Birth/death posterior sampling for variable catalogs with custom likelihoods."""

import numpy as np
from scipy.special import gammaln, logsumexp


def sample_catalog(
    log_likelihood, base_bounds, component_bounds,
    *, maximum_components=3, count_prior_mean=1.0, number_steps=12000,
    burn_steps=4000, thinning=10, proposal_fraction=0.05, seed=None,
    initial_base=None,
):
    """Sample a custom-likelihood catalog with finite uniform parameter priors.

    The count prior is Poisson, normalized over 0..maximum_components. Catalog
    rows are exchangeable. Births draw a row from its normalized prior and insert
    it at a random position; deaths choose a row uniformly. Hastings factors
    include state-dependent move probabilities, row-selection probabilities,
    and the birth density. The augmentation Jacobian is one. Within-dimension
    proposals reflect at prior bounds. This routine is single-chain and does
    not perform evidence estimation or establish convergence.
    """

    base_bounds = np.asarray(base_bounds, dtype=float)
    component_bounds = np.asarray(component_bounds, dtype=float)
    for bounds in (base_bounds, component_bounds):
        if (bounds.ndim != 2 or bounds.shape[1] != 2 or bounds.shape[0] == 0
                or not np.isfinite(bounds).all() or np.any(bounds[:, 1] <= bounds[:, 0])):
            raise ValueError('Parameter bounds must be finite increasing [minimum, maximum] pairs')
    for value in (maximum_components, number_steps, burn_steps, thinning):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError('Counts and iteration settings must be integers')
    if maximum_components < 0 or not 0 <= burn_steps < number_steps or thinning < 1:
        raise ValueError('Provide nonnegative count bounds, burn < steps, and positive thinning')
    if not np.isfinite(count_prior_mean) or count_prior_mean <= 0:
        raise ValueError('count_prior_mean must be finite and positive')
    if not np.isfinite(proposal_fraction) or proposal_fraction <= 0:
        raise ValueError('proposal_fraction must be finite and positive')
    generator = np.random.default_rng(seed)
    base = generator.uniform(base_bounds[:, 0], base_bounds[:, 1]) if initial_base is None else np.asarray(initial_base, dtype=float).copy()
    if (base.shape != (base_bounds.shape[0],) or not np.isfinite(base).all()
            or np.any(base < base_bounds[:, 0]) or np.any(base > base_bounds[:, 1])):
        raise ValueError('Initial base parameters must lie inside their bounds')
    catalog = np.empty((0, component_bounds.shape[0]))
    counts = np.arange(maximum_components + 1)
    log_count_prior = counts * np.log(count_prior_mean) - gammaln(counts + 1)
    log_count_prior -= logsumexp(log_count_prior)
    log_component_volume = np.log(np.diff(component_bounds, axis=1)[:, 0]).sum()
    log_birth_density = -log_component_volume

    def posterior(current_base, current_catalog):
        likelihood = float(log_likelihood(current_base, current_catalog))
        if np.isnan(likelihood) or likelihood == np.inf:
            raise ValueError('Likelihood must return a finite value or negative infinity')
        return likelihood + log_count_prior[len(current_catalog)] - len(current_catalog) * log_component_volume

    def move_probabilities(count):
        weights = np.array([0.5, 0.25 if count < maximum_components else 0.0,
                            0.25 if count > 0 else 0.0])
        return weights / weights.sum()

    def reflected(value, lower, upper):
        width = upper - lower
        phase = (value - lower) % (2.0 * width)
        return lower + min(phase, 2.0 * width - phase)

    current_posterior = posterior(base, catalog)
    if not np.isfinite(current_posterior):
        raise ValueError('Initial empty-catalog state must have a finite likelihood')
    base_samples, catalogs, log_samples, count_trace = [], [], [], []
    proposals = np.zeros(3, dtype=int)
    acceptances = np.zeros(3, dtype=int)
    for step in range(number_steps):
        count = len(catalog)
        probabilities = move_probabilities(count)
        move = generator.choice(3, p=probabilities)
        proposals[move] += 1
        proposed_base, proposed_catalog = base.copy(), catalog.copy()
        log_hastings = 0.0
        if move == 1:
            row = generator.uniform(component_bounds[:, 0], component_bounds[:, 1])
            position = generator.integers(count + 1)
            proposed_catalog = np.insert(catalog, position, row, axis=0)
            log_hastings = np.log(move_probabilities(count + 1)[2]) - np.log(probabilities[1]) - log_birth_density
        elif move == 2:
            position = generator.integers(count)
            proposed_catalog = np.delete(catalog, position, axis=0)
            log_hastings = np.log(move_probabilities(count - 1)[1]) - np.log(probabilities[2]) + log_birth_density
        else:
            coordinate = generator.integers(base.size + catalog.size)
            if coordinate < base.size:
                bounds = base_bounds[coordinate]
                value = base[coordinate] + generator.normal() * proposal_fraction * (bounds[1] - bounds[0])
                proposed_base[coordinate] = reflected(value, *bounds)
            else:
                row, column = divmod(coordinate - base.size, component_bounds.shape[0])
                bounds = component_bounds[column]
                value = catalog[row, column] + generator.normal() * proposal_fraction * (bounds[1] - bounds[0])
                proposed_catalog[row, column] = reflected(value, *bounds)
        proposed_posterior = posterior(proposed_base, proposed_catalog)
        if np.log(generator.random()) < proposed_posterior - current_posterior + log_hastings:
            base, catalog, current_posterior = proposed_base, proposed_catalog, proposed_posterior
            acceptances[move] += 1
        count_trace.append(len(catalog))
        if step >= burn_steps and (step - burn_steps) % thinning == 0:
            base_samples.append(base.copy())
            catalogs.append(catalog.copy())
            log_samples.append(current_posterior)
    return {
        'base_samples': np.asarray(base_samples), 'catalog_samples': catalogs,
        'component_counts': np.asarray([len(sample) for sample in catalogs]),
        'count_trace': np.asarray(count_trace), 'log_posterior': np.asarray(log_samples),
        'proposal_counts': dict(zip(('within', 'birth', 'death'), proposals)),
        'acceptance_counts': dict(zip(('within', 'birth', 'death'), acceptances)),
    }