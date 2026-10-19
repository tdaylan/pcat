import pickle
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from scipy.integrate import trapezoid

from pcat import ProposalContext
from pcat.main import _retr_logproposal_penalty, _retr_proposal_context
from pcat.neural_proposals import NeuralElementProposal


def _state_condition_for_test(gdat, context):
    return np.array([gdat.neural_proposal_features[0], context.element_counts[0]])


def _fit_proposal():
    pytest.importorskip("torch")
    random = np.random.default_rng(12)
    conditions = random.normal(size=(256, 2))
    component = random.integers(0, 2, size=conditions.shape[0])
    means = np.where(component[:, None] == 0, [0.25, 0.75], [0.75, 0.25])
    elements = np.clip(means + random.normal(scale=0.04, size=means.shape), 0., 1.)
    return NeuralElementProposal.fit(
        conditions,
        elements,
        number_components=4,
        hidden_sizes=(16,),
        epochs=20,
        batch_size=64,
        uniform_fraction=0.05,
        seed=3,
    )


def test_neural_proposal_import_does_not_require_torch():
    assert NeuralElementProposal.__name__ == "NeuralElementProposal"


def test_neural_proposal_callbacks_are_worker_serializable():
    cloudpickle = pytest.importorskip("cloudpickle")
    proposal = _fit_proposal()
    callbacks = cloudpickle.loads(cloudpickle.dumps(proposal.callbacks()))
    state = SimpleNamespace(neural_proposal_features=np.array([0.2, -0.1]))
    unit = callbacks["retr_drawpropelem"](state)

    assert np.isfinite(callbacks["retr_lpdfpropelem"](state, unit))


def test_proposal_context_copies_and_freezes_catalog_coordinates():
    source = np.array([[0.2, 0.4]])
    context = ProposalContext(0, (source, np.empty((0, 3))))
    source[0, 0] = 0.9

    np.testing.assert_array_equal(context.unit_elements, [[0.2, 0.4]])
    assert context.element_counts == (1, 0)
    with pytest.raises(ValueError):
        context.unit_elements[0, 0] = 0.5


def test_death_reverse_context_excludes_removed_element():
    model = SimpleNamespace(
        numbpopl=1,
        numbparagenrbase=0,
        numbparagenrelemcuml=[0],
        numbparagenrelemsing=[2],
        indxparagenrelemsing=[np.array([0, 1])],
    )
    unit_parameters = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

    context = _retr_proposal_context(
        model, unit_parameters, [[0, 1, 2]], 0, removed_element=1
    )

    np.testing.assert_array_equal(context.unit_elements, [[0.1, 0.2], [0.5, 0.6]])


@pytest.mark.parametrize(
    ("log_density_forward", "log_density_reverse", "expected"),
    [(-1.7, 0., -1.7), (0., -1.7, 1.7), (-2., -3., 1.)],
)
def test_proposal_penalty_is_forward_minus_reverse(
    log_density_forward, log_density_reverse, expected
):
    assert _retr_logproposal_penalty(
        log_density_forward, log_density_reverse
    ) == pytest.approx(expected)


def test_neural_state_callbacks_use_proposal_context():
    proposal = _fit_proposal()
    callbacks = proposal.state_callbacks(_state_condition_for_test)
    state = SimpleNamespace(neural_proposal_features=np.array([0.2]))
    context = ProposalContext(0, (np.array([[0.1, 0.3]]),))
    unit = callbacks["retr_drawpropelem_state"](state, context)
    callbacks = pickle.loads(pickle.dumps(callbacks))

    assert unit.shape == (2,)
    assert np.all((unit >= 0.) & (unit <= 1.))
    assert np.isfinite(callbacks["retr_lpdfpropelem_state"](state, unit, context))


def test_neural_density_is_normalized_and_callback_density_matches_draw(tmp_path):
    proposal = _fit_proposal()
    coordinate = np.linspace(0., 1., 241)
    grid_first, grid_second = np.meshgrid(coordinate, coordinate, indexing="ij")
    grid = np.column_stack((grid_first.ravel(), grid_second.ravel()))
    log_density = proposal.log_prob_batch(np.array([0.1, -0.2]), grid)
    density = np.exp(log_density).reshape(grid_first.shape)
    integral = trapezoid(trapezoid(density, coordinate, axis=1), coordinate)
    assert integral == pytest.approx(1., abs=0.025)

    state = SimpleNamespace(neural_proposal_features=np.array([0.1, -0.2]))
    callbacks = proposal.callbacks()
    unit = callbacks["retr_drawpropelem"](state)
    assert unit.shape == (2,)
    assert np.all((unit >= 0.) & (unit <= 1.))
    assert callbacks["retr_lpdfpropelem"](state, unit) == pytest.approx(
        proposal.log_prob(state.neural_proposal_features, unit)
    )
    assert proposal.log_prob(state.neural_proposal_features, [1.1, 0.5]) == -np.inf

    path = tmp_path / "proposal.pt"
    proposal.save(path)
    restored = NeuralElementProposal.load(path)
    assert restored.log_prob(state.neural_proposal_features, unit) == pytest.approx(
        proposal.log_prob(state.neural_proposal_features, unit)
    )


def test_neural_proposal_callbacks_drive_pcat_catalog_moves(tmp_path):
    pytest.importorskip("torch")
    from astropy.io import fits
    from pcat import sampling

    random = np.random.default_rng(6)
    proposal = NeuralElementProposal.fit(
        random.normal(size=(64, 1)),
        random.uniform(size=(64, 4)),
        number_components=2,
        hidden_sizes=(8,),
        epochs=4,
        batch_size=32,
        seed=4,
    )
    edges = np.linspace(1.8e5, 1.82e5, 41)  # [m^-1]
    width = np.diff(edges)  # [m^-1]
    centers = 0.5 * (edges[1:] + edges[:-1])  # [m^-1]
    counts = np.rint(1000. + 800. * np.exp(-0.5 * ((centers - 1.81e5) / 60.)**2))  # [count]
    path_input = tmp_path / "neural_catalog" / "data" / "inpt"
    path_input.mkdir(parents=True)
    fits.writeto(path_input / "sbrt.fits", counts[:, None, None, None])
    fits.writeto(path_input / "expo.fits", (1. / width)[:, None, None])

    result = sampling.sample(
        typeexpr="fire",
        typedata="inpt",
        strgexprsbrt="sbrt.fits",
        typeexpo="file",
        strgexpo="expo.fits",
        binsenerfull=edges,
        spectype=["voig"],
        spatdisttype=["line"],
        typeelem=["lghtlinevoig"],
        dictfitt={
            "typeelem": ["lghtlinevoig"],
            "spectype": ["voig"],
            "sbrtbacknorm": [np.full((edges.size - 1, 1, 1), 1000.)],
            "listnamediff": ["back0000"],
        },
        limtparaelem={"flux": (1e3, 1e6), "sigm": (20., 200.), "gamm": (1., 50.)},
        maxmgangdata=1e-4,
        anlytype="spec",
        fittminmnumbelempop0=0,
        fittmaxmnumbelempop0=3,
        inittype="rand",
        typeseed=2,
        numbproc=1,
        probtran=1.,
        probspmr=0.,
        probjump=0.,
        neural_proposal_features=np.array([0.25]),
        numbswep=400,
        numbburn=100,
        numbsamp=300,
        boolcheckconv=False,
        boolmakeplot=False,
        boolmakeplotinit=False,
        booldiag=False,
        typeverb=0,
        pathbase=str(tmp_path / "neural_catalog"),
        strgcnfg="neural_catalog",
        **proposal.callbacks(),
    )

    path_output = tmp_path / "neural_catalog" / "data" / "outp" / "neural_catalog"
    with h5py.File(path_output / "gdatmodi0000post.h5", "r") as output:
        proposal_types = np.asarray(output["listpostindxproptype"][()]).reshape(-1)
        assert np.any(proposal_types == 1)
        assert np.any(proposal_types == 2)
        assert np.isfinite(output["listpostaccplprb"][()]).any()

def test_state_conditioned_neural_callbacks_drive_pcat_birth_death(tmp_path):
    pytest.importorskip("torch")
    from astropy.io import fits
    from pcat import sampling

    random = np.random.default_rng(16)
    proposal = NeuralElementProposal.fit(
        random.normal(size=(64, 2)),
        random.uniform(size=(64, 4)),
        number_components=2,
        hidden_sizes=(8,),
        epochs=4,
        batch_size=32,
        seed=5,
    )
    state_callbacks = proposal.state_callbacks(_state_condition_for_test)

    edges = np.linspace(1.8e5, 1.82e5, 41)
    width = np.diff(edges)
    centers = 0.5 * (edges[1:] + edges[:-1])
    counts = np.rint(1000. + 800. * np.exp(-0.5 * ((centers - 1.81e5) / 60.)**2))
    path_input = tmp_path / "state_neural_catalog" / "data" / "inpt"
    path_input.mkdir(parents=True)
    fits.writeto(path_input / "sbrt.fits", counts[:, None, None, None])
    fits.writeto(path_input / "expo.fits", (1. / width)[:, None, None])

    sampling.sample(
        typeexpr="fire",
        typedata="inpt",
        strgexprsbrt="sbrt.fits",
        typeexpo="file",
        strgexpo="expo.fits",
        binsenerfull=edges,
        spectype=["voig"],
        spatdisttype=["line"],
        typeelem=["lghtlinevoig"],
        dictfitt={
            "typeelem": ["lghtlinevoig"],
            "spectype": ["voig"],
            "sbrtbacknorm": [np.full((edges.size - 1, 1, 1), 1000.)],
            "listnamediff": ["back0000"],
        },
        limtparaelem={"flux": (1e3, 1e6), "sigm": (20., 200.), "gamm": (1., 50.)},
        maxmgangdata=1e-4,
        anlytype="spec",
        fittminmnumbelempop0=0,
        fittmaxmnumbelempop0=3,
        inittype="rand",
        typeseed=8,
        numbproc=1,
        probtran=1.,
        probspmr=0.,
        probjump=0.,
        neural_proposal_features=np.array([0.25]),
        numbswep=400,
        numbburn=100,
        numbsamp=300,
        boolcheckconv=False,
        boolmakeplot=False,
        boolmakeplotinit=False,
        booldiag=False,
        typeverb=0,
        pathbase=str(tmp_path / "state_neural_catalog"),
        strgcnfg="state_neural_catalog",
        **state_callbacks,
    )

    path_output = tmp_path / "state_neural_catalog" / "data" / "outp" / "state_neural_catalog"
    with h5py.File(path_output / "gdatmodi0000post.h5", "r") as output:
        proposal_types = np.asarray(output["listpostindxproptype"][()]).reshape(-1)
        assert np.any(proposal_types == 1)
        assert np.any(proposal_types == 2)
        assert np.isfinite(output["listpostaccplprb"][()]).any()
