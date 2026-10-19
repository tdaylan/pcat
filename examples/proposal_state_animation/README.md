# Split/merge Jacobian example

This short PCAT run makes the split/merge Jacobian visible. It uses a seeded simulated spectrum with two Voigt lines and a one-to-three-line metamodel, then exercises within-model, birth, death, split, and merge proposals over 120 sweeps. The run illustrates the calculation. It is too short to establish posterior convergence.

A split maps a parent line's amplitude-coordinate value $F$ to $rF$ and $(1-r)F$, while the proposal auxiliary variables describe the displacement and split fraction. For this amplitude-fraction transform, PCAT adds $\log|J|=\log F$ to the split log acceptance ratio. A reverse merge uses the negative term. More generally,

$$
\log\alpha = \Delta\log\pi + \log(q_{\rm reverse}/q_{\rm forward}) + \log|J|.
$$

The example records the unclipped log acceptance ratio as well as the clipped probability. Its no-Jacobian values subtract only $\log|J|$ from each saved log ratio, holding the proposal, target-density change, and Hastings term fixed. The comparison uses the same proposed states. The figure also plots the acceptance curves using the median split and merge Jacobians from the run.

The animated comparison below isolates a change-of-variables Jacobian without changing model dimension. Both one-parameter models describe the same positive amplitude: one uses $A$ with a uniform density on $[1,e]$, and one uses $\eta=\log A$ with the transformed density proportional to $e^\eta$. For a linear-to-log jump, the target-density ratio is $A$ and the Jacobian is $1/A$, giving unit acceptance. Omitting the Jacobian makes the reverse acceptance probability $1/A$ and biases the chain toward the logarithmic parameterization. With the default seed, all 40 Jacobian-corrected jumps are accepted versus 31 of 40 without the correction. Other seeds can produce different rates.

Run the example and notebook from the PCAT repository root:

```bash
python examples/proposal_state_animation/proposal_state_animation.py
```

![Retained model, every proposed candidate, and Jacobian acceptance comparison](visuals/jacobian_acceptance_effect.png)

![Animated same-dimensional parameterization jumps, with and without the Jacobian](visuals/parameterization_jacobian_jumps.gif)

PCAT writes `proposal_sequence.gif` for retained states and `proposal_candidates.gif` for every candidate, including rejected proposals. Each candidate frame overlays the current and proposed predictions on the data in one panel and compares their data-minus-model residuals in a second panel. Colors and arrows mark the transition. The split/merge comparison is saved as `jacobian_acceptance_effect.png`. The analytically constructed same-dimensional comparison is saved as `parameterization_jacobian_jumps.gif`.
