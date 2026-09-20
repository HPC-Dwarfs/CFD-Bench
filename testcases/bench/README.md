# Benchmark tiers

Empty on purpose. The small, medium and large tiers land here in their own
change; this directory exists so they have somewhere to go and so the split
between the three kinds of setup is visible from the layout alone.

What belongs here, when it arrives:

  - setups whose grids are sized for scaling runs rather than for physics,
    chosen so the multigrid hierarchy survives the decomposition at the rank
    counts they are meant for;
  - a fixed step count rather than a physical final time, so the work is the
    same whatever the grid;
  - no recorded baselines. These exist to be timed, not compared against a
    stored field. Numerical agreement is the regression setups' job.

See ../flow for the physical setups and ../regression for the shortened
variants those baselines are recorded from.
