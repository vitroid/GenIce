# Physics test

In this section, we test the randomness of the generated hydrogen bond orientations by GenIce package.

## Ring phase statistics

Every ice structure has cyclic paths of the hydrogen bonds, and if we regard the bonds as arrows directing from a donor to an acceptor, there are a variety of the orientations of the arrows.

For example, a closed path of hydrogen bonds in which all the hydrogen bonds are oriented in the same way is called "a homodromic cycle".

![homodromic 8-cycle](https://upload.wikimedia.org/wikipedia/commons/5/50/DC8.png "An example of the homodromic cycle.")

In case of 6-membered rings, the appearance probability of a homodromic ring is 17.5% if the hydrogen bond network is totally randomized. If some correlations remain in the randomized network, the probability of finding a homodromic cycle would deviate from this value. Ideal probability of finding a homodromic cycle (and other type of ring phases) can be estimated quite precisely using the Pauling-like approximation.

The tests `test_ring*.done` prepare various ice structures and evaluate the deviation from the ideal distribution using Kullback-Leibler divergence.  The test shows that GenIce generates the distribution very close to the ideal one.

Ring phase statistics is useful to check the intermediate-range order in the size of rings.


## Long range dipole-dipole correlation

Radial Kirkwood G function is useful to check the (non-)existence of the long-range dipole-dipole correlation, that may remain if the hydrogen bond orientations are not well shuffled.

Radial Kirkwood G function is defined in the following equation.

    g(r) = < u_i Σ u_j >
	
where u is a normalized molecular dipole vector and the summation is taken over the molecules within r from molecule i.

If the structure is well-shuffled and net dipole moment is zero, g(r) smoothly approaches to zero in the long distance limit.  If net dipole is not fixed to zero, it typically converges to about 3 in average.  If the structure is hydrogen-ordered, on the other hand, the g(r) function oscilates with a very large amplitude.  Thus the smoothness of g(r) indicates how the dipoles are correlated.

