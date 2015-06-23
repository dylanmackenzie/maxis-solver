Genetic Algorithm for Weighted MaxIS
====================================

This program uses a parallel genetic algorithm to approximate solutions
to the Weighted Maximum Independent Set problem on large graphs. It
implements both single-threaded and multi-threaded versions of the
algorithm found in Beasley and Chu's [A Genetic Algorithm for the Set
Covering
Problem](http://www.sciencedirect.com/science/article/pii/037722179500159X).
I used the
[BHOSLIB](http://www.nlsde.buaa.edu.cn/~kexu/benchmarks/graph-benchmarks.htm)
MIS benchmark graphs for testing, which are non-weighted. The algorithm
generates optimal solutions for most of the 450 vertex graphs. Finding
optimal search parameters for the larger graphs is a work in progress,
but I suspect that 450 vertices is near the upper limit if an optimal
solution is required in a reasonable amount of time. Like all genetic
algorithms, the solver provides no guarantee that a given solution is
optimal.

This solver is designed mostly to explore the efficacy of genetic
algorithms on highly-constrained NP optimization problems. Emphasis was
put on flexibility over performance. For example, the algorithm can use
a multitude of containers to store individuals in the population,
allowing benchmarks of different container types (e.g. vector&lt;bool&gt; vs
vector&lt;char&gt;). It also allows the user to select from several genetic
operators or implement their own by inheriting from a virtual base
class.

Compilation and Usage
---------------------

The solver must be built by a compiler with c++11 support and depends on
boost. If DEBUG is defined during compilation, the solver will log the
state of the population every few thousand iterations.

Run `bin/graph --help` for usage information.

Algorithm
---------

Familiarity with genetic algorithms is assumed for this section.

The fitness of a cantidate solution is the sum of the vertex weights
contained in that solution. For the unweighted version, every vertex is
given a weight of 1. The solver uses the incremental replacement method
and always replaces the weakest memeber of the population. In doing so,
we ensure that the average fitness of the population increases
monotonically at the cost of a higher likelihood of premature
convergence. This is functionally equivalent to a very high level of
elitism.

The minimum set covering problem is closely related to the minimum
vertex covering problem on a graph, whose solution is the complement to
that graph's minimum independent set. It is therefore reasonable to
expect that approaches which are effective on minimum set covering
problems will also work well for finding the maximum independent set
(MIS). Beasely and Chu describe (as far as I can tell) two novel
improvements on a typical genetic algorithm specifically for the minimum
set covering. Both of these are implemented in this solver.

The first is the heuristic feasibility operator. Given a randomly
generated set, such as one that is created by recombinating or mutating
other sets, it is highly likely that some vertices are now over- or
under-covered. A traditional genetic algorithm would solve this by
implementing a penalty function, but such a function is difficult to get
right, and results in yet another opaque parameter to the algorithm.
Instead, Beasely and Chu's heuristic feasibility operator is applied to
each new cantidate solution after recombination and mutation. The
operator is basically a greedy solver for creating maximal independent
sets from a random set. It first eliminates vertices using a heuristic
to determine which are most likely to appear in the MIS. Then it adds
back vertices where it can without invalidating the independent set
constraint. In this case, the heuristic is the weight of the vertex
divided by the number of neighbors it has. The algorithm runs in
O(|V|^2) time, where |V| is the number of vertices in the graph.

This operator ensures that all solutions created by the genetic process
are independent sets. It also increases the chance of generating
duplicate members of the population, thereby reducing genetic diversity.
To mitigate this, we ensure each newly created solution is unique,
rerunning the genetic operators if it is not. This poses some
difficulties for the parallel version of the algorithm (see below).

The other major improvement is the use of a "generalised fitness-based
crossover operator" for recombination. It uses a parent's fitness to
determine the probability that each bit appears in the child. Although
theoretically we could use any number of parents, the algorithm
currently only selects two for efficiency reasons. According to Beasely
and Chu, this operator creates more diverse solutions "when the two
parent solutions are similar in structure than the one-point or
two-point crossover operators." The heuristic feasibility operator makes
genetically similar parents very likely, as most solutions will contain
vertices that the heuristic considers good ones.

We also use Beasely and Chu's variable rate mutator, although it is more
of a minor optimization than a major one.

Implementation Details
----------------------

The solver is written in C++ and makes heavy use of features introduced
in C++11. Currently, it depends on boost's program options, dynamic
bitset, and iterator facade libraries. Most of the implementation
details for the single-threaded case are fairly obvious. We reorder the
input graph to give "good" vertices lower indices which simplifies the
heuristic feasibility operator. The heuristic feasibility operator is
also specialized for boost's dynamic bitset: it runs about an order of
magnitude faster than the naive implementation due to the availability
of bitwise operators. See src/solver.cpp for more details. This requires
a custom iterator for a dynamic bitset that allows it to be used with
the more general genetic operators. Also, we keep the population sorted by
fitness to facilitate the replacement of the weakest vertex.

The parallel version is a bit more complex. We start by assigning a
subset of the population to each available processor. Each core then
performs a constant number of iterations of the genetic algoritm.  Once
every core has performed its required number of iterations, a separate
thread migrates members of the population at random between threads.
Each processor maintains its own hash table for detecting unique
solutions, but uniqueness is not guaranteed across all processors, as
this would incur a lot of overhead as each thread competes for access to
the global hash table. We also implement a custom synchronization class,
to manage all of these threads.

Limitations and Possible Improvements
-------------------------------------

This was programmed in a little over a week and has many imperfections.
The graph class was implemented in a hurry and in total ignorance of the
far superior boost graph library. It only accepts graphs in the DIMACS
ASCII format, and is probably the least optimized and least organized
part of the code base. As of v0.1, the code attempts to maintain const
correctness, but there are still a few places that would benefit from
the use of smart pointers. These issues could be resolved with a bit of
patience.

The genetic operators are implemented with virtual base
classes, meaning they incur the overhead of a virtual function call at
every iteration of the algorithm. If performance is critical, we could
templatize the genetic solver and remove the virtual base classes, but I
like the clarity of interfaces. This performance penalty may eventually
disappear as gcc's devirtualize optimization matures.

I didn't decide to use boost until late in the development process. As a
result, many functions would benefit from boost's zip iterators. Also, I
intended to write unit tests for the graph class and the genetic
operators, but it is likely that the graph class will be completely
reimplemented. This is currently the highest priority.

The synchronization class is confusing and relies on a combination of
atomic counters and condition variables that act like a semaphore. I
would like to reimplement this with a shared mutex so that valgrind's
helgrind tool can properly analyze the code. The sharing of a single
vector of phenotypes between multiple threads without a lock also
unnerves me, although the vector's size is not modified and each thread
only accesses its own section of the vector. Perhaps there is a way to
'share memory by communicating' instead of 'communicating by shared
memory' in the words of Rob Pike.


