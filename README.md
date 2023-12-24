Download the TSP instances.

Extract them and put them in a folder called `tsp`. 

When installed into your environment (for example using Poetry), you can do:

```
tisp <path to tsp .dat>
```

Alternatively, when using it as a library, you can do:

```python
from tisp import cli

cli()
```

Then run the script. By default it looks for the `gr48` instance.

To run a specific instance, do (the path is interpreted relatively):


```python
from tisp import solve_from_path

solve_from_path("gr96.dat")
```

Or, when you have an instance already loaded, with the costs a simple cost matrix of type `list[list[float]]`:

```python
from tisp import do_branch_and_bound

# inst = ... load instance here

do_branch_and_bound(inst)
```

It only works for fully connected instances. An example is:

```
3
0 20 5
20 0 3
5 3 0
```

"3" is the size of the instance, the rest indicates the cost of edge (i, j), where i is the row and j the column.

Supports Python up to 3.11. Uses mip as the underlying LP solving library. By default it uses included CBC binaries.

### Performance

The main peformance bottleneck is the computation of the min cut to find the constraint to be added to the LP, although for larger instances solving the LP's also takes significant time.

On my 5-year old laptop, `gr48` takes under 10 seconds. Computing it directly as an LP with a solver can take 30-50 seconds. For larger instances, this solver performs even better.