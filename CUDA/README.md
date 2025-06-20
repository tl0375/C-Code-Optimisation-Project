# A simple CFD application solving a lid-cavity problem

This is a simple 2D CFD simulation application solving the lid-driven cavity problem.

The set-up is a square domain (2.0 x 2.0), with Dirichlet boundary conditions on all sides. Three sides are stationary, while one side is moving with a constant velocity tangent to the side.

At the end of execution, a VTK file is produced for visualisation purposes. This can be loaded into VisIt for analysis.

## Building

The application can be built with the provided Makefile. e.g.

```
$ make
```

This will build a lid-cavity binary.

## Running

The application can be ran in its default configuration with:

```
$ ./lid-cavity
```

This will run with a 500 x 500 grid and will output status every 100 iterations. At the end of execution a VTK file will be produced. 

There are numerous other options available. These can be queried with:

```
$ ./lid-cavity --help
```

If you would like to visualise the entire execution, you should enable checkpointing and increase the output frequency. You can do this like so:

```
$ mkdir out
$ ./lid-cavity -x 50 -c -f 50 -o out/my_sim -n 1000
```

This will run a 50 x 50 problem (perfect for manageable visualisation), enable checkpointing every 50 time steps, and will output the VTK files to the out directory. It will run for 1000 steps.