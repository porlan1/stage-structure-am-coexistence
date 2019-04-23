# stage-structure-am-coexistence

## Requirements:
* gsl (https://www.gnu.org/software/gsl/)
* gnuplot (http://www.gnuplot.info/)

## Compiling files to run:
`gcc -lgsl -lm -lblas -o O ODE_system.c -L /usr/local/lib -I /usr/local/include/`


## Running program:
`./O`

## Plotting results
plots are generated using gnuplot (see GNU_PLOTS_MASTER*)

## Example of workflow
```
gcc -lgsl -lm -lblas -o O ODE_system.c -L /usr/local/lib -I /usr/local/include/

./O > CASE0_sp1_alone.dat

gnuplot example.plt
```

I have the terminal set to svg, because other options were not working, the last step could be:

```
gnuplot example.plt > plot.svg
```

you can open the plot.svg file in a web browser
