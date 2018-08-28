# Wyznaczamy środek Polski

This repository consists of the source code used to estimate the centre
of the administrative area of Poland in article
[„Wyznaczamy środek Polski”](http://www.deltami.edu.pl/temat/matematyka/zastosowania/2017/12/30/Wyznaczamy_srodek_Polski/)
(published in *Delta* monthly, January 2018).

The source code may be used to calculate centre of mass of various areas
on the ellipsoid of reference, provided the accurate geodesic data
are available.

## Requirements

* C++ compiler with support for C++ 2011 standard
* installed and configured MPI environment (OpenMPI, MPICH or similar)
* installed GeographicLib and shapefile libraries

## Compilation

`make`

## Running

To reproduce results from the article, “Państwo.*” files included in data
available on the website of Państwowy Rejestr Granic
[http://www.gugik.gov.pl/pzgik/dane-bez-oplat/dane-z-panstwowego-rejestru-granic-i-powierzchni-jednostek-podzialow-terytorialnych-kraju-prg]
have to be downloaded and placed in “data” subdirectory.
Various versions of these data are also available from [http://srodekpolski.pl].

```
USAGE: ./srodek N
```

Program needs one command-line argument *N*, which is the number of grid cells
in one dimension. Total number of grid points (labeled *M* in the article)
will be equal to *M* = *N*<sup>2</sup>.

### Example

(for the original data set used for publication)

```
./srodek 10
    10  312679.529 52.1143376218 19.4236714921 -4083.578 52°06′51.62″ 19°25′25.22″
```

Subsequent values denote:

* chosen value of *N*
* administrative area of Poland (m<sup>2</sup>)
* estimated latitude of the centre (°)
* estimated longitude of the centre (°)
* estimated altitude of the centre relative to the ellipsoid of reference (m)
* latitude (repeated) in degrees-minutes-seconds
* longitude (repeated) in degrees-minutes-seconds
