# Adjoint-state Traveltime Tomography (Trainning Version)

![Repository size](https://img.shields.io/github/repo-size/MIGG-NTU/ATT_Training)

The code is for adjoint-state traveltime tomography (ATT) in the Cartesian coordinates.
This simple version is mainly for training purpose.

The details of the algorithm can be found at
Tong, P. (2021).
Adjoint-state traveltime tomography: Eikonal equation-based methods and application to the Anza area in southern California.
*Journal of Geophysical Research: Solid Earth*, 126, e2021JB021818.
https://doi.org/10.1029/2021JB021818.

## Requirements

In order to run the codes, you have to install [GNU Fortran](https://gcc.gnu.org/fortran/) [[Operating system setup tutorials in Chinese](https://seismo-learn.org/seismology101/computer/setup/)]:

```bash
# Fedora
$ sudo dnf install gcc-gfortran
    
# CentOS
$ sudo yum install gfortran
    
# Ubuntu/Debian
$ sudo apt install gfortran
    
# macOS
$ brew install gfortran
```

[Generic Mapping Tools (GMT)](https://www.generic-mapping-tools.org/) is used to plot figures.
Refer to [[GMT reference book in Chinese](https://docs.gmt-china.org/latest/install/)] for GMT installation and usage.

## Program document

### Part I: Data

```shell
$ cd ./data/
```

Prepare your data in `dataInFormat_D_xyz` following the designed format.
I use `syntheticData.f90` to generate synthetic earthquakes and seismic stations.


```shell
$ cd dataSelection/0_commandCenter/
```

Set the parameters in `parametersData.F90`.
You can select the regions for earthquakes and seismic stations.
The two regions can be different.

```shell
$ ./workflowStep.sh
```

### Part II: CommandCenter

```shell
$ cd ../../../commandCenter
```

You should have the following files ready at `../data/`:
1. `sources`
2. `receivers`
3. `nevtstaray`
4. `evtstaminmax`
5. `traveltimeReceiverGathers`

Set the parameters in `parametersGenerator.F90`. The dimension of the forward grid,
number of multiple grids, iteration number and some others can be defined here.
You can compile and execute `parametersGenerator.F90` at this time.

```shell
$ gfortran -o xpara parametersGenerator.F90
$ ./xpara
$ rm xpara
```

### Part III: Model

```shell
$ cd ../model
```

If this is for a recovery test, you can set up the target velocity model by editing `velocity3d_true.F90`.

```shell
$ gfortran -o xvel velocity3d_true.F90
$ ./xvel
$ rm xvel
```

Otherwise, just define the initial model for seismic tomography by editing `velocity3d.F90`.
```shell
$ gfortran -o xvel velocity3d.F90
$ ./xvel
$ rm xvel
```

### Part IV: Mesh

```shell
$ cd ../mesh
```

Edit `memeshgenerator.F90` to adjust the size of the inversion grid.
The following three parameters should be properly set.

```Fortran
invx = 12
invy = 3
invz = 11
```

In general, we sample one-wavelength anomaly by 5 grid points or the grid interval
is about one fourth of the anomaly wavelength.

```shell
$ gfortran -o xmesh meshgenerator.F90
$ ./xmesh
$ rm xmesh regmesh.mod
```

### Part V: CommandCenter

```shell
$ cd ../commandCenter

# Run it if this is a recovery test
$ ./workflow_obstime.sh

$ ./workflow_inversion.sh
```

The obtained velocity model is located at `../inversion/` as `velocity3d015`.

### Part VI: Figure

You can display the results along the cross-section set by `lineEnds`:

```shell
# Plot result of checkerboard test
$ ./plot-cross-section-checkerboard.sh
# Plot result of real data
$ ./plot-cross-section-real.sh
```

The Bash script will call `zCutVelocity.f90` and the gmt script `vcut.gmt`.

----

## Examples

### Example1: Calculate and visualize individual kernel

```shell
$ cd ./data/
$ cp ./dataInFormat_D_xyz_individual ./dataInFormat_D_xyz
$ cd ../commandCenter/
$ ./workflow_obstime.sh
$ ./workflow_inversion.sh
```
Then we can plot the result, select model index 1 and velocity perturbation bound 0.06.
```shell
$ cd ../figure/
$ ./plot-cross-section-checkerboard.sh
```
In the image, (a) is velocity model after 1 round of iteration, (b) is the individual kernel, (c) is true velocity model. Red stars denote sources, blue triangles denote receivers.

### Example2: Calculate and visualize event kernel

First, clean the output of previous example.
```shell
$ cd ../commandCenter/
$ ./cleanup.sh
```
Then run the following code.
```shell
$ cd ../data/
$ cp ./dataInFormat_D_xyz_event ./dataInFormat_D_xyz
$ cd ../commandCenter/
$ ./workflow_obstime.sh
$ ./workflow_inversion.sh
```
Then we can plot the result,select model index 1 and velocity perturbation bound 0.06.
```shell
$ cd ../figure/
$ ./plot-cross-section-checkerboard2.sh
```
In the image, (a) is velocity model after 1 round of iteration, (b) is the event kernel, (c) is true velocity model. Red stars denote sources, blue triangles denote receivers.

Note: In this example, we show the event kernel. By definition, it is kernel of single event and multiple receivers.
In order to optimize performance, we exchange receivers and event. According to reciprocity principle, the result is equivalent.
In practice, no need to do such optimization, because number of events are almost always much more than receivers.

### Example3: Calculate and visualize misfit kernel

First, clean the output of previous example.
```shell
$ cd ../commandCenter/
$ ./cleanup.sh
```
Then run the following code.
```shell
$ cd ../data/
$ cp ./dataInFormat_D_xyz_full ./dataInFormat_D_xyz
$ cd ../commandCenter/
$ ./workflow_obstime.sh
$ ./workflow_inversion.sh
```
Then we can plot the result,select model index 1 and velocity perturbation bound 0.06.
``` shell
$ cd ../figure/
$ ./plot-cross-section-checkerboard.sh
```
In the image, (a) is velocity model after 1 round of iteration, (b) is the misfit kernel, (c) is true velocity model. Red stars denote sources, blue triangles denote receivers.

This result shows velocity model after first iteration. We can edit `commandCenter/parametersGenerator.F90` to run more iterations.
```Fortran
! Need more time to finish
niter = 15
```
### Example4 Real data from Parkfield
