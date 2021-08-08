# Adjoint-state Traveltime Tomography (Trainning Version)

![Repository size](https://img.shields.io/github/repo-size/MIGG-NTU/ATT_Training)

The code is for adjoint-state traveltime tomography in the Cartesian coordinates.
This version is for trainning purpose only.

## Requirements

In order to run the codes, you have to install the following software:

1.  [GNU Fortran](https://gcc.gnu.org/fortran/) [[Operating system setup tutorials in Chinese](https://seismo-learn.org/seismology101/computer/setup/)]

    ```
    # Fedora
    $ sudo dnf install gcc-gfortran
    
    # CentOS
    $ sudo yum install gfortran
    
    # Ubuntu/Debian
    $ sudo apt install gfortran
    
    # macOS
    $ brew install gfortran
    ```

2.  [Generic Mapping Tools (GMT)](https://www.generic-mapping-tools.org/) [[GMT reference book in Chinese](https://docs.gmt-china.org/latest/install/)]

## Part I: Data

```
$ cd data
```

Prepare your data in `dataInFormat_D_xyz` following the designed format.
I use `syntheticData.f90` to generate synthetic earthquakes and seismic stations.


```
$ cd dataSelection/0_commandCenter/
```

Set the parameters in `parametersData.F90`.
You can select the regions for earthquakes and seismic stations.
The two regions can be different.

```
$ ./workflowStep.sh
```

## Part II: CommandCenter

```
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

```
$ gfortran -o xpara parametersGenerator.F90
$ ./xpara
$ rm xpara
```

## Part III: Model

```
cd ../model
```

If this is for a recovery test, you can set up the target velocity model by editing `velocity3d_true.F90`.

```
$ gfortran -o xvel velocity3d_true.F90
$ ./xvel
$ rm xvel
```

Otherwise, just define the initial model for seismic tomography by editing `velocity3d.F90`.
```
$ gfortran -o xvel velocity3d.F90
$ ./xvel
$ rm xvel
```

## Part IV: Mesh

```
cd ../mesh
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

```
$ gfortran -o xmesh meshgenerator.F90
$ ./xmesh
$ rm xmesh regmesh.mod
```

## Part V: CommandCenter

```
cd ../commandCenter

# Run it if this is a recovery test
$ ./workflow_obstime.sh

$ ./workflow_inversion.sh
```

The obtained velocity model is located at `../inversion/` as `velocity3d015`.

## Part VI: Figure

You can display the results along the cross-section set by `lineEnds`:

```
$ ./plot-cross-section.sh
```

The Bash script will call `zCutVelocity.f90` and the gmt script `vcut.gmt`.
