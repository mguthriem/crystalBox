# crystalBox

Contains a useful class `Box` for holding crystallographic values derived from a CIF file

## Use case

This is intended to be used within scripts run from [mantidworkbench](https://www.mantidproject.org/) a software framework for analysing neutron data. 

## What does Box do?

Once instantiated, via a CIF (Crystallographic Information File), a Box holds useful crystallographic information and derives several list of reflection properties including:

|Attribute | Description |
|----      | ----        |
|nRef      | number of reflections|
|hkl       | list of Miller indices|
|dSpacing  | list of d-spacings |
|fSq       | list of structure factors|
|mult      | list of multiplicities   |
|estInt    | list of estimated intendities|

It also has the following functions:

|Function | Description |
|----      | ----        |
|summary()     | Prints short summary of crystallographic info|
|ticksWS(_yval_) | Creates a mantid workspace containing points at values of d-spacings and heights yVal|
|plot()        | Creates a mantid plot with tick marks at positions of d-spacings | 
|scaleLattice(float)  | Multiplies lattice parameters by constant `float`, recalculates d-spacings and updates plot|
|resetLattice()       | Returns lattice to original values in cif file|

## Setting up 

`crystalBox` uses an assigned folder to store cif files. The location of this folder should be specified with a full path on line 10. A `.csv` called `nickNames.csv` should be created in that folder containing entries like:
```
nickname1, sample1.cif
nickname2, sample2.cif
```
where sample1.cif (etc) are cif files stored in the assigned folder.

After setting up that folder, make sure crystalBox is available on your path, then you can import it into a python script using:
```
import crystalBox as crystal
```
## Instantiating a Box

A Box is instantiated by specifying a CIF file. If only a filename (including extension `.cif`) is given, `crystalBox` will search for this file in both the current directory _and_ the assigned cif folder. 
```
sample1 = crystal.Box('sample1.cif')
```
If a nickname has been set up, you can also use this e.g.
```
sample1 = crystal.Box('nickname1')
```
if you use a nickname, then this is remembered and subsequently used to name the tick workspace if you create it. 
TODO: it should use the filename of a cif file as a nickname

You can change the nickname by editing the attribute `nickName`:
```
sample1.nickName = 'a better nickname'

`crystalBox` has an additional function `showNicknames` to display the list of defined nicknames:
```
crystal.showNicknames()
```

## Using a crystal box

### Access and manipulate some handy crystallographic attributes

Once instantiated any of the attributes listed above are available. e.g. you can print the list of d-spacings
```
for peak in sample1.dSpacing:
   print(peak)
```

You can apply a simple scale factor to lattice parameters using the function `scaleLattice`:

```
sample1.scaleLattice(0.99)
```
will multiply the _a_,_b_ and _c_ lattice parameters by 0.99 and recalculate all d-spacings. If a tick workspace exists it will be updated.

If you want to reset the lattice parameters to the values in the original cif file, this can be done using:
```
sample1.resetLattice()
```
### Ticks workspace

It is often useful to display markers at expected peak positions when inspecting powder diffraction data. This can be done with crystalBox, in mantid workbench, by creating a "ticks" workspace. It contains a set of data points with x-values at the expected d-spacings of Bragg reflections and constant y-values at a user specified position on the y-axis of _yval_.

This is done with the ticksWS method of Box:
```
sample1.ticksWS(yval)
```
Repeating this command with a different _yval_ allows adjustment of the height of the data points. Any existing plots of ticks will be updated.

### Plotting a ticks workspace

Once the ticks workspace exists, you can also create a plot where the positions of Bragg reflections are marked with vertical markers using the `plot()` method. Some simple attributes are available to customise these prior to creating the plot:
```
sample1.markerHeight = 15
sample1.markerColor = 'blue'
sample1.plot()
```
In mantid workbench you can then drag one or more workspaces containing diffraction data onto the plot workspace
