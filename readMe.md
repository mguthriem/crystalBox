# crystalBox

Contains a useful class `Box` for crystalline materials using mantid

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
|ticksWS(yVal) | Creates a mantid workspace containing points at values of d-spacings and heights yVal|
|plot()*        | Creates a mantid plot with tick marks at positions of d-spacings | 
|scaleLattice(float)  | Multiplies lattice parameters by constant `float`, recalculates d-spacings and updates plot|
|resetLattice()       | Returns lattice to original values in cif file|

(*Plot output can be tweaked with attributes `markerHeight` and `markeredgecolor`)

## setting up 

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
## instantiating a crystal box

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
You can change a the nickname by editing the attribute `nickName`:
```
sample1.nickName = 'a better nickname'

