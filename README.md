# mol2surface
Fortran code that convert xyz-file to POSCAR file for VASP calculation <br />
Compile the code:
```markdown
gfortran mol2surface.f90 -o mol2surface.exe
```
For utilisation give the name of the xyz file:
```markdown
./mol2surface.exe file.xyz
```
The code read also the info_vasp file that contains:
Values of the box parameters and angles <br />
a= 26.48795                !a along x-axis <br />
b= 26.3264              !b along y-axis <br />
c= 3.407662            !c along z-axis <br />
alpha= 90              !alpha angle    <br />
beta= 90               !beta angle     <br />
gamma= 63.0            !gamma angle   <br />
Number of type of atoms in the xyz-file <br />
4                         !then give the type of atom<br />
N <br />
C <br />
H <br />
O <br />
Then the code provide a POSCAR file that you can use for VASP calculations.
