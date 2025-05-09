# elphbolt-module-phonon

This is a place for notes on the `phonon.f90` which is located at 
`elphbolt/src/`.


## read_ifc2

This subroutine reads in the 2nd order force constants from the Quantum 
Espresso format. The documentation of QEs PHonon module is 
[here](https://www.quantum-espresso.org/Doc/user_guide_PDF/ph_user_guide.pdf). 
Strangely this is labeled as a user guide but assumes how the phonon module 
of them works.


### Parameters and Variables

The subroutine itself only has two parameters which are the `inout` phonon 
derived type (or class) and an `in` variable of the crystal type.

Quite a number of local variable are declared which I will try to fill up 
with meaning while going through the script

#### 64 bit integers

- `qscell(3)`: qgrid used in the DFPT calculation
- `tipo(crys%numatoms)`: ???
- `t1`: temp `ifc2` tensor location index 5
- `t2`: temp `ifc2` tensor location index 6
- `t3`: temp `ifc2` tensor location index 7
- `i`: dummy index for counting
- `j`: dummy index for counting
- `iat`: temp `ifc2` tensor location index 3 
- `jat`: temp `ifc2` tensor location index 4
- `ibrav`: ???
- `ipol`: temp `ifc2` tensor location index 1
- `jpol`: temp `ifc2` tensor location index 2
- `m1`:
- `m2`:
- `m3`:
- `ntype`: number of atomic species
- `nat`: number of atoms in the basis
- `nfc2`: number of force constants 3*3*`nat`*`nat`

#### 64 bit floating

- `r(crys%numatoms, 3)`: relative coordinates for each atom in basis
- `wscell(3,0:3)`:
- `celldm(6)`: vector of cell dimensions (why 6d?)
- `at(3,3)`: ???
- `mass(crys%numelements)`: mass for each atomic species
- `zeff(crys%numatoms, 3, 3)`: born effective charges
- `eps(3, 3)`: dielectric tensor
- `dnrm2`:
- `massfactor`: species mass conversion factor

#### Characters

- `polar_key`: polar material or not? probably eiter 'T' or 'F'
- `label(crys%numelements)`: label of atomic species

### Routine analysis

Lets go over what the routine does and how it uses the parameters and 
variables.

1. [Allocate 
   arrays](https://fortran-lang.org/learn/best_practices/allocatable_arrays/) 
   `mm` and `rr`. These arrays are declared in the derived type of the 
   `phonon` class. By best practice these *should* be scratch arrays. But I 
   am [reading 
   now](https://fortran-lang.org/learn/quickstart/arrays_strings/) that 
   they are just dynamic arrays.

~~~fortran
allocate(self%mm(crys%numatoms, crys%numatoms))
allocate(self%rr(crys%numatoms, crys%numatoms, 3))
~~~

2. [List-directed 
   I/O](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc5/index.html) 
   is a free-form I/O and it follows the rules given in source on the link. 
   It reads the first row of the ifc2-file with the number of species 
   `ntype`, number of atoms in the basis `nat`, `ibrav` ?? and 
   `celldm`-array

~~~fortran
read(1,*) ntype, nat, ibrav, celldm(1:6)
~~~

3. Special case with short-hand [implied 
   do](https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/io.html). I 
   can only guess that `ibrav` controls the specific lattice or something 
   similar. The implied do is used to read all 9 components of the 
   `at`-array.

~~~fortran
if (ibrav==0) then
   read(1,*) ((at(i,j),i=1,3),j=1,3)
end if
~~~

>[!CAUTION]
> What is `ibrav` and what denotes the special case for `ibrav=0`?

4. Read all species type lines which contain the species `label` and 
   species `mass`

~~~fortran
    do i = 1, ntype
       read(1, *) j, label(i), mass(i)
    end do
    mass = crys%masses/massfactor
~~~

5. For the number of atoms in the basis read `tipo` ?? and the relative 
   coordinates to each other

~~~fortran
do i = 1, nat
   read(1, *) j, tipo(i), r(i, 1:3)
end do
r = transpose(matmul(crys%lattvecs, crys%basis))/bohr2nm
~~~

>[!CAUTION]
> Why does this step override the `r`-array anyways? Why read it in the 
> first place?

6. If it is a polar material read the $\varepsilon$-tensor and for each 
   atom in the basis read the Born-effective charges

~~~fortran
read(1, *) polar_key
if(polar_key == "T") then
   do i = 1, 3
      read(1, *) eps(i, 1:3)
   end do
   do i = 1, nat
      read(1, *)
      do j = 1, 3
         read(1, *) zeff(i, j, 1:3)
      end do
   end do
end if
~~~

7. Which q-grid was used and copy it to `phonon%scell`

~~~fortran
read(1,*) qscell(1:3)
self%scell = qscell
~~~

>[!CAUTION]
> Why copy `qscell` to `scell`?

8. Allocate the `ifc2` tensor of rank 7 and store the number of force 
   constants

~~~fortran
allocate(self%ifc2(3, 3, nat, nat, self%scell(1), self%scell(2), 
self%scell(3)))
nfc2 = 3*3*nat*nat
~~~

9. Read the actual force constants by adresses. Every lines' last value 
   represents one element in the tensor. The rest is just adresses. After 
   this the actual reading of data is completed.
   The tensor has rank 7. 

$$U = U_0 + \frac{1}{2} \sum_{n,n',\alpha,\alpha'} \frac{\partial^2 
U}{\partial \xi^\alpha_n \,\partial\xi^{\alpha'}_{n'}}
\xi^{\alpha}_{n}\xi^{\alpha'}_{n'} + \dots$$


    

~~~fortran
do i = 1, nfc2
   read(1, *) ipol, jpol, iat, jat
   do j = 1, self%scell(1)*self%scell(2)*self%scell(3)
      read(1, *) t1, t2, t3, &
           self%ifc2(ipol, jpol, iat, jat, t1, t2, t3)
   end do
end do
~~~

10. Forcing the conservation of momentum into the `ifc2`-tensor.

~~~fortran
do i = 1, 3
   do j = 1, 3
      do iat = 1, nat
         self%ifc2(i, j, iat, iat, 1, 1, 1) = self%ifc2(i, j, iat, iat, 1, 1, 1) - &
              sum(self%ifc2(i, j, iat, :, :, :, :))
      end do
   end do
end do
~~~


>[!CAUTION]
> The ordering of this loop goes against the column-major nature of Fortran 
> as far as I can tell. It is probably not that import in this part because 
> the subroutine should just read in the file. But still weird.


