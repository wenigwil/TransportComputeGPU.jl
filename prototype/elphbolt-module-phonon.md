# elphbolt-module-phonon

This is a place for notes on the `phonon.f90` which is located at 
`elphbolt/src/`.

## types

### phonon

- **`prefix = 'ph'`** *`string`*
    - prefix for identifying particle type
- **`scell(3)`** *`integer`*
    - number of unitcells in the supercell used for DFPT or finite 
      displacement methods
- **`ifc3(:,:,:,:)`** *`array of floats`*
    - 3rd order force constant tensor with (probably) muxed indices
- **`numrtiplets`** *`integer`*
    - number of triplets in the ifc3-file
    - triplets in the sense analogous to 2nd order force constant atom 
      pairs
- **`R_j(:,:), R_k(:,:)`** *`array of floats`*
    - positions of the 2nd and 3rd unitcell in the supercell for an ifc3 
      triplet
- **`Index_i(:), Index_j(:), Index_k(:)`** *`vector of integers`*
    - labels of primitive cell atoms in the ifc3 triplet
- **`simplex_squared_evals`** *`array of floats`*
    - simplex (???) vertices filled with squared eigenvalues
    - i don't understand any part of this yet
- **`cg_indexlist_irred(:)`** *`vector of integers`*
    - probably a list with group theory contents
- **`cg_indexlist`** *`vector of integers`*
    - probably a list with group theory contents
- **`cgset_indexlist_irred(:)`** *`vector of integers`*
    - probably a list with group theory contents
- **`cgset_indexlist`** *`vector of integers`*
    - probably a list with group theory contents
- **`ifc2(:,:,:,:,:,:,:)`** *`array of floats`*
    - 2nd order force constant tensor
- **`rws(124,0:3)`** *`array of floats`*
    - positions of supercells from stretched lattice vectors of the 
      unitcell
    - stretched by `scell(:)` 
    - supercells are part of a fixed `5x5x5` super-supercell
- **`cell_r(1:3, 0:3)`** *`array of floats`*
    - transposed `crys%lattvecs` (lattvecs are put into inputfile in 
      nanometers) and converted to units of bohr ( by *dividing* with 
      `bohr2nm`)
    - columns of `cell_r` are the cartesian axis and rows number the 
      lattice vectors
    - IMPORTANT: first columns denotes the length of the vector
- **`cell_g(1:3, 0:3)`** *`array of floats`*
    - analogous to `cell_r`
    - transposed of `crys%reclattvecs` which are computed from the 
      `crys%lattvecs`
    - converted from 1/nm to 1/bohr (by multiplying with bohr2nm)
    - first column contains the length of the vector
- **`mm(:,:)`** *`array of floats`*
    - symmetric matrix that contains all masses of all atoms in the 
      unitcell (by the counts of basis vectors *not* species) on the 
      diagonal
    - off-diagonal elements contain the mass factors (sqrt(mi * mj)) of all 
      atom combinations
- **`rr(:,:,:)`** *`array of floats`*
    - tensor that holds the difference of all atom positions in the unit 
      cell
    - can be thought of as a matrix of vectors. Diagonal is the zero-vector
    - off-diagonal elements are the difference in atom positions in the 
      unitcell (for 2 atoms in the unitcell this is of shape *2x2*x3)
- **`ws_cell(:, :)`** *`array of integers`*
    - IMPORTANT: this is *easily* (bad naming) mistaken for `wscell(3)` 
      which contains the numbers of unitcell in the supercell copied from 
      the qe-ifc2-file!!! **This is something else!**
    - takes part in the `dynmat`-symmetry restore process (`wsweights`)
- **`ws_weight(:)`** *`array of floats`*
    - contains the computed `weight` corresponding to each unitcell in the 
      super-supercell
    - length should be `5*5*5*wscell(1)*wscell(2)*wscell(3)` (big numba)

## initialize

This `subroutine` is said to initialize the phonon data type, calculate 
*ground state* phonon properties and read ifc3 data.



## phonon_espresso_precompute

## phonon_espresso



## calculate_phonons


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
- `m1`: iterator for the number of supercells in the super-supercell
- `m2`: iterator for the number of supercells in the super-supercell
- `m3`: iterator for the number of supercells in the super-supercell
- `ntype`: number of atomic species
- `nat`: number of atoms in the basis
- `nfc2`: number of force constants 3*3*`nat`*`nat`

#### 64 bit floating

- `r(crys%numatoms, 3)`: relative coordinates for each atom in basis
- `wscell(3,0:3)`: stretched unitcell lattice vectors. Supercell lattice 
  vectors
- `celldm(6)`: vector of cell dimensions (why 6d?)
- `at(3,3)`: ???
- `mass(crys%numelements)`: mass for each atomic species
- `zeff(crys%numatoms, 3, 3)`: born effective charges
- `eps(3, 3)`: dielectric tensor
- `dnrm2`: LAPACK for euclidean norm of a vector
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
> Why copy `qscell` to `scell`? **ANSWER: `scell` is the number of unit 
> cells, the supercell is made up of! `scell` is used later at step 12**

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


11. Grab reciprocal and direct lattice vectors, convert them to nanometers 
    and save them into arrays. The lattice vectors are stored in columns. 
    Rows goes through the 3 dimensions. The 0-column denotes the length of 
    the lattice vectors

~~~fortran
self%cell_r(:, 1:3) = transpose(crys%lattvecs)/bohr2nm
do i = 1, 3
   self%cell_r(i, 0) = dnrm2(3, self%cell_r(i, 1:3), 1)
end do
self%cell_g(:, 1:3) = transpose(crys%reclattvecs)*bohr2nm
do i = 1, 3
   self%cell_g(i, 0) = dnrm2(3, self%cell_g(i, 1:3), 1)
end do
~~~

12. stretch the unitcell lattice vector by the number of unitcells in each 
    direction which make up the supercell (`scell` contains the number of 
    unit cells in each of the lattice vector directions)

~~~fortran
wscell(1,1:3) = self%cell_r(1,1:3)*self%scell(1)
wscell(2,1:3) = self%cell_r(2,1:3)*self%scell(2)
wscell(3,1:3) = self%cell_r(3,1:3)*self%scell(3)
~~~

13. The supercell is made up of unitcells. Now we build data corresponding 
    to a super-supercell that is made up of 5 supercells in each lattice 
    vector (stretched and contained in `wscell`-array) direction. We are 
    storing all vectors that point to a specific supercell in the 
    super-supercell (columns `1:3` of `self%rws`) AND *half of their 
    squared lengths* (column `0`).

~~~fortran
j = 1
do m1 = -2, 2
   do m2 = -2, 2
      do m3 = -2, 2
         if(all((/m1, m2, m3/).eq.0)) then
            cycle
         end if
         do i = 1, 3
            self%rws(j, i) = wscell(1, i)*m1 + wscell(2, i)*m2 + wscell(3, i)*m3
         end do
         self%rws(j, 0) = 0.5*dot_product(self%rws(j, 1:3), self%rws(j, 1:3))
         j = j + 1
      end do
   end do
end do
~~~

14. We are building symmetric matrices. They hold **a)** Masses and Mass 
    factors (`self%mm`). Diagonal is the mass of atom `i`. All other 
    elements are the mass factors of atoms `i` and `j`. **b)** Distances 
    between atom `i` and `j` (`self%rr`). Diagonal is `0`.

~~~fortran
do i = 1, nat
   self%mm(i, i) = mass(tipo(i))
   self%rr(i, i, :) = 0
   do j = i + 1, nat
      self%mm(i, j) = sqrt(mass(tipo(i))*mass(tipo(j)))
      self%rr(i, j, 1:3) = r(i, 1:3) - r(j, 1:3)
      self%mm(j, i) = self%mm(i, j)
      self%rr(j, i, 1:3) = -self%rr(i, j, 1:3)
   end do
end do
~~~

## read_ifc3

## deallocate_phonon_quantities

## save_xmassvar

## allocate_xmassvar

## clean_xmassvar

