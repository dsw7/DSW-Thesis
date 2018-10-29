* David *

This directory is almost equivalent to:
  ~/DSW-Thesis/these_PDB_to_DFT_part2
  
Except here I am not pulling a Met-aromatic interaction from PDB and
narrowing data down benzene / dimethylsulfide system. I am instead
modelling the translation + rotation of H2S away from a benzene ring:

  C - C          H
 /     \        /
C       C  <-> S
 \     /        \
  C - C          H
  
  C - C            H
 /     \          /
C       C  <- -> S
 \     /          \
  C - C            H
  
  C - C              H
 /     \            /
C       C  <- - -> S
 \     /            \
  C - C              H
  
...

The centroid of the benzene ring is the origin of this frame -> (0, 0, 0)

I am using the homogeneous transformation matrix:

    |cos(a)  0  sin(a)   t|
    |     0  1       0   0|
T = |-sin(a) 0  cos(a)   0|
    |     0  0       0   1|
    
To move the H2S molecule away from the benzene by increment t,
and to rotate the lone pairs away from the pi electron cloud
by angle a.

I then perform DFT on all my rotations + translations:

(norm1, a1, HOMO_E1)
(norm2, a2, HOMO_E2)
(norm3, a3, HOMO_E3)
...
(norm_n, a_n, HOMO_E_n)

I fit this data to a 3-dimensional surface.
