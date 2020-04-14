/+  roessler, lazytrig
!:
::  3D Projection of Rössler Strange Attractor System
::  
::  By ~lagrev-nocfep (N E Davis)
::
::  Mathematical implementation based on
::    https://en.wikipedia.org/wiki/3D_projection
::
|=  [np=@ud]
=<  (gen3d np)
::=<  (genxyzs:roessler np .0.01)
|%
++  gen3d
  |=  [np=@ud]  ^-  (list (list @rs))
  =/  pts  (genxyz-pts:roessler np .0.001)
  =/  cx  .50
  =/  cy  .-50
  =/  cz  .50
  =/  theta-x  .0.0
  =/  theta-y  .0.0
  =/  theta-z  .0.0
  =/  ex  .10
  =/  ey  .10
  =/  ez  .10
  =/  index  0
  =/  acc=(list (list @rs))  ~
  |-  ^-  (list (list @rs))
  ?:  =(index np)  acc
  =/  pt  (snag index pts)
  =/  at  i.-:pt
  =/  ax  i.+<:pt
  =/  ay  i.+>-:pt
  =/  az  i.+>+<:pt
  =/  vt  (render-3d-point ax ay az cx cy cz theta-x theta-y theta-z ex ey ez)
  $(index +(index), acc `(list (list @rs))`[vt acc])
++  render-3d-point
  |=  [ax=@rs ay=@rs az=@rs cx=@rs cy=@rs cz=@rs theta-x=@rs theta-y=@rs theta-z=@rs ex=@rs ey=@rs ez=@rs]  ^-  (list @rs)
  =/  sin-x  (sine:lazytrig theta-x)
  =/  sin-y  (sine:lazytrig theta-y)
  =/  sin-z  (sine:lazytrig theta-z)
  =/  cos-x  (cosine:lazytrig theta-x)
  =/  cos-y  (cosine:lazytrig theta-y)
  =/  cos-z  (cosine:lazytrig theta-z)
  =/  xx  (sub:rs ax cx)
  =/  yy  (sub:rs ay cy)
  =/  zz  (sub:rs az cz)
  =/  dx  (sub:rs (mul:rs cos-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))) (mul:rs sin-y zz))
  =/  dy  (add:rs (mul:rs sin-x (add:rs (mul:rs cos-y zz) (mul:rs sin-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))))) (mul:rs cos-x (add:rs (mul:rs cos-z yy) (mul:rs sin-z xx))))
  =/  dz  (sub:rs (mul:rs cos-x (add:rs (mul:rs cos-y zz) (mul:rs sin-y (add:rs (mul:rs sin-z yy) (mul:rs cos-z xx))))) (mul:rs sin-x (add:rs (mul:rs cos-z yy) (mul:rs sin-z xx))))
  =/  bx  (add:rs (mul:rs ez (div:rs dx dz)) ex)
  =/  by  (add:rs (mul:rs ez (div:rs dy dz)) ey)
  ::~&  "{<ax>} {<ay>} {<az>} -> {<bx>} {<by>}"
  ~[bx by]
--
