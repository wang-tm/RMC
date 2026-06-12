function z = fast_multiply_left(U,Vt,D,y)

z = U*(Vt*y) - D*y;