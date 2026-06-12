function z = fast_multiply_right(V,Ut,Dt,y)

z = V*(Ut*y) - Dt*y;