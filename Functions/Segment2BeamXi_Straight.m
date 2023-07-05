function [xi, J] = Segment2BeamXi_Straight(xi_tilde, S_elem, S_seg)

   
S1 = S_elem(1);
S2 = S_elem(end);

Sa = S_seg(1);
Sb = S_seg(2);

xi = (xi_tilde*(Sb-Sa) + (Sa+Sb) - (S1+S2))/(S2-S1);
J = (Sb-Sa)/2;
    
