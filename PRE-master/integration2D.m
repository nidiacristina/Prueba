function I = integration2D(Np,fvec,Jx,Jy)
% Integrates in the reference 2D element.
% WORKS CURRENTLY ONLY FOR THE SAME NUMBER OF MODES IN BOTH
% DIMENSIONS.

% Inputs:
% Np     : Integration points.
% f      : RHS in physical space, evaluated at QxQ points
% f      : f IS IN LEXICOGRAPHICAL ORDER VECTOR (LOCAL!!)
% Jx, Jy : Jacobians

% Outputs:
% I      : 2D integral 


f = vector_to_matrix_lex( fvec,Np,Np );

%% Grid and bases. 
[~,omgx] = gll(Np);
omgy = omgx;

%% Evaluation Kernel


integralx = zeros(1,Np);
indx = 0;
for j = Np:-1:1
    indx = indx + 1;
    integralx(indx) = f(j,:)*omgx;
end

I = integralx * omgy ;

%%
I = Jx*Jy*I;

end % end of function