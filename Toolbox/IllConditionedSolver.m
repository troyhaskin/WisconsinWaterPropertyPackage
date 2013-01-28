function x = IllConditionedSolver(A,b)
    
    % Grab the parallel row and produce an orthogonal row
    aPara = A(end,:);
    aOrth = GetOrthogonalVector(aPara);
    
    % Form the corrected, full-rank matrix, it's determinant, and  it's adjoint
    Ap    = [A(1:end-1,1:end);aOrth];
    DetAp = det(Ap);
    AdjAp = Adjugate(Ap);
    
    % Grab the RHS of the parallel row
    bN  = b(end);

    % Produce an RHS for the orthogonal row such that InvA*b == InvAp*bp == x
    Sum = 0;
    N   = length(b);
    for j = 1 : N
        AdjApDotB = AdjAp(j,1:N-1) * b(1:N-1);
        Sum       = Sum + aPara(j)*AdjApDotB;
    end
    bNp = (DetAp*bN - Sum)/sum(aPara.*AdjAp(:,N)');
    bp  = [b(1:N-1);bNp];
    
    % Form the inverse
    InvAp = AdjA/DetAp; % Already have the adjoint, might as well use it.
    
    % Solve
    x     = InvAp * bp;
    
end

function Orth = GetOrthogonalVector(ParallelRow)
    
    N    = length(ParallelRow);
    Mask = 1:N-1;
    
    OrthPrime = (2*rand(1,N-1)-1) .* ParallelRow(Mask);
    OrthLast  = -sum(OrthPrime.*ParallelRow(Mask))/ParallelRow(N);
    
    Orth = [OrthPrime,OrthLast];
end
