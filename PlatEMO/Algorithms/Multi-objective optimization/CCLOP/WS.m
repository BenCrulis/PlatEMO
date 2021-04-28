function [U] = WS(V)
    % compute the WS transformation introduced in MOEA/D-AWA
    % function written by Ben Crulis

    U = 1 ./ V;
    S = sum(U, 2);
    U = U ./ repmat(S, [1,3]);
end

