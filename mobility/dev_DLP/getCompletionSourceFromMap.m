function lambda0 = getCompletionSourceFromMap(F, T, Fmap)
%GETCOMPLETIONSOURCEFROMMAP Construct recompletion source for a force map.
%
%   lambda0 = GETCOMPLETIONSOURCEFROMMAP(F,T,Fmap) returns the minimal-norm
%   source vector satisfying Fmap'*lambda0 = [F; T].

D = Fmap' * Fmap;
ab = D \ [F; T];
lambda0 = Fmap * ab;

end
