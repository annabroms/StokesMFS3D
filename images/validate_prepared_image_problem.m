function grid = validate_prepared_image_problem(q,opt,mode)
%VALIDATE_PREPARED_IMAGE_PROBLEM Reject stale image geometry or SVD factors.
prepared = opt.prepared;
if ~isequaln(prepared.q,q)
    error('validate_prepared_image_problem:GeometryMismatch', ...
        'Prepared image centres do not match the requested geometry.');
end
plain = rmfield(opt,'prepared');
if ~isequaln(orderfields(plain),orderfields(prepared.options))
    error('validate_prepared_image_problem:OptionsMismatch', ...
        'Prepared image options do not match the requested solve.');
end
grid = prepared.grid;
if ~isequaln(grid.centres,q)
    error('validate_prepared_image_problem:GeometryMismatch', ...
        'Prepared image grid has different centres.');
end
if mode == "resistance"
    factors = prepared.resistanceFactors;
else
    factors = prepared.mobilityFactors;
end
if numel(factors)~=size(q,1) || isempty(factors)
    error('validate_prepared_image_problem:MissingFactors', ...
        'Prepared image factors are missing for this solver.');
end
for k = 1:size(q,1)
    body = grid.body(k);
    ncols = 3*sum(arrayfun(@(b) size(b.points,1),body.basis));
    nrows = 3*size(body.collocation,1);
    if size(factors(k).Y,1)~=ncols || ...
            size(factors(k).U,1)~=nrows
        error('validate_prepared_image_problem:FactorSizeMismatch', ...
            'Prepared image factor dimensions do not match the grid.');
    end
end
end
