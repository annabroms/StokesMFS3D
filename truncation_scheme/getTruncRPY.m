function R = getTruncRPY(P,q,rin_base,a,Lcut)

M = size(rin_base,1);
R = zeros(P*M*3); 

for k = 1:P
    d = vecnorm(q(k,:)-q,2,2);
    d = d-2;
    ind = find(d<Lcut);

    
   
    rin = [];
    store_ind = [];
    for j = 1:length(ind)
        i = ind(j);
        rin = [rin; rin_base+q(i,:)];
        store_ind = [store_ind; (i-1)*3*M+1:i*3*M];
    end
    
    R(store_ind,store_ind) = generate_blob_matrix_vec(rin,a);

   
       
%     figure()
%     scatter3(rin(:,1),rin(:,2),rin(:,3));
%     hold on
%     scatter3(rout(:,1),rout(:,2),rout(:,3));

end

end