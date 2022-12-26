function func=iterative_mat_gen(K_mat,mesh_num,mesh_size,order_N)

func=eye(mesh_num);
if (order_N~=0)
    for m=1:1:order_N
        tmp01=mesh_size^m;
        tmp02=K_mat^m;
        func=func+tmp01*tmp02;
    end
end