function result=run_synthetic_example()
% Illustrative calculations only; no figure/data file from manuscript.
n=57;dx=3.67712e-3;
[X,Y]=meshgrid((0:n-1)*dx,(0:n-1)*dx);
w=1.2; U=-w*(Y-mean(Y(:)));V=w*(X-mean(X(:)));
[R,Wz,~,~]=pof.planarLiutex(U,V,dx,dx,0.8);
[T,mask,p80]=pof.detectVortices(R,Wz,X,Y,dx,dx); %#ok<ASGLU>
[Pi,valid,info]=pof.horizontalFlux(U,V,dx,dx,20.43e-3);
[r,S3]=pof.structureFunctionS3(U,V,dx,dx,14);
result=struct('nDetected',height(T),'P80',p80,'validPi',nnz(valid),...
 'ellFWHM_m',info.ellFWHM_m,'radialCenters_m',r,'S3',S3,...
 'isManuscriptData',false);
disp(result);
end
