function [Rmb,Rbm] = RotationMatrix(alpha)
Rmb = zeros(length(alpha),2,2); %preallocate memory
Rbm = Rmb;
for i=1:length(alpha)
    Rmb(i,:,:) = [cos(alpha(i)),sin(alpha(i)); -sin(alpha(i)),cos(alpha(i))];%transformation matrix R(m->b)
    Rbm(i,:,:) = squeeze(Rmb(i,:,:))'; %transformation matrix R(b->m)
end
end