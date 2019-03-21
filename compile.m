% mex interface.c
% mex lagrange_remap_1.c 
% mex octovof.c shape.c
% mex PLIC.c
% cmex phitopsi.c
tic;
sError = zeros(4,5,2,2);
vError = zeros(4,5,2,2);
CpTime = zeros(4,5,2,2);
for RM = 1:4
    for M = 1:5
        for am = 1:2
            for sc = 1:2
                %disp([RM M am sc])
                [sError(RM,M,am,sc), vError(RM,M,am,sc), CpTime(RM,M,am,sc)] = advection(RM,M,am,sc);
                disp(CpTime(RM,M,am,sc))
            end
        end
    end
end
toc;
% save images/mat.mat sError vError CpTime
