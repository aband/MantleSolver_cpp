% Test generated linear system with matlab first

fileID = fopen('build/MatrixCheck.dat','r');

A = fscanf(fileID, '%f', [24,Inf]);

fclose(fileID);
