function [] = drawcellscalar(M,N,mark,folder,name)

		  % Get cell center grid
		  filename = strcat(folder, '/build/cellgridx.dat');
		  fileID = fopen(filename, 'r');
		  vertexx = fscanf(fileID, '%f', [1,Inf]);
		  vertexx = reshape(vertexx, M, M);

		  filename = strcat(folder, '/build/cellgridy.dat');
		  fileID = fopen(filename, 'r');
		  vertexy = fscanf(fileID, '%f', [1,Inf]);
		  vertexy = reshape(vertexy, M, M);
