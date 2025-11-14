function [] = drawedgescalar(M, N, mark, folder, name)

        % Get gauss grids
        filename = strcat(folder, '/build/vertgaussgridx.dat');
		  fileID = fopen(filename, 'r');
		  vertgx = fscanf(fileID, '%f', [1,Inf]);
        vertgx = reshape(vertgx', (M+1)*3, N)

        filename = strcat(folder, '/build/vertgaussgridy.dat');
		  fileID = fopen(filename, 'r');
		  vertgy = fscanf(fileID, '%f', [1,Inf])
        vertgy = reshape(vertgy, (M+1)*3, N)

        filename = strcat(folder, '/build/horigaussgridx.dat');
		  fileID = fopen(filename, 'r');
		  horigx = fscanf(fileID, '%f', [1,Inf]);
        horigx = reshape(horigx, M*3, N+1);

        filename = strcat(folder, '/build/horigaussgridy.dat');
		  fileID = fopen(filename, 'r');
		  horigy = fscanf(fileID, '%f', [1,Inf]);
        horigy = reshape(horigy, M*3, N+1);

		  filename = strcat(folder,'/build/',name, string(mark),'.dat');
		  fileID = fopen(filename, 'r');
		  val = fscanf(fileID, '%f', [1,Inf]);

        vertdof = N*(M+1)*3;

        horidof = M*(N+1)*3;

        valvert = val(1:vertdof);
		  valvert = reshape(valvert, (M+1)*3,N)

		  valhori = val(vertdof+1 : vertdof+horidof);
		  valhori = reshape(valhori, M*3,N+1);

		  figure
        surf(vertgx', vertgy', valvert'); 
		  figure
        surf(horigx', horigy', valhori'); 
        
