function [] = drawedgescalar(M, N, mark, folder, name)

        % Get gauss grids
        filename = strcat(folder, '/build/vertgaussgridx.dat');
		  fileID = fopen(filename, 'r');
		  vertgx = fscanf(fileID, '%f', [1,Inf]);
        vertgx = reshape(vertgx, 3, N*(M+1));
        revertgx = vertgx(:,1:M+1);
        for i= 1:N-1 
        revertgx = [revertgx;vertgx(:,i*(M+1)+1:(i+1)*(M+1))];
        end

        filename = strcat(folder, '/build/vertgaussgridy.dat');
		  fileID = fopen(filename, 'r');
		  vertgy = fscanf(fileID, '%f', [1,Inf]);
        vertgy = reshape(vertgy, 3, N*(M+1));
        revertgy = vertgy(:,1:M+1);
        for i= 1:N-1 
        revertgy = [revertgy;vertgy(:,i*(M+1)+1:(i+1)*(M+1))];
        end

        filename = strcat(folder, '/build/horigaussgridx.dat');
		  fileID = fopen(filename, 'r');
		  horigx = fscanf(fileID, '%f', [1,Inf]);
        horigx = reshape(horigx, M*3, N+1);

        filename = strcat(folder, '/build/horigaussgridy.dat');
		  fileID = fopen(filename, 'r');
		  horigy = fscanf(fileID, '%f', [1,Inf]);
        horigy = reshape(horigy, M*3, N+1);

        vertdof = N*(M+1)*3;
        horidof = M*(N+1)*3;

		  filename = strcat(folder,'/build/',name, string(mark),'.dat');
		  fileID = fopen(filename, 'r');
		  val = fscanf(fileID, '%f', [1,Inf]);
        valvert = val(1:vertdof);
		  valvert = reshape(valvert, 3, N*(M+1));
        revalvert = valvert(:,1:M+1);
        for i= 1:N-1 
        revalvert = [revalvert;valvert(:,i*(M+1)+1:(i+1)*(M+1))];
        end

		  valhori = val(vertdof+1 : vertdof+horidof);
		  valhori = reshape(valhori, M*3,N+1);

		  figure
        surf(revertgx, revertgy, revalvert); 
		  title('Vertical Gauss Points');
		  figure
        surf(horigx, horigy, valhori); 
		  title('Horizontal Gauss Points');
       
        figure 
	     surf(revertgx, revertgy, revalvert); 
	     hold on
        surf(horigx, horigy, valhori); 
		  title('All Gauss Points');
        hold off 
