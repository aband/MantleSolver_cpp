

function [] = plot_edge_velocity(M, N, mark, folder, name)

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

        % Compute number of dofs
		  vertdof = N*(M+1)*3;
        horidof = M*(N+1)*3;

        % Get x component
		  %filename = strcat(folder,'/build/',name,'x', string(mark),'.dat');
		  filename = strcat(folder,'/build/',name,'x', int2str(mark),'.dat');

		  fileID = fopen(filename, 'r');
		  val = fscanf(fileID, '%f', [1,Inf]);
        valvert = val(1:vertdof);
		  valvert = reshape(valvert, 3, N*(M+1));
        revalvertx = valvert(:,1:M+1);
        for i= 1:N-1 
        		 revalvertx = [revalvertx;valvert(:,i*(M+1)+1:(i+1)*(M+1))];
        end

        valhorix = val(vertdof+1 : vertdof+horidof);
		  valhorix = reshape(valhorix, M*3,N+1);

        % Get y component
		  %filename = strcat(folder,'/build/',name,'y', string(mark),'.dat');
		  filename = strcat(folder,'/build/',name,'y', int2str(mark),'.dat');

		  fileID = fopen(filename, 'r');
		  val = fscanf(fileID, '%f', [1,Inf]);
        valvert = val(1:vertdof);
		  valvert = reshape(valvert, 3, N*(M+1));
        revalverty = valvert(:,1:M+1);
        for i= 1:N-1 
        		 revalverty = [revalverty;valvert(:,i*(M+1)+1:(i+1)*(M+1))];
        end

        valhoriy = val(vertdof+1 : vertdof+horidof);
		  valhoriy = reshape(valhoriy, M*3,N+1);

		  figure
        quiver(revertgx, revertgy, revalvertx, revalverty); 
		  hold on
        l = streamslice(revertgx,revertgy,revalvertx,revalverty,1);
        %l = streamline(revertgx,revertgy,revalvertx,revalverty);

        set(l,'LineWidth',2);
        set(l,'Color','k')
        hold off
		  title('Vertical Gauss Points');

		  figure
        quiver(horigx, horigy, valhorix, valhoriy); 
		  hold on
        l = streamslice(horigx',horigy',valhorix',valhoriy',1);
        %l = streamline(horigx',horigy',valhorix',valhoriy');

        set(l,'LineWidth',2);
        set(l,'Color','k')
        hold off

		  title('Horizontal Gauss Points');
 
        figure
        quiver(revertgx, revertgy, revalvertx, revalverty); 
		  hold on
        quiver(horigx, horigy, valhorix, valhoriy); 
        hold off	
		  title('Plot Velocities on Gauss Points All together');
