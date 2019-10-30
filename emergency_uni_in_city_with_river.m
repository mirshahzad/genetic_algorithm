function emergencyresponse()
 
%=============================================================================
%Problem: Develop a genetic algorithm for the problem described in Part 1
%         assuming that there is a river that divides the city into two parts, West
%         and East, at x = 5 km. West and East are connected by a bridge located at
%         x = 5 km and y = 5.5 km.
% 
          %Find the optimal locaiton of the emergency response unit and comapre it
%         with the one obtained in Part 1.
%=============================================================================')


num=49;    % Number of generation
 
rand('seed',1.4929e+009);
global city_location
city_location=[0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;
               0.5 1; 1.5 2; 2.5 3; 3.5 4; 4.5 5; 5.5 6; 6.5 7;];

       
weights = [
           3; 4; 1; 2; 1; 3; 8;
           2; 1; 3; 1; 3; 9; 7;
           5; 1; 2; 4; 4; 9; 8;
           4; 2; 1; 1; 2; 5; 9;
           8; 9; 6; 3; 2; 8; 7;
           9; 8; 5; 2; 1; 7; 9;
           8; 9; 6; 1; 1; 8; 9;
           ];

 
city_distance = dist(city_location')
city_distance = weights' .* city_distance
 
nind=100;    % Size of a chromosome population
ngenes=49;  % Number of genes in a chromosome
nvar=2;     % Number of variables
Pc=0.9;      % Crossover probability
Pm=0.001;     % Mutation probability
ngener=49;   % Number of generations
n_show=10;   % Number of generations between showing the progress
xymin=0.5;   % Possible minimum values of parameters "x" and "y"
xymax=6.5;    % Possible maximum values of parameters "x" and "y"


disp(' ')
fprintf(1,' nind=%.0f;    Size of the chromosome population\n',nind);
fprintf(1,' Pc=%.1f;      Crossover probability\n',Pc);
fprintf(1,' Pm=%.3f;    Mutation probability\n',Pm);
fprintf(1,' ngener=%.0f;   Number of generations\n',ngener);
fprintf(1,' n_show=%.0f;   Number of generations between showing the progress\n',n_show);
disp(' ')
 
fprintf(1,'Hit any key to generate a population of %.0f chromosomes.\n',nind);
pause

% chrom = rand([1 49],nind,49)
chrom=[];
 
for k=1:nind
   num=ngenes; city_array=1:num; xxx=[];
   for n=1:num
      a=rand(1);
      for i=1:num
         if a<i/num
            xxx=[xxx city_array(i)];
            break
         end
      end
      city_array(i)=[]; 
      num=num-1;
   end
   chrom(k,:)=xxx;
end
% chrom
rout= chrom;
 
% Calculate the chromosome fitness
ObjV=evalObjFun(rout,city_distance,nind,ngenes);
best=min(ObjV);
ave=mean(ObjV);
[a b]=min(ObjV);
 
%[a, b]=min(ObjV);
for m=1:(ngener/n_show)
   for i=1:n_show
   
      % Fitness evaluation
      fitness=(1./ObjV)';
   
      % Roulette wheel selection
      numsel=round(nind*0.9); 
      % The number of chromosomes to be selected for reproduction
      cumfit=repmat(cumsum(fitness),1,numsel);
      chance=repmat(rand(1,numsel),nind,1)*cumfit(nind,1);
      [selind,j]=find(chance < cumfit & chance >= [zeros(1,numsel);cumfit(1:nind-1,:)]);
      newchrom=chrom(selind,:);
   
      % Crossover Probability
      points=round(rand(floor(numsel/2),1).*(ngenes-1))+1;
      points=[points round(rand(floor(numsel/2),1).*(ngenes-1))+1];
      points=sort((points*(rand(1)<Pc)),2);
      for j=1:length(points(:,1))   
         swap_sect=newchrom(2*j-1:2*j,points(j,1)+1:points(j,2));
         remain_sect=newchrom(2*j-1:2*j,:);
         for k=1:ngenes
            for n=1:length(swap_sect(1,:))
               if newchrom(2*j-1,k)==swap_sect(2,n);
                  remain_sect(1,k)=0;
               end
               if newchrom(2*j,k)==swap_sect(1,n);
                  remain_sect(2,k)=0;
               end
            end
         end
         [a b c1]=find(remain_sect(1,:)); 
         [a b c2]=find(remain_sect(2,:)); 
         remain_sect=[c1; c2];
         newchrom(2*j-1:2*j,:)=[remain_sect(1:2,1:points(j,1)),...
               flipud(newchrom(2*j-1:2*j,points(j,1)+1:points(j,2))),...
               remain_sect(1:2,points(j,1)+1:length(remain_sect(1,:)))];   
      end
 
      % Mutation probability
      for i=1:numsel
         if rand(1)<Pm
            points=sort((round(rand(floor(numsel/2),1).*(ngenes-1))+1)');
            newchrom(i,:)=[newchrom(i,1:points(1)),...
               fliplr(newchrom(i,points(1)+1:points(2))),...
               newchrom(i,points(2)+1:ngenes)];   
         end
      end
 
      % Creating a new population of chromosomes
      if nind-numsel, % Preserving a part of the parent chromosome population
         [ans,Index]=sort(fitness);
         chrom=[chrom(Index(numsel+1:nind),:);newchrom];
      else % Replacing the entire parent chromosome population with a new one
         chrom=newchrom;
      end
 
      % Fitness calculation
      rout=chrom;
      ObjV=evalObjFun(rout,city_distance,nind,ngenes);
   
      best=[best min(ObjV)];
      ave=[ave mean(ObjV)];
      end
 
   [a b]=min(ObjV);
    coordinates = findCoordinates(chrom,city_distance,nind,ngenes,a)
 
end
route = calculateRoute(a,b, chrom,city_distance,nind,ngenes)
 
disp("==================")
disp("minimum cost")
disp(a)
disp("==================")
disp("Best Location")
disp(coordinates(1) + " " + coordinates(2))
disp("==================")
 disp("total cost")
 disp(route)
 
 
 figure
plot(coordinates(:,1),coordinates(:,2),'.r','markersize',25)
grid on
hold on
 
% figure('name','The best rout found');
% plot(city_location(:,1),city_location(:,2),'.r','markersize',25)
% title(['The total distance: ',num2str(a)]);
% grid on
% hold on
 
 
function ObjV=evalObjFun(chrom,city_distance,nind,ngenes)
global city_location
path=0; ObjV=[];
 
for k=1:nind
   for i=1:(ngenes - 1)
     if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
     else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
     else if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),40) + 1.5;
          path=path+city_distance(41,chrom(k,(i+1)));
     else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),41) + 1.5;
          path=path+city_distance(40,chrom(k,(i+1)));
     end
    end
   end
  end
 end
   ObjV(k)=path; path=0;
end
 
function coordinates=findCoordinates(chrom,city_distance,nind,ngenes,a)
global city_location
coordinates = [];path=0;
for k=1:nind
   for i=1:(ngenes - 1)
      if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
      else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
      else if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),40) + 1.5;
          path=path+city_distance(41,chrom(k,(i+1)));
     
      else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),41) + 1.5;
          path=path+city_distance(40,chrom(k,(i+1)));
      end
     end
   end
  end
      if path == a
          coordinates = [city_location(chrom(k,i), 1) city_location(chrom(k,i), 2)]
          break
      end
   end
path=0;
end
 
function ObjV=calculateRoute(xcoord, ycoord, chrom,city_distance,nind,ngenes)
 global city_location
 ObjV = 49 * 49
 label = round(xcoord) + 7 * fix(ycoord); path=0
 for k=1:nind
   for i=1:(ngenes - 1)
    if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
    else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),chrom(k,(i+1)));
    else if mod(chrom(k,i), 7) < 6 && mod(chrom(k,i + 1), 7) >= 6
          path=path+city_distance(chrom(k,i),40) + 1.5;
          path=path+city_distance(41,chrom(k,(i+1)));
    else if mod(chrom(k,i), 7) >= 6 && mod(chrom(k,i + 1), 7) < 6
          path=path+city_distance(chrom(k,i),41) + 1.5;
          path=path+city_distance(40,chrom(k,(i+1)));
    end
   end
  end
 end
   if city_location(chrom(k,i), 1) == xcoord && city_location(chrom(k,i), 2) == ycoord && path <= ObjV
       
       ObjV = path
       break
   end
   end
   path=0;
 end