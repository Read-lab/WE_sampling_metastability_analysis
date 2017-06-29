function General_trans_051117_NANOG_HPC_fullspec_f50_full_parallel(vorloops,loops2,temp_save_location,num_nodes,ReplicasRequired,name,tau,tau_name,final_save_location,vor_init,curr_iteration,species,bash_clump,rstart,full_restart,eigv)
loops2 = str2double(loops2);
num_nodes = str2double(num_nodes);
ReplicasRequired = str2double(ReplicasRequired);
tau = str2double(tau);
species = str2double(species);
bash_clump = str2double(bash_clump);
vorloops = str2double(vorloops);
tic = cputime;
curr_iteration = str2double(curr_iteration);
rstart = str2double(rstart);
full_restart = str2double(full_restart);

save_location = [final_save_location '/trans_' name '_tau' tau_name '.mat'];
temp_dir_transition = [final_save_location '/trans_' name '_tau' tau_name '/'];
if not(exist(temp_dir_transition))
    mkdir(temp_dir_transition);
end

previous_BNG_dir = [temp_save_location '/t' num2str(curr_iteration)];
load(vor_init);
load(eigv);

%%%%%%%%%%%%%%%%%%%%%%%%%%

timestep = tau/2;

matobj_save_info = matfile(save_location,'Writable',true); %saving the updated weights of the regions

%BNG_root = '/home/grad/BioNetGen-2.2.6-stable/bin/run_network'; %BionetGen simulation root
%BNG_root = '/dfs2/elread/rxn-share/BioNetGen-2.2.6-stable/bin/run_network';
% BNG_root = 'C:\Users\Maggie\Documents\BioNetGen-2.2.6-stable\bin\run_network.exe';

BioNetGen_replica_folder = temp_save_location; %Creates folder to store the temporary BioNetGen simulations
if not(exist([BioNetGen_replica_folder '/t0']))
    mkdir([BioNetGen_replica_folder '/t0']);
end


cd(temp_save_location);
if full_restart==1
trans_matrix_prob = zeros(num_nodes);
end
prob  = zeros(num_nodes,1);


    replicas_forwards = load(['curr_out_rep_' num2str(rstart) '.txt']);
    [~,sortind,~] = unique(replicas_forwards(:,end-1));
    replicas_forwards = [replicas_forwards(sortind,2:end)];

replicas_forwards(:,end) = replicas_forwards(:,end)./sum(replicas_forwards(:,end));
    

iteration = curr_iteration + 1;
if iteration > 3
iteration = 0;
end
if rstart+1 < vorloops
for ijk = rstart+1:vorloops
    
    
    
    
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    if ~exist(newBNGdir,'dir')
        mkdir(newBNGdir);
    end
    
    Binsprevx = knnsearch(VoronoiLocs(:,1:species),replicas_forwards(:,1:species));
    dlmwrite([temp_save_location '/curr_rep_ind.txt'],replicas_forwards(:,end-1));
    dlmwrite([temp_save_location '/curr_rep_data.txt'],[[1:1:length(Binsprevx)]',replicas_forwards(:,end)],'precision','%10.5e');
    weights = replicas_forwards(:,end);
unix(['rm -rf ' temp_save_location '/tempout_*.txt']);

    unix([ ' > ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
unix([ ' > ' temp_save_location '/finished.txt']);

    iters = floor(length(Binsprevx)/bash_clump);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%
    temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
    temp{17} = sprintf('OUT_NUM="%d"',ijk);
    temp{18} = sprintf('TMP="%s"',temp_save_location);
    temp{19} = sprintf('TSTEP="%f"',timestep);
    temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
    temp{22} = sprintf('NEW_ID="%d"',iteration);
    temp{24} = sprintf('STOP="%d"',bash_clump);
    temp{32} = sprintf('#$ -t 1-%d',iters);
    
    
    
    
    
    fid = fopen(['sub_BNG_scripts.sh'],'w');
    fprintf(fid,'%s \n',temp{:});
    fclose(fid);
    


    
    if rem(length(Binsprevx),bash_clump) ~=0
        temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
        temp{4}= sprintf('#$ -N temp_fin');
        temp{17} = sprintf('OUT_NUM="%d"',ijk);
        temp{18} = sprintf('TMP="%s"',temp_save_location);
        temp{19} = sprintf('TSTEP="%f"',timestep);
        temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
        temp{22} = sprintf('NEW_ID="%d"',iteration);
        temp{24} = sprintf('STOP="%d"',bash_clump);
        
        temp{32} = sprintf('#temp');
        temp{34} = sprintf('for i in {%d..%d}',iters*bash_clump+1,length(Binsprevx));
        temp{36} = sprintf('	CURRLINE="${i}"');
        temp{46} = sprintf('echo 1 > ${TMP}/finished2.txt');
        
        
        fid = fopen(['sub_BNG_scripts_fin.sh'],'w');
        fprintf(fid,'%s \n',temp{:});
        fclose(fid);
        unix(['qsub sub_BNG_scripts_fin.sh']);
        
    end
    
    unix(['qsub sub_BNG_scripts.sh']);
     [~,Is_done] = unix(['wc -l < finished.txt']);

  while str2num(Is_done) < iters
     pause(10)
     [~,Is_done] = unix(['wc -l < finished.txt']);
     disp(str2num(Is_done));
  end


unix(['rm -rf ' temp_save_location '/finished.txt']);

if  rem(length(Binsprevx),bash_clump) ~=0
  
  while ~exist([temp_save_location '/finished2.txt']);
     pause(10)
     [~,Is_done] = unix(['wc -l < curr_out_rep_' num2str(ijk) '.txt']);
     disp(str2num(Is_done));
  end
unix(['rm -rf ' temp_save_location '/finished2.txt']);

end


unix(['cat ' temp_save_location '/tempout_*.txt >> ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
unix(['rm -rf ' temp_save_location '/finished.txt']);


    replicas_new = load(['curr_out_rep_' num2str(ijk) '.txt']);
    
    [~,sortind,~] = unique(replicas_new(:,end-1));
    replicas_new = [replicas_new(sortind,2:end-1),weights];
    curr_iteration = iteration;
    iteration = iteration+1;
    if iteration > 3
        iteration = 0;
        unix('rm -rf temp.o*');
        unix('rm -rf temp_fin.o*');
        
    end
    
    newVoronoi = zeros(num_nodes,species);
disp(size(replicas_new))
disp(size(newVoronoi))
    new_ind_voronoi = randsample(length(replicas_new(:,1)),1);
disp(size(new_ind_voronoi))
    newVoronoi(1,:) = replicas_new(new_ind_voronoi,1:species);
    for i = 2:num_nodes
        [~,Distances] = knnsearch(newVoronoi(1:i-1,:),replicas_new(:,1:species));
        [~,idxmax] = max(Distances);
        newVoronoi(i,:) = replicas_new(idxmax,1:species);
    end
    
    %
    VoronoiLocs = newVoronoi;
    
    Binsprev = knnsearch(VoronoiLocs,replicas_new(:,1:species));
    
    
    replicas_forwards = WEstep_051117(replicas_new,VoronoiLocs,ReplicasRequired,Binsprev);
    
    toc = cputime-tic;
    display([ijk,toc]);
    
    matobj_save_info.replicas_forwards = replicas_forwards;
    matobj_save_info.VoronoiLocs = VoronoiLocs;
    matobj_save_info.iteration = iteration;
    save([save_location ],'VoronoiLocs','replicas_forwards','iteration','toc','-v7.3')
    
    
    if rem(ijk,2)==0
        dlmwrite([temp_dir_transition 'voronoi' num2str(ijk) '.txt'],VoronoiLocs );
        
    end
    
    
    
    
    
    
    
    
    
    
end
end


for ijk = rstart+1:loops2
    
    
    
    
    newBNGdir = [BioNetGen_replica_folder '/t' num2str(iteration)];
    if ~exist(newBNGdir,'dir')
        mkdir(newBNGdir);
    end
    
    Binsprevx = knnsearch(VoronoiLocs,replicas_forwards(:,1:species));
    dlmwrite([temp_save_location '/curr_rep_ind.txt'],replicas_forwards(:,end-1));
    dlmwrite([temp_save_location '/curr_rep_data.txt'],[[1:1:length(Binsprevx)]',replicas_forwards(:,end)],'precision','%10.5e');
    weights = replicas_forwards(:,end);
    unix([ ' > ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);
    unix([ ' > ' temp_save_location '/finished.txt']);

    iters = floor(length(Binsprevx)/bash_clump);
  unix(['rm -rf ' temp_save_location '/tempout_*.txt']);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
    temp{17} = sprintf('OUT_NUM="%d"',ijk);
    temp{18} = sprintf('TMP="%s"',temp_save_location);
    temp{19} = sprintf('TSTEP="%f"',timestep);
    temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
    temp{22} = sprintf('NEW_ID="%d"',iteration);
    temp{24} = sprintf('STOP="%d"',bash_clump);
    temp{32} = sprintf('#$ -t 1-%d',iters);
    
    
    
    
    
    fid = fopen(['sub_BNG_scripts.sh'],'w');
    fprintf(fid,'%s \n',temp{:});
    fclose(fid);
    

    if rem(length(Binsprevx),bash_clump) ~=0
        temp= regexp(fileread('sub_BNG_batch_template.sh'),'\n','split');
        temp{4}= sprintf('#$ -N temp_fin');
        temp{17} = sprintf('OUT_NUM="%d"',ijk);
        temp{18} = sprintf('TMP="%s"',temp_save_location);
        temp{19} = sprintf('TSTEP="%f"',timestep);
        temp{21} = sprintf('PREV_ID="%d"',curr_iteration);
        temp{22} = sprintf('NEW_ID="%d"',iteration);
        temp{24} = sprintf('STOP="%d"',bash_clump);
        
        temp{32} = sprintf('#temp');
        temp{34} = sprintf('for i in {%d..%d}',iters*bash_clump+1,length(Binsprevx));
        temp{36} = sprintf('	CURRLINE="${i}"');
          temp{46} = sprintf('echo 1 > ${TMP}/finished2.txt');
      
        
        fid = fopen(['sub_BNG_scripts_fin.sh'],'w');
        fprintf(fid,'%s \n',temp{:});
        fclose(fid);
        unix(['qsub sub_BNG_scripts_fin.sh']);
        
    end
    
    unix(['qsub sub_BNG_scripts.sh']);
    
     [~,Is_done] = unix(['wc -l < finished.txt']);

  while str2num(Is_done) < iters
     pause(10)
     [~,Is_done] = unix(['wc -l < finished.txt']);
     disp(str2num(Is_done));
  end


unix(['rm -rf ' temp_save_location '/finished.txt']);
if  rem(length(Binsprevx),bash_clump) ~=0
  
  while ~exist([temp_save_location '/finished2.txt']);
     pause(10)
     [~,Is_done] = unix(['wc -l < curr_out_rep_' num2str(ijk) '.txt']);
     disp(str2num(Is_done));
  end
unix(['rm -rf ' temp_save_location '/finished2.txt']);

end


unix(['rm -rf ' temp_save_location '/finished.txt']);
unix(['cat ' temp_save_location '/tempout_*.txt >> ' temp_save_location '/curr_out_rep_' num2str(ijk) '.txt']);

    replicas_new = load(['curr_out_rep_' num2str(ijk) '.txt']);
    
    
    [~,sortind,~] = unique(replicas_new(:,end-1));
    replicas_new = [replicas_new(sortind,2:end-1),weights];
    curr_iteration = iteration;
    iteration = iteration+1;
    if iteration > 3
        iteration = 0;
        unix('rm -rf temp.o*');
        unix('rm -rf temp_fin.o*');
        
    end
    
    
    
    Binsprev = knnsearch(VoronoiLocs,replicas_new(:,1:species));
    
    for i = 1:(num_nodes)
        prob(i) = prob(i) + sum(replicas_new((Binsprev == i),end));
    end
    
    
    
    count_rows = [Binsprevx,Binsprev ];
    temp_urows = unique(count_rows,'rows');
    for irows = 1:length(temp_urows(:,1))
        outof = temp_urows(irows,1);
        into = temp_urows(irows,2);
        times_visited = (ismember(count_rows,temp_urows(irows,:),'rows'));
        trans_matrix_prob(sub2ind(size(trans_matrix_prob),outof,into)) = trans_matrix_prob(sub2ind(size(trans_matrix_prob),outof,into))+sum(replicas_new(times_visited,end));
        
    end
    
    
    replicas_forwards = WEstep_051117(replicas_new,VoronoiLocs,ReplicasRequired,Binsprev);
    
    toc = cputime-tic;
    display([ijk,toc]);
    
    matobj_save_info.prob = prob;
    matobj_save_info.replicas_forwards = replicas_forwards;
    matobj_save_info.VoronoiLocs = VoronoiLocs;
    matobj_save_info.iteration = iteration;
    matobj_save_info.trans_matrix_prob = trans_matrix_prob;
    save([save_location ],'VoronoiLocs','replicas_forwards','trans_matrix_prob','iteration','toc','-v7.3')
    
    
    if rem(ijk,2)==0
        dlmwrite([temp_dir_transition 'trans_matrix_prob' num2str(ijk) '.txt'],trans_matrix_prob );
        
    end
    
    
    
    
    
    
    
    
    
    
end


toc = cputime-tic;
display(toc);
save([save_location ],'VoronoiLocs','replicas_forwards','trans_matrix_prob','iteration','toc','-v7.3')


end


function [newreps] = WEstep_051117(replicas,Voronoi_List,numreps,bins)
newreps = [];
%This code assumes the columns of the replicas matrix are the species, and that the weights are in the last column
for i = 1:length(Voronoi_List(:,1))
    replicas_in_bini = replicas(bins==i,:);
    if not(isempty(replicas_in_bini))
        newreps_temp = WEstep(replicas_in_bini,numreps);
        if  abs(sum(replicas_in_bini(:,end)) -   sum(newreps_temp(:,end)) > 10^(-12))
            sum(newreps_temp(:,end))
            sum(replicas_in_bini(:,end))
            error('WEstep failure');
        end
        newreps = [newreps; newreps_temp];
        
    end
end


end



function newreps = WEstep(curr_reps,numreps)  %%%%%error in weight calc

blah = length(curr_reps(:,1));
if isempty(blah)
    error('reps_empty!')
end
%%%%%number of repilcas incorred incremented by 1
newreps = curr_reps;

curr_weights = newreps(:,end);
total_weight = sum(curr_weights);
idx_remove = [];
[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);

%%removes smallest weight outlier
while (smallest_weight < total_weight/(3*numreps)) && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    
    min_idx = min_idx+1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(idx_remove,:) = [];
end

w2 = sum(newreps(:,end));

curr_weights =  newreps(:,end);
%%removes largest weight outlier
[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx = 0;
max_factor = 1;

if largest_weight > total_weight/numreps*3
    
    while largest_weight > total_weight/numreps*3
        largest_weight = tmpweight/max_factor;
        max_factor = max_factor+1;
    end
    
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    new_reps_temp(:,end) = repmat(tmpweight/(max_factor+mxx),max_factor+mxx,1);
    newreps = [newreps; new_reps_temp];
    newreps(max_idx,:) = [];
end


curr_weights = newreps(:,end);
curr_len = length(curr_weights);
w3 = sum(newreps(:,end));

%%Merges weights until the number of replicas = numreps

[SortedWeights, SortedInd] = sort(curr_weights);
min_idx = 2;
max_idx = length(SortedInd);
smallest_weight = SortedWeights(1);
idx_remove = [];

while curr_len  > numreps && (min_idx < max_idx)
    smallest_weight = sum(SortedWeights(1:min_idx));
    idx_remove = [ SortedInd(1:min_idx)];
    min_idx = min_idx+1;
    curr_len = curr_len-1;
end
if isempty(idx_remove) == 0
    new_rep_ind = randsample(SortedInd(1:min_idx),1,true,SortedWeights(1:min_idx));
    
    new_reps_temp = newreps(new_rep_ind,:);
    new_reps_temp(:,end) = smallest_weight;
    newreps = [newreps; new_reps_temp];
    newreps(sort(idx_remove),:) = [];
    
end



w4 = sum(newreps(:,end));

%%splits weights until the number of replicas = numreps

curr_weights =  newreps(:,end);
curr_len = length(curr_weights);

[largest_weight, max_idx] = max(curr_weights);
tmpweight = largest_weight;
mxx_fact2 = 0;

max_factor = 1;


if curr_len < numreps
    while curr_len < numreps
        max_factor = max_factor+1;
        curr_len = curr_len+1;
    end
    new_reps_temp = repmat(newreps(max_idx,:),max_factor+mxx,1);
    a = length(new_reps_temp(:,1));
    
    new_reps_temp(:,end) = tmpweight/(max_factor+mxx_fact2).*ones(a,1);
    newreps = [newreps; new_reps_temp];
    shouldbezero = (length(newreps(:,end))-1-numreps);
    
    
    newreps(max_idx,:) = [];
    
    
    shouldbezero2 = (length(newreps(:,end))-numreps);
    
    
    if shouldbezero2~=0
        error(' replica number from WEstep after splitting is wrong');
    end
    
    
end

w5 = sum(newreps(:,end));


%%Checking for errors in the weight calculation.

if (all(abs([total_weight-w2,w2-w3,w3-w4,w4-w5]))) > 10^(-12)
    abs([total_weight-w2,w2-w3,w3-w4,w4-w5])
    error('Error is in one of the three subsections of WE step')
elseif length(newreps(:,end)) ~= numreps
    error('Final replica length is wrong');
end


end





