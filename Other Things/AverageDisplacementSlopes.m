%This script averages the slope of linear fits to the displacements of 
%elements in the middle of shear zone trials for batches of trials. It will
%create a column vector with the average slope of the trials in each row of
%the files array, and the column vector is called avgslope.

%trial names to be analyzed. ROWS of the cell array should contain the
%files that you want averaged
files={'friction_p1_random','friction_p1_random_2','friction_p1_random_3';
    'friction_p2_random','friction_p2_random_2','friction_p2_random_3';
    'friction_p3_random','friction_p3_random_2','friction_p3_random_3';
    'friction_p4_random','friction_p4_random_2','friction_p4_random_3';
    'friction_p5_random','friction_p5_random_2','friction_p5_random_3';
    'friction_p75_random','friction_p75_random_2','friction_p75_random_3';
    'friction_p9_random','friction_p9_random_2','friction_p9_random_3';
    'friction_1_random','friction_1_random_2','friction_1_random_3'};

bandlength=2;
    
P=cell(size(files,1),size(files,2));
slopes=zeros(size(files,1),size(files,2));

for i=1:size(files,1)
    for j=1:size(files,2)
        
        trial=load(files{i,j});

        L=length(trial.Estore);
        N=length(trial.Estore(1).x);

        %find approximate middle vertical point of initial array
        middle=(max(trial.Estore(1).y)+min(trial.Estore(1).y))/2;

        %define the size of the band of included elements
        verticalband=[middle+mean(trial.r)*bandlength,...
            middle-mean(trial.r)*bandlength];

        %find elements initially within the band
        idx=zeros(1,length(trial.Estore(1).x));
        count=1;
        for k=1:N
            if(trial.Estore(1).y(k)<verticalband(1) &&...
                    trial.Estore(1).y(k)>verticalband(2))
                idx(count)=k;
                count=count+1;
            end
        end
        idx(idx==0)=[];

        %bubble sort element indices by horizontal locations
        l=length(idx);
        done=0;
        while(~done)
            done=1;
            for k=1:l-1
                if(trial.Estore(1).x(idx(k))>trial.Estore(1).x(idx(k+1)))
                    temp=idx(k);
                    idx(k)=idx(k+1);
                    idx(k+1)=temp;
                    done=0;
                end
            end
        end
    
        l=trial.Estore(1).x(idx);
        displacement=(trial.Estore(L).y(idx)-trial.Estore(1).y(idx));

        P{i,j}=polyfit(l,displacement,1);
        slopes(i,j)=P{i,j}(1);
        
    end
end

avgslope=zeros(size(slopes,1),1);
for i=1:size(slopes,1)
    avgslope(i)=mean(slopes(i,:));
end