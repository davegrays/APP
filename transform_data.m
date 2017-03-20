function new_array_3D = transform_data(subject_array_3D,logtransform,ranknormalize)

[sx,~,~]=size(subject_array_3D);

%data transformations
if logtransform == 1 && ranknormalize == 0
    disp('log transforming data');
    subject_array_3D=log(subject_array_3D);
elseif ranknormalize == 1
    disp('rank-normal transforming data');
    for i=1:sx
        for j=i+1:sx
            outcome=subject_array_3D(i,j,:);outcome=outcome(:);
            rank = tiedrank(outcome);
            p = rank / ( length(rank) + 1 ); %# +1 to avoid Inf for the max point
            subject_array_3D(i,j,:) = norminv( p, 0, 1 ) * std(outcome) + mean(outcome);
        end
    end
end

new_array_3D = subject_array_3D;